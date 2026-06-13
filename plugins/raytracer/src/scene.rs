//! Scene raytracing orchestration.
//!
//! Collects primitives from the viewer, builds a BVH, resolves settings,
//! extracts camera matrices, and dispatches GPU raytracing.

use std::sync::OnceLock;

use patinae_color::NamedPalette;
use patinae_plugin::prelude::*;
use patinae_render::{
    DisplayedGeometry, DisplayedMaterial, DisplayedMeshVertex, DisplayedPrimitive,
    GeometryExportOptions,
};
use patinae_scene::{normalize_matrix, Object};
use patinae_settings::ResolvedSettings;

/// Reflect-scale adjustment for classic shading so multiple positional lights
/// do not pile up.
fn compute_reflect_scale(light_count: i32, light_dirs: &[[f32; 3]]) -> f32 {
    let num_pos_lights = (light_count - 1).max(0) as usize;
    if num_pos_lights == 0 {
        return 1.0;
    }
    let mut sum = 0.0f32;
    for dir in light_dirs.iter().take(num_pos_lights).copied() {
        let len = (dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]).sqrt();
        if len > 0.0 {
            sum += 1.0 - dir[2] / len;
        }
    }
    sum *= 0.5;
    if sum > 0.0 {
        1.0 / sum
    } else {
        1.0
    }
}

use crate::bvh::Bvh;
use crate::collect::{collect_from_molecule, CollectOptions, RayColorResolver};
use crate::gpu::{raytrace, RaytraceParams};
use crate::primitive::{GpuCylinder, GpuSphere, GpuTriangle, PrimitiveCollector, Primitives};
use crate::settings::{RaytraceSettings, ResolvedRaySettings};

const RAY_LINE_RADIUS: f32 = 0.035;
const RAY_POINT_RADIUS: f32 = 0.035;

static RAY_GPU_RUNTIME: OnceLock<RayGpuRuntime> = OnceLock::new();

// ---------------------------------------------------------------------------
// Error
// ---------------------------------------------------------------------------

#[derive(Debug, thiserror::Error)]
pub(crate) enum RaytraceSceneError {
    #[error("No primitives to raytrace")]
    NoPrimitives,
    #[error("BVH construction failed: {0}")]
    BvhFailed(String),
    #[error("Displayed geometry export failed: {0}")]
    GeometryExport(String),
    #[error("Raytracing failed: {0}")]
    RaytraceFailed(String),
}

struct RayGpuRuntime {
    device: wgpu::Device,
    queue: wgpu::Queue,
}

// ---------------------------------------------------------------------------
// Shared helpers
// ---------------------------------------------------------------------------

/// Collect all raytracing primitives from the object registry and build a BVH.
fn collect_scene_primitives(
    viewer: &mut dyn ViewerLike,
) -> Result<(Primitives, Bvh), RaytraceSceneError> {
    let mut export_error = None;
    let primitives = match viewer.export_displayed_geometry(&GeometryExportOptions::default()) {
        Ok(displayed) => primitives_from_displayed_geometry(&displayed),
        Err(err) => {
            export_error = Some(err);
            collect_scene_primitives_from_molecules(viewer)
        }
    };

    if primitives.is_empty() {
        if let Some(err) = export_error {
            return Err(RaytraceSceneError::GeometryExport(format!(
                "{err}; molecular fallback found no primitives"
            )));
        }
        return Err(RaytraceSceneError::NoPrimitives);
    }

    let bvh = Bvh::build(&primitives).map_err(|e| RaytraceSceneError::BvhFailed(e.to_string()))?;

    Ok((primitives, bvh))
}

fn collect_scene_primitives_from_molecules(viewer: &dyn ViewerLike) -> Primitives {
    let mut collector = PrimitiveCollector::new();
    let colors = RayColorResolver::new(viewer.named_palette(), &viewer.session().palette);

    for name in viewer.objects().names() {
        let Some(molecule_object) = viewer.objects().get_molecule(name) else {
            continue;
        };
        if !molecule_object.state().enabled {
            continue;
        }
        let Some(coord_set) = molecule_object.display_coord_set() else {
            continue;
        };

        let resolved = ResolvedSettings::resolve(viewer.settings(), molecule_object.overrides());
        let options = CollectOptions {
            sphere_scale: resolved.sphere.scale,
            stick_radius: resolved.stick.radius,
            sphere_transparency: resolved.sphere.transparency,
            stick_transparency: resolved.stick.transparency,
            collect_spheres: molecule_object
                .state()
                .visible_reps
                .is_visible(patinae_mol::RepMask::SPHERES),
            collect_sticks: molecule_object
                .state()
                .visible_reps
                .is_visible(patinae_mol::RepMask::STICKS),
            collect_cartoon: false,
            collect_surface: false,
        };
        let primitives =
            collect_from_molecule(molecule_object.molecule(), coord_set, &colors, &options);
        collector.add_spheres(primitives.spheres);
        collector.add_cylinders(primitives.cylinders);
        collector.add_triangles(primitives.triangles);
    }

    collector.build()
}

fn ray_gpu_runtime() -> Result<&'static RayGpuRuntime, RaytraceSceneError> {
    if let Some(runtime) = RAY_GPU_RUNTIME.get() {
        return Ok(runtime);
    }

    let runtime = create_ray_gpu_runtime()?;
    let _ = RAY_GPU_RUNTIME.set(runtime);
    RAY_GPU_RUNTIME.get().ok_or_else(|| {
        RaytraceSceneError::RaytraceFailed("plugin GPU runtime was not initialized".into())
    })
}

fn create_ray_gpu_runtime() -> Result<RayGpuRuntime, RaytraceSceneError> {
    let instance = wgpu::Instance::new(&wgpu::InstanceDescriptor::default());
    let adapter = pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
        power_preference: wgpu::PowerPreference::HighPerformance,
        force_fallback_adapter: false,
        compatible_surface: None,
    }))
    .map_err(|err| RaytraceSceneError::RaytraceFailed(format!("GPU adapter unavailable: {err}")))?;

    let (device, queue) = pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor {
        label: Some("raytracer.plugin.device"),
        required_features: wgpu::Features::empty(),
        required_limits: adapter.limits(),
        memory_hints: wgpu::MemoryHints::Performance,
        experimental_features: wgpu::ExperimentalFeatures::default(),
        trace: wgpu::Trace::Off,
    }))
    .map_err(|err| RaytraceSceneError::RaytraceFailed(format!("GPU device unavailable: {err}")))?;

    Ok(RayGpuRuntime { device, queue })
}

fn primitives_from_displayed_geometry(displayed: &DisplayedGeometry) -> Primitives {
    let mut collector = PrimitiveCollector::new();
    for object in &displayed.objects {
        for primitive in &object.primitives {
            match primitive {
                DisplayedPrimitive::Mesh { mesh, .. } => {
                    for tri in mesh.vertices.chunks_exact(3) {
                        collector
                            .add_triangle(ray_triangle_from_vertices(&tri[0], &tri[1], &tri[2]));
                    }
                }
                DisplayedPrimitive::Sphere {
                    center,
                    radius,
                    material,
                    ..
                } => {
                    collector.add_sphere(GpuSphere::new(
                        *center,
                        *radius,
                        material.rgba,
                        material.transparency,
                    ));
                }
                DisplayedPrimitive::Cylinder {
                    start,
                    end,
                    radius,
                    material_start,
                    material_end,
                    ..
                } => {
                    collector.add_cylinder(GpuCylinder::new(
                        *start,
                        *end,
                        *radius,
                        material_start.rgba,
                        material_end.rgba,
                        material_start.transparency.max(material_end.transparency),
                    ));
                }
                DisplayedPrimitive::LineSegment {
                    start,
                    end,
                    material_start,
                    material_end,
                    ..
                } => {
                    // Screen-space lines have no physical radius. The ray
                    // adapter renders them as a thin world-space cylinder so
                    // line/mesh displays still produce a visible ray image.
                    collector.add_cylinder(GpuCylinder::new(
                        *start,
                        *end,
                        RAY_LINE_RADIUS,
                        material_start.rgba,
                        material_end.rgba,
                        material_start.transparency.max(material_end.transparency),
                    ));
                }
                DisplayedPrimitive::PointSample {
                    position, material, ..
                } => {
                    // Screen-space points are semantic samples; ray uses a
                    // small physical sphere approximation.
                    collector.add_sphere(GpuSphere::new(
                        *position,
                        RAY_POINT_RADIUS,
                        material.rgba,
                        material.transparency,
                    ));
                }
            }
        }
    }
    collector.build()
}

fn ray_triangle_from_vertices(
    a: &DisplayedMeshVertex,
    b: &DisplayedMeshVertex,
    c: &DisplayedMeshVertex,
) -> GpuTriangle {
    let material = triangle_material(a.material, b.material, c.material);
    GpuTriangle::new(
        a.position,
        b.position,
        c.position,
        a.normal,
        b.normal,
        c.normal,
        material.rgba,
        material.transparency,
    )
}

fn triangle_material(
    a: DisplayedMaterial,
    b: DisplayedMaterial,
    c: DisplayedMaterial,
) -> DisplayedMaterial {
    let rgba = [
        (a.rgba[0] + b.rgba[0] + c.rgba[0]) / 3.0,
        (a.rgba[1] + b.rgba[1] + c.rgba[1]) / 3.0,
        (a.rgba[2] + b.rgba[2] + c.rgba[2]) / 3.0,
        (a.rgba[3] + b.rgba[3] + c.rgba[3]) / 3.0,
    ];
    DisplayedMaterial {
        base_rgba: rgba,
        rep_rgba: rgba,
        rgba,
        transparency: a.transparency.max(b.transparency).max(c.transparency),
    }
}

/// Resolve all raytracing settings from the app settings and plugin settings.
fn resolve_raytrace_settings(
    settings: &patinae_settings::Settings,
    named_palette: &NamedPalette,
    clear_color: [f32; 3],
    ray: &ResolvedRaySettings,
) -> RaytraceSettings {
    let classic = &settings.shading.classic;

    let (ambient, direct, reflect, specular, shininess) = if ray.use_custom {
        (
            ray.custom_ambient,
            ray.custom_direct,
            ray.custom_reflect,
            ray.custom_specular,
            ray.custom_shininess,
        )
    } else {
        (
            classic.ambient,
            classic.direct,
            classic.reflect,
            classic.specular,
            classic.shininess,
        )
    };

    let ray_trace_color = if ray.color >= 0 {
        named_palette
            .get_by_index(ray.color as u32)
            .map(|c| c.to_rgba(1.0))
            .unwrap_or([0.0, 0.0, 0.0, 1.0])
    } else {
        [0.0, 0.0, 0.0, 1.0]
    };

    let light_count = match settings.shading.mode {
        patinae_settings::ShadingMode::Skripkin => 1,
        patinae_settings::ShadingMode::Classic | patinae_settings::ShadingMode::Full => {
            classic.light_count
        }
    };

    let light_dirs_3: Vec<[f32; 3]> = vec![
        classic.light,
        classic.light2,
        classic.light3,
        classic.light4,
        classic.light5,
        classic.light6,
        classic.light7,
        classic.light8,
        classic.light9,
    ];

    let mut light_dirs = [[0.0f32; 4]; 9];
    for (i, dir) in light_dirs_3.iter().enumerate() {
        light_dirs[i] = [dir[0], dir[1], dir[2], 0.0];
    }

    let reflect_scale = compute_reflect_scale(light_count, &light_dirs_3);
    let reflect_scaled = reflect * reflect_scale;
    let (direct_adjusted, reflect_final) = if light_count < 2 {
        ((direct + reflect_scaled).min(1.0), 0.0)
    } else {
        (direct, reflect_scaled)
    };

    RaytraceSettings {
        light_dirs,
        light_count,
        spec_count: classic.spec_count,
        ambient,
        direct: direct_adjusted,
        reflect: reflect_final,
        specular,
        shininess,
        bg_color: [clear_color[0], clear_color[1], clear_color[2], 1.0],
        fog_start: settings.shading.common.fog_start,
        fog_end: 1.0,
        fog_density: settings.shading.common.fog,
        fog_color: [clear_color[0], clear_color[1], clear_color[2], 1.0],
        ray_shadow: ray.shadow,
        ray_max_passes: ray.max_passes as u32,
        ray_trace_fog: ray.fog > 0.0,
        ray_transparency_shadows: ray.transparency_shadows,
        ray_trace_mode: ray.mode,
        ray_trace_color,
        ray_opaque_background: match ray.opaque_background {
            -1 => settings.ui.opaque_background,
            0 => false,
            _ => true,
        },
        transparency_mode: settings.surface.transparency_mode,
        ray_trace_slope_factor: ray.slope_factor,
        ray_trace_depth_factor: ray.depth_factor,
        ray_trace_disco_factor: ray.disco_factor,
        ray_trace_gain: ray.gain,
        silhouette_thickness: settings.shading.common.silhouette_width,
        silhouette_depth_jump: settings.shading.common.silhouette_depth_jump,
    }
}

/// Camera matrices needed by the raytracing pipeline.
struct CameraMatrices {
    view: [[f32; 4]; 4],
    proj: [[f32; 4]; 4],
    view_inv: [[f32; 4]; 4],
    proj_inv: [[f32; 4]; 4],
    camera_pos: [f32; 4],
}

/// Compute camera matrices from a view matrix and a projection matrix.
fn extract_camera_matrices(
    view_matrix: lin_alg::f32::Mat4,
    proj_matrix: lin_alg::f32::Mat4,
) -> CameraMatrices {
    let mut view_inv = view_matrix
        .inverse()
        .unwrap_or(lin_alg::f32::Mat4::new_identity());
    normalize_matrix(&mut view_inv);
    let camera_pos = [view_inv.data[12], view_inv.data[13], view_inv.data[14], 1.0];

    let mut proj_inv = proj_matrix
        .inverse()
        .unwrap_or(lin_alg::f32::Mat4::new_identity());
    normalize_matrix(&mut proj_inv);

    CameraMatrices {
        view: matrix_to_array(view_matrix),
        proj: matrix_to_array(proj_matrix),
        view_inv: matrix_to_array(view_inv),
        proj_inv: matrix_to_array(proj_inv),
        camera_pos,
    }
}

// ---------------------------------------------------------------------------
// Scene raytracing (mutable viewer — sets camera aspect temporarily)
// ---------------------------------------------------------------------------

pub(crate) fn raytrace_scene(
    viewer: &mut dyn ViewerLike,
    ray: &ResolvedRaySettings,
    width: Option<u32>,
    height: Option<u32>,
    antialias: u32,
) -> Result<(Vec<u8>, u32, u32), RaytraceSceneError> {
    let final_width = width.unwrap_or(1024);
    let final_height = height.unwrap_or(768);

    let (primitives, bvh) = collect_scene_primitives(viewer)?;

    let rt_settings = resolve_raytrace_settings(
        viewer.settings(),
        viewer.named_palette(),
        viewer.clear_color(),
        ray,
    );

    // Camera matrices — temporarily set aspect ratio for the output dimensions
    let camera = viewer.camera_mut();
    let original_aspect = camera.aspect();
    camera.set_aspect(final_width as f32 / final_height as f32);
    let cam = extract_camera_matrices(camera.view_matrix(), camera.projection_matrix());
    camera.set_aspect(original_aspect);

    let gpu = ray_gpu_runtime()?;

    let params = RaytraceParams::new(final_width, final_height)
        .with_antialias(antialias)
        .with_camera(
            cam.view,
            cam.proj,
            cam.view_inv,
            cam.proj_inv,
            cam.camera_pos,
        )
        .with_settings(rt_settings);

    let image_data = raytrace(&gpu.device, &gpu.queue, &primitives, &bvh, &params)
        .map_err(|e| RaytraceSceneError::RaytraceFailed(e.to_string()))?;

    Ok((image_data, final_width, final_height))
}

// ---------------------------------------------------------------------------
// Helper
// ---------------------------------------------------------------------------

#[inline]
fn matrix_to_array(m: lin_alg::f32::Mat4) -> [[f32; 4]; 4] {
    <[[f32; 4]; 4]>::from(m)
}

#[cfg(test)]
mod tests {
    use super::*;
    use patinae_render::{DisplayedMesh, DisplayedObjectGeometry, ObjectId, RepKind};

    #[test]
    fn displayed_cartoon_mesh_converts_to_ray_triangle() {
        let material = DisplayedMaterial::from_rgba([0.8, 0.1, 0.2, 1.0]);
        let displayed = DisplayedGeometry {
            objects: vec![DisplayedObjectGeometry {
                object_id: ObjectId(1),
                primitives: vec![DisplayedPrimitive::Mesh {
                    rep: RepKind::Cartoon,
                    mesh: DisplayedMesh {
                        vertices: vec![
                            DisplayedMeshVertex {
                                position: [0.0, 0.0, 0.0],
                                normal: [0.0, 0.0, 1.0],
                                owner_atom_id: 0,
                                material,
                                flags: 0,
                            },
                            DisplayedMeshVertex {
                                position: [1.0, 0.0, 0.0],
                                normal: [0.0, 0.0, 1.0],
                                owner_atom_id: 0,
                                material,
                                flags: 0,
                            },
                            DisplayedMeshVertex {
                                position: [0.0, 1.0, 0.0],
                                normal: [0.0, 0.0, 1.0],
                                owner_atom_id: 0,
                                material,
                                flags: 0,
                            },
                        ],
                    },
                }],
            }],
        };

        let primitives = primitives_from_displayed_geometry(&displayed);

        assert_eq!(primitives.total_count(), 1);
        assert_eq!(primitives.triangles.len(), 1);
    }

    #[test]
    #[ignore = "GPU smoke; run explicitly when validating ray runtime"]
    fn ray_smoke_writes_nonempty_png_from_displayed_geometry() {
        let gpu = ray_gpu_runtime().expect("ray smoke runtime");

        let material = DisplayedMaterial::from_rgba([0.2, 0.8, 0.4, 1.0]);
        let displayed = DisplayedGeometry {
            objects: vec![DisplayedObjectGeometry {
                object_id: ObjectId(1),
                primitives: vec![DisplayedPrimitive::Sphere {
                    rep: RepKind::Sphere,
                    owner_atom_id: 0,
                    center: [0.0, 0.0, -5.0],
                    radius: 1.2,
                    material,
                }],
            }],
        };
        let primitives = primitives_from_displayed_geometry(&displayed);
        let bvh = Bvh::build(&primitives).expect("ray smoke bvh");
        let params = RaytraceParams::new(64, 64);

        let image_data = raytrace(&gpu.device, &gpu.queue, &primitives, &bvh, &params)
            .expect("ray smoke render");
        assert!(image_data.iter().any(|&byte| byte != 0));

        let path = std::path::Path::new("/private/tmp/patinae-ray-export-smoke.png");
        image::save_buffer(path, &image_data, 64, 64, image::ColorType::Rgba8)
            .expect("save ray smoke png");
        assert!(std::fs::metadata(path).unwrap().len() > 0);
    }
}
