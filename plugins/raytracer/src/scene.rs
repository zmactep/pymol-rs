//! Scene raytracing orchestration.
//!
//! Collects primitives from the viewer, builds a BVH, resolves settings,
//! extracts camera matrices, and dispatches GPU raytracing.

use pymol_color::{ChainColors, ElementColors, NamedColors};
use pymol_plugin::prelude::*;
use pymol_render::MeshVertex;
use pymol_scene::{compute_reflect_scale, normalize_matrix, MoleculeObject, Object, ObjectRegistry};

use crate::bvh::Bvh;
use crate::collect::{collect_from_molecule, CollectOptions, RayColorResolver};
use crate::gpu::{raytrace, RaytraceParams};
use crate::primitive::{GpuCylinder, GpuTriangle, Primitives, PrimitiveCollector};
use crate::settings::{RaytraceSettings, ResolvedRaySettings};

// ---------------------------------------------------------------------------
// Error
// ---------------------------------------------------------------------------

#[derive(Debug, thiserror::Error)]
pub(crate) enum RaytraceSceneError {
    #[error("No primitives to raytrace")]
    NoPrimitives,
    #[error("BVH construction failed: {0}")]
    BvhFailed(String),
    #[error("Raytracing failed: {0}")]
    RaytraceFailed(String),
}

// ---------------------------------------------------------------------------
// Shared helpers
// ---------------------------------------------------------------------------

/// Collect all raytracing primitives from the object registry and build a BVH.
fn collect_scene_primitives(
    registry: &ObjectRegistry,
    named_colors: &NamedColors,
    settings: &pymol_settings::Settings,
) -> Result<(Primitives, Bvh), RaytraceSceneError> {
    let element_colors = ElementColors::default();
    let chain_colors = ChainColors;
    let color_resolver = RayColorResolver::new(named_colors, &element_colors, &chain_colors);

    let options = CollectOptions {
        sphere_scale: settings.sphere.scale,
        stick_radius: settings.stick.radius,
        ..Default::default()
    };

    let mut collector = PrimitiveCollector::new();
    let names: Vec<String> = registry.names().map(|s| s.to_string()).collect();

    for name in &names {
        if let Some(mol_obj) = registry.get_molecule(name) {
            if !mol_obj.is_enabled() {
                continue;
            }

            let molecule = mol_obj.molecule();
            if let Some(coord_set) = molecule.get_coord_set(mol_obj.display_state()) {
                let prims =
                    collect_from_molecule(molecule, coord_set, &color_resolver, &options);
                collector.add_spheres(prims.spheres);
                collector.add_cylinders(prims.cylinders);
                collector.add_triangles(prims.triangles);
            }

            collector.add_triangles(collect_triangles_from_mesh_reps(mol_obj));

            if let Some(edges) = mol_obj.get_mesh_edges() {
                let mesh_radius = 0.01;
                let cylinders: Vec<GpuCylinder> = edges
                    .chunks(2)
                    .filter(|pair| pair.len() == 2)
                    .map(|pair| {
                        GpuCylinder::new(
                            pair[0].position,
                            pair[1].position,
                            mesh_radius,
                            pair[0].color,
                            pair[1].color,
                            0.0,
                        )
                    })
                    .collect();
                collector.add_cylinders(cylinders);
            }
        }
    }

    let primitives = collector.build();
    if primitives.is_empty() {
        return Err(RaytraceSceneError::NoPrimitives);
    }

    let bvh =
        Bvh::build(&primitives).map_err(|e| RaytraceSceneError::BvhFailed(e.to_string()))?;

    Ok((primitives, bvh))
}

/// Resolve all raytracing settings from the app settings and plugin settings.
fn resolve_raytrace_settings(
    settings: &pymol_settings::Settings,
    named_colors: &NamedColors,
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
        named_colors
            .get_by_index(ray.color as u32)
            .map(|c| c.to_rgba(1.0))
            .unwrap_or([0.0, 0.0, 0.0, 1.0])
    } else {
        [0.0, 0.0, 0.0, 1.0]
    };

    let light_count = match settings.shading.mode {
        pymol_settings::ShadingMode::Skripkin => 1,
        pymol_settings::ShadingMode::Classic | pymol_settings::ShadingMode::Full => {
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

    let (primitives, bvh) =
        collect_scene_primitives(viewer.objects(), viewer.named_colors(), viewer.settings())?;

    let rt_settings = resolve_raytrace_settings(
        viewer.settings(),
        viewer.named_colors(),
        viewer.clear_color(),
        ray,
    );

    // Camera matrices — temporarily set aspect ratio for the output dimensions
    let camera = viewer.camera_mut();
    let original_aspect = camera.aspect();
    camera.set_aspect(final_width as f32 / final_height as f32);
    let cam = extract_camera_matrices(camera.view_matrix(), camera.projection_matrix());
    camera.set_aspect(original_aspect);

    let device = viewer
        .gpu_device()
        .ok_or_else(|| RaytraceSceneError::RaytraceFailed("No GPU device".into()))?;
    let queue = viewer
        .gpu_queue()
        .ok_or_else(|| RaytraceSceneError::RaytraceFailed("No GPU queue".into()))?;

    let params = RaytraceParams::new(final_width, final_height)
        .with_antialias(antialias)
        .with_camera(cam.view, cam.proj, cam.view_inv, cam.proj_inv, cam.camera_pos)
        .with_settings(rt_settings);

    let image_data = raytrace(device, queue, &primitives, &bvh, &params)
        .map_err(|e| RaytraceSceneError::RaytraceFailed(e.to_string()))?;

    Ok((image_data, final_width, final_height))
}

// ---------------------------------------------------------------------------
// Preview raytracing (immutable — no viewer mutation required)
// ---------------------------------------------------------------------------

/// Raytrace a preview image using only immutable references.
///
/// Unlike [`raytrace_scene`] this never mutates the camera; instead it
/// computes the projection matrix for the requested aspect ratio directly.
#[allow(clippy::too_many_arguments)]
pub(crate) fn raytrace_preview(
    registry: &ObjectRegistry,
    camera: &pymol_scene::Camera,
    settings: &pymol_settings::Settings,
    named_colors: &NamedColors,
    clear_color: [f32; 3],
    ray: &ResolvedRaySettings,
    device: &wgpu::Device,
    queue: &wgpu::Queue,
    width: u32,
    height: u32,
) -> Result<(Vec<u8>, u32, u32), RaytraceSceneError> {
    let (primitives, bvh) = collect_scene_primitives(registry, named_colors, settings)?;

    let rt_settings = resolve_raytrace_settings(settings, named_colors, clear_color, ray);

    let preview_aspect = width as f32 / height as f32;
    let cam = extract_camera_matrices(
        camera.view_matrix(),
        camera.projection_matrix_for_aspect(preview_aspect),
    );

    let params = RaytraceParams::new(width, height)
        .with_antialias(1)
        .with_camera(cam.view, cam.proj, cam.view_inv, cam.proj_inv, cam.camera_pos)
        .with_settings(rt_settings);

    let image_data = raytrace(device, queue, &primitives, &bvh, &params)
        .map_err(|e| RaytraceSceneError::RaytraceFailed(e.to_string()))?;

    Ok((image_data, width, height))
}

// ---------------------------------------------------------------------------
// Collect triangles from mesh representations
// ---------------------------------------------------------------------------

fn collect_triangles_from_mesh_reps(mol_obj: &MoleculeObject) -> Vec<GpuTriangle> {
    let mut triangles = Vec::new();

    let convert_mesh = |vertices: &[MeshVertex], indices: &[u32]| -> Vec<GpuTriangle> {
        let mut tris = Vec::with_capacity(indices.len() / 3);
        for chunk in indices.chunks(3) {
            if chunk.len() == 3 {
                let i0 = chunk[0] as usize;
                let i1 = chunk[1] as usize;
                let i2 = chunk[2] as usize;
                if i0 < vertices.len() && i1 < vertices.len() && i2 < vertices.len() {
                    let v0 = &vertices[i0];
                    let v1 = &vertices[i1];
                    let v2 = &vertices[i2];
                    let avg_alpha = (v0.color[3] + v1.color[3] + v2.color[3]) / 3.0;
                    let color = [
                        (v0.color[0] + v1.color[0] + v2.color[0]) / 3.0,
                        (v0.color[1] + v1.color[1] + v2.color[1]) / 3.0,
                        (v0.color[2] + v1.color[2] + v2.color[2]) / 3.0,
                        avg_alpha,
                    ];
                    let transparency = 1.0 - avg_alpha;
                    tris.push(GpuTriangle::new(
                        v0.position,
                        v1.position,
                        v2.position,
                        v0.normal,
                        v1.normal,
                        v2.normal,
                        color,
                        transparency,
                    ));
                }
            }
        }
        tris
    };

    if let Some((vertices, indices)) = mol_obj.get_cartoon_mesh() {
        triangles.extend(convert_mesh(vertices, indices));
    }
    if let Some((vertices, indices)) = mol_obj.get_surface_mesh() {
        triangles.extend(convert_mesh(vertices, indices));
    }
    if let Some((vertices, indices)) = mol_obj.get_ribbon_mesh() {
        triangles.extend(convert_mesh(vertices, indices));
    }

    triangles
}

// ---------------------------------------------------------------------------
// Helper
// ---------------------------------------------------------------------------

#[inline]
fn matrix_to_array(m: lin_alg::f32::Mat4) -> [[f32; 4]; 4] {
    <[[f32; 4]; 4]>::from(m)
}
