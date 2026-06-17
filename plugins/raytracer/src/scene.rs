//! Scene raytracing orchestration.
//!
//! Collects primitives from the viewer, builds a BVH, resolves settings,
//! extracts camera matrices, and dispatches GPU raytracing.

use patinae_color::NamedPalette;
use patinae_plugin::prelude::*;
#[cfg(test)]
use patinae_render::{DisplayedGeometry, TraceGeometryChunk};
use patinae_scene::normalize_matrix;

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

use crate::gpu::RaytraceParams;
#[cfg(test)]
use crate::primitive::{GpuCylinder, GpuSphere};
#[cfg(test)]
use crate::primitive::{GpuTriangle, PrimitiveCollector, Primitives};
use crate::settings::{RaytraceSettings, ResolvedRaySettings};

#[cfg(test)]
const RAY_LINE_RADIUS: f32 = 0.035;
#[cfg(test)]
const RAY_POINT_RADIUS: f32 = 0.035;
// ---------------------------------------------------------------------------
// Error
// ---------------------------------------------------------------------------

#[derive(Debug, thiserror::Error)]
pub(crate) enum RaytraceSceneError {
    #[error("Render artifact ray path failed: {0}")]
    RenderArtifacts(String),
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum RaytraceSceneTarget {
    Export,
    Viewport,
}

pub(crate) enum RaytraceSceneOutput {
    CpuImage {
        data: Vec<u8>,
        width: u32,
        height: u32,
        profile_lines: Vec<String>,
    },
    GpuViewport {
        width: u32,
        height: u32,
        profile_lines: Vec<String>,
    },
}

// ---------------------------------------------------------------------------
// Test/support helpers
// ---------------------------------------------------------------------------

#[cfg(test)]
fn primitives_from_displayed_geometry(displayed: &DisplayedGeometry) -> Primitives {
    let mut collector = PrimitiveCollector::new();
    add_trace_geometry_primitives(
        &TraceGeometryChunk::from_displayed(displayed),
        &mut collector,
    );
    collector.build()
}

#[cfg(test)]
fn add_trace_geometry_primitives(chunk: &TraceGeometryChunk, collector: &mut PrimitiveCollector) {
    for sphere in &chunk.spheres {
        collector.add_sphere(GpuSphere::new(
            sphere.center,
            sphere.radius,
            sphere.material.rgba,
            sphere.material.transparency,
        ));
    }

    for cylinder in &chunk.cylinders {
        collector.add_cylinder(GpuCylinder::new(
            cylinder.start,
            cylinder.end,
            cylinder.radius,
            cylinder.material_start.rgba,
            cylinder.material_end.rgba,
            cylinder
                .material_start
                .transparency
                .max(cylinder.material_end.transparency),
        ));
    }

    for triangle in &chunk.triangles {
        collector.add_triangle(GpuTriangle::new(
            triangle.positions[0],
            triangle.positions[1],
            triangle.positions[2],
            triangle.normals[0],
            triangle.normals[1],
            triangle.normals[2],
            triangle.material.rgba,
            triangle.material.transparency,
        ));
    }

    for line in &chunk.line_segments {
        // Screen-space lines have no physical radius. The ray adapter renders
        // them as a thin world-space cylinder so line/mesh displays remain
        // visible in the ray image.
        collector.add_cylinder(GpuCylinder::new(
            line.start,
            line.end,
            RAY_LINE_RADIUS,
            line.material_start.rgba,
            line.material_end.rgba,
            line.material_start
                .transparency
                .max(line.material_end.transparency),
        ));
    }

    for point in &chunk.point_samples {
        // Screen-space points are semantic samples; ray uses a small physical
        // sphere approximation.
        collector.add_sphere(GpuSphere::new(
            point.position,
            RAY_POINT_RADIUS,
            point.material.rgba,
            point.material.transparency,
        ));
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
    target: RaytraceSceneTarget,
) -> Result<RaytraceSceneOutput, RaytraceSceneError> {
    let final_width = width.unwrap_or(1024);
    let final_height = height.unwrap_or(768);

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

    let snapshot = viewer
        .open_render_artifact_snapshot()
        .map_err(RaytraceSceneError::RenderArtifacts)?;
    let snapshot_id = snapshot.snapshot_id;
    let render_result = crate::artifact_gpu::raytrace_artifacts(
        viewer,
        &snapshot,
        &params,
        match target {
            RaytraceSceneTarget::Export => crate::artifact_gpu::RaytraceArtifactTarget::CpuReadback,
            RaytraceSceneTarget::Viewport => {
                crate::artifact_gpu::RaytraceArtifactTarget::ViewportGpu
            }
        },
    )
    .map_err(RaytraceSceneError::RenderArtifacts);
    let close_result = viewer.close_render_artifact_snapshot(snapshot_id);
    if let Err(err) = close_result {
        return Err(RaytraceSceneError::RenderArtifacts(format!(
            "failed to close render artifact snapshot: {err}"
        )));
    }
    match render_result? {
        crate::artifact_gpu::RaytraceArtifactOutput::CpuImage {
            data,
            profile_lines,
        } => Ok(RaytraceSceneOutput::CpuImage {
            data,
            width: final_width,
            height: final_height,
            profile_lines,
        }),
        crate::artifact_gpu::RaytraceArtifactOutput::ViewportGpu { profile_lines } => {
            Ok(RaytraceSceneOutput::GpuViewport {
                width: final_width,
                height: final_height,
                profile_lines,
            })
        }
    }
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
    use patinae_render::{
        DisplayedMaterial, DisplayedMesh, DisplayedMeshVertex, DisplayedObjectGeometry,
        DisplayedPrimitive, ObjectId, RepKind, TraceLineSegment, TraceMaterial, TracePointSample,
    };

    fn primitives_from_trace_geometry_chunk(chunk: &TraceGeometryChunk) -> Primitives {
        let mut collector = PrimitiveCollector::new();
        add_trace_geometry_primitives(chunk, &mut collector);
        collector.build()
    }

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
    fn trace_triangle_matches_displayed_triangle_count() {
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
        let displayed_primitives = primitives_from_displayed_geometry(&displayed);
        let trace = TraceGeometryChunk::from_displayed(&displayed);
        let trace_primitives = primitives_from_trace_geometry_chunk(&trace);

        assert_eq!(
            trace_primitives.total_count(),
            displayed_primitives.total_count()
        );
        assert_eq!(trace_primitives.triangles.len(), 1);
    }

    #[test]
    fn trace_line_and_point_samples_map_to_cylinder_and_sphere() {
        let material = TraceMaterial {
            rgba: [0.2, 0.8, 0.4, 1.0],
            transparency: 0.0,
        };
        let trace = TraceGeometryChunk {
            line_segments: vec![TraceLineSegment {
                start: [0.0, 0.0, 0.0],
                end: [1.0, 0.0, 0.0],
                width_px: 1.0,
                material_start: material,
                material_end: material,
            }],
            point_samples: vec![TracePointSample {
                position: [0.0, 1.0, 0.0],
                radius_px: 1.0,
                material,
            }],
            ..TraceGeometryChunk::default()
        };

        let primitives = primitives_from_trace_geometry_chunk(&trace);

        assert_eq!(primitives.cylinders.len(), 1);
        assert_eq!(primitives.spheres.len(), 1);
        assert_eq!(primitives.total_count(), 2);
    }
}
