//! GPU raytracing plugin for PyMOL-RS
//!
//! Provides the `ray` command and ray tracing settings. All raytracing logic
//! (BVH, compute shaders, primitive collection) lives in this plugin.

use pymol_plugin::prelude::*;
use pymol_plugin::{define_plugin_settings, pymol_plugin};

use pymol_color::{ChainColors, ElementColors};
use pymol_render::MeshVertex;
use pymol_scene::{
    compute_reflect_scale, normalize_matrix, MoleculeObject, Object, ViewportImage,
};

// Raytracer engine modules (moved from crates/pymol-raytracer)
pub mod bvh;
pub mod collect;
pub mod context;
pub mod edge_pipeline;
pub mod error;
pub mod pipeline;
pub mod primitive;
pub mod render;

use bvh::Bvh;
use collect::{collect_from_molecule, CollectOptions, RayColorResolver};
use context::raytrace;
use primitive::{GpuCylinder, GpuTriangle, PrimitiveCollector};
use render::{RaytraceParams, RaytraceSettings};

// ---------------------------------------------------------------------------
// Plugin settings
// ---------------------------------------------------------------------------

define_plugin_settings! {
    RaySettings {
        max_passes: i32 = 25, name = "ray_max_passes",
            min = 1.0, max = 100.0;
        shadow: bool = true, name = "ray_shadow";
        trace_mode: i32 = 0, name = "ray_trace_mode";
        trace_fog: f32 = -1.0, name = "ray_trace_fog";
        trace_color: i32 = -6, name = "ray_trace_color";
        opaque_background: i32 = -1, name = "ray_opaque_background";
        transparency_shadows: bool = true, name = "ray_transparency_shadows";
        trace_depth_factor: f32 = 0.1, name = "ray_trace_depth_factor",
            min = 0.0, max = 1.0;
        trace_slope_factor: f32 = 0.6, name = "ray_trace_slope_factor",
            min = 0.0, max = 2.0;
        trace_disco_factor: f32 = 0.05, name = "ray_trace_disco_factor",
            min = 0.0, max = 1.0;
        trace_gain: f32 = 0.12, name = "ray_trace_gain",
            min = 0.0, max = 1.0;
    }
}

// ---------------------------------------------------------------------------
// Plugin declaration
// ---------------------------------------------------------------------------

pymol_plugin! {
    name: "raytracer",
    version: "0.2.3",
    description: "GPU compute shader raytracing — provides the 'ray' command",
    commands: [RayCommand],
    settings: [RaySettings],
}

// ---------------------------------------------------------------------------
// RayCommand
// ---------------------------------------------------------------------------

struct RayCommand;

impl Command for RayCommand {
    fn name(&self) -> &str {
        "ray"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "ray" performs ray-tracing and saves the resulting image to a file.
    Ray tracing produces high-quality images with proper shadows, lighting,
    and transparency effects.

USAGE

    ray [ width [, height [, antialias [, filename [, quiet ]]]]]

ARGUMENTS

    width = integer: width in pixels (default: current window width)
    height = integer: height in pixels (default: current window height)
    antialias = integer: antialiasing level 1-4 (default: antialias setting)
        1 = no antialiasing
        2 = 2x2 supersampling
        3 = 3x3 supersampling
        4 = 4x4 supersampling
    filename = string: output file path (default: no file saved, returns to display)
    quiet = 0/1: suppress feedback (default: 0)

EXAMPLES

    ray                          # Raytrace at current resolution
    ray 1920, 1080               # Raytrace at 1080p
    ray 1920, 1080, 2            # Raytrace at 1080p with 2x2 AA
    ray width=1920, height=1080, filename=output.png

SEE ALSO

    png
"#
    }

    fn execute<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        args: &ParsedCommand,
    ) -> CmdResult {
        use std::time::Instant;

        let width = args
            .get_int(0)
            .or_else(|| args.get_named_int("width"))
            .map(|v| v as u32);

        let height = args
            .get_int(1)
            .or_else(|| args.get_named_int("height"))
            .map(|v| v as u32);

        let antialias = args
            .get_int(2)
            .or_else(|| args.get_named_int("antialias"))
            .unwrap_or_else(|| ctx.viewer.settings().ui.antialias as i64) as u32;

        let filename = args
            .get_str(3)
            .or_else(|| args.get_named_str("filename"));

        let quiet = args
            .get_bool(4)
            .or_else(|| args.get_named_bool("quiet"))
            .unwrap_or(false);

        // Ensure representations are built (mutable borrow, released immediately)
        ctx.viewer.prepare_render();

        // Read ray settings from plugin's dynamic settings
        let ray_settings = read_ray_settings(|name| {
            ctx.dynamic_setting(name)
                .and_then(|e| e.store.read().ok())
                .and_then(|s| s.get(name).cloned())
        });

        // Determine output dimensions
        let viewport = ctx.viewer.viewport_size();
        let final_width = width.unwrap_or(viewport.0.max(1024));
        let final_height = height.unwrap_or(viewport.1.max(768));

        let start = Instant::now();

        // Perform raytracing (handles borrow-checker dance internally)
        let (image_data, w, h) = raytrace_scene(
            ctx.viewer,
            &ray_settings,
            Some(final_width),
            Some(final_height),
            antialias,
        )
        .map_err(|e| CmdError::execution(format!("Ray tracing failed: {e}")))?;

        let elapsed = start.elapsed();

        if let Some(path) = filename {
            let path = expand_path(path);
            let path = if path.extension().map(|e| e.to_ascii_lowercase()) != Some("png".into()) {
                path.with_extension("png")
            } else {
                path.to_path_buf()
            };

            image::save_buffer(&path, &image_data, w, h, image::ColorType::Rgba8)
                .map_err(|e| CmdError::execution(format!("Failed to save PNG: {e}")))?;

            if !quiet {
                ctx.print(&format!(
                    " Ray: render time {:02}:{:02}:{:02}.{:02}  ({}x{})",
                    elapsed.as_secs() / 3600,
                    (elapsed.as_secs() % 3600) / 60,
                    elapsed.as_secs() % 60,
                    elapsed.subsec_millis() / 10,
                    w,
                    h
                ));
                ctx.print(&format!(" Saved \"{}\"", path.display()));
            }
        } else {
            ctx.viewer.set_viewport_image(Some(ViewportImage {
                data: image_data,
                width: w,
                height: h,
            }));

            if !quiet {
                ctx.print(&format!(
                    " Ray: render time {:02}:{:02}:{:02}.{:02}  ({}x{})",
                    elapsed.as_secs() / 3600,
                    (elapsed.as_secs() % 3600) / 60,
                    elapsed.as_secs() % 60,
                    elapsed.subsec_millis() / 10,
                    w,
                    h
                ));
                ctx.print(
                    " Ray trace complete. Use 'png filename' to save, or interact to dismiss.",
                );
            }
        }

        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Helper: expand ~ in paths
// ---------------------------------------------------------------------------

fn expand_path(path: &str) -> std::path::PathBuf {
    if path.starts_with('~') {
        if let Some(home) = dirs_hint() {
            return std::path::PathBuf::from(path.replacen('~', &home, 1));
        }
    }
    std::path::PathBuf::from(path)
}

fn dirs_hint() -> Option<String> {
    std::env::var("HOME").ok()
}

// ---------------------------------------------------------------------------
// Dynamic settings reader
// ---------------------------------------------------------------------------

struct ResolvedRaySettings {
    shadow: bool,
    max_passes: i32,
    trace_mode: i32,
    trace_fog: f32,
    trace_color: i32,
    opaque_background: i32,
    transparency_shadows: bool,
    trace_depth_factor: f32,
    trace_slope_factor: f32,
    trace_disco_factor: f32,
    trace_gain: f32,
}

fn read_ray_settings(get_setting: impl Fn(&str) -> Option<pymol_settings::SettingValue>) -> ResolvedRaySettings {
    use pymol_settings::SettingValue;

    let get_bool = |name: &str, default: bool| -> bool {
        match get_setting(name) {
            Some(SettingValue::Bool(b)) => b,
            _ => default,
        }
    };

    let get_i32 = |name: &str, default: i32| -> i32 {
        match get_setting(name) {
            Some(SettingValue::Int(i)) => i as i32,
            _ => default,
        }
    };

    let get_f32 = |name: &str, default: f32| -> f32 {
        match get_setting(name) {
            Some(SettingValue::Float(f)) => f,
            _ => default,
        }
    };

    ResolvedRaySettings {
        shadow: get_bool("ray_shadow", true),
        max_passes: get_i32("ray_max_passes", 25),
        trace_mode: get_i32("ray_trace_mode", 0),
        trace_fog: get_f32("ray_trace_fog", -1.0),
        trace_color: get_i32("ray_trace_color", -6),
        opaque_background: get_i32("ray_opaque_background", -1),
        transparency_shadows: get_bool("ray_transparency_shadows", true),
        trace_depth_factor: get_f32("ray_trace_depth_factor", 0.1),
        trace_slope_factor: get_f32("ray_trace_slope_factor", 0.6),
        trace_disco_factor: get_f32("ray_trace_disco_factor", 0.05),
        trace_gain: get_f32("ray_trace_gain", 0.12),
    }
}

// ---------------------------------------------------------------------------
// Collect triangles from mesh reps
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
// Convert a lin_alg Mat4 to [[f32; 4]; 4]
// ---------------------------------------------------------------------------

#[inline]
fn matrix_to_array(m: lin_alg::f32::Mat4) -> [[f32; 4]; 4] {
    <[[f32; 4]; 4]>::from(m)
}

// ---------------------------------------------------------------------------
// Scene raytracing
// ---------------------------------------------------------------------------

#[derive(Debug, thiserror::Error)]
enum RaytraceSceneError {
    #[error("No primitives to raytrace")]
    NoPrimitives,
    #[error("BVH construction failed: {0}")]
    BvhFailed(String),
    #[error("Raytracing failed: {0}")]
    RaytraceFailed(String),
}

fn raytrace_scene(
    viewer: &mut dyn ViewerLike,
    ray: &ResolvedRaySettings,
    width: Option<u32>,
    height: Option<u32>,
    antialias: u32,
) -> Result<(Vec<u8>, u32, u32), RaytraceSceneError> {
    let final_width = width.unwrap_or(1024);
    let final_height = height.unwrap_or(768);

    // --- Phase 1: collect primitives (immutable borrow of viewer) ---
    let mut all_primitives = PrimitiveCollector::new();
    {
        let named_colors = viewer.named_colors();
        let element_colors = ElementColors::default();
        let chain_colors = ChainColors;
        let color_resolver = RayColorResolver::new(named_colors, &element_colors, &chain_colors);

        let settings = viewer.settings();
        let sphere_scale = settings.sphere.scale;
        let stick_radius = settings.stick.radius;
        let options = CollectOptions {
            sphere_scale,
            stick_radius,
            ..Default::default()
        };

        let registry = viewer.objects();
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
                    all_primitives.add_spheres(prims.spheres);
                    all_primitives.add_cylinders(prims.cylinders);
                    all_primitives.add_triangles(prims.triangles);
                }

                let triangles = collect_triangles_from_mesh_reps(mol_obj);
                all_primitives.add_triangles(triangles);

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
                    all_primitives.add_cylinders(cylinders);
                }
            }
        }
    }

    let primitives = all_primitives.build();
    if primitives.is_empty() {
        return Err(RaytraceSceneError::NoPrimitives);
    }

    let bvh =
        Bvh::build(&primitives).map_err(|e| RaytraceSceneError::BvhFailed(e.to_string()))?;

    // --- Phase 2: read settings into owned values (immutable borrow of viewer) ---
    let rt_settings = {
        let settings = viewer.settings();
        let named_colors = viewer.named_colors();
        let clear_color = viewer.clear_color();
        let classic = &settings.shading.classic;

        let ambient = classic.ambient;
        let direct = classic.direct;
        let reflect = classic.reflect;
        let specular = classic.specular;
        let shininess = classic.shininess;

        let transparency_mode = settings.surface.transparency_mode;

        let ray_trace_color_idx = ray.trace_color;
        let ray_trace_color = if ray_trace_color_idx >= 0 {
            named_colors
                .get_by_index(ray_trace_color_idx as u32)
                .map(|c| c.to_rgba(1.0))
                .unwrap_or([0.0, 0.0, 0.0, 1.0])
        } else {
            [0.0, 0.0, 0.0, 1.0]
        };

        let shading_mode = settings.shading.mode;
        let light_count = match shading_mode {
            pymol_settings::ShadingMode::Skripkin => 1,
            pymol_settings::ShadingMode::Classic | pymol_settings::ShadingMode::Full => {
                classic.light_count
            }
        };
        let spec_count = classic.spec_count;

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

        let rt_settings = RaytraceSettings {
            light_dirs,
            light_count,
            spec_count,
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
            ray_trace_fog: ray.trace_fog > 0.0,
            ray_transparency_shadows: ray.transparency_shadows,
            ray_trace_mode: ray.trace_mode,
            ray_trace_color,
            ray_opaque_background: ray.opaque_background,
            transparency_mode,
            ray_trace_slope_factor: ray.trace_slope_factor,
            ray_trace_depth_factor: ray.trace_depth_factor,
            ray_trace_disco_factor: ray.trace_disco_factor,
            ray_trace_gain: ray.trace_gain,
            silhouette_thickness: settings.shading.common.silhouette_width,
            silhouette_depth_jump: settings.shading.common.silhouette_depth_jump,
            ..Default::default()
        };

        rt_settings
    };

    // --- Phase 3: camera matrices (mutable borrow of viewer) ---
    let camera = viewer.camera_mut();
    let original_aspect = camera.aspect();
    camera.set_aspect(final_width as f32 / final_height as f32);

    let view_matrix = matrix_to_array(camera.view_matrix());
    let proj_matrix = matrix_to_array(camera.projection_matrix());

    let mut view_inv = camera
        .view_matrix()
        .inverse()
        .unwrap_or(lin_alg::f32::Mat4::new_identity());
    normalize_matrix(&mut view_inv);

    let camera_pos =
        lin_alg::f32::Vec3::new(view_inv.data[12], view_inv.data[13], view_inv.data[14]);
    let view_inv_matrix = matrix_to_array(view_inv);

    let mut proj_inv = camera
        .projection_matrix()
        .inverse()
        .unwrap_or(lin_alg::f32::Mat4::new_identity());
    normalize_matrix(&mut proj_inv);
    let proj_inv_matrix = matrix_to_array(proj_inv);

    camera.set_aspect(original_aspect);

    // --- Phase 4: GPU raytrace (immutable borrow for device/queue) ---
    let device = viewer
        .gpu_device()
        .ok_or_else(|| RaytraceSceneError::RaytraceFailed("No GPU device".into()))?;
    let queue = viewer
        .gpu_queue()
        .ok_or_else(|| RaytraceSceneError::RaytraceFailed("No GPU queue".into()))?;

    let params = RaytraceParams::new(final_width, final_height)
        .with_antialias(antialias)
        .with_camera(
            view_matrix,
            proj_matrix,
            view_inv_matrix,
            proj_inv_matrix,
            [camera_pos.x, camera_pos.y, camera_pos.z, 1.0],
        )
        .with_settings(rt_settings);

    let image_data = raytrace(device, queue, &primitives, &bvh, &params)
        .map_err(|e| RaytraceSceneError::RaytraceFailed(e.to_string()))?;

    Ok((image_data, final_width, final_height))
}
