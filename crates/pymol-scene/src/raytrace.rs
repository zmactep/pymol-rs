//! Raytracing support for scene rendering
//!
//! This module provides the unified raytracing functionality used by both
//! the standalone Viewer and the GUI's ViewerAdapter.

use pymol_color::{ChainColors, ElementColors, NamedColors};
use pymol_raytracer::{
    collect_from_molecule, matrix_to_array, raytrace, Bvh, CollectOptions, GpuTriangle,
    PrimitiveCollector, RayColorResolver, RaytraceParams, RaytraceSettings,
};
use pymol_render::MeshVertex;
use pymol_settings::GlobalSettings;

use crate::object::{MoleculeObject, Object, ObjectRegistry};
use crate::Camera;

/// Error type for raytracing operations
#[derive(Debug, thiserror::Error)]
pub enum RaytraceError {
    #[error("No render context available")]
    NoRenderContext,
    #[error("No primitives to raytrace")]
    NoPrimitives,
    #[error("BVH construction failed: {0}")]
    BvhFailed(String),
    #[error("Raytracing failed: {0}")]
    RaytraceFailed(String),
    #[error("Image save failed: {0}")]
    ImageSaveFailed(String),
}

/// Normalize a matrix by dividing by data[15] if needed.
///
/// lin_alg's inverse can produce scaled results where the homogeneous
/// coordinate (w) is not 1.0. This function normalizes the matrix to
/// ensure proper affine transform properties.
fn normalize_inverse_matrix(mat: lin_alg::f32::Mat4) -> lin_alg::f32::Mat4 {
    let w = mat.data[15];
    if w.abs() > 1e-6 && (w - 1.0).abs() > 1e-6 {
        let mut normalized = mat.clone();
        for i in 0..16 {
            normalized.data[i] /= w;
        }
        normalized
    } else {
        mat
    }
}

/// Input parameters for raytracing a scene
///
/// This struct groups all the references needed to perform raytracing,
/// avoiding the need to pass many individual parameters.
pub struct RaytraceInput<'a> {
    /// wgpu device
    pub device: &'a wgpu::Device,
    /// wgpu queue
    pub queue: &'a wgpu::Queue,
    /// Camera for view matrices
    pub camera: &'a mut Camera,
    /// Object registry containing molecules
    pub registry: &'a ObjectRegistry,
    /// Global settings
    pub settings: &'a GlobalSettings,
    /// Named colors
    pub named_colors: &'a NamedColors,
    /// Element colors
    pub element_colors: &'a ElementColors,
    /// Chain colors
    pub chain_colors: &'a ChainColors,
    /// Background/clear color
    pub clear_color: [f32; 3],
    /// Default size (width, height) if not specified
    pub default_size: (u32, u32),
}

/// Collect triangles from mesh-based representations (cartoon, surface, mesh)
pub fn collect_triangles_from_mesh_reps(mol_obj: &MoleculeObject) -> Vec<GpuTriangle> {
    let mut triangles = Vec::new();

    // Helper to convert mesh vertices/indices to GpuTriangles
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
                    // Use average color for the triangle
                    let color = [
                        (v0.color[0] + v1.color[0] + v2.color[0]) / 3.0,
                        (v0.color[1] + v1.color[1] + v2.color[1]) / 3.0,
                        (v0.color[2] + v1.color[2] + v2.color[2]) / 3.0,
                        (v0.color[3] + v1.color[3] + v2.color[3]) / 3.0,
                    ];
                    tris.push(GpuTriangle::new(
                        v0.position,
                        v1.position,
                        v2.position,
                        v0.normal,
                        v1.normal,
                        v2.normal,
                        color,
                        0.0, // transparency
                    ));
                }
            }
        }
        tris
    };

    // Collect from cartoon representation
    if let Some((vertices, indices)) = mol_obj.get_cartoon_mesh() {
        triangles.extend(convert_mesh(vertices, indices));
    }

    // Collect from surface representation
    if let Some((vertices, indices)) = mol_obj.get_surface_mesh() {
        triangles.extend(convert_mesh(vertices, indices));
    }

    // Collect from mesh representation
    if let Some((vertices, indices)) = mol_obj.get_mesh_data() {
        triangles.extend(convert_mesh(vertices, indices));
    }

    triangles
}

/// Perform raytracing on a scene
///
/// This is the main entry point for raytracing. It collects primitives from
/// all enabled molecules, builds a BVH, and performs GPU raytracing.
///
/// # Arguments
/// * `input` - All required input parameters grouped in RaytraceInput
/// * `width` - Optional output width (uses default_size if None)
/// * `height` - Optional output height (uses default_size if None)
/// * `antialias` - Antialiasing level (1 = no AA, 2+ = supersampling)
///
/// # Returns
/// * `Ok((image_data, width, height))` - RGBA image data and dimensions
/// * `Err(RaytraceError)` - If raytracing fails
pub fn raytrace_scene(
    input: &mut RaytraceInput,
    width: Option<u32>,
    height: Option<u32>,
    antialias: u32,
) -> Result<(Vec<u8>, u32, u32), RaytraceError> {
    let final_width = width.unwrap_or(input.default_size.0);
    let final_height = height.unwrap_or(input.default_size.1);

    // Collect primitives from all enabled molecule objects
    let mut all_primitives = PrimitiveCollector::new();
    let color_resolver = RayColorResolver::new(
        input.named_colors,
        input.element_colors,
        input.chain_colors,
    );

    // Get sphere/stick settings
    let sphere_scale = input.settings.get_float(pymol_settings::id::sphere_scale);
    let stick_radius = input.settings.get_float(pymol_settings::id::stick_radius);
    let options = CollectOptions {
        sphere_scale,
        stick_radius,
        ..Default::default()
    };

    for name in input.registry.names().map(|s| s.to_string()).collect::<Vec<_>>() {
        if let Some(mol_obj) = input.registry.get_molecule(&name) {
            if !mol_obj.is_enabled() {
                continue;
            }

            let molecule = mol_obj.molecule();
            if let Some(coord_set) = molecule.get_coord_set(mol_obj.display_state()) {
                let prims = collect_from_molecule(molecule, coord_set, &color_resolver, &options);
                all_primitives.add_spheres(prims.spheres);
                all_primitives.add_cylinders(prims.cylinders);
                all_primitives.add_triangles(prims.triangles);
            }

            // Collect triangles from mesh-based representations (cartoon, surface, mesh)
            let triangles = collect_triangles_from_mesh_reps(mol_obj);
            all_primitives.add_triangles(triangles);
        }
    }

    let primitives = all_primitives.build();
    if primitives.is_empty() {
        return Err(RaytraceError::NoPrimitives);
    }

    // Build BVH
    let bvh = Bvh::build(&primitives)
        .map_err(|e| RaytraceError::BvhFailed(e.to_string()))?;

    // Update camera aspect ratio for raytracing
    let original_aspect = input.camera.aspect();
    input.camera.set_aspect(final_width as f32 / final_height as f32);

    // Get camera matrices
    let view_matrix = matrix_to_array(input.camera.view_matrix());
    let proj_matrix = matrix_to_array(input.camera.projection_matrix());

    // Get and normalize inverse matrices
    let view_inv = input.camera.view_matrix().inverse().unwrap_or(lin_alg::f32::Mat4::new_identity());
    let view_inv_normalized = normalize_inverse_matrix(view_inv);

    // Extract camera position from normalized inverse view matrix
    let camera_pos = lin_alg::f32::Vec3::new(
        view_inv_normalized.data[12],
        view_inv_normalized.data[13],
        view_inv_normalized.data[14],
    );

    let view_inv_matrix = matrix_to_array(view_inv_normalized);

    let proj_inv = input.camera.projection_matrix().inverse().unwrap_or(lin_alg::f32::Mat4::new_identity());
    let proj_inv_normalized = normalize_inverse_matrix(proj_inv);
    let proj_inv_matrix = matrix_to_array(proj_inv_normalized);

    // Restore camera aspect
    input.camera.set_aspect(original_aspect);

    // Build raytracing settings from global settings
    let ray_shadow = input.settings.get_bool(pymol_settings::id::ray_shadow);
    let ray_trace_fog = input.settings.get_float(pymol_settings::id::ray_trace_fog) > 0.0;
    let ray_max_passes = input.settings.get_int(pymol_settings::id::ray_max_passes) as u32;
    let ray_transparency_shadows = input.settings.get_bool(pymol_settings::id::ray_transparency_shadows);

    let ambient = input.settings.get_float(pymol_settings::id::ambient);
    let direct = input.settings.get_float(pymol_settings::id::direct);
    let reflect = input.settings.get_float(pymol_settings::id::reflect);
    let specular = input.settings.get_float(pymol_settings::id::specular);
    let shininess = input.settings.get_float(pymol_settings::id::shininess);

    // Ray trace mode settings
    let ray_trace_mode = input.settings.get_int(pymol_settings::id::ray_trace_mode);
    let ray_opaque_background = input.settings.get_int(pymol_settings::id::ray_opaque_background);

    // Resolve ray_trace_color (default -6 = black in PyMOL)
    let ray_trace_color_idx = input.settings.get_color(pymol_settings::id::ray_trace_color);
    let ray_trace_color = if ray_trace_color_idx >= 0 {
        // Positive index: look up in named colors
        input.named_colors
            .get_by_index(ray_trace_color_idx as u32)
            .map(|c| c.to_rgba(1.0))
            .unwrap_or([0.0, 0.0, 0.0, 1.0])
    } else {
        // Negative index: PyMOL special colors, -6 = black
        [0.0, 0.0, 0.0, 1.0]
    };

    // Edge detection parameters (PyMOL gradient-of-gradient algorithm)
    let ray_trace_slope_factor = input.settings.get_float(pymol_settings::id::ray_trace_slope_factor);
    let ray_trace_depth_factor = input.settings.get_float(pymol_settings::id::ray_trace_depth_factor);
    let ray_trace_disco_factor = input.settings.get_float(pymol_settings::id::ray_trace_disco_factor);
    let ray_trace_gain = input.settings.get_float(pymol_settings::id::ray_trace_gain);

    let settings = RaytraceSettings {
        // PyMOL default light direction (from upper-left-front)
        light_dir: [-0.4, -0.4, -1.0, 0.0],
        ambient,
        direct,
        reflect,
        specular,
        shininess,
        bg_color: [input.clear_color[0], input.clear_color[1], input.clear_color[2], 1.0],
        fog_start: input.settings.get_float(pymol_settings::id::fog_start),
        fog_end: 1.0,
        fog_density: input.settings.get_float(pymol_settings::id::fog),
        fog_color: [input.clear_color[0], input.clear_color[1], input.clear_color[2], 1.0],
        ray_shadow,
        ray_max_passes,
        ray_trace_fog,
        ray_transparency_shadows,
        ray_trace_mode,
        ray_trace_color,
        ray_opaque_background,
        ray_trace_slope_factor,
        ray_trace_depth_factor,
        ray_trace_disco_factor,
        ray_trace_gain,
        ..Default::default()
    };

    let params = RaytraceParams::new(final_width, final_height)
        .with_antialias(antialias)
        .with_camera(
            view_matrix,
            proj_matrix,
            view_inv_matrix,
            proj_inv_matrix,
            [camera_pos.x, camera_pos.y, camera_pos.z, 1.0],
        )
        .with_settings(settings);

    // Perform raytracing
    let image_data = raytrace(input.device, input.queue, &primitives, &bvh, &params)
        .map_err(|e| RaytraceError::RaytraceFailed(e.to_string()))?;

    Ok((image_data, final_width, final_height))
}

/// Perform raytracing and save to a PNG file
///
/// Convenience function that calls `raytrace_scene` and saves the result.
pub fn raytrace_to_file(
    input: &mut RaytraceInput,
    path: impl AsRef<std::path::Path>,
    width: Option<u32>,
    height: Option<u32>,
    antialias: u32,
) -> Result<(u32, u32), RaytraceError> {
    let (image_data, final_width, final_height) = raytrace_scene(input, width, height, antialias)?;

    // Save to PNG
    image::save_buffer(
        path,
        &image_data,
        final_width,
        final_height,
        image::ColorType::Rgba8,
    )
    .map_err(|e| RaytraceError::ImageSaveFailed(e.to_string()))?;

    Ok((final_width, final_height))
}
