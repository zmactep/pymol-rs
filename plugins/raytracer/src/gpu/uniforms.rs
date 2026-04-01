//! GPU uniform buffer layout for raytracing.

use crate::bvh::Bvh;
use crate::primitive::Primitives;

use super::RaytraceParams;

/// GPU uniforms for raytracing (maps 1:1 to the WGSL `Uniforms` struct).
#[repr(C)]
#[derive(Copy, Clone, Debug, bytemuck::Pod, bytemuck::Zeroable)]
pub struct RaytraceUniforms {
    // Camera
    pub view_matrix: [[f32; 4]; 4],
    pub proj_matrix: [[f32; 4]; 4],
    pub view_inv_matrix: [[f32; 4]; 4],
    pub proj_inv_matrix: [[f32; 4]; 4],
    pub camera_pos: [f32; 4],

    // Viewport
    pub viewport: [f32; 4], // width, height, 1/width, 1/height

    // Multi-light support
    pub light_dirs: [[f32; 4]; 9],
    pub light_count: i32,
    pub spec_count: i32,
    pub _pad_light_0: i32,
    pub _pad_light_1: i32,

    // Lighting parameters
    pub ambient: f32,
    pub direct: f32,
    pub reflect: f32,
    pub specular: f32,
    pub shininess: f32,
    pub _pad_light1: f32,
    pub _pad_light2: f32,
    pub _pad_light3: f32,

    // Background
    pub bg_color: [f32; 4],

    // Fog
    pub fog_start: f32,
    pub fog_end: f32,
    pub fog_density: f32,
    pub _pad0: f32,
    pub fog_color: [f32; 4],

    // Primitive counts
    pub sphere_count: u32,
    pub cylinder_count: u32,
    pub triangle_count: u32,
    pub bvh_node_count: u32,

    // Ray settings
    pub ray_shadow: u32,
    pub ray_max_passes: u32,
    pub ray_trace_fog: u32,
    pub ray_transparency_shadows: u32,

    // Ray trace mode settings
    pub ray_trace_mode: u32,
    pub ray_opaque_background: u32,
    pub transparency_mode: u32,
    pub _pad2: u32,
    pub ray_trace_color: [f32; 4],
}

impl RaytraceUniforms {
    /// Create uniforms from raytracing parameters.
    pub fn from_params(params: &RaytraceParams, primitives: &Primitives, bvh: &Bvh) -> Self {
        let settings = &params.settings;
        let supersample = params.antialias.max(1);
        let width = (params.width * supersample) as f32;
        let height = (params.height * supersample) as f32;

        Self {
            view_matrix: params.view_matrix,
            proj_matrix: params.proj_matrix,
            view_inv_matrix: params.view_inv_matrix,
            proj_inv_matrix: params.proj_inv_matrix,
            camera_pos: params.camera_pos,
            viewport: [width, height, 1.0 / width, 1.0 / height],
            light_dirs: settings.light_dirs,
            light_count: settings.light_count,
            spec_count: settings.spec_count,
            _pad_light_0: 0,
            _pad_light_1: 0,
            ambient: settings.ambient,
            direct: settings.direct,
            reflect: settings.reflect,
            specular: settings.specular,
            shininess: settings.shininess,
            _pad_light1: 0.0,
            _pad_light2: 0.0,
            _pad_light3: 0.0,
            bg_color: settings.bg_color,
            fog_start: settings.fog_start,
            fog_end: settings.fog_end,
            fog_density: settings.fog_density,
            _pad0: 0.0,
            fog_color: settings.fog_color,
            sphere_count: primitives.spheres.len() as u32,
            cylinder_count: primitives.cylinders.len() as u32,
            triangle_count: primitives.triangles.len() as u32,
            bvh_node_count: bvh.nodes.len() as u32,
            ray_shadow: if settings.ray_shadow { 1 } else { 0 },
            ray_max_passes: settings.ray_max_passes,
            ray_trace_fog: if settings.ray_trace_fog { 1 } else { 0 },
            ray_transparency_shadows: if settings.ray_transparency_shadows { 1 } else { 0 },
            ray_trace_mode: settings.ray_trace_mode as u32,
            ray_opaque_background: if settings.ray_opaque_background { 1 } else { 0 },
            transparency_mode: settings.transparency_mode as u32,
            _pad2: 0,
            ray_trace_color: settings.ray_trace_color,
        }
    }
}
