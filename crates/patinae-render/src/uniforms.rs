//! Frame-wide uniforms — bind group 0.
//!
//! Layout MUST mirror `FrameUniforms` in `shaders/common/frame.wgsl`. The
//! struct is `bytemuck::Pod` and is uploaded directly to a uniform buffer.

use bytemuck::{Pod, Zeroable};

/// Maximum number of directional lights accepted by the shader. The host
/// clamps light count to 1..=10, where 1 means ambient only and larger values
/// add directional lights.
pub const MAX_LIGHTS: usize = 9;

/// Default direction vectors for `light` .. `light9`. Used when the host has
/// not yet pushed shading settings and as the seed for resets.
const DEFAULT_LIGHT_DIRS: [[f32; 4]; MAX_LIGHTS] = [
    [-0.4, -0.4, -1.0, 0.0],
    [-0.55, -0.7, 0.15, 0.0],
    [0.3, -0.6, -0.2, 0.0],
    [-1.2, 0.3, -0.2, 0.0],
    [0.3, 0.6, -0.75, 0.0],
    [-0.3, 0.5, 0.0, 0.0],
    [0.9, -0.1, -0.15, 0.0],
    [1.3, 2.0, 0.8, 0.0],
    [-1.7, -0.5, 1.2, 0.0],
];

/// Per-frame uniform block. Bound at group 0, binding 0.
///
/// All matrices are column-major (wgsl convention). Vec4 fields are used as
/// packed parameter tuples — see field docs for what each lane carries.
#[repr(C)]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
pub struct FrameUniforms {
    pub view_proj: [[f32; 4]; 4],
    pub view: [[f32; 4]; 4],
    pub proj: [[f32; 4]; 4],
    pub view_inv: [[f32; 4]; 4],
    /// Inverse projection — used by SSAO to reconstruct view-space
    /// position from a depth sample and screen-space UV.
    pub proj_inv: [[f32; 4]; 4],
    /// Up to `MAX_LIGHTS` directional lights in **view space**. Each `xyz`
    /// is the light direction vector, toward the light; shaders negate it to
    /// produce the from-fragment-to-light vector. `w` is reserved.
    pub light_dirs: [[f32; 4]; MAX_LIGHTS],
    /// `(ambient, direct, reflect, specular)` for classic multi-light shading.
    pub light_intensity: [f32; 4],
    /// `(shininess, spec_direct, spec_direct_power, depth_cue_factor)`.
    pub light_spec: [f32; 4],
    /// `(light_count_as_f32, spec_count_as_f32, _, _)`. The shader casts
    /// to `i32` — uniform i32 is supported but mixing in `Pod` adds
    /// padding noise; floats are simpler.
    pub light_counts: [f32; 4],
    /// `(start, end, density, _)` — fog params in view-space depth.
    pub fog: [f32; 4],
    /// `(r, g, b, _)` — fog tint (host typically passes scene `clear_color`).
    pub fog_color: [f32; 4],
    /// `(near, far, picking_scale, scene_max_depth)`. `scene_max_depth`
    /// rescales the WBOIT weight function so it stays scale-invariant
    /// across small (1UBQ) and huge (3J3Q) structures.
    pub clip: [f32; 4],
    /// `(width, height, 1/width, 1/height)` in physical pixels.
    pub viewport: [f32; 4],
    /// `(time_seconds, dt_seconds, frame_index_as_f32, _)`.
    pub time: [f32; 4],
}

impl Default for FrameUniforms {
    fn default() -> Self {
        let identity = [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        Self {
            view_proj: identity,
            view: identity,
            proj: identity,
            view_inv: identity,
            proj_inv: identity,
            light_dirs: DEFAULT_LIGHT_DIRS,
            // Classic multi-light defaults.
            light_intensity: [0.14, 0.45, 0.45, 1.0],
            light_spec: [55.0, 0.0, 55.0, 0.0],
            // light_count = 2 (ambient + 1 directional); spec_count = -1 (all).
            light_counts: [2.0, -1.0, 0.0, 0.0],
            fog: [0.0, 1.0, 0.0, 0.0],
            fog_color: [0.0, 0.0, 0.0, 0.0],
            clip: [0.1, 1000.0, 0.5, 200.0],
            viewport: [1.0, 1.0, 1.0, 1.0],
            time: [0.0, 0.0, 0.0, 0.0],
        }
    }
}

impl FrameUniforms {
    pub const SIZE: u64 = std::mem::size_of::<FrameUniforms>() as u64;

    /// Update the viewport lane to match a physical-pixel size.
    pub fn set_viewport(&mut self, width: u32, height: u32) {
        let w = width.max(1) as f32;
        let h = height.max(1) as f32;
        self.viewport = [w, h, 1.0 / w, 1.0 / h];
    }

    /// Set the scene max depth used by the WBOIT weight. Larger structures
    /// need a larger value so the weight stays in a sensible range.
    pub fn set_scene_max_depth(&mut self, depth: f32) {
        self.clip[3] = depth.max(1.0);
    }
}
