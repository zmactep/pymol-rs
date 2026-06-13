// Frame uniforms — bind group 0, binding 0.
//
// Single source of truth for view/projection, lighting, fog, clipping, and
// viewport state. Every pipeline that draws geometry binds this group.
// Layout MUST match `FrameUniforms` in src/uniforms.rs (bytemuck::Pod).

const MAX_LIGHTS: u32 = 9u;

struct FrameUniforms {
    view_proj:        mat4x4<f32>,
    view:             mat4x4<f32>,
    proj:             mat4x4<f32>,
    view_inv:         mat4x4<f32>,
    // Inverse projection — SSAO reconstructs view-space position
    // from depth + screen UV via `proj_inv * (ndc.x, ndc.y, z, 1)`.
    proj_inv:         mat4x4<f32>,
    // Up to MAX_LIGHTS directional lights in **view space** (camera-attached).
    // Each .xyz is the light direction vector toward the light.
    light_dirs:       array<vec4<f32>, MAX_LIGHTS>,
    // (ambient, direct, reflect, specular)
    light_intensity:  vec4<f32>,
    // (shininess, spec_direct, spec_direct_power, depth_cue_factor)
    light_spec:       vec4<f32>,
    // (light_count_as_f32, spec_count_as_f32, _, _) — cast to i32 in shader.
    light_counts:     vec4<f32>,
    // (start, end, density, _)
    fog:              vec4<f32>,
    // (r, g, b, _)
    fog_color:        vec4<f32>,
    // (near, far, picking_scale, scene_max_depth) — scene_max_depth rescales WBOIT.
    clip:             vec4<f32>,
    // (width, height, 1/width, 1/height) in physical pixels.
    viewport:         vec4<f32>,
    // (time_seconds, dt_seconds, frame_index, _).
    time:             vec4<f32>,
};

@group(0) @binding(0)
var<uniform> frame: FrameUniforms;
