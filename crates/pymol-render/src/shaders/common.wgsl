// Shared uniforms, bind groups, and utility functions for all shaders.
// Concatenated before each per-shader file via concat!(include_str!(...)).

const MAX_LIGHTS: u32 = 9u;

struct GlobalUniforms {
    view_proj: mat4x4<f32>,
    view: mat4x4<f32>,
    view_inv: mat4x4<f32>,
    proj: mat4x4<f32>,
    camera_pos: vec4<f32>,
    light_dirs: array<vec4<f32>, 9>,
    light_count: i32,
    spec_count: i32,
    _pad_light_0: i32,
    _pad_light_1: i32,
    ambient: f32,
    direct: f32,
    spec_direct: f32,
    spec_direct_power: f32,
    reflect: f32,
    specular: f32,
    shininess: f32,
    _pad0: f32,
    fog_start: f32,
    fog_end: f32,
    fog_density: f32,
    depth_cue: f32,
    fog_color: vec4<f32>,
    bg_color: vec4<f32>,
    viewport: vec4<f32>,
    clip_planes: vec4<f32>,
}

@group(0) @binding(0)
var<uniform> uniforms: GlobalUniforms;

// Shadow sampling (group 1)
struct ShadowParams {
    shadow_count: u32,
    grid_size: u32,
    bias: f32,
    intensity: f32,
    mode: u32,       // 0=disabled, 1=AO (skripkin), 2=directional (full)
    pcf_samples: u32, // PCF kernel size (N×N, 1 = no PCF)
    _pad2: u32,
    _pad3: u32,
}

@group(1) @binding(0) var shadow_atlas: texture_2d<f32>;
@group(1) @binding(1) var shadow_sampler: sampler;
@group(1) @binding(2) var<uniform> shadow_params: ShadowParams;
@group(1) @binding(3) var<storage, read> shadow_matrices: array<mat4x4<f32>>;

fn apply_fog(color: vec3<f32>, depth: f32) -> vec3<f32> {
    if uniforms.fog_density <= 0.0 {
        return color;
    }
    let fog_factor = clamp((uniforms.fog_end - depth) / (uniforms.fog_end - uniforms.fog_start), 0.0, 1.0);
    return mix(uniforms.fog_color.rgb, color, fog_factor);
}

fn apply_depth_cue(color: vec3<f32>, depth: f32) -> vec3<f32> {
    if uniforms.depth_cue <= 0.0 {
        return color;
    }
    let near = uniforms.clip_planes.x;
    let far = uniforms.clip_planes.y;
    let normalized_depth = clamp((depth - near) / (far - near), 0.0, 1.0);
    let cue = 1.0 - uniforms.depth_cue * normalized_depth * 0.5;
    return color * cue;
}

// Sample shadow for a single light (directional shadow mode) with N×N PCF
fn sample_light_shadow(world_pos: vec3<f32>, light_index: u32) -> f32 {
    let shadow_pos = shadow_matrices[light_index] * vec4<f32>(world_pos, 1.0);
    let ndc = shadow_pos.xyz / shadow_pos.w;
    let uv = vec2<f32>(ndc.x * 0.5 + 0.5, 1.0 - (ndc.y * 0.5 + 0.5));
    if uv.x < 0.0 || uv.x > 1.0 || uv.y < 0.0 || uv.y > 1.0 { return 1.0; }
    let grid = shadow_params.grid_size;
    let col = light_index % grid;
    let row = light_index / grid;
    let base_tile_uv = vec2<f32>(
        (f32(col) + uv.x) / f32(grid),
        (f32(row) + uv.y) / f32(grid),
    );
    let n = i32(shadow_params.pcf_samples);
    let texel = vec2<f32>(1.0) / vec2<f32>(textureDimensions(shadow_atlas));
    let half_offset = f32(n - 1) * 0.5;
    var lit = 0.0;
    for (var dy = 0i; dy < n; dy++) {
        for (var dx = 0i; dx < n; dx++) {
            let offset = (vec2<f32>(f32(dx), f32(dy)) - half_offset) * texel;
            let d = textureSampleLevel(shadow_atlas, shadow_sampler, base_tile_uv + offset, 0.0).r;
            if ndc.z <= d + shadow_params.bias { lit += 1.0; }
        }
    }
    lit /= f32(n * n);
    return mix(1.0 - shadow_params.intensity, 1.0, lit);
}

// Compute ambient occlusion from multi-directional shadow maps
fn compute_shadow_ao(world_pos: vec3<f32>) -> f32 {
    if shadow_params.shadow_count == 0u {
        return 1.0;
    }

    var lit_count = 0u;
    let count = shadow_params.shadow_count;
    let grid = shadow_params.grid_size;

    for (var i = 0u; i < count; i++) {
        let shadow_pos = shadow_matrices[i] * vec4<f32>(world_pos, 1.0);
        let ndc = shadow_pos.xyz / shadow_pos.w;
        let uv = vec2<f32>(ndc.x * 0.5 + 0.5, 1.0 - (ndc.y * 0.5 + 0.5));

        if uv.x < 0.0 || uv.x > 1.0 || uv.y < 0.0 || uv.y > 1.0 {
            lit_count += 1u;
            continue;
        }

        let col = i % grid;
        let row = i / grid;
        let tile_uv = vec2<f32>(
            (f32(col) + uv.x) / f32(grid),
            (f32(row) + uv.y) / f32(grid),
        );

        let shadow_depth = textureSampleLevel(shadow_atlas, shadow_sampler, tile_uv, 0.0).r;
        let fragment_depth = ndc.z;

        if fragment_depth <= shadow_depth + shadow_params.bias {
            lit_count += 1u;
        }
    }

    let fraction_lit = f32(lit_count) / f32(count);
    return mix(1.0, fraction_lit, shadow_params.intensity);
}
