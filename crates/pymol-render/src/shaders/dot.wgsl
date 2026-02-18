// Dot shader for dot surface representation
// Renders dots as small billboard circles
// Note: Dots don't use full lighting, but uniforms must match other shaders

struct GlobalUniforms {
    view_proj: mat4x4<f32>,
    view: mat4x4<f32>,
    view_inv: mat4x4<f32>,
    proj: mat4x4<f32>,
    camera_pos: vec4<f32>,
    // Multi-light support (not used by dots, but must match layout)
    light_dirs: array<vec4<f32>, 9>,
    light_count: i32,
    spec_count: i32,  // Number of lights contributing specular (-1 = all)
    _pad_light_0: i32,
    _pad_light_1: i32,
    // Headlight (camera light) parameters
    ambient: f32,
    direct: f32,
    spec_direct: f32,
    spec_direct_power: f32,
    // Positional light parameters
    reflect: f32,
    specular: f32,
    shininess: f32,
    _pad0: f32,
    // Fog parameters
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

// Shadow sampling (group 1) — multi-directional shadow map AO
struct ShadowParams {
    shadow_count: u32,
    grid_size: u32,
    bias: f32,
    intensity: f32,
}

@group(1) @binding(0) var shadow_atlas: texture_2d<f32>;
@group(1) @binding(1) var shadow_sampler: sampler;
@group(1) @binding(2) var<uniform> shadow_params: ShadowParams;
@group(1) @binding(3) var<storage, read> shadow_matrices: array<mat4x4<f32>>;

// Billboard vertex
struct BillboardVertex {
    @location(10) offset: vec2<f32>,
}

// Instance data (per dot)
struct DotInstance {
    @location(0) position: vec3<f32>,
    @location(1) size: f32,
    @location(2) color: vec4<f32>,
}

struct VertexOutput {
    @builtin(position) clip_position: vec4<f32>,
    @location(0) color: vec4<f32>,
    @location(1) uv: vec2<f32>,
    @location(2) view_depth: f32,
    @location(3) world_pos: vec3<f32>,
}

@vertex
fn vs_main(billboard: BillboardVertex, instance: DotInstance) -> VertexOutput {
    var output: VertexOutput;

    // Transform dot position to view space
    let pos_world = vec4<f32>(instance.position, 1.0);
    let pos_view = (uniforms.view * pos_world).xyz;

    // Create billboard quad in view space
    // Size is in screen pixels, convert to view space units
    let pixel_size = instance.size * 0.01; // Approximate scaling
    let billboard_pos = pos_view + vec3<f32>(billboard.offset * pixel_size, 0.0);

    output.clip_position = uniforms.proj * vec4<f32>(billboard_pos, 1.0);
    output.color = instance.color;
    output.uv = billboard.offset;
    output.view_depth = -pos_view.z;
    output.world_pos = instance.position;

    return output;
}

// Apply fog
fn apply_fog(color: vec3<f32>, depth: f32) -> vec3<f32> {
    if uniforms.fog_density <= 0.0 {
        return color;
    }
    let fog_factor = clamp((uniforms.fog_end - depth) / (uniforms.fog_end - uniforms.fog_start), 0.0, 1.0);
    return mix(uniforms.fog_color.rgb, color, fog_factor);
}

// Apply depth cue
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

@fragment
fn fs_main(input: VertexOutput) -> @location(0) vec4<f32> {
    // Discard pixels outside the circle
    let dist = length(input.uv);
    if dist > 1.0 {
        discard;
    }

    // Soft edge for anti-aliasing
    let alpha = 1.0 - smoothstep(0.8, 1.0, dist);

    var color = input.color.rgb;

    // Apply multi-directional shadow AO
    let ao = compute_shadow_ao(input.world_pos);
    color *= ao;

    // Apply depth cue and fog
    color = apply_depth_cue(color, input.view_depth);
    color = apply_fog(color, input.view_depth);

    return vec4<f32>(color, input.color.a * alpha);
}

// Multi-directional shadow AO: for each shadow direction, project fragment
// into shadow map tile and compare depths to determine occlusion.
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

        // NDC to UV (0-1 range)
        let uv = vec2<f32>(ndc.x * 0.5 + 0.5, 1.0 - (ndc.y * 0.5 + 0.5));

        // Skip if outside shadow map
        if uv.x < 0.0 || uv.x > 1.0 || uv.y < 0.0 || uv.y > 1.0 {
            lit_count += 1u;
            continue;
        }

        // Compute atlas UV for tile i
        let col = i % grid;
        let row = i / grid;
        let tile_uv = vec2<f32>(
            (f32(col) + uv.x) / f32(grid),
            (f32(row) + uv.y) / f32(grid),
        );

        let shadow_depth = textureSample(shadow_atlas, shadow_sampler, tile_uv).r;
        let fragment_depth = ndc.z;

        if fragment_depth <= shadow_depth + shadow_params.bias {
            lit_count += 1u;
        }
    }

    let fraction_lit = f32(lit_count) / f32(count);
    return mix(1.0, fraction_lit, shadow_params.intensity);
}
