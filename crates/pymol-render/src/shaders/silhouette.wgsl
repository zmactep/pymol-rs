// ChimeraX-style depth-only silhouette edge detection shader
//
// Renders black outlines at depth discontinuities by sampling the depth buffer
// in 4 directions. Thickness is controlled by sampling distance, not dilation.

struct SilhouetteParams {
    // (1/width, 1/height, thickness, depth_jump)
    step_and_params: vec4<f32>,
    color: vec4<f32>,
}

@group(0) @binding(0) var depth_tex: texture_2d<f32>;
@group(0) @binding(1) var depth_sampler: sampler;
@group(0) @binding(2) var<uniform> params: SilhouetteParams;

struct VertexOutput {
    @builtin(position) position: vec4<f32>,
    @location(0) uv: vec2<f32>,
}

// Fullscreen triangle from vertex index (3 vertices, no buffer needed)
@vertex
fn vs_main(@builtin(vertex_index) vi: u32) -> VertexOutput {
    var out: VertexOutput;
    // Positions: (-1,-1), (3,-1), (-1,3) — covers entire screen
    let x = f32(i32(vi & 1u) * 4 - 1);
    let y = f32(i32(vi >> 1u) * 4 - 1);
    out.position = vec4<f32>(x, y, 0.0, 1.0);
    out.uv = vec2<f32>(x * 0.5 + 0.5, 1.0 - (y * 0.5 + 0.5));
    return out;
}

@fragment
fn fs_main(in: VertexOutput) -> @location(0) vec4<f32> {
    let center_depth = textureSample(depth_tex, depth_sampler, in.uv).r;

    // Background pixels: no edge
    if center_depth >= 0.999 {
        discard;
    }

    let thickness = params.step_and_params.z;
    let depth_jump = params.step_and_params.w;
    let dx = params.step_and_params.x * thickness;
    let dy = params.step_and_params.y * thickness;

    // Sample 4 neighbors
    let d_left  = textureSample(depth_tex, depth_sampler, in.uv + vec2<f32>(-dx, 0.0)).r;
    let d_right = textureSample(depth_tex, depth_sampler, in.uv + vec2<f32>( dx, 0.0)).r;
    let d_up    = textureSample(depth_tex, depth_sampler, in.uv + vec2<f32>(0.0, -dy)).r;
    let d_down  = textureSample(depth_tex, depth_sampler, in.uv + vec2<f32>(0.0,  dy)).r;

    // Perspective correction: depth differences are smaller for farther objects
    let threshold = depth_jump * (1.0 - center_depth * 0.9);

    // Edge if a neighbor is significantly farther (larger depth) or is background
    let is_edge = (d_left - center_depth > threshold) ||
                  (d_right - center_depth > threshold) ||
                  (d_up - center_depth > threshold) ||
                  (d_down - center_depth > threshold) ||
                  (d_left >= 0.999 && center_depth < 0.99) ||
                  (d_right >= 0.999 && center_depth < 0.99) ||
                  (d_up >= 0.999 && center_depth < 0.99) ||
                  (d_down >= 0.999 && center_depth < 0.99);

    if !is_edge {
        discard;
    }

    return params.color;
}
