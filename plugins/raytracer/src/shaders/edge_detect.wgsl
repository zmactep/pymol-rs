// Silhouette-style depth-only edge detection shader for ray_trace_mode 1, 2, 3
//
// Uses the same algorithm as the real-time silhouette pipeline:
// - Sample 4 depth neighbors at configurable distance (thickness)
// - Perspective-corrected depth threshold
// - Edge if neighbor is significantly farther or is background
//
// Parameters:
// - thickness: Sampling distance in pixels (higher = thicker outlines)
// - depth_jump: Depth discontinuity threshold (lower = more sensitive)

struct EdgeParams {
    viewport: vec2<f32>,
    thickness: f32,
    depth_jump: f32,
}

@group(0) @binding(0) var depth_tex: texture_2d<f32>;
@group(0) @binding(1) var edge_output: texture_storage_2d<r32float, write>;
@group(0) @binding(2) var<uniform> params: EdgeParams;

// Sample depth with bounds checking
fn sample_depth(coord: vec2<i32>) -> f32 {
    let clamped = clamp(coord, vec2<i32>(0), vec2<i32>(params.viewport) - vec2<i32>(1));
    return textureLoad(depth_tex, clamped, 0).r;
}

// Check if a depth value represents background
fn is_background(d: f32) -> bool {
    return d <= 0.0 || d >= 0.999;
}

@compute @workgroup_size(8, 8, 1)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let coord = vec2<i32>(gid.xy);

    // Bounds check
    if f32(coord.x) >= params.viewport.x || f32(coord.y) >= params.viewport.y {
        return;
    }

    let center_depth = sample_depth(coord);

    // Background pixels: no edge
    if is_background(center_depth) {
        textureStore(edge_output, coord, vec4<f32>(0.0, 0.0, 0.0, 1.0));
        return;
    }

    let step = i32(max(params.thickness, 1.0));

    // Sample 4 neighbors at thickness distance
    let d_left  = sample_depth(coord + vec2<i32>(-step, 0));
    let d_right = sample_depth(coord + vec2<i32>( step, 0));
    let d_up    = sample_depth(coord + vec2<i32>(0, -step));
    let d_down  = sample_depth(coord + vec2<i32>(0,  step));

    // Perspective correction: depth differences are smaller for farther objects
    let threshold = params.depth_jump * (1.0 - center_depth * 0.9);

    // Edge if a neighbor is significantly farther (larger depth) or is background
    let is_edge = (d_left - center_depth > threshold) ||
                  (d_right - center_depth > threshold) ||
                  (d_up - center_depth > threshold) ||
                  (d_down - center_depth > threshold) ||
                  (is_background(d_left) && center_depth < 0.99) ||
                  (is_background(d_right) && center_depth < 0.99) ||
                  (is_background(d_up) && center_depth < 0.99) ||
                  (is_background(d_down) && center_depth < 0.99);

    var edge = 0.0;
    if is_edge {
        edge = 1.0;
    }

    textureStore(edge_output, coord, vec4<f32>(edge, 0.0, 0.0, 1.0));
}
