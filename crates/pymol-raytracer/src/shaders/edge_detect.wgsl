// PyMOL-accurate edge detection shader
// Pass 2 of multi-pass edge detection for ray_trace_mode 1, 2, 3
//
// This implements PyMOL's algorithm from Ray.cpp:
// 1. Compute depth gradients per pixel
// 2. Compare gradient changes with neighbors (second derivative)
// 3. Edge detected based on slope_factor, depth_factor thresholds
// 4. Apply edge coherence filter to remove noise

struct EdgeParams {
    viewport: vec2<f32>,
    slope_factor: f32,      // Gradient magnitude difference threshold
    depth_factor: f32,      // Gradient direction difference threshold  
    disco_factor: f32,      // Gradient discontinuity threshold
    gain: f32,              // Controls sensitivity AND line thickness
    _pad0: f32,
    _pad1: f32,
}

@group(0) @binding(0) var depth_tex: texture_2d<f32>;
@group(0) @binding(1) var normal_tex: texture_2d<f32>;
@group(0) @binding(2) var edge_output: texture_storage_2d<r32float, write>;
@group(0) @binding(3) var<uniform> params: EdgeParams;

// Sample depth with bounds checking
fn sample_depth(coord: vec2<i32>) -> f32 {
    let clamped = clamp(coord, vec2<i32>(0), vec2<i32>(params.viewport) - vec2<i32>(1));
    return textureLoad(depth_tex, clamped, 0).r;
}

// Check if a depth value represents background
fn is_background(d: f32) -> bool {
    return d <= 0.0 || d >= 0.999;
}

// Gaussian-weighted depth sampling for smoother gradients
fn sample_depth_smooth(coord: vec2<i32>) -> f32 {
    var sum = 0.0;
    var weight = 0.0;
    
    // 3x3 Gaussian kernel weights
    let w = array<f32, 9>(
        1.0, 2.0, 1.0,
        2.0, 4.0, 2.0,
        1.0, 2.0, 1.0
    );
    
    for (var dy = -1; dy <= 1; dy++) {
        for (var dx = -1; dx <= 1; dx++) {
            let d = sample_depth(coord + vec2<i32>(dx, dy));
            if !is_background(d) {
                let idx = (dy + 1) * 3 + (dx + 1);
                sum += d * w[idx];
                weight += w[idx];
            }
        }
    }
    
    if weight > 0.0 {
        return sum / weight;
    }
    return sample_depth(coord);
}

// Compute depth gradient using Sobel operator on smoothed depth
fn compute_gradient_sobel(coord: vec2<i32>) -> vec3<f32> {
    let d = sample_depth(coord);
    
    if is_background(d) {
        return vec3<f32>(0.0, 0.0, 0.0);
    }
    
    // Sample 3x3 neighborhood with smoothing
    let d_tl = sample_depth_smooth(coord + vec2<i32>(-1, -1));
    let d_t  = sample_depth_smooth(coord + vec2<i32>(0, -1));
    let d_tr = sample_depth_smooth(coord + vec2<i32>(1, -1));
    let d_l  = sample_depth_smooth(coord + vec2<i32>(-1, 0));
    let d_r  = sample_depth_smooth(coord + vec2<i32>(1, 0));
    let d_bl = sample_depth_smooth(coord + vec2<i32>(-1, 1));
    let d_b  = sample_depth_smooth(coord + vec2<i32>(0, 1));
    let d_br = sample_depth_smooth(coord + vec2<i32>(1, 1));
    
    // Sobel X: detect vertical edges
    var gx = 0.0;
    var wx = 0.0;
    if !is_background(d_tl) { gx -= d_tl; wx += 1.0; }
    if !is_background(d_l)  { gx -= 2.0 * d_l; wx += 2.0; }
    if !is_background(d_bl) { gx -= d_bl; wx += 1.0; }
    if !is_background(d_tr) { gx += d_tr; wx += 1.0; }
    if !is_background(d_r)  { gx += 2.0 * d_r; wx += 2.0; }
    if !is_background(d_br) { gx += d_br; wx += 1.0; }
    if wx > 0.0 { gx /= (wx * 0.25); }
    
    // Sobel Y: detect horizontal edges
    var gy = 0.0;
    var wy = 0.0;
    if !is_background(d_tl) { gy -= d_tl; wy += 1.0; }
    if !is_background(d_t)  { gy -= 2.0 * d_t; wy += 2.0; }
    if !is_background(d_tr) { gy -= d_tr; wy += 1.0; }
    if !is_background(d_bl) { gy += d_bl; wy += 1.0; }
    if !is_background(d_b)  { gy += 2.0 * d_b; wy += 2.0; }
    if !is_background(d_br) { gy += d_br; wy += 1.0; }
    if wy > 0.0 { gy /= (wy * 0.25); }
    
    // Apply gain
    gx *= params.gain;
    gy *= params.gain;
    
    let mag = sqrt(gx * gx + gy * gy);
    return vec3<f32>(gx, gy, mag);
}

// Detect silhouette edges (where surface meets background)
fn detect_silhouette(coord: vec2<i32>, radius: i32) -> f32 {
    let center_depth = sample_depth(coord);
    
    if is_background(center_depth) {
        return 0.0;
    }
    
    // Check neighbors within radius for background
    for (var dy = -radius; dy <= radius; dy++) {
        for (var dx = -radius; dx <= radius; dx++) {
            if dx == 0 && dy == 0 { continue; }
            let neighbor_depth = sample_depth(coord + vec2<i32>(dx, dy));
            if is_background(neighbor_depth) {
                return 1.0;
            }
        }
    }
    
    return 0.0;
}

// PyMOL-style gradient comparison edge detection
fn detect_edge_gradient(coord: vec2<i32>) -> f32 {
    let center_depth = sample_depth(coord);
    
    if is_background(center_depth) {
        return 0.0;
    }
    
    let g = compute_gradient_sobel(coord);
    
    var max_slope = 0.0;
    var max_depth = 0.0;
    var min_dot = 1.0;
    
    let eps = 0.0001;
    
    // 8-connected neighbors
    let offsets = array<vec2<i32>, 8>(
        vec2<i32>(-1, -1), vec2<i32>(0, -1), vec2<i32>(1, -1),
        vec2<i32>(-1, 0),                     vec2<i32>(1, 0),
        vec2<i32>(-1, 1),  vec2<i32>(0, 1),  vec2<i32>(1, 1)
    );
    
    for (var i = 0u; i < 8u; i++) {
        let neighbor_coord = coord + offsets[i];
        let neighbor_depth = sample_depth(neighbor_coord);
        
        if is_background(neighbor_depth) {
            continue;
        }
        
        let pg = compute_gradient_sobel(neighbor_coord);
        
        // Gradient magnitude difference
        let slope_diff = abs(g.z - pg.z);
        max_slope = max(max_slope, slope_diff);
        
        // Gradient direction difference (squared)
        let ddx = g.x - pg.x;
        let ddy = g.y - pg.y;
        let depth_diff = ddx * ddx + ddy * ddy;
        max_depth = max(max_depth, depth_diff);
        
        // Dot product of normalized gradients
        if g.z > eps && pg.z > eps {
            let gn = vec2<f32>(g.x / g.z, g.y / g.z);
            let pn = vec2<f32>(pg.x / pg.z, pg.y / pg.z);
            let d = gn.x * pn.x + gn.y * pn.y;
            min_dot = min(min_dot, d);
        }
    }
    
    var edge = 0.0;
    
    // Gradient magnitude discontinuity
    if max_slope > params.slope_factor {
        edge = max(edge, smoothstep(params.slope_factor, params.slope_factor * 2.0, max_slope));
    }
    
    // Gradient direction discontinuity
    if max_depth > params.depth_factor {
        edge = max(edge, smoothstep(params.depth_factor, params.depth_factor * 2.0, max_depth));
    }
    
    // Gradient reversal
    if min_dot < -params.disco_factor {
        edge = max(edge, smoothstep(-params.disco_factor, -params.disco_factor * 2.0, -min_dot));
    }
    
    return edge;
}

// Count edge neighbors for coherence filtering
// Returns the number of neighboring pixels that also have edges
fn count_edge_neighbors(coord: vec2<i32>, threshold: f32) -> i32 {
    var count = 0;
    
    // Check 8 neighbors
    let offsets = array<vec2<i32>, 8>(
        vec2<i32>(-1, -1), vec2<i32>(0, -1), vec2<i32>(1, -1),
        vec2<i32>(-1, 0),                     vec2<i32>(1, 0),
        vec2<i32>(-1, 1),  vec2<i32>(0, 1),  vec2<i32>(1, 1)
    );
    
    for (var i = 0u; i < 8u; i++) {
        let neighbor_coord = coord + offsets[i];
        let neighbor_depth = sample_depth(neighbor_coord);
        
        if !is_background(neighbor_depth) {
            let e = detect_edge_gradient(neighbor_coord);
            if e > threshold {
                count += 1;
            }
        }
    }
    
    return count;
}

@compute @workgroup_size(8, 8, 1)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let coord = vec2<i32>(gid.xy);
    
    if f32(coord.x) >= params.viewport.x || f32(coord.y) >= params.viewport.y {
        return;
    }
    
    let center_depth = sample_depth(coord);
    
    if is_background(center_depth) {
        textureStore(edge_output, coord, vec4<f32>(0.0, 0.0, 0.0, 1.0));
        return;
    }
    
    // Fixed dilation radius of 1 for consistent line thickness
    let dilation_radius = 1;
    
    // Detect edges using gradient-based method
    var edge = detect_edge_gradient(coord);
    
    // Only dilate weak edges (not already strong) - prevents over-thickening
    if edge > 0.1 && edge < 0.5 {
        // Check immediate neighbors only for slight thickening
        let offsets = array<vec2<i32>, 4>(
            vec2<i32>(0, -1), vec2<i32>(-1, 0), vec2<i32>(1, 0), vec2<i32>(0, 1)
        );
        for (var i = 0u; i < 4u; i++) {
            let sample_coord = coord + offsets[i];
            let sample_d = sample_depth(sample_coord);
            if !is_background(sample_d) {
                let e = detect_edge_gradient(sample_coord);
                if e > 0.5 {
                    edge = max(edge, 0.6); // Boost weak edges near strong edges
                }
            }
        }
    }
    
    // Add silhouette edges (object boundary) with minimal dilation
    edge = max(edge, detect_silhouette(coord, dilation_radius));
    
    // Apply coherence filter: only keep edges that have at least 1 neighboring edge
    // This removes completely isolated noise pixels while preserving contour lines
    if edge > 0.3 && edge < 0.8 {
        let neighbor_count = count_edge_neighbors(coord, 0.25);
        if neighbor_count < 1 {
            edge = 0.0; // Remove isolated edge pixels (noise)
        }
    }
    
    // Final threshold to make lines solid
    if edge > 0.4 {
        edge = 1.0;
    } else {
        edge = 0.0;
    }
    
    textureStore(edge_output, coord, vec4<f32>(edge, 0.0, 0.0, 1.0));
}
