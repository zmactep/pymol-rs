// Normal-based edge detection shader for ray_trace_mode 1, 2, 3
// 
// This implements PyMOL's edge detection algorithm using:
// 1. Normal discontinuities - detect surface creases where normals change sharply
// 2. Depth discontinuities - detect silhouettes where depth changes sharply
// 3. Background edges - detect where surface meets background
//
// Parameters:
// - slope_factor: Normal dot product threshold (higher = less sensitive, range 0-1)
// - depth_factor: Depth difference threshold (higher = less sensitive)
// - gain: Edge line thickness control via dilation radius

struct EdgeParams {
    viewport: vec2<f32>,
    slope_factor: f32,      // Normal discontinuity threshold (dot product)
    depth_factor: f32,      // Depth discontinuity threshold
    disco_factor: f32,      // Reserved for future use
    gain: f32,              // Controls edge line thickness
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

// Sample normal with bounds checking
// Normals are stored as [0,1] range, need to convert to [-1,1]
fn sample_normal(coord: vec2<i32>) -> vec3<f32> {
    let clamped = clamp(coord, vec2<i32>(0), vec2<i32>(params.viewport) - vec2<i32>(1));
    let n = textureLoad(normal_tex, clamped, 0).rgb;
    // Convert from [0,1] to [-1,1] and normalize
    return normalize(n * 2.0 - 1.0);
}

// Check if a depth value represents background
fn is_background(d: f32) -> bool {
    return d <= 0.0 || d >= 0.999;
}

// Detect edges at a single pixel using normal and depth discontinuities
fn detect_edge(coord: vec2<i32>) -> f32 {
    let center_depth = sample_depth(coord);
    
    // Background pixels have no edges
    if is_background(center_depth) {
        return 0.0;
    }
    
    let center_normal = sample_normal(coord);
    
    // Check 4-connected neighbors (produces cleaner lines than 8-connected)
    let offsets = array<vec2<i32>, 4>(
        vec2<i32>(0, -1),   // top
        vec2<i32>(-1, 0),   // left
        vec2<i32>(1, 0),    // right
        vec2<i32>(0, 1)     // bottom
    );
    
    for (var i = 0u; i < 4u; i++) {
        let neighbor_coord = coord + offsets[i];
        let neighbor_depth = sample_depth(neighbor_coord);
        
        // Silhouette edge: foreground pixel adjacent to background
        if is_background(neighbor_depth) {
            return 1.0;
        }
        
        // Normal discontinuity edge: surface crease
        // When normals differ significantly, we're at a crease/edge
        let neighbor_normal = sample_normal(neighbor_coord);
        let ndot = dot(center_normal, neighbor_normal);
        if ndot < params.slope_factor {
            return 1.0;
        }
        
        // Depth discontinuity edge: major depth jump
        // This catches edges that normals might miss (e.g., overlapping surfaces)
        let depth_diff = abs(center_depth - neighbor_depth);
        if depth_diff > params.depth_factor {
            return 1.0;
        }
    }
    
    return 0.0;
}

// Dilate edges to control line thickness
// radius controls the dilation amount
fn detect_edge_dilated(coord: vec2<i32>, radius: i32) -> f32 {
    // First check center pixel
    if detect_edge(coord) > 0.5 {
        return 1.0;
    }
    
    // If radius > 0, check if any nearby pixel is an edge
    if radius > 0 {
        for (var dy = -radius; dy <= radius; dy++) {
            for (var dx = -radius; dx <= radius; dx++) {
                if dx == 0 && dy == 0 { continue; }
                
                let sample_coord = coord + vec2<i32>(dx, dy);
                let sample_depth = sample_depth(sample_coord);
                
                // Only dilate from foreground edges
                if !is_background(sample_depth) {
                    if detect_edge(sample_coord) > 0.5 {
                        return 1.0;
                    }
                }
            }
        }
    }
    
    return 0.0;
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
    
    // Calculate dilation radius from gain parameter
    // gain 0.12 (PyMOL default) -> radius 0 (no dilation, 1px lines)
    // gain 0.06 -> radius 1 (slight dilation, ~2px lines)
    // gain 0.03 -> radius 2 (more dilation, ~3px lines)
    let dilation_radius = i32(clamp(floor(0.12 / max(params.gain, 0.01)), 0.0, 2.0));
    
    // Detect edges with optional dilation for line thickness
    let edge = detect_edge_dilated(coord, dilation_radius);
    
    textureStore(edge_output, coord, vec4<f32>(edge, 0.0, 0.0, 1.0));
}
