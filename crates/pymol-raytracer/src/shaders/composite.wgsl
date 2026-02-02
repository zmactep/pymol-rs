// Composite shader for ray trace modes
// Pass 3 of multi-pass edge detection - combines color, edges based on mode

// Uniforms for compositing
struct CompositeParams {
    viewport: vec2<f32>,
    mode: u32,               // 0=normal, 1=color+outline, 2=outline only, 3=quantized+outline
    use_transparent_bg: u32, // 0=opaque, 1=transparent
    edge_color: vec3<f32>,   // Color for edges (usually black)
    quantize_levels: f32,    // Number of quantization levels for mode 3
    bg_color: vec3<f32>,     // Background color
    _pad: f32,
}

@group(0) @binding(0) var color_tex: texture_2d<f32>;
@group(0) @binding(1) var edge_tex: texture_2d<f32>;
@group(0) @binding(2) var depth_tex: texture_2d<f32>;
@group(0) @binding(3) var output: texture_storage_2d<rgba8unorm, write>;
@group(0) @binding(4) var<uniform> params: CompositeParams;

// Quantize color for posterized look (mode 3)
fn quantize_color(color: vec3<f32>, levels: f32) -> vec3<f32> {
    return floor(color * levels + 0.5) / levels;
}

@compute @workgroup_size(8, 8, 1)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let coord = vec2<i32>(gid.xy);
    
    // Bounds check
    if f32(coord.x) >= params.viewport.x || f32(coord.y) >= params.viewport.y {
        return;
    }
    
    let color = textureLoad(color_tex, coord, 0);
    let edge = textureLoad(edge_tex, coord, 0).r;
    let depth = textureLoad(depth_tex, coord, 0).r;
    
    var final_color: vec4<f32>;
    
    // Check if this is a background pixel (no surface hit)
    let is_background = depth <= 0.0 || depth >= 0.999;
    
    if is_background {
        // Background handling
        if params.use_transparent_bg == 1u {
            final_color = vec4<f32>(0.0, 0.0, 0.0, 0.0);
        } else {
            final_color = vec4<f32>(params.bg_color, 1.0);
        }
    } else {
        // Surface pixel - apply mode-specific compositing
        if params.mode == 0u {
            // Mode 0: Normal rendering (no edge effects)
            final_color = color;
        } else if params.mode == 1u {
            // Mode 1: Normal color + outline
            let blended = mix(color.rgb, params.edge_color, edge);
            final_color = vec4<f32>(blended, color.a);
        } else if params.mode == 2u {
            // Mode 2: Outline only
            // PyMOL behavior: edges get edge_color, non-edges get white (or transparent)
            if edge > 0.5 {
                // Edge pixel: draw with edge color
                final_color = vec4<f32>(params.edge_color, 1.0);
            } else if params.use_transparent_bg == 1u {
                // Non-edge surface with transparent bg: make transparent
                final_color = vec4<f32>(0.0, 0.0, 0.0, 0.0);
            } else {
                // Non-edge surface with opaque bg: white
                final_color = vec4<f32>(1.0, 1.0, 1.0, 1.0);
            }
        } else if params.mode == 3u {
            // Mode 3: Quantized colors + outline
            let quantized = quantize_color(color.rgb, params.quantize_levels);
            let blended = mix(quantized, params.edge_color, edge);
            final_color = vec4<f32>(blended, color.a);
        } else {
            // Fallback: just use color
            final_color = color;
        }
    }
    
    // Clamp output
    final_color = clamp(final_color, vec4<f32>(0.0), vec4<f32>(1.0));
    
    textureStore(output, coord, final_color);
}
