// Common shader utilities for PyMOL-RS rendering
// This file is included in other shaders via copy-paste for now
// (WGSL doesn't have a native include mechanism)
//
// PyMOL Multi-Light Model:
// - light_count = 1: Ambient only (no directional lights)
// - light_count = 2: Ambient + 1 directional light
// - light_count = 3-10: Ambient + (N-1) directional lights
// - Headlight (direct): Always from camera direction, ensures front-facing surfaces are lit
// - Positional lights (reflect): World-space directional lights for depth/shadow cues

const MAX_LIGHTS: u32 = 9u;

// Global uniforms shared by all shaders
struct GlobalUniforms {
    view_proj: mat4x4<f32>,
    view: mat4x4<f32>,
    view_inv: mat4x4<f32>,
    proj: mat4x4<f32>,
    camera_pos: vec4<f32>,
    // Multi-light support
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

// PyMOL multi-light model (reference implementation)
// - Headlight (from camera direction) controlled by 'direct' setting
// - Positional lights controlled by light_count and 'reflect'/'specular' settings
// - spec_count controls how many positional lights contribute specular (-1 = all)
//
// Arguments:
// - normal: Surface normal in view space
// - view_dir: Direction from fragment to camera in view space
// - base_color: Surface base color
// - view_matrix: View matrix for transforming light directions
// - light_dirs: Array of light directions in world space
// - light_count: Number of active lights (1 = ambient only, 2+ = with positional)
// - spec_count: Number of lights contributing specular (-1 = all positional lights)
// - ambient, direct, reflect, specular, shininess, spec_direct, spec_direct_power: Lighting params
fn pymol_lighting_multi(
    normal: vec3<f32>,
    view_dir: vec3<f32>,
    base_color: vec3<f32>,
    view_matrix: mat4x4<f32>,
    light_dirs: array<vec4<f32>, 9>,
    light_count: i32,
    spec_count: i32,
    ambient: f32,
    direct: f32,
    reflect: f32,
    specular: f32,
    shininess: f32,
    spec_direct: f32,
    spec_direct_power: f32
) -> vec3<f32> {
    let n = normalize(normal);
    let v = normalize(view_dir);
    
    // Ambient contribution
    var color = base_color * ambient;
    
    // === HEADLIGHT (camera light) ===
    // Light comes from camera direction - always active
    let headlight_ndotl = max(dot(n, v), 0.0);
    color += base_color * headlight_ndotl * direct;
    
    // Headlight specular (Blinn-Phong)
    if spec_direct > 0.0 && headlight_ndotl > 0.0 {
        let spec = pow(headlight_ndotl, spec_direct_power);
        color += vec3<f32>(1.0) * spec * spec_direct;
    }
    
    // === POSITIONAL LIGHTS ===
    // light_count = 1: ambient only (no positional lights)
    // light_count = 2: 1 positional light
    // light_count = N: N-1 positional lights
    let num_pos_lights = max(light_count - 1, 0);
    
    // Resolve spec_count: -1 means all positional lights contribute specular
    let effective_spec_count = select(spec_count, num_pos_lights, spec_count < 0);
    
    for (var i = 0; i < num_pos_lights; i++) {
        // Transform light direction to view space and negate (light points toward surface)
        let light_dir_world = light_dirs[i].xyz;
        let light_dir_view = normalize((view_matrix * vec4<f32>(light_dir_world, 0.0)).xyz);
        let l = normalize(-light_dir_view);
        
        // Diffuse contribution (always applied)
        let ndotl = max(dot(n, l), 0.0);
        color += base_color * ndotl * reflect;
        
        // Specular contribution (Blinn-Phong) - only if i < effective_spec_count
        if specular > 0.0 && ndotl > 0.0 && i < effective_spec_count {
            let h = normalize(l + v);
            let ndoth = max(dot(n, h), 0.0);
            let spec = pow(ndoth, shininess);
            color += vec3<f32>(1.0) * spec * specular;
        }
    }
    
    return color;
}

// Apply fog based on depth
fn apply_fog(color: vec3<f32>, depth: f32, fog_start: f32, fog_end: f32, fog_color: vec3<f32>) -> vec3<f32> {
    let fog_factor = clamp((fog_end - depth) / (fog_end - fog_start), 0.0, 1.0);
    return mix(fog_color, color, fog_factor);
}

// Apply depth cue (PyMOL-style darkening with depth)
fn apply_depth_cue(color: vec3<f32>, depth: f32, factor: f32) -> vec3<f32> {
    let cue = 1.0 - factor * (1.0 - depth);
    return color * cue;
}

// Convert linear depth to normalized [0, 1] range
fn linearize_depth(depth: f32, near: f32, far: f32) -> f32 {
    return (2.0 * near) / (far + near - depth * (far - near));
}
