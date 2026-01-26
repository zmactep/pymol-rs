// Common shader utilities for PyMOL-RS rendering
// This file is included in other shaders via copy-paste for now
// (WGSL doesn't have a native include mechanism)
//
// PyMOL Lighting Model:
// - Headlight (direct): Always from camera direction, ensures front-facing surfaces are lit
// - Positional light (reflect): World-space directional light for depth/shadow cues

// Global uniforms shared by all shaders
struct GlobalUniforms {
    view_proj: mat4x4<f32>,
    view: mat4x4<f32>,
    view_inv: mat4x4<f32>,
    proj: mat4x4<f32>,
    camera_pos: vec4<f32>,
    light_dir: vec4<f32>,
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

// PyMOL dual-light model
// Combines headlight (camera light) with positional directional light
//
// Arguments:
// - normal: Surface normal in view space
// - view_dir: Direction from fragment to camera in view space
// - light_dir_view: Positional light direction transformed to view space
// - base_color: Surface base color
// - uniforms: Global uniform buffer reference (for lighting parameters)
fn pymol_lighting(
    normal: vec3<f32>,
    view_dir: vec3<f32>,
    light_dir_view: vec3<f32>,
    base_color: vec3<f32>,
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
    let ambient_color = base_color * ambient;
    
    // === CAMERA FILL LIGHT ===
    // Constant illumination for all visible surfaces - doesn't depend on normal
    // This provides consistent base lighting regardless of surface orientation
    let fill_factor = 0.5; // 50% of headlight is constant fill
    let camera_fill = base_color * direct * fill_factor;
    
    // === HEADLIGHT (camera light, normal-dependent) ===
    // Use actual view direction as headlight - varies with surface angle for depth
    let headlight_dir = v;
    let headlight_ndotl = max(dot(n, headlight_dir), 0.0);
    let headlight_diffuse = base_color * headlight_ndotl * direct * (1.0 - fill_factor);
    
    // Headlight specular (Blinn-Phong)
    var headlight_specular = vec3<f32>(0.0);
    if spec_direct > 0.0 && headlight_ndotl > 0.0 {
        let headlight_ndoth = max(dot(n, v), 0.0);
        let headlight_spec = pow(headlight_ndoth, spec_direct_power);
        headlight_specular = vec3<f32>(1.0) * headlight_spec * spec_direct;
    }
    
    // === POSITIONAL LIGHT (world-space directional light) ===
    let pos_light_dir = normalize(-light_dir_view);
    
    // Positional light diffuse
    let pos_ndotl = max(dot(n, pos_light_dir), 0.0);
    let pos_diffuse = base_color * pos_ndotl * reflect;
    
    // Positional light specular (Blinn-Phong)
    var pos_specular = vec3<f32>(0.0);
    if specular > 0.0 && pos_ndotl > 0.0 {
        let pos_h = normalize(pos_light_dir + v);
        let pos_ndoth = max(dot(n, pos_h), 0.0);
        let pos_spec = pow(pos_ndoth, shininess);
        pos_specular = vec3<f32>(1.0) * pos_spec * specular;
    }
    
    // Combine all contributions
    return ambient_color + camera_fill + headlight_diffuse + headlight_specular + pos_diffuse + pos_specular;
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
