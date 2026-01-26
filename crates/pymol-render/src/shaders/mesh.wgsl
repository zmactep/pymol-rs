// Mesh shader with PyMOL dual-light model for surface and mesh rendering
//
// PyMOL Lighting Model:
// - Headlight (direct): Always from camera direction, ensures front-facing surfaces are lit
// - Positional light (reflect): World-space directional light for depth/shadow cues

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

@group(0) @binding(0)
var<uniform> uniforms: GlobalUniforms;

struct VertexInput {
    @location(0) position: vec3<f32>,
    @location(1) normal: vec3<f32>,
    @location(2) color: vec4<f32>,
}

struct VertexOutput {
    @builtin(position) clip_position: vec4<f32>,
    @location(0) color: vec4<f32>,
    @location(1) normal_view: vec3<f32>,
    @location(2) position_view: vec3<f32>,
    @location(3) view_depth: f32,
}

@vertex
fn vs_main(input: VertexInput) -> VertexOutput {
    var output: VertexOutput;
    
    let world_pos = vec4<f32>(input.position, 1.0);
    let view_pos = uniforms.view * world_pos;
    
    output.clip_position = uniforms.view_proj * world_pos;
    output.color = input.color;
    
    // Transform normal to view space (using normal matrix = inverse transpose of view)
    // For now, assume no non-uniform scaling, so we can use view matrix directly
    output.normal_view = (uniforms.view * vec4<f32>(input.normal, 0.0)).xyz;
    output.position_view = view_pos.xyz;
    output.view_depth = -view_pos.z;
    
    return output;
}

// Simplified camera-centric lighting model
// Uses only headlight (from camera) for rotation-independent illumination
fn pymol_lighting(normal: vec3<f32>, view_dir: vec3<f32>, light_dir_view: vec3<f32>, base_color: vec3<f32>) -> vec3<f32> {
    let n = normalize(normal);
    let v = normalize(view_dir);
    
    // Higher ambient for better base illumination
    let ambient_color = base_color * (uniforms.ambient + 0.15);
    
    // === HEADLIGHT (camera light) ===
    // Light comes from camera direction - provides 3D depth while being rotation-independent
    // Combine direct and reflect into a single camera light for consistency
    let headlight_intensity = uniforms.direct + uniforms.reflect;
    let headlight_ndotl = max(dot(n, v), 0.0);
    let headlight_diffuse = base_color * headlight_ndotl * headlight_intensity;
    
    // Headlight specular (Blinn-Phong)
    // For headlight, half vector H = normalize(L + V) = normalize(V + V) = V
    var headlight_specular = vec3<f32>(0.0);
    let spec_intensity = uniforms.spec_direct + uniforms.specular;
    if spec_intensity > 0.0 && headlight_ndotl > 0.0 {
        let spec_power = max(uniforms.spec_direct_power, uniforms.shininess);
        let spec = pow(headlight_ndotl, spec_power);
        headlight_specular = vec3<f32>(1.0) * spec * spec_intensity;
    }
    
    // Combine contributions
    return ambient_color + headlight_diffuse + headlight_specular;
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
fn fs_main(input: VertexOutput, @builtin(front_facing) front_facing: bool) -> @location(0) vec4<f32> {
    // Get the interpolated normal in view space
    var normal = normalize(input.normal_view);
    
    // View direction (from fragment to camera, camera at origin in view space)
    let view_dir = normalize(-input.position_view);
    
    // Two-sided lighting: ensure normal faces the camera
    // If dot(normal, view_dir) < 0, the normal points away from the camera
    // This handles both back-facing triangles and inconsistent vertex normals
    if dot(normal, view_dir) < 0.0 {
        normal = -normal;
    }
    
    // Transform light direction to view space
    let light_dir_view = (uniforms.view * vec4<f32>(uniforms.light_dir.xyz, 0.0)).xyz;
    
    // Calculate lighting using PyMOL dual-light model
    var color = pymol_lighting(normal, view_dir, light_dir_view, input.color.rgb);
    
    // Apply depth cue and fog
    color = apply_depth_cue(color, input.view_depth);
    color = apply_fog(color, input.view_depth);
    
    return vec4<f32>(color, input.color.a);
}
