// Mesh shader with PyMOL multi-light model for surface and mesh rendering
//
// PyMOL Lighting Model:
// - light_count = 1: Ambient only (no directional lights)
// - light_count = 2: Ambient + 1 directional light
// - light_count = 3-10: Ambient + (N-1) directional lights
// - Headlight (direct): Always from camera direction, ensures front-facing surfaces are lit
// - Positional lights (reflect): World-space directional lights for depth/shadow cues

const MAX_LIGHTS: u32 = 9u;

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

// PyMOL multi-light model
// - Headlight (from camera direction) controlled by 'direct' setting
// - Positional lights controlled by light_count and 'reflect'/'specular' settings
// - spec_count controls how many positional lights contribute specular (-1 = all)
fn pymol_lighting(normal: vec3<f32>, view_dir: vec3<f32>, base_color: vec3<f32>) -> vec3<f32> {
    let n = normalize(normal);
    let v = normalize(view_dir);
    
    // Ambient contribution
    var color = base_color * uniforms.ambient;
    
    // === HEADLIGHT (camera light) ===
    // Light comes from camera direction - always active
    let headlight_ndotl = max(dot(n, v), 0.0);
    color += base_color * headlight_ndotl * uniforms.direct;
    
    // Headlight specular (Blinn-Phong)
    if uniforms.spec_direct > 0.0 && headlight_ndotl > 0.0 {
        let spec = pow(headlight_ndotl, uniforms.spec_direct_power);
        color += vec3<f32>(1.0) * spec * uniforms.spec_direct;
    }
    
    // === POSITIONAL LIGHTS ===
    // light_count = 1: ambient only (no positional lights)
    // light_count = 2: 1 positional light
    // light_count = N: N-1 positional lights
    let num_pos_lights = max(uniforms.light_count - 1, 0);
    
    // Resolve spec_count: -1 means all positional lights contribute specular
    let effective_spec_count = select(uniforms.spec_count, num_pos_lights, uniforms.spec_count < 0);
    
    for (var i = 0; i < num_pos_lights; i++) {
        // Light directions are already in view space (PyMOL convention)
        // Negate because PyMOL defines direction FROM light TO scene
        let l = normalize(-uniforms.light_dirs[i].xyz);
        
        // Diffuse contribution (always applied)
        let ndotl = max(dot(n, l), 0.0);
        color += base_color * ndotl * uniforms.reflect;
        
        // Specular contribution (Blinn-Phong) - only if i < effective_spec_count
        if uniforms.specular > 0.0 && ndotl > 0.0 && i < effective_spec_count {
            let h = normalize(l + v);
            let ndoth = max(dot(n, h), 0.0);
            let spec = pow(ndoth, uniforms.shininess);
            color += vec3<f32>(1.0) * spec * uniforms.specular;
        }
    }
    
    return color;
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
    
    // Calculate lighting using PyMOL multi-light model
    var color = pymol_lighting(normal, view_dir, input.color.rgb);
    
    // Apply depth cue and fog
    color = apply_depth_cue(color, input.view_depth);
    color = apply_fog(color, input.view_depth);
    
    return vec4<f32>(color, input.color.a);
}
