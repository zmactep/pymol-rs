// Line shader for bond and wireframe rendering
// Note: Lines don't use lighting, but uniforms must match other shaders

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
    @location(1) color: vec4<f32>,
}

struct VertexOutput {
    @builtin(position) clip_position: vec4<f32>,
    @location(0) color: vec4<f32>,
    @location(1) view_depth: f32,
}

@vertex
fn vs_main(input: VertexInput) -> VertexOutput {
    var output: VertexOutput;
    
    let world_pos = vec4<f32>(input.position, 1.0);
    let view_pos = uniforms.view * world_pos;
    
    output.clip_position = uniforms.view_proj * world_pos;
    output.color = input.color;
    output.view_depth = -view_pos.z; // Positive depth into screen
    
    return output;
}

// Apply fog based on depth
fn apply_fog(color: vec3<f32>, depth: f32, fog_start: f32, fog_end: f32, fog_color: vec3<f32>) -> vec3<f32> {
    if uniforms.fog_density <= 0.0 {
        return color;
    }
    let fog_factor = clamp((fog_end - depth) / (fog_end - fog_start), 0.0, 1.0);
    return mix(fog_color, color, fog_factor);
}

// Apply depth cue
fn apply_depth_cue(color: vec3<f32>, depth: f32, factor: f32) -> vec3<f32> {
    if factor <= 0.0 {
        return color;
    }
    let near = uniforms.clip_planes.x;
    let far = uniforms.clip_planes.y;
    let normalized_depth = clamp((depth - near) / (far - near), 0.0, 1.0);
    let cue = 1.0 - factor * normalized_depth * 0.5;
    return color * cue;
}

@fragment
fn fs_main(input: VertexOutput) -> @location(0) vec4<f32> {
    var color = input.color.rgb;
    
    // Apply depth cue
    color = apply_depth_cue(color, input.view_depth, uniforms.depth_cue);
    
    // Apply fog
    color = apply_fog(
        color,
        input.view_depth,
        uniforms.fog_start,
        uniforms.fog_end,
        uniforms.fog_color.rgb
    );
    
    return vec4<f32>(color, input.color.a);
}
