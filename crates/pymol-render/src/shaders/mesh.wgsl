// Mesh shader for surface and mesh rendering.
// Requires: common.wgsl + lighting.wgsl prepended via concat!(include_str!(...))

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
    @location(4) world_pos: vec3<f32>,
}

@vertex
fn vs_main(input: VertexInput) -> VertexOutput {
    var output: VertexOutput;

    let world_pos = vec4<f32>(input.position, 1.0);
    let view_pos = uniforms.view * world_pos;

    output.clip_position = uniforms.view_proj * world_pos;
    output.color = input.color;
    output.normal_view = (uniforms.view * vec4<f32>(input.normal, 0.0)).xyz;
    output.position_view = view_pos.xyz;
    output.view_depth = -view_pos.z;
    output.world_pos = input.position;

    return output;
}

@fragment
fn fs_main(input: VertexOutput, @builtin(front_facing) front_facing: bool) -> @location(0) vec4<f32> {
    var normal = normalize(input.normal_view);
    let view_dir = normalize(-input.position_view);

    if dot(normal, view_dir) < 0.0 {
        normal = -normal;
    }

    // Calculate lighting — branch on shadow mode
    var color: vec3<f32>;
    if shadow_params.mode == 2u {
        color = pymol_lighting_shadowed(normal, view_dir, input.color.rgb, input.world_pos);
    } else {
        color = pymol_lighting(normal, view_dir, input.color.rgb);
        let ao = compute_shadow_ao(input.world_pos);
        color *= ao;
    }

    color = apply_depth_cue(color, input.view_depth);
    color = apply_fog(color, input.view_depth);

    return vec4<f32>(color, input.color.a);
}
