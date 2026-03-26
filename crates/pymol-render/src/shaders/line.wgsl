// Line shader for bond and wireframe rendering.
// Requires: common.wgsl prepended via concat!(include_str!(...))

struct VertexInput {
    @location(0) position: vec3<f32>,
    @location(1) color: vec4<f32>,
}

struct VertexOutput {
    @builtin(position) clip_position: vec4<f32>,
    @location(0) color: vec4<f32>,
    @location(1) view_depth: f32,
    @location(2) world_pos: vec3<f32>,
}

@vertex
fn vs_main(input: VertexInput) -> VertexOutput {
    var output: VertexOutput;

    let world_pos = vec4<f32>(input.position, 1.0);
    let view_pos = uniforms.view * world_pos;

    output.clip_position = uniforms.view_proj * world_pos;
    output.color = input.color;
    output.view_depth = -view_pos.z;
    output.world_pos = input.position;

    return output;
}

@fragment
fn fs_main(input: VertexOutput) -> @location(0) vec4<f32> {
    var color = input.color.rgb;

    // Apply shadow — mode-aware
    if shadow_params.mode == 2u {
        color *= sample_light_shadow(input.world_pos, 0u);
    } else {
        let ao = compute_shadow_ao(input.world_pos);
        color *= ao;
    }

    color = apply_depth_cue(color, input.view_depth);
    color = apply_fog(color, input.view_depth);

    return vec4<f32>(color, input.color.a);
}
