// Dot shader for dot surface representation — small billboard circles.
// Requires: common.wgsl prepended via concat!(include_str!(...))

struct BillboardVertex {
    @location(10) offset: vec2<f32>,
}

struct DotInstance {
    @location(0) position: vec3<f32>,
    @location(1) size: f32,
    @location(2) color: vec4<f32>,
}

struct VertexOutput {
    @builtin(position) clip_position: vec4<f32>,
    @location(0) color: vec4<f32>,
    @location(1) uv: vec2<f32>,
    @location(2) view_depth: f32,
    @location(3) world_pos: vec3<f32>,
}

@vertex
fn vs_main(billboard: BillboardVertex, instance: DotInstance) -> VertexOutput {
    var output: VertexOutput;

    let pos_world = vec4<f32>(instance.position, 1.0);
    let pos_view = (uniforms.view * pos_world).xyz;

    let pixel_size = instance.size * 0.01;
    let billboard_pos = pos_view + vec3<f32>(billboard.offset * pixel_size, 0.0);

    output.clip_position = uniforms.proj * vec4<f32>(billboard_pos, 1.0);
    output.color = instance.color;
    output.uv = billboard.offset;
    output.view_depth = -pos_view.z;
    output.world_pos = instance.position;

    return output;
}

@fragment
fn fs_main(input: VertexOutput) -> @location(0) vec4<f32> {
    let dist = length(input.uv);
    if dist > 1.0 {
        discard;
    }

    let alpha = 1.0 - smoothstep(0.8, 1.0, dist);

    var color = input.color.rgb;

    color = apply_depth_cue(color, input.view_depth);
    color = apply_fog(color, input.view_depth);

    return vec4<f32>(color, input.color.a * alpha);
}
