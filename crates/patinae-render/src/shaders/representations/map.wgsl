// Generic map contour rendering. Geometry is renderer-neutral contour output
// from patinae-algos; no SceneStore or atom ownership is involved.

// {{INCLUDE_FRAME}}
// {{INCLUDE_WBOIT}}
// {{INCLUDE_LIGHTING}}

struct MapParams {
    color: vec4<f32>,
    model: mat4x4<f32>,
};

@group(2) @binding(0) var<uniform> map_params: MapParams;

struct MapVertex {
    @location(0) position: vec3<f32>,
    @location(1) normal: vec3<f32>,
};

struct VsOut {
    @builtin(position) clip_position: vec4<f32>,
    @location(0) view_pos: vec3<f32>,
    @location(1) view_normal: vec3<f32>,
    @location(2) world_pos: vec3<f32>,
};

@vertex
fn vs_main(v: MapVertex) -> VsOut {
    var out: VsOut;
    let world = map_params.model * vec4<f32>(v.position, 1.0);
    let view_pos = (frame.view * world).xyz;
    let world_n = normalize((map_params.model * vec4<f32>(v.normal, 0.0)).xyz);
    out.clip_position = frame.proj * vec4<f32>(view_pos, 1.0);
    out.view_pos = view_pos;
    out.view_normal = normalize((frame.view * vec4<f32>(world_n, 0.0)).xyz);
    out.world_pos = world.xyz;
    return out;
}

@fragment
fn fs_main(input: VsOut) -> TranslucentOut {
    let n = normalize(input.view_normal);
    let view_dir = -input.view_pos;
    var lit = classic_lighting(n, view_dir, map_params.color.rgb);
    lit = apply_shadow(lit, input.world_pos);
    lit = apply_fog(lit, -input.view_pos.z);
    lit = apply_depth_cue(lit, -input.view_pos.z);
    return write_translucent(lit, map_params.color.a, -input.view_pos.z, frame.clip.w);
}

@fragment
fn fs_opaque(input: VsOut) -> @location(0) vec4<f32> {
    let n = normalize(input.view_normal);
    let view_dir = -input.view_pos;
    var lit = classic_lighting(n, view_dir, map_params.color.rgb);
    lit = apply_shadow(lit, input.world_pos);
    lit = apply_fog(lit, -input.view_pos.z);
    lit = apply_depth_cue(lit, -input.view_pos.z);
    return vec4<f32>(lit, 1.0);
}
