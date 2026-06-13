// WBOIT composite — full-screen resolve of accum/reveal into final color.
//
// final.rgb = accum.rgb / max(accum.a, ε) * (1 - reveal) + dst.rgb * reveal
// final.a   = 1 - reveal
//
// Pipeline writes to the swap-chain target with blend disabled (it's an
// override write per pixel, no incremental blend with previous frame).

@group(0) @binding(0) var accum_tex:   texture_2d<f32>;
@group(0) @binding(1) var reveal_tex:  texture_2d<f32>;
@group(0) @binding(2) var samp:        sampler;

struct VsOut {
    @builtin(position) clip_pos: vec4<f32>,
    @location(0)       uv:       vec2<f32>,
};

@vertex
fn vs_main(@builtin(vertex_index) vid: u32) -> VsOut {
    // Full-screen triangle, no buffer required.
    let p = array<vec2<f32>, 3>(
        vec2<f32>(-1.0, -1.0),
        vec2<f32>( 3.0, -1.0),
        vec2<f32>(-1.0,  3.0),
    );
    let uv = array<vec2<f32>, 3>(
        vec2<f32>(0.0, 1.0),
        vec2<f32>(2.0, 1.0),
        vec2<f32>(0.0, -1.0),
    );
    var out: VsOut;
    out.clip_pos = vec4<f32>(p[vid], 0.0, 1.0);
    out.uv = uv[vid];
    return out;
}

@fragment
fn fs_main(in: VsOut) -> @location(0) vec4<f32> {
    let accum  = textureSample(accum_tex,  samp, in.uv);
    let reveal = textureSample(reveal_tex, samp, in.uv).r;

    let avg = accum.rgb / max(accum.a, 1e-5);
    let coverage = 1.0 - reveal;
    // Composite onto a transparent destination — caller can blend on top.
    return vec4<f32>(avg * coverage, coverage);
}
