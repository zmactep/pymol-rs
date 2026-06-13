// Composite the SSAO factor into the host colour target via a
// multiplicative blend. The fragment shader outputs only the AO
// darkening factor; the pipeline's blend state computes
// `output = dst * src` (i.e. `host_color *= ao_factor`).
//
// FS reads `ssao_tex` (R8Unorm, blurred occlusion in [0, 1]).
// `params.intensity` is the user-facing strength multiplier.

struct ComposeParams {
    intensity: f32,
    _pad0:     f32,
    _pad1:     f32,
    _pad2:     f32,
};

@group(0) @binding(0) var<uniform> params: ComposeParams;
@group(0) @binding(1) var ssao_tex:  texture_2d<f32>;

struct VsOut {
    @builtin(position) clip:  vec4<f32>,
    @location(0)       uv:    vec2<f32>,
};

@vertex
fn vs_main(@builtin(vertex_index) vid: u32) -> VsOut {
    // Full-screen triangle.
    var pos = vec2<f32>(-1.0, -1.0);
    if (vid == 1u) { pos = vec2<f32>(3.0, -1.0); }
    if (vid == 2u) { pos = vec2<f32>(-1.0, 3.0); }
    var out: VsOut;
    out.clip = vec4<f32>(pos, 0.0, 1.0);
    // UV in [0, 1], y-flipped to match texture-space convention used
    // elsewhere in patinae-render.
    out.uv = vec2<f32>(pos.x * 0.5 + 0.5, 1.0 - (pos.y * 0.5 + 0.5));
    return out;
}

@fragment
fn fs_main(input: VsOut) -> @location(0) vec4<f32> {
    let dims = textureDimensions(ssao_tex);
    let px = vec2<i32>(input.uv * vec2<f32>(f32(dims.x), f32(dims.y)));
    let ao = textureLoad(ssao_tex, px, 0).r;
    let factor = clamp(1.0 - params.intensity * (1.0 - ao), 0.0, 1.0);
    // Pipeline blend is (src=Zero, dst=Src) → host_color *= rgb. Alpha
    // unaffected.
    return vec4<f32>(factor, factor, factor, 1.0);
}
