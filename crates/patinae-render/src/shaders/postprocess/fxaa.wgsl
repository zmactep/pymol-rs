// Simplified Lottes-style FXAA — fast-approximate anti-aliasing.
// Computes a 5-tap luma cross around the pixel, detects edges, and
// blends toward the perpendicular gradient via bilinear sampling.
//
// Trade-off vs full FXAA 3.11: skips the iterative edge-endpoint
// walk. Good enough on impostor-rep silhouettes (smooth curves where
// most of the AA win comes from cheap blending), saves ~150 lines.

struct FxaaParams {
    /// `(1/width, 1/height)` of the source colour target.
    inv_texel_size: vec2<f32>,
    /// Edges with luma range below this are skipped (anti-grain). Set
    /// from the user-facing `fxaa_edge_min` setting once exposed; for
    /// now a constant 0.0312 (Lottes default).
    edge_min:    f32,
    /// Relative-to-max edge threshold. 0.125 = Lottes default.
    edge_max:    f32,
};

@group(0) @binding(0) var<uniform> params: FxaaParams;
@group(0) @binding(1) var src_color: texture_2d<f32>;
@group(0) @binding(2) var src_sampler: sampler;

fn rgb2luma(rgb: vec3<f32>) -> f32 {
    return sqrt(dot(rgb, vec3<f32>(0.299, 0.587, 0.114)));
}

struct VsOut {
    @builtin(position) clip: vec4<f32>,
    @location(0)       uv:   vec2<f32>,
};

@vertex
fn vs_main(@builtin(vertex_index) vid: u32) -> VsOut {
    var pos = vec2<f32>(-1.0, -1.0);
    if (vid == 1u) { pos = vec2<f32>(3.0, -1.0); }
    if (vid == 2u) { pos = vec2<f32>(-1.0, 3.0); }
    var out: VsOut;
    out.clip = vec4<f32>(pos, 0.0, 1.0);
    out.uv = vec2<f32>(pos.x * 0.5 + 0.5, 1.0 - (pos.y * 0.5 + 0.5));
    return out;
}

@fragment
fn fs_main(input: VsOut) -> @location(0) vec4<f32> {
    let uv = input.uv;
    let inv = params.inv_texel_size;
    let center = textureSampleLevel(src_color, src_sampler, uv, 0.0);

    let luma_c = rgb2luma(center.rgb);
    let luma_n = rgb2luma(textureSampleLevel(src_color, src_sampler, uv + vec2<f32>( 0.0, -inv.y), 0.0).rgb);
    let luma_s = rgb2luma(textureSampleLevel(src_color, src_sampler, uv + vec2<f32>( 0.0,  inv.y), 0.0).rgb);
    let luma_w = rgb2luma(textureSampleLevel(src_color, src_sampler, uv + vec2<f32>(-inv.x, 0.0), 0.0).rgb);
    let luma_e = rgb2luma(textureSampleLevel(src_color, src_sampler, uv + vec2<f32>( inv.x, 0.0), 0.0).rgb);

    let luma_min = min(luma_c, min(min(luma_n, luma_s), min(luma_w, luma_e)));
    let luma_max = max(luma_c, max(max(luma_n, luma_s), max(luma_w, luma_e)));
    let range = luma_max - luma_min;

    // Below the threshold the edge is too soft / too dark to be worth
    // smoothing. Pass through.
    if (range < max(params.edge_min, luma_max * params.edge_max)) {
        return center;
    }

    // Edge orientation: stronger horizontal (left/right luma diff)
    // means the edge is mostly vertical-running ⇒ blend along Y.
    let edge_h = abs(luma_w + luma_e - 2.0 * luma_c);
    let edge_v = abs(luma_n + luma_s - 2.0 * luma_c);
    let is_horizontal = edge_h >= edge_v;

    // Pick the step direction perpendicular to the edge.
    let step_dir = select(
        vec2<f32>(0.0, inv.y),
        vec2<f32>(inv.x, 0.0),
        is_horizontal,
    );

    // Decide which side to bias toward (steeper-gradient neighbour).
    let l1 = select(luma_w, luma_n, is_horizontal);
    let l2 = select(luma_e, luma_s, is_horizontal);
    let g1 = abs(l1 - luma_c);
    let g2 = abs(l2 - luma_c);
    let bias_sign = select(1.0, -1.0, g1 > g2);

    // Sample at ±0.5 texel offset using bilinear — this is the
    // FXAA-lite trick: a single bilinear tap at a non-integer offset
    // does a half-and-half blend for free.
    let sample_uv = uv + step_dir * (0.5 * bias_sign);
    let blended = textureSampleLevel(src_color, src_sampler, sample_uv, 0.0);

    // Mix amount scales with edge strength relative to local luma.
    let blend_amount = clamp(range / (luma_max + 1e-4), 0.0, 0.75);
    let rgb = mix(center.rgb, blended.rgb, blend_amount);
    return vec4<f32>(rgb, center.a);
}
