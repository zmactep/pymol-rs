// Edge-aware bilateral blur for the SSAO output. Run twice
// (horizontal + vertical) with `dir = (1, 0)` / `(0, 1)` selected by
// `params.dir_x / dir_y`.
//
// Per pixel:
//   Read center AO + center depth from `frame.depth`.
//   For each of 4 tap offsets (±1, ±2 pixels in `dir`), read tap AO +
//   tap depth. Weight = gaussian(spatial) × gaussian(depth_diff).
//   Output weighted-average AO.

struct BlurParams {
    /// Step in pixels per tap. `(1, 0)` = horizontal pass, `(0, 1)` = vertical.
    dir:           vec2<i32>,
    /// Depth-aware falloff. Larger → blur crosses more depth bands.
    depth_sigma:   f32,
    _pad:          f32,
};

@group(0) @binding(0) var<uniform> params: BlurParams;
@group(0) @binding(1) var depth_tex: texture_depth_2d;
@group(0) @binding(2) var src_ssao:  texture_2d<f32>;
@group(0) @binding(3) var out_ssao:  texture_storage_2d<r32float, write>;

@compute @workgroup_size(8, 8, 1)
fn cs_blur(@builtin(global_invocation_id) gid: vec3<u32>) {
    let dims = textureDimensions(out_ssao);
    if (gid.x >= dims.x || gid.y >= dims.y) {
        return;
    }
    let px = vec2<i32>(gid.xy);
    let center_depth = textureLoad(depth_tex, px, 0);
    if (center_depth >= 1.0) {
        // Background — pass through.
        let v = textureLoad(src_ssao, px, 0).r;
        textureStore(out_ssao, px, vec4<f32>(v, 0.0, 0.0, 1.0));
        return;
    }
    let center_ao = textureLoad(src_ssao, px, 0).r;

    // Spatial gaussian weights for offsets {-2, -1, 0, 1, 2} pixels
    // (σ ≈ 1.4 → roughly 0.054 / 0.244 / 0.404 / 0.244 / 0.054).
    let weights = array<f32, 5>(0.054, 0.244, 0.404, 0.244, 0.054);
    let offsets = array<i32, 5>(-2, -1, 0, 1, 2);

    var sum: f32 = 0.0;
    var w_total: f32 = 0.0;
    for (var i: i32 = 0; i < 5; i = i + 1) {
        let off = offsets[i];
        let tap = px + params.dir * off;
        let tx = clamp(tap.x, 0, i32(dims.x) - 1);
        let ty = clamp(tap.y, 0, i32(dims.y) - 1);
        let tap_clamped = vec2<i32>(tx, ty);
        let tap_depth = textureLoad(depth_tex, tap_clamped, 0);
        if (tap_depth >= 1.0) {
            continue;
        }
        let tap_ao = textureLoad(src_ssao, tap_clamped, 0).r;
        let depth_diff = abs(tap_depth - center_depth);
        let depth_weight = exp(-depth_diff * depth_diff / (params.depth_sigma * params.depth_sigma + 1e-6));
        let w = weights[i] * depth_weight;
        sum = sum + tap_ao * w;
        w_total = w_total + w;
    }
    let result = select(center_ao, sum / w_total, w_total > 1e-6);
    textureStore(out_ssao, px, vec4<f32>(clamp(result, 0.0, 1.0), 0.0, 0.0, 1.0));
}
