// Weighted Blended Order-Independent Transparency (McGuire & Bavoil 2013).
//
// Every fragment that the geometry pipeline emits passes through
// `write_translucent`. α=1 fragments are a degenerate case of the same
// formula; the opaque fast-path is a separate pipeline choice.
//
// Required render targets when this helper is used:
//   @location(0) accum  : Rgba16Float, blend = (One, One)             — additive
//   @location(1) reveal : R16Float,    blend = (Zero, OneMinusSrcColor) — multiplicative
//
// Compose pass divides accum.rgb by accum.a and blends with reveal.

struct TranslucentOut {
    @location(0) accum:  vec4<f32>,
    @location(1) reveal: f32,
};

// `z_view` is positive view-space depth (i.e. -view.z for a right-handed
// camera looking down -Z). `scene_max_depth` (frame.clip.w) normalises
// the depth bias so the weight is scale-invariant across scenes.
fn wboit_weight(z_view: f32, alpha: f32, scene_max_depth: f32) -> f32 {
    let d = z_view / max(scene_max_depth, 1.0);
    let denom = 1e-5 + pow(d, 4.0);
    return alpha * clamp(0.03 / denom, 1e-2, 3e3);
}

fn write_translucent(color: vec3<f32>, alpha: f32, z_view: f32, scene_max_depth: f32) -> TranslucentOut {
    let w = wboit_weight(z_view, alpha, scene_max_depth);
    var out: TranslucentOut;
    out.accum  = vec4<f32>(color * alpha * w, alpha * w);
    out.reveal = alpha;
    return out;
}
