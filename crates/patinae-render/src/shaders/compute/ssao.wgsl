// Crytek-style screen-space ambient occlusion (depth-only).
//
// Per pixel:
//   1. Sample center depth from `frame.depth`. Background (depth=1.0)
//      → AO=1.0 (no darkening).
//   2. Reconstruct view-space position P via `frame.proj_inv`.
//   3. Reconstruct view-space normal N from screen-space derivatives of
//      depth (sample +x and +y neighbours, reconstruct their view-space
//      positions, take cross product). This is the "depth gradient
//      normal" — accurate on interior surfaces, slightly noisy at
//      silhouettes (acceptable for AO).
//   4. Generate a per-pixel rotation from a hash of pixel coord.
//   5. For each of N_SAMPLES hemisphere samples (passed in via uniform):
//      transform sample to view space (TBN basis with N), project back
//      to screen UV, read neighbour depth. If neighbour is in front of
//      the sample AND within `params.radius`, count as occluded.
//   6. Output AO = 1.0 - occlusion_count / N_SAMPLES, biased to fit
//      `R8Unorm` storage.

// {{INCLUDE_FRAME}}

const N_SAMPLES: u32 = 32u;
const PI: f32 = 3.14159265359;

struct SsaoParams {
    /// Sphere radius in view-space units (Å).
    radius:    f32,
    /// Bias applied to the sample depth before comparing — suppresses
    /// self-occlusion on near-flat surfaces.
    bias:      f32,
    /// Time-modulated phase so consecutive frames jitter their sample
    /// pattern (helps temporal blur).
    frame_phase: f32,
    _pad:      f32,
    /// Precomputed hemisphere samples. Z-axis is the surface normal
    /// (TBN-aligned). Length-biased toward the surface center so
    /// shorter samples dominate the integral. CPU side rebuilds these
    /// once on `SsaoCompute::new`.
    samples:   array<vec4<f32>, 32>,
};

// `frame: FrameUniforms` is declared at @group(0) @binding(0) by
// INCLUDE_FRAME above.
@group(0) @binding(1) var<uniform> params: SsaoParams;
@group(0) @binding(2) var depth_tex: texture_depth_2d;
@group(0) @binding(3) var out_ssao:  texture_storage_2d<r32float, write>;

/// Reconstruct a view-space position from screen-space UV (in [0,1])
/// and clip-space z. `frame.proj_inv * vec4(ndc_xy, z, 1)`.
fn view_pos_from_uv(uv: vec2<f32>, depth: f32) -> vec3<f32> {
    let ndc = vec3<f32>(uv.x * 2.0 - 1.0, 1.0 - uv.y * 2.0, depth);
    let h = frame.proj_inv * vec4<f32>(ndc, 1.0);
    return h.xyz / h.w;
}

/// Cheap hash for per-pixel random rotation. `vec2 → vec3` so we get a
/// rotation axis + angle. Output components in `[0, 1)`.
fn hash3(seed: vec2<f32>) -> vec3<f32> {
    let p = vec3<f32>(
        dot(seed, vec2<f32>(127.1, 311.7)),
        dot(seed, vec2<f32>(269.5, 183.3)),
        dot(seed, vec2<f32>(419.2, 371.9)),
    );
    return fract(sin(p) * 43758.5453);
}

@compute @workgroup_size(8, 8, 1)
fn cs_ssao(@builtin(global_invocation_id) gid: vec3<u32>) {
    let dims = textureDimensions(out_ssao);
    if (gid.x >= dims.x || gid.y >= dims.y) {
        return;
    }
    let px = vec2<i32>(gid.xy);
    let depth = textureLoad(depth_tex, px, 0);
    if (depth >= 1.0) {
        // Background — no occlusion.
        textureStore(out_ssao, px, vec4<f32>(1.0, 0.0, 0.0, 1.0));
        return;
    }
    let dim_f = vec2<f32>(f32(dims.x), f32(dims.y));
    let uv = (vec2<f32>(gid.xy) + vec2<f32>(0.5)) / dim_f;
    let p = view_pos_from_uv(uv, depth);

    // View-space normal from depth gradient. Use +1px in screen-X and
    // +1px in screen-Y, reconstruct their positions, take cross.
    // Falls back to the "near" cross at silhouettes (cross has length
    // < 1e-6) so the AO doesn't blow up.
    let dx_uv = uv + vec2<f32>(1.0 / dim_f.x, 0.0);
    let dy_uv = uv + vec2<f32>(0.0, 1.0 / dim_f.y);
    let dx_px = vec2<i32>(min(px.x + 1, i32(dims.x) - 1), px.y);
    let dy_px = vec2<i32>(px.x, min(px.y + 1, i32(dims.y) - 1));
    let dx_depth = textureLoad(depth_tex, dx_px, 0);
    let dy_depth = textureLoad(depth_tex, dy_px, 0);
    let p_dx = view_pos_from_uv(dx_uv, dx_depth);
    let p_dy = view_pos_from_uv(dy_uv, dy_depth);
    var normal = normalize(cross(p_dx - p, p_dy - p));
    // Flip if normal faces away from the camera (view-space camera is
    // at origin, looking down -Z; valid normals point toward +Z).
    if (normal.z < 0.0) {
        normal = -normal;
    }

    // Per-pixel rotation: rebuild an arbitrary tangent perpendicular to
    // the normal, rotated by a hashed angle. Wrap the hemisphere with a
    // TBN basis (N=normal, T=tangent, B=N×T).
    let h = hash3(vec2<f32>(gid.xy) + vec2<f32>(params.frame_phase));
    let raw_t = vec3<f32>(h.x - 0.5, h.y - 0.5, 0.0);
    var tangent = raw_t - normal * dot(normal, raw_t);
    let tlen = length(tangent);
    if (tlen < 1e-3) {
        // Pathological alignment — pick world-X fallback.
        tangent = vec3<f32>(1.0, 0.0, 0.0);
    } else {
        tangent = tangent / tlen;
    }
    let bitangent = cross(normal, tangent);

    var occlusion: f32 = 0.0;
    for (var i: u32 = 0u; i < N_SAMPLES; i = i + 1u) {
        let s = params.samples[i].xyz;
        // Sample in TBN basis, scaled by user radius.
        let view_offset = (tangent * s.x + bitangent * s.y + normal * s.z) * params.radius;
        let sample_view = p + view_offset;
        // Project back to clip → NDC → UV.
        let sample_clip = frame.proj * vec4<f32>(sample_view, 1.0);
        if (sample_clip.w <= 0.0) {
            continue;
        }
        let sample_ndc = sample_clip.xyz / sample_clip.w;
        let sample_uv = vec2<f32>(sample_ndc.x * 0.5 + 0.5, 1.0 - (sample_ndc.y * 0.5 + 0.5));
        if (sample_uv.x < 0.0 || sample_uv.x > 1.0 ||
            sample_uv.y < 0.0 || sample_uv.y > 1.0) {
            continue;
        }
        let sample_px = vec2<i32>(sample_uv * dim_f);
        let scene_depth = textureLoad(depth_tex, sample_px, 0);
        if (scene_depth >= 1.0) {
            continue;
        }
        // Reconstruct the actual scene point at that UV; compare its
        // view-space z to the sample's. Distance attenuation prevents
        // far occluders from contributing.
        let scene_view = view_pos_from_uv(sample_uv, scene_depth);
        // View space: camera at origin looking down -Z. `z` is negative
        // for points in front of camera; less-negative = closer.
        // Sample is occluded iff the scene at its UV is between sample
        // and camera, i.e. `scene_view.z > sample_view.z` (less
        // negative). `dz = sample.z - scene.z < -bias` enforces that
        // with the bias dead-zone.
        let dz = sample_view.z - scene_view.z;
        let range = smoothstep(0.0, 1.0, params.radius / abs(p.z - scene_view.z + 1e-4));
        if (dz < -params.bias) {
            occlusion = occlusion + range;
        }
    }
    let ao = 1.0 - occlusion / f32(N_SAMPLES);
    textureStore(out_ssao, px, vec4<f32>(clamp(ao, 0.0, 1.0), 0.0, 0.0, 1.0));
}
