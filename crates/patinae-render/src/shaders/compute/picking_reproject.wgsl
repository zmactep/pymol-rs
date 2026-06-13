// Picking-texture reprojection — warps the previous frame's picking
// pixels into the current frame's view without re-rasterizing geometry.
//
// Two passes:
//   1) `cs_min_depth` — every source pixel atomicMin's its post-projection
//                       depth into `best_depth[target_idx]`, claiming the
//                       nearest source for that target.
//   2) `cs_write`     — every source pixel re-reads `best_depth`; if its
//                       depth matches, it writes its packed picking value
//                       to `picking_out` via `textureStore`. This races
//                       only on ties, which produce identical outputs.
//
// `picking_in` / `depth_in` come from the previous frame; `picking_out` is
// a storage view of the *current* frame's picking texture. Pixels whose
// reprojected ndc lies outside [-1, 1]³ are dropped silently — they would
// have been "no atom" anyway. `view_proj_prev_inv` is the inverse of the
// view-proj matrix used to render the source picking texture;
// `view_proj_cur` is this frame's matrix. The host fills both via
// `ReprojectParams`.

struct ReprojectParams {
    view_proj_prev_inv: mat4x4<f32>,
    view_proj_cur:      mat4x4<f32>,
};

@group(0) @binding(0) var<uniform>      params:       ReprojectParams;
@group(0) @binding(1) var               picking_in:   texture_2d<u32>;
@group(0) @binding(2) var               depth_in:     texture_2d<f32>;
@group(0) @binding(3) var<storage, read_write> best_depth: array<atomic<u32>>;
@group(0) @binding(4) var               picking_out:  texture_storage_2d<rg32uint, write>;

const INF_DEPTH: u32 = 0xFFFFFFFFu;

// Reprojects pixel `(x_src, y_src)` from the previous frame into the
// current frame. Returns `(target_xy, depth_u32_for_atomic, packed_id)`
// or `None` (via a flag) if the result lies outside the current frustum
// or the source pixel was the "no atom" sentinel.
struct ReprojectHit {
    valid:       bool,
    target_xy:   vec2<u32>,
    depth_u32:   u32,
    packed:      vec2<u32>,
};

fn reproject(src_xy: vec2<u32>, dims: vec2<u32>) -> ReprojectHit {
    var out: ReprojectHit;
    out.valid = false;

    let coord = vec2<i32>(src_xy);
    let packed = textureLoad(picking_in, coord, 0).xy;
    // Sentinel — empty pixel in the previous frame.
    if packed.x == 0u && packed.y == 0u {
        return out;
    }

    let depth_prev = textureLoad(depth_in, coord, 0).r;
    // Unproject (xy_ndc, depth_ndc) through inv(view_proj_prev).
    let uv = (vec2<f32>(src_xy) + vec2<f32>(0.5)) / vec2<f32>(dims);
    let ndc_prev = vec4<f32>(
        uv.x * 2.0 - 1.0,
        1.0 - uv.y * 2.0,  // flip Y so picking texel (0,0) is top-left
        depth_prev,
        1.0,
    );
    let world_h = params.view_proj_prev_inv * ndc_prev;
    if abs(world_h.w) < 1e-6 {
        return out;
    }
    let world = world_h.xyz / world_h.w;

    // Project through view_proj_cur.
    let ndc_cur_h = params.view_proj_cur * vec4<f32>(world, 1.0);
    if ndc_cur_h.w <= 0.0 {
        return out;
    }
    let ndc_cur = ndc_cur_h.xyz / ndc_cur_h.w;
    if ndc_cur.x < -1.0 || ndc_cur.x > 1.0
        || ndc_cur.y < -1.0 || ndc_cur.y > 1.0
        || ndc_cur.z < 0.0 || ndc_cur.z > 1.0
    {
        return out;
    }

    let uv_cur = vec2<f32>(ndc_cur.x * 0.5 + 0.5, 0.5 - ndc_cur.y * 0.5);
    let target_f = uv_cur * vec2<f32>(dims);
    let target_xy = vec2<u32>(
        u32(clamp(target_f.x, 0.0, f32(dims.x) - 1.0)),
        u32(clamp(target_f.y, 0.0, f32(dims.y) - 1.0)),
    );
    // `depth_u32` is the ndc-z mapped to u32, monotonically — perfect for
    // atomicMin (smaller raw bits ⇒ closer to camera in [0, 1] range).
    let depth_clamped = clamp(ndc_cur.z, 0.0, 1.0);
    let depth_u32 = u32(depth_clamped * f32(0xFFFFFF00u));

    out.valid = true;
    out.target_xy = target_xy;
    out.depth_u32 = depth_u32;
    out.packed = packed;
    return out;
}

@compute @workgroup_size(8, 8, 1)
fn cs_min_depth(@builtin(global_invocation_id) gid: vec3<u32>) {
    let dims = textureDimensions(picking_in);
    if gid.x >= dims.x || gid.y >= dims.y {
        return;
    }
    let hit = reproject(gid.xy, dims);
    if !hit.valid {
        return;
    }
    let idx = hit.target_xy.y * dims.x + hit.target_xy.x;
    atomicMin(&best_depth[idx], hit.depth_u32);
}

@compute @workgroup_size(8, 8, 1)
fn cs_write(@builtin(global_invocation_id) gid: vec3<u32>) {
    let dims = textureDimensions(picking_in);
    if gid.x >= dims.x || gid.y >= dims.y {
        return;
    }
    let hit = reproject(gid.xy, dims);
    if !hit.valid {
        return;
    }
    let idx = hit.target_xy.y * dims.x + hit.target_xy.x;
    // Only the thread whose depth matches `best_depth` writes. Ties write
    // identical values from the same source pixel layout — harmless.
    let best = atomicLoad(&best_depth[idx]);
    if hit.depth_u32 == best {
        textureStore(
            picking_out,
            vec2<i32>(hit.target_xy),
            vec4<u32>(hit.packed.x, hit.packed.y, 0u, 0u),
        );
    }
}
