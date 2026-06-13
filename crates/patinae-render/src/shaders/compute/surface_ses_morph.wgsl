// Connolly SES probe-erosion morph.
//
// Reads the analytical vdW SDF written by `surface_vdw_sdf.wgsl` and produces
// the SES scalar field consumed by `surface_mc.wgsl` (with `invert_inside=1`,
// `iso=0`).
//
// ## Connolly SES, in one paragraph
//
// The solvent-accessible volume is `SAV = { x : dist(x, A) <= probe }` where
// `A = { p : f_vdw(p) >= probe }` is the set of allowed probe positions. The
// molecular-surface (Connolly SES) region is `C = R^3 \ SAV` — the set of
// points the probe surface can never reach. Algebraically:
//
//     x in C  <=>  forall y in ball(x, probe): f_vdw(y) < probe
//             <=>  max_{y in ball(x, probe)} f_vdw(y) < probe
//
// So the field
//
//     g_ses(x) = max_{y in ball(x, probe)} f_vdw(y) - probe
//
// is negative inside the molecular volume, zero on the SES surface, and
// positive in solvent-accessible space. Marching cubes consumes it with
// `inside <=> value <= iso=0` (the SES branch of `corner_inside`).
//
// ## Implementation
//
// One thread per voxel, brute-force max over a ball stencil of radius
// `stencil_radius_voxels`. The stencil extent is computed CPU-side from
// `ceil(probe_radius / voxel_size)`. For typical probe = 1.4 Å and voxel
// ~0.5 Å this is a 7^3 = 343 sample loop per voxel — not free, but the
// dispatch is dwarfed by the MC pass that follows. A separable
// approximation isn't safe here (the ball is not the box), and would
// over-erode in concave pockets; brute force keeps the geometry honest.
//
// ## Bind group layout (group 0)
//
// | Binding | Resource                                    | Access  |
// |--------:|---------------------------------------------|---------|
// |       0 | `MorphParams` uniform                       | uniform |
// |       1 | `sdf_in` r32float sampled texture (load)    | sampled |
// |       2 | `field_out` rgba16float storage texture     | write   |
//
// The vdW SDF is bound as a regular `texture_3d<f32>` (read via `textureLoad`,
// no sampler). `texture_storage_3d<r32float, read>` is not baseline-supported
// on WebGPU 1.0 — read-only storage on r32 formats requires an opt-in
// adapter feature that we don't take. The plain sampled-texture path is
// universally available and just as fast for our access pattern (no
// filtering, point loads).

struct MorphParams {
    dims: vec3<u32>,
    voxel_size: f32,
    probe_radius: f32,
    stencil_radius_voxels: u32,
    _pad0: u32,
    _pad1: u32,
};

@group(0) @binding(0) var<uniform> params: MorphParams;
@group(0) @binding(1) var sdf_in: texture_3d<f32>;
@group(0) @binding(2) var field_out: texture_storage_3d<rgba16float, write>;

fn sample_clamped(i: i32, j: i32, k: i32) -> f32 {
    let dx = i32(params.dims.x);
    let dy = i32(params.dims.y);
    let dz = i32(params.dims.z);
    let ci = clamp(i, 0, dx - 1);
    let cj = clamp(j, 0, dy - 1);
    let ck = clamp(k, 0, dz - 1);
    return textureLoad(sdf_in, vec3<i32>(ci, cj, ck), 0).r;
}

@compute @workgroup_size(4, 4, 4)
fn cs_main(@builtin(global_invocation_id) gid: vec3<u32>) {
    if (gid.x >= params.dims.x || gid.y >= params.dims.y || gid.z >= params.dims.z) {
        return;
    }
    let bx = i32(gid.x);
    let by = i32(gid.y);
    let bz = i32(gid.z);
    let r = i32(params.stencil_radius_voxels);
    // Compare squared offsets in voxel units against the squared stencil
    // radius. Inclusive check (<=) keeps the boundary voxels in the ball.
    let r2 = f32(r * r);

    var max_sdf: f32 = -1.0e30;
    for (var dk: i32 = -r; dk <= r; dk = dk + 1) {
        for (var dj: i32 = -r; dj <= r; dj = dj + 1) {
            for (var di: i32 = -r; di <= r; di = di + 1) {
                let d2 = f32(di * di + dj * dj + dk * dk);
                if (d2 > r2) {
                    continue;
                }
                let v = sample_clamped(bx + di, by + dj, bz + dk);
                max_sdf = max(max_sdf, v);
            }
        }
    }

    let g_ses = max_sdf - params.probe_radius;
    // .xyz is unused — MC computes its own gradient from the .w channel via
    // central differences.
    textureStore(field_out, vec3<i32>(bx, by, bz), vec4<f32>(0.0, 0.0, 0.0, g_ses));
}
