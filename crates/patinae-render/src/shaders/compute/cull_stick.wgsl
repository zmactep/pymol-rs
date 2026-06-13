// Per-rep frustum cull — stick variant.
//
// StickInstance carries two endpoints + a bond radius in `p0_radius.w`.
// Bounding sphere = midpoint, half-length + radius.

// {{INCLUDE_CULL_COMMON}}

struct StickInstance {
    p0_radius: vec4<f32>,   // xyz = p0, w = bond radius
    p1_pad:    vec4<f32>,   // xyz = p1
    groups:    vec2<u32>,
    _pad:      vec2<u32>,
};

@group(0) @binding(1) var<storage, read>       raw_inst: array<StickInstance>;
@group(0) @binding(3) var<storage, read_write> out_inst: array<StickInstance>;

@compute @workgroup_size(64)
fn cs_cull_stick(
    @builtin(global_invocation_id) gid: vec3<u32>,
    @builtin(num_workgroups)       nwg: vec3<u32>,
) {
    let i = gid.x + gid.y * nwg.x * 64u;
    if (i >= raw_count[0]) {
        return;
    }
    let inst = raw_inst[i];
    let p0 = inst.p0_radius.xyz;
    let p1 = inst.p1_pad.xyz;
    let center = 0.5 * (p0 + p1);
    let half_len = 0.5 * length(p1 - p0);
    let radius = half_len + inst.p0_radius.w;
    if (should_cull(center, radius)) {
        return;
    }
    let slot = atomicAdd(&indirect[1], 1u);
    if (slot < params.raw_capacity) {
        out_inst[slot] = inst;
    }
}
