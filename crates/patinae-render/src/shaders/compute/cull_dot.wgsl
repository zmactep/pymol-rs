// Per-rep frustum cull — dot variant.
//
// DotAtomInstance carries the atom center + vdW radius. Frustum culling runs
// once per atom, not once per generated surface sample.

// {{INCLUDE_CULL_COMMON}}

struct DotAtomInstance {
    center: vec3<f32>,
    group_id: u32,
    vdw_radius: f32,
    _pad0: u32,
    _pad1: u32,
    _pad2: u32,
};

@group(0) @binding(1) var<storage, read>       raw_inst: array<DotAtomInstance>;
@group(0) @binding(3) var<storage, read_write> out_inst: array<DotAtomInstance>;

@compute @workgroup_size(64)
fn cs_cull_dot(
    @builtin(global_invocation_id) gid: vec3<u32>,
    @builtin(num_workgroups)       nwg: vec3<u32>,
) {
    let i = gid.x + gid.y * nwg.x * 64u;
    if (i >= raw_count[0] || i >= params.raw_capacity) {
        return;
    }
    let inst = raw_inst[i];
    if (should_cull(inst.center, inst.vdw_radius + params.kind_radius)) {
        return;
    }
    let slot = atomicAdd(&indirect[1], 1u);
    if (slot < params.raw_capacity) {
        out_inst[slot] = inst;
    }
}
