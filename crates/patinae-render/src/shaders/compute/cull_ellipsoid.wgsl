// Per-rep frustum cull — ellipsoid variant.
//
// EllipsoidInstance carries center + three axis vectors. `max_extent`
// is pre-baked on host as the longest axis length (= worst-case
// bounding sphere radius).

// {{INCLUDE_CULL_COMMON}}

struct EllipsoidInstance {
    center:     vec3<f32>,
    max_extent: f32,
    axis0:      vec3<f32>,
    group_id:   u32,
    axis1:      vec3<f32>,
    _pad0:      u32,
    axis2:      vec3<f32>,
    _pad1:      u32,
};

@group(0) @binding(1) var<storage, read>       raw_inst: array<EllipsoidInstance>;
@group(0) @binding(3) var<storage, read_write> out_inst: array<EllipsoidInstance>;

@compute @workgroup_size(64)
fn cs_cull_ellipsoid(
    @builtin(global_invocation_id) gid: vec3<u32>,
    @builtin(num_workgroups)       nwg: vec3<u32>,
) {
    let i = gid.x + gid.y * nwg.x * 64u;
    if (i >= raw_count[0]) {
        return;
    }
    let inst = raw_inst[i];
    if (should_cull(inst.center, inst.max_extent)) {
        return;
    }
    let slot = atomicAdd(&indirect[1], 1u);
    if (slot < params.raw_capacity) {
        out_inst[slot] = inst;
    }
}
