// Per-rep frustum cull — sphere variant.
//
// See `cull_common.wgsl` for shared structs and helpers. Per-kind
// shaders add only the instance struct + bindings 1/3 + entry point.

// {{INCLUDE_CULL_COMMON}}

struct SphereInstance {
    center:   vec3<f32>,
    radius:   f32,
    group_id: u32,
    _pad:     u32,
};

@group(0) @binding(1) var<storage, read>       raw_inst: array<SphereInstance>;
@group(0) @binding(3) var<storage, read_write> out_inst: array<SphereInstance>;

@compute @workgroup_size(64)
fn cs_cull_sphere(
    @builtin(global_invocation_id) gid: vec3<u32>,
    @builtin(num_workgroups)       nwg: vec3<u32>,
) {
    let i = gid.x + gid.y * nwg.x * 64u;
    if (i >= raw_count[0]) {
        return;
    }
    let inst = raw_inst[i];
    if (should_cull(inst.center, inst.radius)) {
        return;
    }
    let slot = atomicAdd(&indirect[1], 1u);
    if (slot < params.raw_capacity) {
        out_inst[slot] = inst;
    }
}
