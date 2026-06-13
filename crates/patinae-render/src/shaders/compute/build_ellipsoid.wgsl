// Ellipsoid instance builder. One thread per `EllipsoidBuildEntry`
// (sparse build-input pattern — only atoms whose ellipsoid axes were
// pre-computed CPU-side via Jacobi eigendecomp + isotropic fallback).
//
// Reads the shared atom centre from `SceneStore.coords`; emits one
// `EllipsoidInstance` per entry into a compacted buffer with the
// instance count atomically incremented inside the indirect draw args.
//
// Group 0 — `SceneStore` (compute side; mirrors render's `INCLUDE_SCENE`
//           layout but at @group(0)).
// Group 1 — per-rep build state.

struct ObjectEntry {
    atom_offset: u32,
    atom_count:  u32,
    bond_offset: u32,
    bond_count:  u32,
    object_id:   u32,
    flags:       u32,
    _pad0_a:     u32,
    _pad0_b:     u32,
    model_matrix: mat4x4<f32>,
};

@group(0) @binding(0) var<uniform>           obj    : ObjectEntry;
@group(0) @binding(2) var<storage, read>     coords : array<vec4<f32>>;

struct EllipsoidBuildEntry {
    axis0:      vec3<f32>,
    atom_local: u32,
    axis1:      vec3<f32>,
    _pad0:      u32,
    axis2:      vec3<f32>,
    _pad1:      u32,
};

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

struct EllipsoidBuildParams {
    entry_count: u32,
    _pad0:       u32,
    _pad1:       u32,
    _pad2:       u32,
};

@group(1) @binding(0) var<uniform>             params:      EllipsoidBuildParams;
@group(1) @binding(1) var<storage, read>       build_input: array<EllipsoidBuildEntry>;
@group(1) @binding(2) var<storage, read_write> instances:   array<EllipsoidInstance>;
@group(1) @binding(3) var<storage, read_write> raw_count:   array<atomic<u32>, 1>;

@compute @workgroup_size(64)
fn cs_main(
    @builtin(global_invocation_id) gid_in: vec3<u32>,
    @builtin(num_workgroups) nwg: vec3<u32>,
) {
    // 2D dispatch handles entry_count > 65535*64.
    let i = gid_in.x + gid_in.y * nwg.x * 64u;
    if i >= params.entry_count { return; }

    let entry = build_input[i];
    let center = coords[obj.atom_offset + entry.atom_local].xyz;

    let l0 = length(entry.axis0);
    let l1 = length(entry.axis1);
    let l2 = length(entry.axis2);
    let max_extent = max(l0, max(l1, l2));
    if max_extent < 1e-4 { return; }

    let dst = atomicAdd(&raw_count[0], 1u);
    var ins: EllipsoidInstance;
    ins.center     = center;
    ins.max_extent = max_extent;
    ins.axis0      = entry.axis0;
    ins.group_id   = entry.atom_local;
    ins.axis1      = entry.axis1;
    ins._pad0      = 0u;
    ins.axis2      = entry.axis2;
    ins._pad1      = 0u;
    instances[dst] = ins;
}
