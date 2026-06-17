struct PrimitiveMetadata {
    sphere_count: u32,
    cylinder_count: u32,
    capsule_count: u32,
    triangle_count: u32,
    primitive_count: u32,
    triangle_capacity: u32,
    visible_triangle_count: u32,
    overflow: u32,
}

struct FinalizeParams {
    direct_triangle_count: u32,
    visible_triangle_capacity: u32,
    counter_index: u32,
    bvh_leaf_size: u32,
    max_workgroups_per_dimension: u32,
    _pad0: u32,
    _pad1: u32,
    _pad2: u32,
}

@group(0) @binding(0) var<storage, read_write> metadata: PrimitiveMetadata;
@group(0) @binding(1) var<storage, read_write> visible_counts: array<atomic<u32>>;
@group(0) @binding(2) var<uniform> params: FinalizeParams;
@group(0) @binding(3) var<storage, read_write> leaf_dispatch_args: array<u32>;

fn workgroups_for_items(item_count: u32) -> vec3<u32> {
    let leaf_size = max(params.bvh_leaf_size, 1u);
    let leaf_count = max(1u, (item_count + leaf_size - 1u) / leaf_size);
    let total_workgroups = max(1u, (leaf_count + 127u) / 128u);
    let max_workgroups = max(params.max_workgroups_per_dimension, 1u);
    let workgroups_x = min(total_workgroups, max_workgroups);
    let workgroups_y = (total_workgroups + workgroups_x - 1u) / workgroups_x;
    return vec3<u32>(workgroups_x, workgroups_y, 1u);
}

@compute @workgroup_size(1)
fn finalize_streaming_metadata() {
    let visible_count = atomicLoad(&visible_counts[params.counter_index]);
    let stored_visible_count = min(visible_count, params.visible_triangle_capacity);
    let triangle_count = params.direct_triangle_count + stored_visible_count;

    metadata.triangle_count = triangle_count;
    metadata.visible_triangle_count = stored_visible_count;
    metadata.primitive_count =
        metadata.sphere_count + metadata.cylinder_count + metadata.capsule_count + triangle_count;
    metadata.overflow = select(0u, 1u, visible_count > params.visible_triangle_capacity);

    let dispatch_args = workgroups_for_items(metadata.primitive_count);
    leaf_dispatch_args[0] = dispatch_args.x;
    leaf_dispatch_args[1] = dispatch_args.y;
    leaf_dispatch_args[2] = dispatch_args.z;
}
