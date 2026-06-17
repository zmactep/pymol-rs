struct BvhParams {
    leaf_slots: u32,
    leaf_start: u32,
    level_start: u32,
    level_count: u32,
    dispatch_width: u32,
    _pad0: u32,
    _pad1: u32,
}

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

struct Sphere {
    center: vec3<f32>,
    radius: f32,
    color: vec4<f32>,
    transparency: f32,
    _pad_a: f32, _pad_b: f32, _pad_c: f32,
}

struct Cylinder {
    start: vec3<f32>,
    radius: f32,
    end: vec3<f32>,
    _pad0: f32,
    color1: vec4<f32>,
    color2: vec4<f32>,
    transparency: f32,
    _pad_a: f32, _pad_b: f32, _pad_c: f32,
}

struct Capsule {
    start: vec3<f32>,
    radius: f32,
    end: vec3<f32>,
    _pad0: f32,
    color1: vec4<f32>,
    color2: vec4<f32>,
    transparency: f32,
    _pad_a: f32, _pad_b: f32, _pad_c: f32,
}

// {{INCLUDE_ARTIFACT_TRIANGLE}}

struct BvhNode {
    min: vec3<f32>,
    left_or_first: u32,
    max: vec3<f32>,
    count: u32,
}

struct PrimitiveBounds {
    encoded: u32,
    valid: bool,
    min: vec3<f32>,
    max: vec3<f32>,
}

@group(0) @binding(0) var<storage, read> spheres: array<Sphere>;
@group(0) @binding(1) var<storage, read> cylinders: array<Cylinder>;
@group(0) @binding(2) var<storage, read> capsules: array<Capsule>;
@group(0) @binding(3) var<storage, read> triangles: array<Triangle>;
@group(0) @binding(4) var<storage, read_write> bvh_nodes: array<BvhNode>;
@group(0) @binding(5) var<storage, read_write> bvh_indices: array<u32>;
@group(0) @binding(6) var<uniform> params: BvhParams;
@group(0) @binding(7) var<storage, read> metadata: PrimitiveMetadata;

const LEAF_SIZE: u32 = 4u;
const INF: f32 = 1.0e30;
const EMPTY_NODE: u32 = 0xffffffffu;
const EPSILON: f32 = 0.0001;

fn empty_node() -> BvhNode {
    return BvhNode(vec3<f32>(0.0), EMPTY_NODE, vec3<f32>(0.0), 0u);
}

fn invalid_bounds() -> PrimitiveBounds {
    return PrimitiveBounds(0u, false, vec3<f32>(0.0), vec3<f32>(0.0));
}

fn is_empty_node(node: BvhNode) -> bool {
    return node.left_or_first == EMPTY_NODE && node.count == 0u;
}

fn sphere_bounds(index: u32) -> PrimitiveBounds {
    let sphere = spheres[index];
    if sphere.radius <= 0.0 || sphere.transparency >= 1.0 || sphere.color.a <= 0.0 {
        return invalid_bounds();
    }
    let r = vec3<f32>(sphere.radius);
    return PrimitiveBounds(index, true, sphere.center - r, sphere.center + r);
}

fn cylinder_bounds(index: u32) -> PrimitiveBounds {
    let cylinder = cylinders[index];
    if cylinder.radius <= 0.0 || cylinder.transparency >= 1.0 {
        return invalid_bounds();
    }
    if length(cylinder.end - cylinder.start) <= EPSILON {
        return invalid_bounds();
    }
    let r = vec3<f32>(cylinder.radius);
    return PrimitiveBounds(
        (1u << 30u) | index,
        true,
        min(cylinder.start, cylinder.end) - r,
        max(cylinder.start, cylinder.end) + r,
    );
}

fn capsule_bounds(index: u32) -> PrimitiveBounds {
    let capsule = capsules[index];
    if capsule.radius <= 0.0 || capsule.transparency >= 1.0 {
        return invalid_bounds();
    }
    if length(capsule.end - capsule.start) <= EPSILON {
        return invalid_bounds();
    }
    let r = vec3<f32>(capsule.radius);
    return PrimitiveBounds(
        (3u << 30u) | index,
        true,
        min(capsule.start, capsule.end) - r,
        max(capsule.start, capsule.end) + r,
    );
}

fn triangle_bounds(index: u32) -> PrimitiveBounds {
    let tri = triangles[index];
    if tri.color.a <= 0.0 {
        return invalid_bounds();
    }
    let v0 = tri.v0.xyz;
    let v1 = tri.v1.xyz;
    let v2 = tri.v2.xyz;
    let area_normal = cross(v1 - v0, v2 - v0);
    if length(area_normal) <= EPSILON {
        return invalid_bounds();
    }
    return PrimitiveBounds(
        (2u << 30u) | index,
        true,
        min(v0, min(v1, v2)),
        max(v0, max(v1, v2)),
    );
}

fn primitive_bounds(primitive_index: u32) -> PrimitiveBounds {
    if primitive_index < metadata.sphere_count {
        return sphere_bounds(primitive_index);
    }
    let cylinder_start = metadata.sphere_count;
    let capsule_start = metadata.sphere_count + metadata.cylinder_count;
    let triangle_start = capsule_start + metadata.capsule_count;
    if primitive_index < capsule_start {
        return cylinder_bounds(primitive_index - cylinder_start);
    }
    if primitive_index < triangle_start {
        return capsule_bounds(primitive_index - capsule_start);
    }
    if primitive_index < metadata.primitive_count {
        return triangle_bounds(primitive_index - triangle_start);
    }
    return invalid_bounds();
}

fn active_leaf_count() -> u32 {
    return max(1u, (metadata.primitive_count + LEAF_SIZE - 1u) / LEAF_SIZE);
}

fn child_node_or_empty(node_index: u32, child_leaf_start: u32) -> BvhNode {
    if child_leaf_start >= active_leaf_count() {
        return empty_node();
    }
    return bvh_nodes[node_index];
}

@compute @workgroup_size(128)
fn build_leaves(@builtin(global_invocation_id) gid: vec3<u32>) {
    let leaf_id = gid.y * params.dispatch_width + gid.x;
    if leaf_id >= params.leaf_slots {
        return;
    }

    let node_index = params.leaf_start + leaf_id;
    let first = leaf_id * LEAF_SIZE;
    if first >= metadata.primitive_count {
        bvh_nodes[node_index] = empty_node();
        return;
    }

    var valid_count = 0u;
    var bmin = vec3<f32>(INF);
    var bmax = vec3<f32>(-INF);
    for (var i = 0u; i < LEAF_SIZE; i = i + 1u) {
        let primitive_index = first + i;
        if primitive_index < metadata.primitive_count {
            let bounds = primitive_bounds(primitive_index);
            if bounds.valid {
                bvh_indices[first + valid_count] = bounds.encoded;
                bmin = min(bmin, bounds.min);
                bmax = max(bmax, bounds.max);
                valid_count = valid_count + 1u;
            }
        }
    }
    if valid_count == 0u {
        bvh_nodes[node_index] = empty_node();
        return;
    }
    bvh_nodes[node_index] = BvhNode(bmin, first, bmax, valid_count);
}

@compute @workgroup_size(128)
fn build_internal(@builtin(global_invocation_id) gid: vec3<u32>) {
    let local_id = gid.y * params.dispatch_width + gid.x;
    if local_id >= params.level_count {
        return;
    }
    let parent_start = params.level_start - params.level_count;
    let parent_index = parent_start + local_id;
    let left_index = params.level_start + local_id * 2u;
    let right_index = left_index + 1u;
    let child_level_count = params.level_count * 2u;
    let child_leaf_span = max(1u, params.leaf_slots / child_level_count);
    let left_leaf_start = local_id * 2u * child_leaf_span;
    let right_leaf_start = left_leaf_start + child_leaf_span;
    let left = child_node_or_empty(left_index, left_leaf_start);
    let right = child_node_or_empty(right_index, right_leaf_start);
    if is_empty_node(left) && is_empty_node(right) {
        bvh_nodes[parent_index] = empty_node();
        return;
    }
    if is_empty_node(left) {
        bvh_nodes[parent_index] = right;
        return;
    }
    if is_empty_node(right) {
        bvh_nodes[parent_index] = left;
        return;
    }
    let bmin = min(left.min, right.min);
    let bmax = max(left.max, right.max);
    bvh_nodes[parent_index] = BvhNode(bmin, left_index, bmax, 0u);
}
