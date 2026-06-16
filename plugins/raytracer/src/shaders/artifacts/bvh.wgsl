struct BvhParams {
    primitive_count: u32,
    sphere_count: u32,
    cylinder_count: u32,
    triangle_count: u32,
    leaf_slots: u32,
    leaf_start: u32,
    level_start: u32,
    level_count: u32,
    dispatch_width: u32,
    _pad0: u32,
    _pad1: u32,
    _pad2: u32,
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

struct Triangle {
    v0: vec3<f32>, _pad0: f32,
    v1: vec3<f32>, _pad1: f32,
    v2: vec3<f32>, _pad2: f32,
    n0: vec3<f32>, _pad3: f32,
    n1: vec3<f32>, _pad4: f32,
    n2: vec3<f32>, _pad5: f32,
    color: vec4<f32>,
    transparency: f32,
    _pad_a: f32, _pad_b: f32, _pad_c: f32,
}

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
@group(0) @binding(2) var<storage, read> triangles: array<Triangle>;
@group(0) @binding(3) var<storage, read_write> bvh_nodes: array<BvhNode>;
@group(0) @binding(4) var<storage, read_write> bvh_indices: array<u32>;
@group(0) @binding(5) var<uniform> params: BvhParams;

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

fn triangle_bounds(index: u32) -> PrimitiveBounds {
    let tri = triangles[index];
    if tri.transparency >= 1.0 || tri.color.a <= 0.0 {
        return invalid_bounds();
    }
    let area_normal = cross(tri.v1 - tri.v0, tri.v2 - tri.v0);
    if length(area_normal) <= EPSILON {
        return invalid_bounds();
    }
    return PrimitiveBounds(
        (2u << 30u) | index,
        true,
        min(tri.v0, min(tri.v1, tri.v2)),
        max(tri.v0, max(tri.v1, tri.v2)),
    );
}

fn primitive_bounds(primitive_index: u32) -> PrimitiveBounds {
    if primitive_index < params.sphere_count {
        return sphere_bounds(primitive_index);
    }
    let cylinder_start = params.sphere_count;
    let triangle_start = params.sphere_count + params.cylinder_count;
    if primitive_index < triangle_start {
        return cylinder_bounds(primitive_index - cylinder_start);
    }
    if primitive_index < params.primitive_count {
        return triangle_bounds(primitive_index - triangle_start);
    }
    return invalid_bounds();
}

@compute @workgroup_size(128)
fn build_leaves(@builtin(global_invocation_id) gid: vec3<u32>) {
    let leaf_id = gid.y * params.dispatch_width + gid.x;
    if leaf_id >= params.leaf_slots {
        return;
    }

    let node_index = params.leaf_start + leaf_id;
    let first = leaf_id * LEAF_SIZE;
    if first >= params.primitive_count {
        bvh_nodes[node_index] = empty_node();
        return;
    }

    var valid_count = 0u;
    var bmin = vec3<f32>(INF);
    var bmax = vec3<f32>(-INF);
    for (var i = 0u; i < LEAF_SIZE; i = i + 1u) {
        let primitive_index = first + i;
        if primitive_index < params.primitive_count {
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
    let left = bvh_nodes[left_index];
    let right = bvh_nodes[right_index];
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
