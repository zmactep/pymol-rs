struct ConvertParams {
    view_matrix: mat4x4<f32>,
    proj_matrix: mat4x4<f32>,
    source_vertex_count: u32,
    source_triangle_start: u32,
    source_triangle_count: u32,
    triangle_offset: u32,
    output_triangle_capacity: u32,
    atom_offset: u32,
    rep_slot: u32,
    transparency: f32,
    dispatch_width: u32,
    counter_index: u32,
    _pad0: u32,
    _pad1: u32,
}

// {{INCLUDE_ARTIFACT_TRIANGLE}}
// {{INCLUDE_ARTIFACT_ATOM_ALPHA}}
// {{INCLUDE_ARTIFACT_COLOR}}

@group(0) @binding(0) var<storage, read> std_vertex_words: array<u32>;
@group(0) @binding(1) var<storage, read> color_lut: array<ColorLutEntry>;
@group(0) @binding(2) var<storage, read_write> triangles: array<Triangle>;
@group(0) @binding(3) var<storage, read> scene_atoms: array<AtomGpu>;
@group(0) @binding(4) var<storage, read> draw_args: array<u32>;
@group(0) @binding(5) var<storage, read_write> visible_counts: array<atomic<u32>>;
@group(0) @binding(6) var<uniform> params: ConvertParams;

const CLIP_EPS: f32 = 0.0001;

// {{INCLUDE_ARTIFACT_STD_VERTEX}}

fn clip_position(position: vec3<f32>) -> vec4<f32> {
    return params.proj_matrix * (params.view_matrix * vec4<f32>(position, 1.0));
}

fn clip_margin(clip: vec4<f32>) -> f32 {
    return CLIP_EPS * max(abs(clip.w), 1.0);
}

fn all3(a: bool, b: bool, c: bool) -> bool {
    return a && b && c;
}

fn triangle_intersects_output(v0: vec3<f32>, v1: vec3<f32>, v2: vec3<f32>) -> bool {
    let c0 = clip_position(v0);
    let c1 = clip_position(v1);
    let c2 = clip_position(v2);
    let m0 = clip_margin(c0);
    let m1 = clip_margin(c1);
    let m2 = clip_margin(c2);

    if all3(c0.w <= 0.0, c1.w <= 0.0, c2.w <= 0.0) {
        return false;
    }
    if all3(c0.x < -c0.w - m0, c1.x < -c1.w - m1, c2.x < -c2.w - m2) {
        return false;
    }
    if all3(c0.x > c0.w + m0, c1.x > c1.w + m1, c2.x > c2.w + m2) {
        return false;
    }
    if all3(c0.y < -c0.w - m0, c1.y < -c1.w - m1, c2.y < -c2.w - m2) {
        return false;
    }
    if all3(c0.y > c0.w + m0, c1.y > c1.w + m1, c2.y > c2.w + m2) {
        return false;
    }
    if all3(c0.z < -m0, c1.z < -m1, c2.z < -m2) {
        return false;
    }
    if all3(c0.z > c0.w + m0, c1.z > c1.w + m1, c2.z > c2.w + m2) {
        return false;
    }
    return true;
}

@compute @workgroup_size(128)
fn build_visible_triangles(@builtin(global_invocation_id) gid: vec3<u32>) {
    let triangle_index = gid.y * params.dispatch_width + gid.x;
    if triangle_index >= params.source_triangle_count {
        return;
    }

    let active_vertex_count = min(draw_args[0], params.source_vertex_count) / 3u * 3u;
    let source_triangle_index = params.source_triangle_start + triangle_index;
    if source_triangle_index >= active_vertex_count / 3u {
        return;
    }

    let base_vertex = source_triangle_index * 3u;
    let v0 = base_vertex;
    let v1 = base_vertex + 1u;
    let v2 = base_vertex + 2u;
    if !triangle_intersects_output(
        vertex_position(v0),
        vertex_position(v1),
        vertex_position(v2),
    ) {
        return;
    }

    let visible_index = atomicAdd(&visible_counts[params.counter_index], 1u);
    if visible_index >= params.output_triangle_capacity {
        return;
    }
    let out_index = params.triangle_offset + visible_index;
    let material = material_for_group(vertex_group_id(v0));

    triangles[out_index] = Triangle(
        vec4<f32>(vertex_position(v0), 0.0),
        vec4<f32>(vertex_position(v1), 0.0),
        vec4<f32>(vertex_position(v2), 0.0),
        material,
        vec4<u32>(vertex_normal_oct(v0), vertex_normal_oct(v1), vertex_normal_oct(v2), 0u),
    );
}
