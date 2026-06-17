struct CullParams {
    view_matrix: mat4x4<f32>,
    proj_matrix: mat4x4<f32>,
    source_vertex_count: u32,
    triangle_offset: u32,
    output_triangle_count: u32,
    atom_offset: u32,
    rep_slot: u32,
    transparency: f32,
    dispatch_width: u32,
    counter_index: u32,
}

@group(0) @binding(0) var<storage, read> std_vertex_words: array<u32>;
@group(0) @binding(1) var<storage, read_write> visible_counts: array<atomic<u32>>;
@group(0) @binding(2) var<uniform> params: CullParams;

const CLIP_EPS: f32 = 0.0001;

fn vertex_word(vertex_index: u32, lane: u32) -> u32 {
    return std_vertex_words[vertex_index * 6u + lane];
}

fn vertex_position(vertex_index: u32) -> vec3<f32> {
    return vec3<f32>(
        bitcast<f32>(vertex_word(vertex_index, 0u)),
        bitcast<f32>(vertex_word(vertex_index, 1u)),
        bitcast<f32>(vertex_word(vertex_index, 2u)),
    );
}

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
fn count_visible_triangles(@builtin(global_invocation_id) gid: vec3<u32>) {
    let triangle_index = gid.y * params.dispatch_width + gid.x;
    if triangle_index >= params.source_vertex_count / 3u {
        return;
    }

    let base_vertex = triangle_index * 3u;
    if triangle_intersects_output(
        vertex_position(base_vertex),
        vertex_position(base_vertex + 1u),
        vertex_position(base_vertex + 2u),
    ) {
        _ = atomicAdd(&visible_counts[params.counter_index], 1u);
    }
}
