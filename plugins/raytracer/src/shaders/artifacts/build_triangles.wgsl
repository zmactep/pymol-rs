struct ConvertParams {
    vertex_capacity: u32,
    triangle_offset: u32,
    source_triangle_start: u32,
    output_triangle_count: u32,
    atom_offset: u32,
    rep_slot: u32,
    transparency: f32,
    dispatch_width: u32,
}

// {{INCLUDE_ARTIFACT_TRIANGLE}}
// {{INCLUDE_ARTIFACT_ATOM_ALPHA}}
// {{INCLUDE_ARTIFACT_COLOR}}

@group(0) @binding(0) var<storage, read> std_vertex_words: array<u32>;
@group(0) @binding(1) var<storage, read> color_lut: array<ColorLutEntry>;
@group(0) @binding(2) var<storage, read_write> triangles: array<Triangle>;
@group(0) @binding(3) var<storage, read> scene_atoms: array<AtomGpu>;
@group(0) @binding(4) var<storage, read> draw_args: array<u32>;
@group(0) @binding(5) var<uniform> params: ConvertParams;

// {{INCLUDE_ARTIFACT_STD_VERTEX}}

@compute @workgroup_size(128)
fn build_triangles(@builtin(global_invocation_id) gid: vec3<u32>) {
    let triangle_index = gid.y * params.dispatch_width + gid.x;
    if triangle_index >= params.output_triangle_count {
        return;
    }

    let active_vertex_count = min(draw_args[0], params.vertex_capacity) / 3u * 3u;
    let triangle_capacity = active_vertex_count / 3u;
    let source_triangle_index = params.source_triangle_start + triangle_index;
    if source_triangle_index >= triangle_capacity {
        return;
    }

    let out_index = params.triangle_offset + triangle_index;
    let base_vertex = source_triangle_index * 3u;
    let v0 = base_vertex;
    let v1 = base_vertex + 1u;
    let v2 = base_vertex + 2u;
    let material = material_for_group(vertex_group_id(v0));

    triangles[out_index] = Triangle(
        vec4<f32>(vertex_position(v0), 0.0),
        vec4<f32>(vertex_position(v1), 0.0),
        vec4<f32>(vertex_position(v2), 0.0),
        material,
        vec4<u32>(vertex_normal_oct(v0), vertex_normal_oct(v1), vertex_normal_oct(v2), 0u),
    );
}
