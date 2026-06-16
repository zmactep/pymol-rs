struct ConvertParams {
    vertex_capacity: u32,
    triangle_offset: u32,
    atom_offset: u32,
    rep_slot: u32,
    transparency: f32,
    dispatch_width: u32,
    _pad0: u32,
    _pad1: u32,
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

// {{INCLUDE_ARTIFACT_COLOR}}

@group(0) @binding(0) var<storage, read> std_vertex_words: array<u32>;
@group(0) @binding(1) var<storage, read> color_lut: array<ColorLutEntry>;
@group(0) @binding(2) var<storage, read_write> triangles: array<Triangle>;
@group(0) @binding(3) var<storage, read> draw_args: array<u32>;
@group(0) @binding(4) var<uniform> params: ConvertParams;

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

fn vertex_normal_oct(vertex_index: u32) -> u32 {
    return vertex_word(vertex_index, 3u);
}

fn vertex_group_id(vertex_index: u32) -> u32 {
    return vertex_word(vertex_index, 4u);
}

fn sign_extend_16(value: u32) -> i32 {
    let low = i32(value & 0xffffu);
    if low >= 32768 {
        return low - 65536;
    }
    return low;
}

fn oct_decode(packed: u32) -> vec3<f32> {
    var x = max(f32(sign_extend_16(packed)) / 32767.0, -1.0);
    var y = max(f32(sign_extend_16(packed >> 16u)) / 32767.0, -1.0);
    let z = 1.0 - abs(x) - abs(y);
    var out: vec3<f32>;
    if z < 0.0 {
        let ox = (1.0 - abs(y)) * select(-1.0, 1.0, x >= 0.0);
        let oy = (1.0 - abs(x)) * select(-1.0, 1.0, y >= 0.0);
        out = vec3<f32>(ox, oy, z);
    } else {
        out = vec3<f32>(x, y, z);
    }
    return normalize(out);
}

fn material_for_group(group_id: u32) -> vec4<f32> {
    let entry = color_lut[params.atom_offset + group_id];
    let packed = rep_packed_color(entry.reps, params.rep_slot);
    var rgba = select(unpack_rep_rgb8(packed), entry.base, packed == REP_COLOR_INHERIT);
    rgba.a = clamp(rgba.a * (1.0 - params.transparency), 0.0, 1.0);
    return rgba;
}

fn invalid_triangle() -> Triangle {
    return Triangle(
        vec3<f32>(0.0), 0.0,
        vec3<f32>(0.0), 0.0,
        vec3<f32>(0.0), 0.0,
        vec3<f32>(0.0), 0.0,
        vec3<f32>(0.0), 0.0,
        vec3<f32>(0.0), 0.0,
        vec4<f32>(0.0),
        1.0,
        0.0, 0.0, 0.0,
    );
}

@compute @workgroup_size(128)
fn build_triangles(@builtin(global_invocation_id) gid: vec3<u32>) {
    let triangle_index = gid.y * params.dispatch_width + gid.x;
    let triangle_capacity = params.vertex_capacity / 3u;
    if triangle_index >= triangle_capacity {
        return;
    }

    let out_index = params.triangle_offset + triangle_index;
    let active_vertex_count = min(draw_args[0], params.vertex_capacity);
    if triangle_index >= active_vertex_count / 3u {
        triangles[out_index] = invalid_triangle();
        return;
    }

    let base_vertex = triangle_index * 3u;
    let v0 = base_vertex;
    let v1 = base_vertex + 1u;
    let v2 = base_vertex + 2u;
    let material = material_for_group(vertex_group_id(v0));

    triangles[out_index] = Triangle(
        vertex_position(v0), 0.0,
        vertex_position(v1), 0.0,
        vertex_position(v2), 0.0,
        oct_decode(vertex_normal_oct(v0)), 0.0,
        oct_decode(vertex_normal_oct(v1)), 0.0,
        oct_decode(vertex_normal_oct(v2)), 0.0,
        material,
        1.0 - material.a,
        0.0, 0.0, 0.0,
    );
}
