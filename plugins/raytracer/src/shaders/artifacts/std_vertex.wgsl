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
