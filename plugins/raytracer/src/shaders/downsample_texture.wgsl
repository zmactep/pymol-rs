struct DownsampleParams {
    src_width: u32,
    dst_width: u32,
    dst_height: u32,
    factor: u32,
}

@group(0) @binding(0) var src_texture: texture_2d<f32>;
@group(0) @binding(1) var dst_texture: texture_storage_2d<rgba8unorm, write>;
@group(0) @binding(2) var<uniform> params: DownsampleParams;

@compute @workgroup_size(8, 8, 1)
fn main(@builtin(global_invocation_id) global_id: vec3<u32>) {
    if global_id.x >= params.dst_width || global_id.y >= params.dst_height {
        return;
    }

    var rgba = vec4<f32>(0.0);
    for (var sy = 0u; sy < params.factor; sy = sy + 1u) {
        for (var sx = 0u; sx < params.factor; sx = sx + 1u) {
            let src_x = global_id.x * params.factor + sx;
            let src_y = global_id.y * params.factor + sy;
            rgba = rgba + textureLoad(src_texture, vec2<i32>(i32(src_x), i32(src_y)), 0);
        }
    }

    let scale = 1.0 / f32(params.factor * params.factor);
    textureStore(dst_texture, vec2<i32>(i32(global_id.x), i32(global_id.y)), rgba * scale);
}
