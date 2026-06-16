struct DownsampleParams {
    src_width: u32,
    dst_width: u32,
    dst_height: u32,
    factor: u32,
}

@group(0) @binding(0) var<storage, read> src_pixels: array<u32>;
@group(0) @binding(1) var<storage, read_write> dst_pixels: array<u32>;
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
            rgba = rgba + unpack4x8unorm(src_pixels[src_y * params.src_width + src_x]);
        }
    }

    let scale = 1.0 / f32(params.factor * params.factor);
    let dst_index = global_id.y * params.dst_width + global_id.x;
    dst_pixels[dst_index] = pack4x8unorm(rgba * scale);
}
