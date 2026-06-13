@group(0) @binding(0) var id_tex: texture_2d<u32>;
@group(0) @binding(1) var<storage, read> object_atom_offsets: array<u32>;
@group(0) @binding(2) var<storage, read> marker_lut: array<u32>;

const MARKER_SELECTED: u32 = 1u << 0u;
const MARKER_HOVER: u32 = 1u << 1u;

@vertex
fn vs_main(@builtin(vertex_index) vid: u32) -> @builtin(position) vec4<f32> {
    let x = f32((vid << 1u) & 2u);
    let y = f32(vid & 2u);
    return vec4<f32>(x * 2.0 - 1.0, 1.0 - y * 2.0, 0.0, 1.0);
}

fn unpack_object_id(r: u32) -> u32 {
    return (r >> 16u) & 0xFFFu;
}

fn unpack_atom_id(r: u32, g: u32) -> u32 {
    return ((g >> 16u) << 16u) | (r & 0xFFFFu);
}

@fragment
fn fs_main(@builtin(position) pos: vec4<f32>) -> @location(0) vec4<f32> {
    let dims = textureDimensions(id_tex);
    let coord = vec2<i32>(floor(pos.xy));
    if coord.x < 0 || coord.y < 0 || coord.x >= i32(dims.x) || coord.y >= i32(dims.y) {
        discard;
    }

    let packed = textureLoad(id_tex, coord, 0).rg;
    if packed.x == 0u && packed.y == 0u {
        discard;
    }

    let object_id = unpack_object_id(packed.x);
    if object_id >= 4096u {
        discard;
    }

    let global_atom = object_atom_offsets[object_id] + unpack_atom_id(packed.x, packed.y);
    if global_atom >= arrayLength(&marker_lut) {
        discard;
    }
    let marker = marker_lut[global_atom];
    let selected = select(0.0, 1.0, (marker & MARKER_SELECTED) != 0u);
    let hover = select(0.0, 1.0, (marker & MARKER_HOVER) != 0u);
    if selected == 0.0 && hover == 0.0 {
        discard;
    }
    return vec4<f32>(selected, hover, 0.0, 1.0);
}
