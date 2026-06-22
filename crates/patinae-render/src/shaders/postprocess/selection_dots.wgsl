// Lite selected atom dots. Draws compact selected-index buffers, not
// every atom in the scene.

// {{INCLUDE_FRAME}}
// {{INCLUDE_SCENE}}

struct SelectionDotsParams {
    radius_px: f32,
    view_bias: f32,
    _pad0: f32,
    _pad1: f32,
};

@group(3) @binding(0) var<uniform> selection_dots: SelectionDotsParams;
@group(3) @binding(1) var<storage, read> selected_indices: array<u32>;

const MARKER_SELECTED: u32 = 1u << 0u;
const SELECTED_FILL: vec4<f32> = vec4<f32>(1.0, 0.0, 0.7843, 0.82);
const SELECTED_RIM: vec4<f32> = vec4<f32>(1.0, 0.82, 1.0, 0.96);

struct VsOut {
    @builtin(position) clip_position: vec4<f32>,
    @location(0) uv: vec2<f32>,
    @location(1) @interpolate(flat) enabled: u32,
};

fn quad_corner(vertex_index: u32) -> vec2<f32> {
    let v = vertex_index % 6u;
    if v == 0u { return vec2<f32>(-1.0, -1.0); }
    if v == 1u { return vec2<f32>( 1.0, -1.0); }
    if v == 2u { return vec2<f32>( 1.0,  1.0); }
    if v == 3u { return vec2<f32>(-1.0, -1.0); }
    if v == 4u { return vec2<f32>( 1.0,  1.0); }
    return vec2<f32>(-1.0, 1.0);
}

fn inactive_out() -> VsOut {
    var out: VsOut;
    out.clip_position = vec4<f32>(2.0, 2.0, 1.0, 1.0);
    out.uv = vec2<f32>(2.0, 2.0);
    out.enabled = 0u;
    return out;
}

@vertex
fn vs_main(
    @builtin(vertex_index) vertex_index: u32,
    @builtin(instance_index) instance_index: u32,
) -> VsOut {
    if instance_index >= arrayLength(&selected_indices) {
        return inactive_out();
    }

    let local_id = selected_indices[instance_index];
    if local_id >= obj.atom_count {
        return inactive_out();
    }

    let global_id = obj.atom_offset + local_id;
    if !scene_visible(global_id) {
        return inactive_out();
    }
    if (scene_marker(global_id) & MARKER_SELECTED) == 0u {
        return inactive_out();
    }

    var view_pos = (frame.view * vec4<f32>(scene_coord(global_id), 1.0)).xyz;
    view_pos.z = view_pos.z + selection_dots.view_bias;

    var clip = frame.proj * vec4<f32>(view_pos, 1.0);
    let corner = quad_corner(vertex_index);
    let pixel_to_clip = vec2<f32>(frame.viewport.z * 2.0, frame.viewport.w * 2.0);
    let offset = corner * selection_dots.radius_px * pixel_to_clip * clip.w;
    clip = vec4<f32>(clip.x + offset.x, clip.y + offset.y, clip.z, clip.w);

    var out: VsOut;
    out.clip_position = clip;
    out.uv = corner;
    out.enabled = 1u;
    return out;
}

@fragment
fn fs_main(input: VsOut) -> @location(0) vec4<f32> {
    if input.enabled == 0u {
        discard;
    }
    let d2 = dot(input.uv, input.uv);
    if d2 > 1.0 {
        discard;
    }
    let rim = smoothstep(0.58, 0.88, d2);
    return mix(SELECTED_FILL, SELECTED_RIM, rim);
}
