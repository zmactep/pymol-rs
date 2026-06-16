struct RepColorLutEntry {
    sphere: u32,
    stick: u32,
    line: u32,
    dot: u32,
    cartoon: u32,
    ribbon: u32,
    surface: u32,
    mesh: u32,
    ellipsoid: u32,
    _pad0: u32,
    _pad1: u32,
    _pad2: u32,
}

struct ColorLutEntry {
    base: vec4<f32>,
    reps: RepColorLutEntry,
}

const REP_COLOR_INHERIT: u32 = 0xffffffffu;

fn rep_packed_color(reps: RepColorLutEntry, slot: u32) -> u32 {
    switch slot {
        case 0u: { return reps.sphere; }
        case 1u: { return reps.stick; }
        case 2u: { return reps.line; }
        case 3u: { return reps.dot; }
        case 4u: { return reps.cartoon; }
        case 5u: { return reps.ribbon; }
        case 6u: { return reps.surface; }
        case 7u: { return reps.mesh; }
        case 8u: { return reps.ellipsoid; }
        default: { return REP_COLOR_INHERIT; }
    }
}

fn unpack_rep_rgb8(packed: u32) -> vec4<f32> {
    return vec4<f32>(
        f32(packed & 0xffu) / 255.0,
        f32((packed >> 8u) & 0xffu) / 255.0,
        f32((packed >> 16u) & 0xffu) / 255.0,
        1.0,
    );
}
