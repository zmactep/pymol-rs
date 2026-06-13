// Octahedral normal decoding — mirrors `oct_encode` in src/representations/mesh.rs.
// 16 bits per component, packed into a u32. See the Cigolle et al. 2014 paper.

fn oct_decode_snorm16(packed: u32) -> vec2<f32> {
    let lo = i32(packed & 0xFFFFu);
    let hi = i32((packed >> 16u) & 0xFFFFu);
    let x_int = select(lo, lo - 0x10000, lo >= 0x8000);
    let y_int = select(hi, hi - 0x10000, hi >= 0x8000);
    let x = max(f32(x_int) / 32767.0, -1.0);
    let y = max(f32(y_int) / 32767.0, -1.0);
    return vec2<f32>(x, y);
}

fn oct_decode(packed: u32) -> vec3<f32> {
    let e = oct_decode_snorm16(packed);
    var n = vec3<f32>(e.x, e.y, 1.0 - abs(e.x) - abs(e.y));
    if n.z < 0.0 {
        let t = (vec2<f32>(1.0) - abs(n.yx)) * select(vec2<f32>(-1.0), vec2<f32>(1.0), n.xy >= vec2<f32>(0.0));
        n.x = t.x;
        n.y = t.y;
    }
    return normalize(n);
}

fn oct_to_snorm16(v: f32) -> u32 {
    let s = i32(round(clamp(v, -1.0, 1.0) * 32767.0));
    return u32(s) & 0xFFFFu;
}

// Mirror of `oct_encode` in src/representations/mesh.rs. Used by GPU passes
// that need to emit `StdVertex.normal_oct` (e.g. surface marching-cubes).
fn oct_encode(n_in: vec3<f32>) -> u32 {
    let len = max(length(n_in), 1e-12);
    var nx = n_in.x / len;
    var ny = n_in.y / len;
    let nz = n_in.z / len;
    let inv = 1.0 / (abs(nx) + abs(ny) + abs(nz) + 1e-12);
    nx = nx * inv;
    ny = ny * inv;
    if nz < 0.0 {
        let sx = select(-1.0, 1.0, nx >= 0.0);
        let sy = select(-1.0, 1.0, ny >= 0.0);
        let tx = (1.0 - abs(ny)) * sx;
        let ty = (1.0 - abs(nx)) * sy;
        nx = tx;
        ny = ty;
    }
    return oct_to_snorm16(nx) | (oct_to_snorm16(ny) << 16u);
}
