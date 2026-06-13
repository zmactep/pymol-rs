//! Shared mesh-vertex types for `StdVertex`-based representations
//! (cartoon, surface, mesh-wireframe).
//!
//! `MeshRep` as a separate representation no longer exists — atom-driven
//! mesh wireframe is handled by `SurfaceRep` in `SurfaceMode::Mesh`
//! (mirrors CartoonRep's Cartoon/Ribbon mode pattern). The shared
//! `surface_mc.wgsl` compute kernel emits LineList output when
//! `McParams.emit_lines = 1`. See `representations/surface/mod.rs`.

use bytemuck::{Pod, Zeroable};

/// 24-byte standard mesh vertex. `normal_oct` packs an
/// octahedral-encoded unit normal (16 bits per component); `flags` is reserved
/// for cap / sheet-arrow style toggles.
#[repr(C)]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
pub struct StdVertex {
    pub position: [f32; 3],
    pub normal_oct: u32,
    pub group_id: u32,
    pub flags: u32,
}

impl StdVertex {
    pub const SIZE: u64 = std::mem::size_of::<Self>() as u64;

    pub fn vertex_layout() -> wgpu::VertexBufferLayout<'static> {
        wgpu::VertexBufferLayout {
            array_stride: Self::SIZE,
            step_mode: wgpu::VertexStepMode::Vertex,
            attributes: &[
                wgpu::VertexAttribute {
                    format: wgpu::VertexFormat::Float32x3,
                    offset: 0,
                    shader_location: 0,
                },
                wgpu::VertexAttribute {
                    format: wgpu::VertexFormat::Uint32,
                    offset: 12,
                    shader_location: 1,
                },
                wgpu::VertexAttribute {
                    format: wgpu::VertexFormat::Uint32,
                    offset: 16,
                    shader_location: 2,
                },
                wgpu::VertexAttribute {
                    format: wgpu::VertexFormat::Uint32,
                    offset: 20,
                    shader_location: 3,
                },
            ],
        }
    }
}

/// Encode a unit-length normal into a single u32 (16-bit snorm per axis).
/// Mirrors `oct_decode` in `shaders/common/octahedral.wgsl`.
pub fn oct_encode(n: [f32; 3]) -> u32 {
    let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt().max(1e-12);
    let mut nx = n[0] / len;
    let mut ny = n[1] / len;
    let nz = n[2] / len;
    let inv = 1.0 / (nx.abs() + ny.abs() + nz.abs() + 1e-12);
    nx *= inv;
    ny *= inv;
    if nz < 0.0 {
        let tx = (1.0 - ny.abs()) * if nx >= 0.0 { 1.0 } else { -1.0 };
        let ty = (1.0 - nx.abs()) * if ny >= 0.0 { 1.0 } else { -1.0 };
        nx = tx;
        ny = ty;
    }
    let to_snorm16 = |v: f32| -> u32 {
        let s = (v.clamp(-1.0, 1.0) * 32767.0).round() as i32;
        (s as u32) & 0xFFFF
    };
    to_snorm16(nx) | (to_snorm16(ny) << 16)
}

#[cfg(test)]
mod tests {
    use super::oct_encode;

    fn round_trip(n: [f32; 3]) -> [f32; 3] {
        // CPU mirror of the WGSL decoder, for unit testing without GPU.
        let packed = oct_encode(n);
        let lo = (packed & 0xFFFF) as i16;
        let hi = ((packed >> 16) & 0xFFFF) as i16;
        let x = (lo as f32 / 32767.0).max(-1.0);
        let y = (hi as f32 / 32767.0).max(-1.0);
        let mut nx = x;
        let mut ny = y;
        let nz = 1.0 - nx.abs() - ny.abs();
        if nz < 0.0 {
            let tx = (1.0 - ny.abs()) * if nx >= 0.0 { 1.0 } else { -1.0 };
            let ty = (1.0 - nx.abs()) * if ny >= 0.0 { 1.0 } else { -1.0 };
            nx = tx;
            ny = ty;
        }
        let l = (nx * nx + ny * ny + nz * nz).sqrt().max(1e-12);
        [nx / l, ny / l, nz / l]
    }

    fn close(a: [f32; 3], b: [f32; 3], eps: f32) -> bool {
        (a[0] - b[0]).abs() < eps && (a[1] - b[1]).abs() < eps && (a[2] - b[2]).abs() < eps
    }

    #[test]
    fn oct_round_trip_axes() {
        for &n in &[
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
        ] {
            let r = round_trip(n);
            assert!(close(n, r, 1e-3), "axis {:?} round-tripped to {:?}", n, r);
        }
    }

    #[test]
    fn oct_round_trip_diagonals() {
        let s = (1.0_f32 / 3.0).sqrt();
        for &sign in &[[s, s, s], [-s, s, -s], [s, -s, -s]] {
            let r = round_trip(sign);
            assert!(
                close(sign, r, 1e-3),
                "diag {:?} round-tripped to {:?}",
                sign,
                r
            );
        }
    }
}
