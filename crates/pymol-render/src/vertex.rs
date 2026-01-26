//! Vertex formats for molecular representations
//!
//! These structures define the GPU-friendly vertex layouts for each
//! representation type. All structs derive `Pod` and `Zeroable` for
//! safe conversion to byte slices.

use bytemuck::{Pod, Zeroable};

/// Vertex data for sphere impostor rendering
///
/// Each sphere is rendered as a billboard quad with ray-sphere intersection
/// computed in the fragment shader. This vertex contains per-instance data.
#[repr(C)]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
pub struct SphereVertex {
    /// Sphere center position in world space
    pub center: [f32; 3],
    /// Sphere radius
    pub radius: f32,
    /// RGBA color
    pub color: [f32; 4],
}

impl SphereVertex {
    /// Vertex buffer layout for sphere instances
    pub fn layout() -> wgpu::VertexBufferLayout<'static> {
        wgpu::VertexBufferLayout {
            array_stride: std::mem::size_of::<Self>() as wgpu::BufferAddress,
            step_mode: wgpu::VertexStepMode::Instance,
            attributes: &[
                // center
                wgpu::VertexAttribute {
                    format: wgpu::VertexFormat::Float32x3,
                    offset: 0,
                    shader_location: 0,
                },
                // radius
                wgpu::VertexAttribute {
                    format: wgpu::VertexFormat::Float32,
                    offset: 12,
                    shader_location: 1,
                },
                // color
                wgpu::VertexAttribute {
                    format: wgpu::VertexFormat::Float32x4,
                    offset: 16,
                    shader_location: 2,
                },
            ],
        }
    }
}

/// Vertex data for cylinder impostor rendering (sticks/bonds)
///
/// Each cylinder is rendered as an oriented quad with ray-cylinder intersection
/// computed in the fragment shader. Supports half-bond coloring.
#[repr(C)]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
pub struct CylinderVertex {
    /// Start position of cylinder axis
    pub start: [f32; 3],
    /// Cylinder radius
    pub radius: f32,
    /// End position of cylinder axis
    pub end: [f32; 3],
    /// Flags: bit 0 = cap1, bit 1 = cap2
    pub flags: u32,
    /// RGBA color at start
    pub color1: [f32; 4],
    /// RGBA color at end
    pub color2: [f32; 4],
}

impl CylinderVertex {
    /// Flag for rendering cap at start
    pub const CAP_START: u32 = 1;
    /// Flag for rendering cap at end
    pub const CAP_END: u32 = 2;
    /// Flag for both caps
    pub const CAP_BOTH: u32 = Self::CAP_START | Self::CAP_END;

    /// Vertex buffer layout for cylinder instances
    pub fn layout() -> wgpu::VertexBufferLayout<'static> {
        wgpu::VertexBufferLayout {
            array_stride: std::mem::size_of::<Self>() as wgpu::BufferAddress,
            step_mode: wgpu::VertexStepMode::Instance,
            attributes: &[
                // start
                wgpu::VertexAttribute {
                    format: wgpu::VertexFormat::Float32x3,
                    offset: 0,
                    shader_location: 0,
                },
                // radius
                wgpu::VertexAttribute {
                    format: wgpu::VertexFormat::Float32,
                    offset: 12,
                    shader_location: 1,
                },
                // end
                wgpu::VertexAttribute {
                    format: wgpu::VertexFormat::Float32x3,
                    offset: 16,
                    shader_location: 2,
                },
                // flags
                wgpu::VertexAttribute {
                    format: wgpu::VertexFormat::Uint32,
                    offset: 28,
                    shader_location: 3,
                },
                // color1
                wgpu::VertexAttribute {
                    format: wgpu::VertexFormat::Float32x4,
                    offset: 32,
                    shader_location: 4,
                },
                // color2
                wgpu::VertexAttribute {
                    format: wgpu::VertexFormat::Float32x4,
                    offset: 48,
                    shader_location: 5,
                },
            ],
        }
    }
}

/// Vertex data for line rendering
///
/// Simple line representation connecting bonded atoms.
#[repr(C)]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
pub struct LineVertex {
    /// Position in world space
    pub position: [f32; 3],
    /// RGBA color
    pub color: [f32; 4],
}

impl LineVertex {
    /// Vertex buffer layout for line vertices
    pub fn layout() -> wgpu::VertexBufferLayout<'static> {
        wgpu::VertexBufferLayout {
            array_stride: std::mem::size_of::<Self>() as wgpu::BufferAddress,
            step_mode: wgpu::VertexStepMode::Vertex,
            attributes: &[
                // position
                wgpu::VertexAttribute {
                    format: wgpu::VertexFormat::Float32x3,
                    offset: 0,
                    shader_location: 0,
                },
                // color
                wgpu::VertexAttribute {
                    format: wgpu::VertexFormat::Float32x4,
                    offset: 12,
                    shader_location: 1,
                },
            ],
        }
    }
}

/// Vertex data for dot rendering
///
/// Dots are rendered as point sprites or small quads.
#[repr(C)]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
pub struct DotVertex {
    /// Position in world space
    pub position: [f32; 3],
    /// Point size (for point sprites)
    pub size: f32,
    /// RGBA color
    pub color: [f32; 4],
}

impl DotVertex {
    /// Vertex buffer layout for dot vertices
    pub fn layout() -> wgpu::VertexBufferLayout<'static> {
        wgpu::VertexBufferLayout {
            array_stride: std::mem::size_of::<Self>() as wgpu::BufferAddress,
            step_mode: wgpu::VertexStepMode::Instance,
            attributes: &[
                // position
                wgpu::VertexAttribute {
                    format: wgpu::VertexFormat::Float32x3,
                    offset: 0,
                    shader_location: 0,
                },
                // size
                wgpu::VertexAttribute {
                    format: wgpu::VertexFormat::Float32,
                    offset: 12,
                    shader_location: 1,
                },
                // color
                wgpu::VertexAttribute {
                    format: wgpu::VertexFormat::Float32x4,
                    offset: 16,
                    shader_location: 2,
                },
            ],
        }
    }
}

/// Vertex data for mesh rendering
///
/// Standard mesh vertex with position, normal, and color for lit rendering.
#[repr(C)]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
pub struct MeshVertex {
    /// Position in world space
    pub position: [f32; 3],
    /// Normal vector for lighting
    pub normal: [f32; 3],
    /// RGBA color
    pub color: [f32; 4],
}

impl MeshVertex {
    /// Vertex buffer layout for mesh vertices
    pub fn layout() -> wgpu::VertexBufferLayout<'static> {
        wgpu::VertexBufferLayout {
            array_stride: std::mem::size_of::<Self>() as wgpu::BufferAddress,
            step_mode: wgpu::VertexStepMode::Vertex,
            attributes: &[
                // position
                wgpu::VertexAttribute {
                    format: wgpu::VertexFormat::Float32x3,
                    offset: 0,
                    shader_location: 0,
                },
                // normal
                wgpu::VertexAttribute {
                    format: wgpu::VertexFormat::Float32x3,
                    offset: 12,
                    shader_location: 1,
                },
                // color
                wgpu::VertexAttribute {
                    format: wgpu::VertexFormat::Float32x4,
                    offset: 24,
                    shader_location: 2,
                },
            ],
        }
    }
}

/// Billboard quad vertices for impostor rendering
///
/// These are the corner offsets for a unit quad centered at origin.
/// Used by sphere and cylinder impostors.
pub const BILLBOARD_VERTICES: [[f32; 2]; 4] = [
    [-1.0, -1.0], // bottom-left
    [1.0, -1.0],  // bottom-right
    [1.0, 1.0],   // top-right
    [-1.0, 1.0],  // top-left
];

/// Indices for a quad made of two triangles (u16 for wgpu::IndexFormat::Uint16)
pub const QUAD_INDICES: [u16; 6] = [0, 1, 2, 0, 2, 3];

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sphere_vertex_size() {
        assert_eq!(std::mem::size_of::<SphereVertex>(), 32);
    }

    #[test]
    fn test_cylinder_vertex_size() {
        assert_eq!(std::mem::size_of::<CylinderVertex>(), 64);
    }

    #[test]
    fn test_line_vertex_size() {
        assert_eq!(std::mem::size_of::<LineVertex>(), 28);
    }

    #[test]
    fn test_dot_vertex_size() {
        assert_eq!(std::mem::size_of::<DotVertex>(), 32);
    }

    #[test]
    fn test_mesh_vertex_size() {
        assert_eq!(std::mem::size_of::<MeshVertex>(), 40);
    }
}
