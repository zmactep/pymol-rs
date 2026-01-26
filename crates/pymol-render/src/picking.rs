//! GPU-based picking system
//!
//! Provides color-coded rendering for object/atom selection via mouse picking.
//! Objects and atoms are rendered with unique colors that encode their IDs,
//! then the pixel at the mouse position is read back to determine what was hit.

use crate::vertex::MeshVertex;
use bytemuck::{Pod, Zeroable};

/// Result of a picking operation
#[derive(Debug, Clone)]
pub struct PickResult {
    /// Object ID (or index in registry)
    pub object_id: u32,
    /// Atom index within the object (u32::MAX if not an atom)
    pub atom_index: u32,
    /// Depth value (0.0 = near, 1.0 = far)
    pub depth: f32,
    /// World-space position of the hit point (if computed)
    pub position: Option<[f32; 3]>,
}

impl PickResult {
    /// Create a new pick result
    pub fn new(object_id: u32, atom_index: u32, depth: f32) -> Self {
        Self {
            object_id,
            atom_index,
            depth,
            position: None,
        }
    }

    /// Check if this hit an atom (vs just an object)
    pub fn is_atom_hit(&self) -> bool {
        self.atom_index != u32::MAX
    }

    /// Create a "no hit" result
    pub fn none() -> Self {
        Self {
            object_id: u32::MAX,
            atom_index: u32::MAX,
            depth: 1.0,
            position: None,
        }
    }

    /// Check if anything was hit
    pub fn is_hit(&self) -> bool {
        self.object_id != u32::MAX
    }
}

/// Picking vertex format
///
/// Encodes object and atom IDs as colors for GPU picking.
#[repr(C)]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
pub struct PickingVertex {
    /// Position in world space
    pub position: [f32; 3],
    /// Object ID encoded in R channel
    pub object_id: u32,
    /// Atom ID encoded in G,B channels (low 16 bits in G, high 16 bits in B)
    pub atom_id: u32,
}

impl PickingVertex {
    /// Vertex buffer layout for picking
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
                // object_id
                wgpu::VertexAttribute {
                    format: wgpu::VertexFormat::Uint32,
                    offset: 12,
                    shader_location: 1,
                },
                // atom_id
                wgpu::VertexAttribute {
                    format: wgpu::VertexFormat::Uint32,
                    offset: 16,
                    shader_location: 2,
                },
            ],
        }
    }
}

/// Picking ID encoder/decoder
///
/// Encodes object and atom IDs into RGBA colors and decodes them back.
#[derive(Debug, Clone, Copy)]
pub struct PickingId {
    /// Object ID (0-255 for simple cases)
    pub object_id: u32,
    /// Atom index within the object
    pub atom_index: u32,
}

impl PickingId {
    /// Create a new picking ID
    pub fn new(object_id: u32, atom_index: u32) -> Self {
        Self { object_id, atom_index }
    }

    /// Create a picking ID for just an object (no specific atom)
    pub fn object(object_id: u32) -> Self {
        Self {
            object_id,
            atom_index: u32::MAX,
        }
    }

    /// Encode as RGBA color (u8 per channel)
    ///
    /// Layout:
    /// - R: object_id low byte
    /// - G: atom_index bits 0-7
    /// - B: atom_index bits 8-15
    /// - A: atom_index bits 16-23
    pub fn to_rgba(&self) -> [u8; 4] {
        [
            (self.object_id & 0xFF) as u8,
            (self.atom_index & 0xFF) as u8,
            ((self.atom_index >> 8) & 0xFF) as u8,
            ((self.atom_index >> 16) & 0xFF) as u8,
        ]
    }

    /// Decode from RGBA color
    pub fn from_rgba(rgba: [u8; 4]) -> Self {
        let object_id = rgba[0] as u32;
        let atom_index = rgba[1] as u32 | ((rgba[2] as u32) << 8) | ((rgba[3] as u32) << 16);
        Self { object_id, atom_index }
    }

    /// Encode as float RGBA (for shader use)
    pub fn to_rgba_f32(&self) -> [f32; 4] {
        let rgba = self.to_rgba();
        [
            rgba[0] as f32 / 255.0,
            rgba[1] as f32 / 255.0,
            rgba[2] as f32 / 255.0,
            rgba[3] as f32 / 255.0,
        ]
    }

    /// Check if this represents "no hit" (background)
    pub fn is_none(&self) -> bool {
        // Background is rendered as all zeros
        self.object_id == 0 && self.atom_index == 0
    }
}

/// Configuration for picking operations
#[derive(Debug, Clone)]
pub struct PickingConfig {
    /// Size of the picking region in pixels (1 = single pixel)
    pub region_size: u32,
    /// Whether to pick the closest hit or any hit
    pub depth_test: bool,
    /// Whether to include labels in picking
    pub pick_labels: bool,
}

impl Default for PickingConfig {
    fn default() -> Self {
        Self {
            region_size: 1,
            depth_test: true,
            pick_labels: true,
        }
    }
}

/// Convert mesh vertices to picking vertices
pub fn mesh_to_picking_vertices(
    mesh: &[MeshVertex],
    object_id: u32,
    atom_indices: Option<&[u32]>,
) -> Vec<PickingVertex> {
    mesh.iter()
        .enumerate()
        .map(|(i, v)| {
            let atom_id = atom_indices
                .and_then(|indices| indices.get(i).copied())
                .unwrap_or(u32::MAX);
            PickingVertex {
                position: v.position,
                object_id,
                atom_id,
            }
        })
        .collect()
}

/// Simple ray for picking calculations
#[derive(Debug, Clone, Copy)]
pub struct Ray {
    /// Ray origin
    pub origin: [f32; 3],
    /// Ray direction (normalized)
    pub direction: [f32; 3],
}

impl Ray {
    /// Create a new ray
    pub fn new(origin: [f32; 3], direction: [f32; 3]) -> Self {
        // Normalize direction
        let len = (direction[0] * direction[0]
            + direction[1] * direction[1]
            + direction[2] * direction[2])
            .sqrt();
        let direction = if len > 0.0 {
            [direction[0] / len, direction[1] / len, direction[2] / len]
        } else {
            [0.0, 0.0, -1.0]
        };
        Self { origin, direction }
    }

    /// Create a ray from screen coordinates and view/projection matrices
    ///
    /// # Arguments
    /// * `x`, `y` - Screen coordinates (pixels)
    /// * `width`, `height` - Screen dimensions
    /// * `view_inv` - Inverse view matrix
    /// * `proj_inv` - Inverse projection matrix
    pub fn from_screen(
        x: f32,
        y: f32,
        width: f32,
        height: f32,
        view_inv: &[[f32; 4]; 4],
        proj_inv: &[[f32; 4]; 4],
    ) -> Self {
        // Convert screen coords to NDC (-1 to 1)
        let ndc_x = (2.0 * x / width) - 1.0;
        let ndc_y = 1.0 - (2.0 * y / height); // Flip Y

        // Near and far points in NDC
        let near_ndc = [ndc_x, ndc_y, -1.0, 1.0];
        let far_ndc = [ndc_x, ndc_y, 1.0, 1.0];

        // Transform through inverse projection
        let near_clip = Self::mat4_mul_vec4(proj_inv, near_ndc);
        let far_clip = Self::mat4_mul_vec4(proj_inv, far_ndc);

        // Perspective divide
        let near_view = [
            near_clip[0] / near_clip[3],
            near_clip[1] / near_clip[3],
            near_clip[2] / near_clip[3],
        ];
        let far_view = [
            far_clip[0] / far_clip[3],
            far_clip[1] / far_clip[3],
            far_clip[2] / far_clip[3],
        ];

        // Transform to world space
        let near_world = Self::mat4_mul_point(view_inv, near_view);
        let far_world = Self::mat4_mul_point(view_inv, far_view);

        // Ray direction
        let direction = [
            far_world[0] - near_world[0],
            far_world[1] - near_world[1],
            far_world[2] - near_world[2],
        ];

        Self::new(near_world, direction)
    }

    fn mat4_mul_vec4(m: &[[f32; 4]; 4], v: [f32; 4]) -> [f32; 4] {
        [
            m[0][0] * v[0] + m[1][0] * v[1] + m[2][0] * v[2] + m[3][0] * v[3],
            m[0][1] * v[0] + m[1][1] * v[1] + m[2][1] * v[2] + m[3][1] * v[3],
            m[0][2] * v[0] + m[1][2] * v[1] + m[2][2] * v[2] + m[3][2] * v[3],
            m[0][3] * v[0] + m[1][3] * v[1] + m[2][3] * v[2] + m[3][3] * v[3],
        ]
    }

    fn mat4_mul_point(m: &[[f32; 4]; 4], p: [f32; 3]) -> [f32; 3] {
        let w = m[0][3] * p[0] + m[1][3] * p[1] + m[2][3] * p[2] + m[3][3];
        let inv_w = if w.abs() > 1e-6 { 1.0 / w } else { 1.0 };
        [
            (m[0][0] * p[0] + m[1][0] * p[1] + m[2][0] * p[2] + m[3][0]) * inv_w,
            (m[0][1] * p[0] + m[1][1] * p[1] + m[2][1] * p[2] + m[3][1]) * inv_w,
            (m[0][2] * p[0] + m[1][2] * p[1] + m[2][2] * p[2] + m[3][2]) * inv_w,
        ]
    }

    /// Test intersection with a sphere
    ///
    /// Returns the distance along the ray to the intersection point, or None if no hit.
    pub fn intersect_sphere(&self, center: [f32; 3], radius: f32) -> Option<f32> {
        let oc = [
            self.origin[0] - center[0],
            self.origin[1] - center[1],
            self.origin[2] - center[2],
        ];

        let a = self.direction[0] * self.direction[0]
            + self.direction[1] * self.direction[1]
            + self.direction[2] * self.direction[2];
        let b = 2.0
            * (oc[0] * self.direction[0]
                + oc[1] * self.direction[1]
                + oc[2] * self.direction[2]);
        let c = oc[0] * oc[0] + oc[1] * oc[1] + oc[2] * oc[2] - radius * radius;

        let discriminant = b * b - 4.0 * a * c;
        if discriminant < 0.0 {
            return None;
        }

        let sqrt_d = discriminant.sqrt();
        let t1 = (-b - sqrt_d) / (2.0 * a);
        let t2 = (-b + sqrt_d) / (2.0 * a);

        if t1 > 0.0 {
            Some(t1)
        } else if t2 > 0.0 {
            Some(t2)
        } else {
            None
        }
    }

    /// Get a point along the ray
    pub fn at(&self, t: f32) -> [f32; 3] {
        [
            self.origin[0] + t * self.direction[0],
            self.origin[1] + t * self.direction[1],
            self.origin[2] + t * self.direction[2],
        ]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_picking_id_encode_decode() {
        let id = PickingId::new(42, 12345);
        let rgba = id.to_rgba();
        let decoded = PickingId::from_rgba(rgba);

        assert_eq!(decoded.object_id, 42);
        assert_eq!(decoded.atom_index, 12345);
    }

    #[test]
    fn test_picking_id_object_only() {
        let id = PickingId::object(5);
        assert_eq!(id.object_id, 5);
        assert_eq!(id.atom_index, u32::MAX);
    }

    #[test]
    fn test_picking_id_large_atom_index() {
        let id = PickingId::new(255, 16_777_215); // Max values for encoding
        let rgba = id.to_rgba();
        let decoded = PickingId::from_rgba(rgba);

        assert_eq!(decoded.object_id, 255);
        assert_eq!(decoded.atom_index, 16_777_215);
    }

    #[test]
    fn test_ray_sphere_intersection() {
        let ray = Ray::new([0.0, 0.0, 0.0], [0.0, 0.0, -1.0]);
        let center = [0.0, 0.0, -5.0];
        let radius = 1.0;

        let t = ray.intersect_sphere(center, radius);
        assert!(t.is_some());
        assert!((t.unwrap() - 4.0).abs() < 0.001); // Should hit at z = -4
    }

    #[test]
    fn test_ray_sphere_miss() {
        let ray = Ray::new([0.0, 0.0, 0.0], [0.0, 0.0, -1.0]);
        let center = [10.0, 0.0, -5.0]; // Sphere is off to the side
        let radius = 1.0;

        let t = ray.intersect_sphere(center, radius);
        assert!(t.is_none());
    }

    #[test]
    fn test_pick_result() {
        let result = PickResult::new(1, 42, 0.5);
        assert!(result.is_hit());
        assert!(result.is_atom_hit());

        let none = PickResult::none();
        assert!(!none.is_hit());
        assert!(!none.is_atom_hit());
    }
}
