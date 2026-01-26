//! Mesh representation for surface and mesh rendering
//!
//! Renders triangle meshes with Phong lighting.

use crate::buffer::GrowableBuffer;
use crate::color_resolver::ColorResolver;
use crate::representation::Representation;
use crate::vertex::MeshVertex;

use pymol_mol::{CoordSet, ObjectMolecule, RepMask};
use pymol_settings::SettingResolver;

/// Mesh representation
///
/// Renders triangle meshes with Phong lighting. This is a general-purpose
/// mesh renderer that can be used for surfaces, isomeshes, and other
/// triangle-based geometry.
pub struct MeshRep {
    /// Vertex data (CPU side)
    vertices: Vec<MeshVertex>,
    /// Index data (CPU side)
    indices: Vec<u32>,
    /// GPU vertex buffer
    vertex_buffer: GrowableBuffer,
    /// GPU index buffer
    index_buffer: GrowableBuffer,
    /// Pipeline for rendering
    pipeline: Option<wgpu::RenderPipeline>,
    /// Uniform bind group
    bind_group: Option<wgpu::BindGroup>,
    /// Whether the representation needs to be rebuilt
    dirty: bool,
    /// Number of indices to render
    index_count: u32,
}

impl MeshRep {
    /// Create a new mesh representation
    pub fn new() -> Self {
        Self {
            vertices: Vec::new(),
            indices: Vec::new(),
            vertex_buffer: GrowableBuffer::new("Mesh Vertices", wgpu::BufferUsages::VERTEX),
            index_buffer: GrowableBuffer::new("Mesh Indices", wgpu::BufferUsages::INDEX),
            pipeline: None,
            bind_group: None,
            dirty: true,
            index_count: 0,
        }
    }

    /// Set the render pipeline and bind group
    pub fn set_pipeline(&mut self, pipeline: wgpu::RenderPipeline, bind_group: wgpu::BindGroup) {
        self.pipeline = Some(pipeline);
        self.bind_group = Some(bind_group);
    }

    /// Set mesh data directly
    ///
    /// This allows setting pre-computed mesh data (e.g., from surface generation).
    pub fn set_mesh(&mut self, vertices: Vec<MeshVertex>, indices: Vec<u32>) {
        self.vertices = vertices;
        self.indices = indices;
        self.index_count = self.indices.len() as u32;
        self.dirty = true;
    }

    /// Add a triangle to the mesh
    pub fn add_triangle(
        &mut self,
        p1: [f32; 3],
        p2: [f32; 3],
        p3: [f32; 3],
        color: [f32; 4],
    ) {
        // Calculate face normal
        let v1 = [p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]];
        let v2 = [p3[0] - p1[0], p3[1] - p1[1], p3[2] - p1[2]];
        let normal = Self::normalize([
            v1[1] * v2[2] - v1[2] * v2[1],
            v1[2] * v2[0] - v1[0] * v2[2],
            v1[0] * v2[1] - v1[1] * v2[0],
        ]);

        let base_idx = self.vertices.len() as u32;

        self.vertices.push(MeshVertex {
            position: p1,
            normal,
            color,
        });
        self.vertices.push(MeshVertex {
            position: p2,
            normal,
            color,
        });
        self.vertices.push(MeshVertex {
            position: p3,
            normal,
            color,
        });

        self.indices.push(base_idx);
        self.indices.push(base_idx + 1);
        self.indices.push(base_idx + 2);
    }

    /// Clear the mesh
    pub fn clear(&mut self) {
        self.vertices.clear();
        self.indices.clear();
        self.index_count = 0;
    }

    fn normalize(v: [f32; 3]) -> [f32; 3] {
        let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        if len > 0.0 {
            [v[0] / len, v[1] / len, v[2] / len]
        } else {
            [0.0, 0.0, 1.0]
        }
    }
}

impl Default for MeshRep {
    fn default() -> Self {
        Self::new()
    }
}

impl Representation for MeshRep {
    fn build(
        &mut self,
        molecule: &ObjectMolecule,
        coord_set: &CoordSet,
        colors: &ColorResolver,
        _settings: &SettingResolver,
    ) {
        // For mesh representation, we typically don't build from molecular data
        // directly - the mesh data is usually set via set_mesh() from a surface
        // generation algorithm. However, we can generate a simple representation
        // here for atoms that have MESH visible.

        self.clear();

        // Simple icosahedron-based mesh for atoms with MESH visible
        for (atom_idx, coord) in coord_set.iter_with_atoms() {
            let atom = match molecule.get_atom(atom_idx) {
                Some(a) => a,
                None => continue,
            };

            // Check if mesh representation is visible for this atom
            if !atom.visible_reps.is_visible(RepMask::MESH) {
                continue;
            }

            let radius = atom.effective_vdw();
            let color = colors.resolve_atom(atom, molecule);
            let center = [coord.x, coord.y, coord.z];

            // Generate a simple icosahedron mesh
            self.add_icosahedron(center, radius, color);
        }

        self.index_count = self.indices.len() as u32;
        self.dirty = false;
    }

    fn upload(&mut self, device: &wgpu::Device, queue: &wgpu::Queue) {
        if !self.vertices.is_empty() {
            self.vertex_buffer.write(device, queue, &self.vertices);
        }
        if !self.indices.is_empty() {
            self.index_buffer.write(device, queue, &self.indices);
        }
    }

    fn render<'a>(&'a self, render_pass: &mut wgpu::RenderPass<'a>) {
        if self.index_count == 0 {
            return;
        }

        // The caller is responsible for setting pipeline and bind_group.
        // We set our own vertex and index buffers.
        if let (Some(vertices), Some(indices)) = (
            self.vertex_buffer.buffer(),
            self.index_buffer.buffer(),
        ) {
            render_pass.set_vertex_buffer(0, vertices.slice(..));
            render_pass.set_index_buffer(indices.slice(..), wgpu::IndexFormat::Uint32);
            render_pass.draw_indexed(0..self.index_count, 0, 0..1);
        }
    }

    fn is_dirty(&self) -> bool {
        self.dirty
    }

    fn set_dirty(&mut self) {
        self.dirty = true;
    }

    fn primitive_count(&self) -> usize {
        self.index_count as usize / 3 // Triangles
    }
}

impl MeshRep {
    /// Add an icosahedron mesh centered at a point
    fn add_icosahedron(&mut self, center: [f32; 3], radius: f32, color: [f32; 4]) {
        // Icosahedron vertices
        let t = (1.0 + 5.0_f32.sqrt()) / 2.0;
        let vertices = [
            [-1.0, t, 0.0],
            [1.0, t, 0.0],
            [-1.0, -t, 0.0],
            [1.0, -t, 0.0],
            [0.0, -1.0, t],
            [0.0, 1.0, t],
            [0.0, -1.0, -t],
            [0.0, 1.0, -t],
            [t, 0.0, -1.0],
            [t, 0.0, 1.0],
            [-t, 0.0, -1.0],
            [-t, 0.0, 1.0],
        ];

        // Normalize and scale vertices
        let scale = radius / (1.0 + t * t).sqrt();
        let scaled: Vec<[f32; 3]> = vertices
            .iter()
            .map(|v| {
                [
                    center[0] + v[0] * scale,
                    center[1] + v[1] * scale,
                    center[2] + v[2] * scale,
                ]
            })
            .collect();

        // Icosahedron faces (20 triangles)
        let faces = [
            [0, 11, 5],
            [0, 5, 1],
            [0, 1, 7],
            [0, 7, 10],
            [0, 10, 11],
            [1, 5, 9],
            [5, 11, 4],
            [11, 10, 2],
            [10, 7, 6],
            [7, 1, 8],
            [3, 9, 4],
            [3, 4, 2],
            [3, 2, 6],
            [3, 6, 8],
            [3, 8, 9],
            [4, 9, 5],
            [2, 4, 11],
            [6, 2, 10],
            [8, 6, 7],
            [9, 8, 1],
        ];

        for face in &faces {
            self.add_triangle(scaled[face[0]], scaled[face[1]], scaled[face[2]], color);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mesh_rep_new() {
        let rep = MeshRep::new();
        assert!(rep.is_dirty());
        assert!(rep.is_empty());
    }

    #[test]
    fn test_add_triangle() {
        let mut rep = MeshRep::new();
        rep.add_triangle(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 0.0, 0.0, 1.0],
        );

        assert_eq!(rep.vertices.len(), 3);
        assert_eq!(rep.indices.len(), 3);
    }
}
