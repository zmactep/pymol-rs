//! Line representation for bond visualization
//!
//! Renders bonds as simple lines connecting atom centers.

use crate::buffer::GrowableBuffer;
use crate::color_resolver::ColorResolver;
use crate::representation::Representation;
use crate::vertex::LineVertex;

use pymol_mol::{CoordSet, ObjectMolecule, RepMask};
use pymol_settings::SettingResolver;

/// Line representation for bonds
///
/// Renders bonds as simple lines. Each bond is drawn as two line segments
/// meeting at the midpoint, with each half colored by its respective atom.
pub struct LineRep {
    /// Vertex data (CPU side)
    vertices: Vec<LineVertex>,
    /// GPU vertex buffer
    vertex_buffer: GrowableBuffer,
    /// Pipeline for rendering
    pipeline: Option<wgpu::RenderPipeline>,
    /// Whether the representation needs to be rebuilt
    dirty: bool,
    /// Number of vertices to render
    vertex_count: u32,
}

impl LineRep {
    /// Create a new line representation
    pub fn new() -> Self {
        Self {
            vertices: Vec::new(),
            vertex_buffer: GrowableBuffer::new("Line Vertices", wgpu::BufferUsages::VERTEX),
            pipeline: None,
            dirty: true,
            vertex_count: 0,
        }
    }

    /// Set the render pipeline
    pub fn set_pipeline(&mut self, pipeline: wgpu::RenderPipeline) {
        self.pipeline = Some(pipeline);
    }
}

impl Default for LineRep {
    fn default() -> Self {
        Self::new()
    }
}

impl Representation for LineRep {
    fn build(
        &mut self,
        molecule: &ObjectMolecule,
        coord_set: &CoordSet,
        colors: &ColorResolver,
        _settings: &SettingResolver,
    ) {
        self.vertices.clear();

        // Iterate over bonds and create line segments
        for bond in molecule.bonds() {
            let atom1_idx = bond.atom1;
            let atom2_idx = bond.atom2;

            // Get atom data
            let atom1 = match molecule.get_atom(atom1_idx) {
                Some(a) => a,
                None => continue,
            };
            let atom2 = match molecule.get_atom(atom2_idx) {
                Some(a) => a,
                None => continue,
            };

            // Check if atoms have lines representation visible
            if !atom1.visible_reps.is_visible(RepMask::LINES)
                && !atom2.visible_reps.is_visible(RepMask::LINES)
            {
                continue;
            }

            // Get coordinates
            let pos1 = match coord_set.get_atom_coord(atom1_idx) {
                Some(p) => p,
                None => continue,
            };
            let pos2 = match coord_set.get_atom_coord(atom2_idx) {
                Some(p) => p,
                None => continue,
            };

            // Get colors
            let color1 = colors.resolve_atom(atom1, molecule);
            let color2 = colors.resolve_atom(atom2, molecule);

            // Calculate midpoint for half-bond coloring
            let midpoint = [
                (pos1.x + pos2.x) * 0.5,
                (pos1.y + pos2.y) * 0.5,
                (pos1.z + pos2.z) * 0.5,
            ];

            // First half: atom1 to midpoint
            self.vertices.push(LineVertex {
                position: [pos1.x, pos1.y, pos1.z],
                color: color1,
            });
            self.vertices.push(LineVertex {
                position: midpoint,
                color: color1,
            });

            // Second half: midpoint to atom2
            self.vertices.push(LineVertex {
                position: midpoint,
                color: color2,
            });
            self.vertices.push(LineVertex {
                position: [pos2.x, pos2.y, pos2.z],
                color: color2,
            });
        }

        self.vertex_count = self.vertices.len() as u32;
        self.dirty = false;
    }

    fn upload(&mut self, device: &wgpu::Device, queue: &wgpu::Queue) {
        if !self.vertices.is_empty() {
            self.vertex_buffer.write(device, queue, &self.vertices);
        }
    }

    fn render<'a>(&'a self, render_pass: &mut wgpu::RenderPass<'a>) {
        if self.vertex_count == 0 {
            return;
        }

        // The caller is responsible for setting the pipeline
        if let Some(buffer) = self.vertex_buffer.buffer() {
            render_pass.set_vertex_buffer(0, buffer.slice(..));
            render_pass.draw(0..self.vertex_count, 0..1);
        }
    }

    fn is_dirty(&self) -> bool {
        self.dirty
    }

    fn set_dirty(&mut self) {
        self.dirty = true;
    }

    fn primitive_count(&self) -> usize {
        self.vertex_count as usize / 2 // Lines have 2 vertices each
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_line_rep_new() {
        let rep = LineRep::new();
        assert!(rep.is_dirty());
        assert!(rep.is_empty());
    }
}
