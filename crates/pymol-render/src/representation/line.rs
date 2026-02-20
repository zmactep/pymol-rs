//! Line representation for bond visualization
//!
//! Renders bonds as simple lines connecting atom centers.
//! Supports multiple bond visualization (double, triple, aromatic) when valence display is enabled.

use crate::buffer::GrowableBuffer;
use crate::color_resolver::ColorResolver;
use crate::representation::Representation;
use crate::vertex::LineVertex;

use pymol_mol::{CoordSet, ObjectMolecule, RepMask};
use pymol_settings::SettingResolver;

use super::bond_utils::{
    calculate_perpendicular_with_neighbor, find_neighbor_position, get_bond_offsets, normalize,
};

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
        settings: &SettingResolver,
    ) {
        self.vertices.clear();

        let valence_enabled = settings.get_bool_if_defined(pymol_settings::id::valence).unwrap_or(true);
        let valence_size = settings.get_float_if_defined(pymol_settings::id::valence_size).unwrap_or(0.06);
        let line_color = settings.get_color(pymol_settings::id::line_color);

        // Scale valence_size for visible line separation.
        // The raw valence_size (default 0.06 Å) is too small for lines.
        // Similar to sticks which use stick_radius * 2.5, we scale up to get
        // a visible separation (~0.15 Å for double bonds with default settings).
        let line_valence_offset = valence_size * 2.5;

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

            // Skip bond if either atom has lines hidden
            if !atom1.repr.visible_reps.is_visible(RepMask::LINES)
                || !atom2.repr.visible_reps.is_visible(RepMask::LINES)
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

            // Get colors (3-level fallback: per-atom → settings → base)
            let color1 = colors.resolve_rep_color(atom1, atom1.repr.colors.line, line_color);
            let color2 = colors.resolve_rep_color(atom2, atom2.repr.colors.line, line_color);

            // Get offsets for multiple bonds (or single offset of 0.0 for single bonds)
            let offsets = if valence_enabled {
                get_bond_offsets(bond.order)
            } else {
                &[0.0_f32] as &[f32]
            };

            // Calculate perpendicular direction for offsetting multiple bonds
            let bond_dir = normalize([
                pos2.x - pos1.x,
                pos2.y - pos1.y,
                pos2.z - pos1.z,
            ]);

            // Find a neighbor atom to determine the molecular plane
            // This ensures double bonds stay in-plane (e.g., in aromatic rings)
            let neighbor_pos = find_neighbor_position(
                molecule,
                coord_set,
                atom1_idx,
                atom2_idx,
            );
            let perp = calculate_perpendicular_with_neighbor(
                bond_dir,
                [pos1.x, pos1.y, pos1.z],
                [pos2.x, pos2.y, pos2.z],
                neighbor_pos,
            );

            // Generate line(s) for this bond
            for &offset_factor in offsets {
                let offset = offset_factor * line_valence_offset;
                let off = [perp[0] * offset, perp[1] * offset, perp[2] * offset];

                // Offset positions
                let p1 = [pos1.x + off[0], pos1.y + off[1], pos1.z + off[2]];
                let p2 = [pos2.x + off[0], pos2.y + off[1], pos2.z + off[2]];

                // Calculate midpoint for half-bond coloring
                let midpoint = [
                    (p1[0] + p2[0]) * 0.5,
                    (p1[1] + p2[1]) * 0.5,
                    (p1[2] + p2[2]) * 0.5,
                ];

                // First half: atom1 to midpoint
                self.vertices.push(LineVertex {
                    position: p1,
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
                    position: p2,
                    color: color2,
                });
            }
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

    fn clear_cpu_data(&mut self) {
        self.vertices = Vec::new();
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
