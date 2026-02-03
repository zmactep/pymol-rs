//! Line representation for bond visualization
//!
//! Renders bonds as simple lines connecting atom centers.
//! Supports multiple bond visualization (double, triple, aromatic) when valence display is enabled.

use crate::buffer::GrowableBuffer;
use crate::color_resolver::ColorResolver;
use crate::representation::Representation;
use crate::vertex::LineVertex;

use pymol_mol::{AtomIndex, BondOrder, CoordSet, ObjectMolecule, RepMask};
use pymol_settings::SettingResolver;

/// Calculate perpendicular vector for bond offset using neighbor atom
///
/// Returns a normalized vector perpendicular to the bond direction that lies
/// in the plane defined by the bond and a neighboring atom. This ensures
/// double/triple bond lines stay within the molecular plane (e.g., in a ring).
///
/// If no neighbor is available, falls back to using a fixed reference axis.
fn calculate_perpendicular_with_neighbor(
    bond_dir: [f32; 3],
    pos1: [f32; 3],
    pos2: [f32; 3],
    neighbor_pos: Option<[f32; 3]>,
) -> [f32; 3] {
    if let Some(neighbor) = neighbor_pos {
        // Calculate midpoint of the bond
        let mid = [
            (pos1[0] + pos2[0]) * 0.5,
            (pos1[1] + pos2[1]) * 0.5,
            (pos1[2] + pos2[2]) * 0.5,
        ];

        // Vector from midpoint to neighbor (lies in molecular plane)
        let to_neighbor = [
            neighbor[0] - mid[0],
            neighbor[1] - mid[1],
            neighbor[2] - mid[2],
        ];

        // Vector rejection: remove the component of to_neighbor parallel to bond_dir
        // This gives us a vector IN the molecular plane that is perpendicular to bond_dir
        let dot = bond_dir[0] * to_neighbor[0] + bond_dir[1] * to_neighbor[1] + bond_dir[2] * to_neighbor[2];
        let perp = [
            to_neighbor[0] - bond_dir[0] * dot,
            to_neighbor[1] - bond_dir[1] * dot,
            to_neighbor[2] - bond_dir[2] * dot,
        ];
        let len_sq = perp[0] * perp[0] + perp[1] * perp[1] + perp[2] * perp[2];

        if len_sq > 0.0001 {
            return normalize(perp);
        }
    }

    // Fallback: use fixed reference axis
    calculate_perpendicular_fallback(bond_dir)
}

/// Fallback perpendicular calculation using fixed reference axes
fn calculate_perpendicular_fallback(bond_dir: [f32; 3]) -> [f32; 3] {
    // Cross with Y-axis (up vector)
    let up = [0.0_f32, 1.0, 0.0];
    let mut perp = cross(bond_dir, up);

    // If bond is nearly parallel to Y-axis, use X-axis instead
    let len_sq = perp[0] * perp[0] + perp[1] * perp[1] + perp[2] * perp[2];
    if len_sq < 0.0001 {
        let right = [1.0_f32, 0.0, 0.0];
        perp = cross(bond_dir, right);
    }

    normalize(perp)
}

/// Cross product of two 3D vectors
#[inline]
fn cross(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Normalize a 3D vector
#[inline]
fn normalize(v: [f32; 3]) -> [f32; 3] {
    let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if len > 0.0001 {
        [v[0] / len, v[1] / len, v[2] / len]
    } else {
        [0.0, 1.0, 0.0] // Fallback
    }
}

/// Get the offset factors for multiple bond lines based on bond order
///
/// Returns factors that will be multiplied by valence_size to get actual offsets.
fn get_bond_offsets(order: BondOrder) -> &'static [f32] {
    match order {
        BondOrder::Double => &[-0.5, 0.5],
        BondOrder::Triple => &[-1.0, 0.0, 1.0],
        BondOrder::Aromatic => &[-0.5, 0.5],
        _ => &[0.0], // Single, Unknown
    }
}

/// Find the position of a neighbor atom for determining the molecular plane
///
/// For a bond A-B, finds another atom bonded to A (not B) or bonded to B (not A)
/// and returns its position. This neighbor defines the local molecular plane.
fn find_neighbor_position(
    molecule: &ObjectMolecule,
    coord_set: &CoordSet,
    atom1_idx: AtomIndex,
    atom2_idx: AtomIndex,
) -> Option<[f32; 3]> {
    // Try to find a neighbor of atom1 (excluding atom2)
    for neighbor_idx in molecule.bonded_atoms(atom1_idx) {
        if neighbor_idx != atom2_idx {
            if let Some(pos) = coord_set.get_atom_coord(neighbor_idx) {
                return Some([pos.x, pos.y, pos.z]);
            }
        }
    }

    // If no neighbor found for atom1, try atom2's neighbors
    for neighbor_idx in molecule.bonded_atoms(atom2_idx) {
        if neighbor_idx != atom1_idx {
            if let Some(pos) = coord_set.get_atom_coord(neighbor_idx) {
                return Some([pos.x, pos.y, pos.z]);
            }
        }
    }

    None
}

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

        // Get valence display settings
        // Setting ID 64 = valence (bool), ID 135 = valence_size (float)
        let valence_enabled = settings.get_bool_if_defined(64).unwrap_or(true);
        let valence_size = settings.get_float_if_defined(135).unwrap_or(0.06);

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
