//! Wire surface representation for mesh visualization
//!
//! Generates a molecular surface (same as surface representation) but renders
//! it as wireframe lines instead of filled triangles, matching PyMOL's mesh
//! representation behavior.

use std::collections::HashSet;

use crate::buffer::GrowableBuffer;
use crate::color_resolver::ColorResolver;
use crate::representation::Representation;
use crate::vertex::LineVertex;

use super::surface::coloring::{self, AtomColor};
use super::surface::distance_field::{self, SurfaceAtom, SurfaceType};
use super::surface::grid::Grid3D;
use super::surface::marching_cubes;

use pymol_mol::{CoordSet, ObjectMolecule, RepMask};
use pymol_settings::SettingResolver;

/// Wire surface representation (mesh wireframe)
///
/// Generates a molecular surface and renders it as wireframe lines.
/// This matches PyMOL's "mesh" representation which shows the surface
/// as a wireframe grid rather than filled triangles.
pub struct WireSurfaceRep {
    /// Vertex data (CPU side) - pairs of LineVertex for each edge
    vertices: Vec<LineVertex>,
    /// GPU vertex buffer
    vertex_buffer: GrowableBuffer,
    /// Whether the representation needs to be rebuilt
    dirty: bool,
    /// Number of vertices to render
    vertex_count: u32,
    /// Surface type
    surface_type: SurfaceType,
    /// Solvent probe radius (Angstroms)
    probe_radius: f32,
    /// Surface quality (-4 to 4)
    quality: i32,
}

impl WireSurfaceRep {
    /// Create a new wire surface representation
    pub fn new() -> Self {
        Self {
            vertices: Vec::new(),
            vertex_buffer: GrowableBuffer::new(
                "Wire Surface Vertices",
                wgpu::BufferUsages::VERTEX,
            ),
            dirty: true,
            vertex_count: 0,
            surface_type: SurfaceType::default(),
            probe_radius: 1.4,
            quality: 0,
        }
    }

    /// Set surface quality (-4 to 4)
    pub fn set_quality(&mut self, quality: i32) {
        let quality = quality.clamp(-4, 4);
        if self.quality != quality {
            self.quality = quality;
            self.dirty = true;
        }
    }

    /// Get grid spacing from quality setting (same as SurfaceRep)
    fn quality_to_spacing(quality: i32) -> f32 {
        match quality {
            -4 => 2.0,
            -3 => 1.5,
            -2 => 1.0,
            -1 => 0.75,
            0 => 0.5,
            1 => 0.35,
            2 => 0.25,
            3 => 0.18,
            _ => 0.125,
        }
    }

    /// Get the vertices (for raytracing)
    pub fn vertices(&self) -> &[LineVertex] {
        &self.vertices
    }

    /// Clear the wire surface
    pub fn clear(&mut self) {
        self.vertices.clear();
        self.vertex_count = 0;
    }

    /// Generate wireframe line vertices from atoms using marching cubes.
    fn generate_wire_vertices(
        atoms: &[SurfaceAtom],
        atom_colors: &[AtomColor],
        surface_type: SurfaceType,
        probe_radius: f32,
        quality: i32,
    ) -> Vec<LineVertex> {
        if atoms.is_empty() {
            return Vec::new();
        }

        let padding = probe_radius + 2.0;
        let (min, max) = distance_field::compute_bounds_with_radii(atoms, padding);

        let box_size = [max[0] - min[0], max[1] - min[1], max[2] - min[2]];

        let mut spacing = Self::quality_to_spacing(quality);

        let max_grid_points: f32 = match quality {
            4 => 1_000_000.0,
            3 => 700_000.0,
            2 => 500_000.0,
            1 => 350_000.0,
            0 => 250_000.0,
            -1 => 175_000.0,
            -2 => 125_000.0,
            -3 => 80_000.0,
            _ => 50_000.0,
        };

        let max_dim = box_size[0].max(box_size[1]).max(box_size[2]);
        let min_spacing_for_limit = max_dim / max_grid_points.cbrt();
        if spacing < min_spacing_for_limit {
            spacing = min_spacing_for_limit;
        }
        spacing = spacing.max(0.5);

        let mut sdf_grid = Grid3D::from_bounds(min, max, spacing, 0.0);

        distance_field::compute_distance_field(
            &mut sdf_grid,
            atoms,
            surface_type,
            probe_radius,
        );

        let result = marching_cubes::extract_isosurface_smooth(&sdf_grid, 0.0);

        if result.positions.is_empty() {
            return Vec::new();
        }

        let result = marching_cubes::weld_vertices(&result, spacing * 0.5);

        let vertex_colors = coloring::color_vertices(&result.positions, atom_colors);

        let mut edges: HashSet<(u32, u32)> = HashSet::new();
        for tri in result.indices.chunks(3) {
            if tri.len() < 3 {
                continue;
            }
            let (a, b, c) = (tri[0], tri[1], tri[2]);
            edges.insert((a.min(b), a.max(b)));
            edges.insert((b.min(c), b.max(c)));
            edges.insert((c.min(a), c.max(a)));
        }

        let mut vertices = Vec::with_capacity(edges.len() * 2);
        for (i0, i1) in edges {
            vertices.push(LineVertex {
                position: result.positions[i0 as usize],
                color: vertex_colors[i0 as usize],
            });
            vertices.push(LineVertex {
                position: result.positions[i1 as usize],
                color: vertex_colors[i1 as usize],
            });
        }

        vertices
    }
}

impl Default for WireSurfaceRep {
    fn default() -> Self {
        Self::new()
    }
}

impl Representation for WireSurfaceRep {
    fn build(
        &mut self,
        molecule: &ObjectMolecule,
        coord_set: &CoordSet,
        colors: &ColorResolver,
        settings: &SettingResolver,
    ) {
        self.clear();

        let surface_type_setting = settings.get_int_if_defined(pymol_settings::id::surface_type).unwrap_or(0);
        let surface_solvent = settings.get_bool_if_defined(pymol_settings::id::surface_solvent).unwrap_or(false);
        self.surface_type = if surface_solvent {
            SurfaceType::SolventAccessible
        } else {
            SurfaceType::from_setting(surface_type_setting)
        };

        let individual_chains = settings.get_bool(pymol_settings::id::surface_individual_chains);
        let mesh_color = settings.get_color(pymol_settings::id::mesh_color);

        // Collect atoms with MESH visibility
        let mut atoms: Vec<SurfaceAtom> = Vec::new();
        let mut atom_colors: Vec<AtomColor> = Vec::new();
        let mut chain_ids: Vec<String> = Vec::new();

        for (atom_idx, coord) in coord_set.iter_with_atoms() {
            let atom = match molecule.get_atom(atom_idx) {
                Some(a) => a,
                None => continue,
            };

            if !atom.repr.visible_reps.is_visible(RepMask::MESH) {
                continue;
            }

            let position = [coord.x, coord.y, coord.z];
            let radius = atom.effective_vdw();
            let color = colors.resolve_rep_color(atom, atom.repr.colors.mesh, mesh_color);

            atoms.push(SurfaceAtom {
                position,
                radius,
                atom_index: atom_idx.into(),
            });

            atom_colors.push(AtomColor {
                position,
                color,
            });

            if individual_chains {
                chain_ids.push(atom.residue.chain.clone());
            }
        }

        if atoms.is_empty() {
            self.dirty = false;
            return;
        }

        if individual_chains {
            use std::collections::BTreeMap;

            let mut groups: BTreeMap<&str, (Vec<SurfaceAtom>, Vec<AtomColor>)> = BTreeMap::new();
            for i in 0..atoms.len() {
                let entry = groups.entry(chain_ids[i].as_str()).or_default();
                entry.0.push(atoms[i].clone());
                entry.1.push(atom_colors[i].clone());
            }

            for (chain_atoms, chain_colors) in groups.values() {
                let verts = Self::generate_wire_vertices(
                    chain_atoms,
                    chain_colors,
                    self.surface_type,
                    self.probe_radius,
                    self.quality,
                );
                self.vertices.extend(verts);
            }
        } else {
            self.vertices = Self::generate_wire_vertices(
                &atoms,
                &atom_colors,
                self.surface_type,
                self.probe_radius,
                self.quality,
            );
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
    fn test_wire_surface_rep_new() {
        let rep = WireSurfaceRep::new();
        assert!(rep.is_dirty());
        assert!(rep.is_empty());
    }

    #[test]
    fn test_quality_to_spacing() {
        assert!((WireSurfaceRep::quality_to_spacing(-4) - 2.0).abs() < 0.001);
        assert!((WireSurfaceRep::quality_to_spacing(0) - 0.5).abs() < 0.001);
        assert!((WireSurfaceRep::quality_to_spacing(4) - 0.125).abs() < 0.001);
    }
}
