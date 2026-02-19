//! Ribbon representation for protein backbone visualization
//!
//! This module provides ribbon rendering of protein backbones as a uniform smooth tube.
//! Unlike the cartoon representation which uses different profiles for helices, sheets,
//! and loops, ribbon renders the entire backbone with the same circular cross-section.
//!
//! # Differences from Cartoon
//!
//! - Uniform circular tube profile for ALL secondary structure types
//! - Same width throughout the entire protein backbone
//! - No arrow visualization on beta sheets
//! - No differentiation between helices, sheets, and loops
//! - Uses ribbon-specific settings (ribbon_radius, ribbon_sampling, etc.)
//!
//! # Example
//!
//! ```ignore
//! use pymol_render::{RibbonRep, RenderContext, ColorResolver};
//! use pymol_mol::ObjectMolecule;
//!
//! let mut ribbon = RibbonRep::new();
//! ribbon.build(&molecule, &coord_set, &color_resolver, &settings);
//! ribbon.upload(&device, &queue);
//! ribbon.render(&mut render_pass);
//! ```

use crate::buffer::GrowableBuffer;
use crate::color_resolver::ColorResolver;
use crate::representation::Representation;
use crate::vertex::MeshVertex;

use pymol_mol::{CoordSet, ObjectMolecule};
use pymol_settings::SettingResolver;

use super::cartoon::backbone::CartoonSmoothSettings;
use super::cartoon::build_cartoon_geometry;
use super::cartoon::geometry::CartoonGeometrySettings;
use super::cartoon::spline::InterpolationSettings;

/// Ribbon representation for protein secondary structure
///
/// Renders protein backbones as smooth ribbons without arrow heads on beta sheets.
/// Uses mesh-based rendering with the standard mesh shader for proper lighting.
pub struct RibbonRep {
    /// Vertex data (CPU side)
    vertices: Vec<MeshVertex>,
    /// Index data (CPU side)
    indices: Vec<u32>,
    /// GPU vertex buffer
    vertex_buffer: GrowableBuffer,
    /// GPU index buffer
    index_buffer: GrowableBuffer,
    /// Render pipeline
    pipeline: Option<wgpu::RenderPipeline>,
    /// Uniform bind group
    bind_group: Option<wgpu::BindGroup>,
    /// Whether the representation needs rebuild
    dirty: bool,
    /// Number of indices to render
    index_count: u32,
}

impl RibbonRep {
    /// Create a new ribbon representation
    pub fn new() -> Self {
        Self {
            vertices: Vec::new(),
            indices: Vec::new(),
            vertex_buffer: GrowableBuffer::new("Ribbon Vertices", wgpu::BufferUsages::VERTEX),
            index_buffer: GrowableBuffer::new("Ribbon Indices", wgpu::BufferUsages::INDEX),
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

    /// Get the vertices (for debugging/testing)
    pub fn vertices(&self) -> &[MeshVertex] {
        &self.vertices
    }

    /// Get the indices (for debugging/testing)
    pub fn indices(&self) -> &[u32] {
        &self.indices
    }
}

impl Default for RibbonRep {
    fn default() -> Self {
        Self::new()
    }
}

impl Representation for RibbonRep {
    fn build(
        &mut self,
        molecule: &ObjectMolecule,
        coord_set: &CoordSet,
        colors: &ColorResolver,
        settings: &SettingResolver,
    ) {
        // Ribbon-specific setting IDs
        const RIBBON_SAMPLING: u16 = 19;
        const RIBBON_POWER: u16 = 17;
        const RIBBON_POWER_B: u16 = 18;
        const RIBBON_THROW: u16 = 121;
        const RIBBON_RADIUS: u16 = 20;

        // Smoothing settings shared with cartoon
        const CARTOON_GAP_CUTOFF: u16 = 750;
        const CARTOON_SMOOTH_CYCLES: u16 = 259;
        const CARTOON_FLAT_CYCLES: u16 = 260;
        const CARTOON_SMOOTH_FIRST: u16 = 257;
        const CARTOON_SMOOTH_LAST: u16 = 258;
        const CARTOON_SMOOTH_LOOPS: u16 = 114;
        const CARTOON_REFINE_NORMALS: u16 = 112;

        let gap_cutoff = settings.get_int_if_defined(CARTOON_GAP_CUTOFF).unwrap_or(10);
        let sampling = settings.get_int_if_defined(RIBBON_SAMPLING).unwrap_or(1);
        let power = settings.get_float_if_defined(RIBBON_POWER).unwrap_or(2.0);
        let power_b = settings.get_float_if_defined(RIBBON_POWER_B).unwrap_or(0.5);
        let throw = settings.get_float_if_defined(RIBBON_THROW).unwrap_or(1.35);
        let ribbon_radius = settings.get_float_if_defined(RIBBON_RADIUS).unwrap_or(0.2);

        let smooth_settings = CartoonSmoothSettings {
            smooth_cycles: settings.get_int_if_defined(CARTOON_SMOOTH_CYCLES).unwrap_or(2) as u32,
            flat_cycles: settings.get_int_if_defined(CARTOON_FLAT_CYCLES).unwrap_or(4) as u32,
            smooth_loops: settings.get_bool_if_defined(CARTOON_SMOOTH_LOOPS).unwrap_or(false),
            smooth_first: settings.get_int_if_defined(CARTOON_SMOOTH_FIRST).unwrap_or(1) as u32,
            smooth_last: settings.get_int_if_defined(CARTOON_SMOOTH_LAST).unwrap_or(1) as u32,
            refine_normals: settings.get_bool_if_defined(CARTOON_REFINE_NORMALS).unwrap_or(true),
        };

        let subdivisions = if sampling < 0 { 10u32 } else { (sampling as u32).max(7) };

        let interp_settings = InterpolationSettings {
            power_a: power,
            power_b,
            throw_factor: throw,
        };

        // Ribbon uses a uniform tube for all secondary structures
        let tube_radius = if ribbon_radius > 0.0 { ribbon_radius } else { 0.3 };
        let geom_settings = CartoonGeometrySettings {
            helix_width: tube_radius,
            helix_height: tube_radius,
            sheet_width: tube_radius,
            sheet_height: tube_radius,
            loop_radius: tube_radius,
            quality: 32,
            round_helices: true,
            fancy_sheets: false,
            fancy_helices: false,
            dumbbell_length: 1.6,
            dumbbell_width: 0.17,
            dumbbell_radius: 0.16,
            arrow_tip_scale: 1.0,
            arrow_length: 0,
            arrow_residues: 0,
            uniform_tube: true,
        };

        let (vertices, indices) = build_cartoon_geometry(
            molecule,
            coord_set,
            colors,
            &smooth_settings,
            &interp_settings,
            &geom_settings,
            subdivisions,
            gap_cutoff,
            pymol_mol::RepMask::RIBBON,
        );

        self.vertices = vertices;
        self.indices = indices;
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
        self.index_count as usize / 3 // triangles
    }

    fn clear_cpu_data(&mut self) {
        self.vertices = Vec::new();
        self.indices = Vec::new();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use lin_alg::f32::Vec3;
    use pymol_mol::{Atom, CoordSet, Element, ObjectMolecule};
    use pymol_settings::GlobalSettings;

    #[test]
    fn test_ribbon_rep_new() {
        let rep = RibbonRep::new();
        assert!(rep.is_dirty());
        assert!(rep.is_empty());
    }

    fn create_test_peptide() -> ObjectMolecule {
        let mut mol = ObjectMolecule::new("peptide");

        let atoms_data = vec![
            ("N", Element::Nitrogen, "ALA", 1, (0.0, 0.0, 0.0)),
            ("CA", Element::Carbon, "ALA", 1, (1.45, 0.0, 0.0)),
            ("C", Element::Carbon, "ALA", 1, (2.0, 1.4, 0.0)),
            ("O", Element::Oxygen, "ALA", 1, (1.4, 2.4, 0.0)),
            ("N", Element::Nitrogen, "GLY", 2, (3.3, 1.5, 0.0)),
            ("CA", Element::Carbon, "GLY", 2, (4.0, 2.8, 0.0)),
            ("C", Element::Carbon, "GLY", 2, (5.5, 2.7, 0.0)),
            ("O", Element::Oxygen, "GLY", 2, (6.1, 1.6, 0.0)),
            ("N", Element::Nitrogen, "ALA", 3, (6.1, 3.9, 0.0)),
            ("CA", Element::Carbon, "ALA", 3, (7.5, 4.1, 0.0)),
            ("C", Element::Carbon, "ALA", 3, (8.0, 5.5, 0.0)),
            ("O", Element::Oxygen, "ALA", 3, (7.2, 6.4, 0.0)),
        ];

        let mut coords = Vec::new();
        for (name, element, resn, resv, coord) in &atoms_data {
            let mut atom = Atom::new(*name, *element);
            atom.set_residue(*resn, *resv, "A");
            atom.repr.visible_reps.set_visible(pymol_mol::RepMask::RIBBON);
            mol.add_atom(atom);
            coords.push(Vec3::new(coord.0, coord.1, coord.2));
        }

        mol.add_coord_set(CoordSet::from_vec3(&coords));
        mol.classify_atoms();
        mol
    }

    #[test]
    fn test_ribbon_build_simple_peptide() {
        use pymol_color::{ChainColors, ElementColors, NamedColors};

        let mol = create_test_peptide();
        let coord_set = mol.get_coord_set(0).unwrap();
        let settings = GlobalSettings::new();
        let settings_resolver = pymol_settings::SettingResolver::global(&settings);

        // Create color tables
        let named_colors = NamedColors::new();
        let element_colors = ElementColors::new();
        let chain_colors = ChainColors;
        let color_resolver = ColorResolver::new(&named_colors, &element_colors, &chain_colors);

        let mut ribbon = RibbonRep::new();
        ribbon.build(&mol, coord_set, &color_resolver, &settings_resolver);

        assert!(!ribbon.vertices.is_empty(), "Ribbon should have vertices");
        assert!(!ribbon.indices.is_empty(), "Ribbon should have indices");
        assert!(!ribbon.is_empty(), "Ribbon should not be empty");
    }

    #[test]
    fn test_ribbon_with_secondary_structure() {
        use pymol_color::{ChainColors, ElementColors, NamedColors};

        // Create a peptide with secondary structure assigned
        let mut mol = ObjectMolecule::new("ss_peptide");

        // Create 6 residues - 2 helix, 2 sheet, 2 loop
        let ss_sequence = [
            (pymol_mol::SecondaryStructure::Helix, "ALA", 1),
            (pymol_mol::SecondaryStructure::Helix, "GLY", 2),
            (pymol_mol::SecondaryStructure::Sheet, "VAL", 3),
            (pymol_mol::SecondaryStructure::Sheet, "LEU", 4),
            (pymol_mol::SecondaryStructure::Loop, "PRO", 5),
            (pymol_mol::SecondaryStructure::Loop, "SER", 6),
        ];

        let mut coords = Vec::new();
        let base_x = 0.0f32;

        for (i, (ss, resn, resv)) in ss_sequence.iter().enumerate() {
            let x = base_x + (i as f32) * 3.8; // ~CA-CA distance

            // N atom
            let mut n = Atom::new("N", Element::Nitrogen);
            n.set_residue(*resn, *resv, "A");
            n.ss_type = *ss;
            n.repr.visible_reps.set_visible(pymol_mol::RepMask::RIBBON);
            mol.add_atom(n);
            coords.push(Vec3::new(x - 1.0, 0.0, 0.0));

            // CA atom
            let mut ca = Atom::new("CA", Element::Carbon);
            ca.set_residue(*resn, *resv, "A");
            ca.ss_type = *ss;
            ca.repr.visible_reps.set_visible(pymol_mol::RepMask::RIBBON);
            mol.add_atom(ca);
            coords.push(Vec3::new(x, 0.0, 0.0));

            // C atom
            let mut c = Atom::new("C", Element::Carbon);
            c.set_residue(*resn, *resv, "A");
            c.ss_type = *ss;
            c.repr.visible_reps.set_visible(pymol_mol::RepMask::RIBBON);
            mol.add_atom(c);
            coords.push(Vec3::new(x + 0.5, 0.5, 0.0));

            // O atom
            let mut o = Atom::new("O", Element::Oxygen);
            o.set_residue(*resn, *resv, "A");
            o.ss_type = *ss;
            o.repr.visible_reps.set_visible(pymol_mol::RepMask::RIBBON);
            mol.add_atom(o);
            coords.push(Vec3::new(x + 0.5, 1.5, 0.0));
        }

        mol.add_coord_set(CoordSet::from_vec3(&coords));
        mol.classify_atoms();

        // Verify protein residues were classified
        let protein_count = mol.residues().filter(|r| r.is_protein()).count();
        assert_eq!(protein_count, 6, "Should have 6 protein residues");

        let coord_set = mol.get_coord_set(0).unwrap();
        let settings = GlobalSettings::new();
        let settings_resolver = pymol_settings::SettingResolver::global(&settings);

        let named_colors = NamedColors::new();
        let element_colors = ElementColors::new();
        let chain_colors = ChainColors;
        let color_resolver = ColorResolver::new(&named_colors, &element_colors, &chain_colors);

        // Build ribbon
        let mut ribbon = RibbonRep::new();
        ribbon.build(&mol, coord_set, &color_resolver, &settings_resolver);

        assert!(!ribbon.vertices.is_empty(), "Ribbon should have vertices");
        assert!(!ribbon.indices.is_empty(), "Ribbon should have indices");

        // Print some diagnostic info
        eprintln!(
            "Ribbon vertices: {}, indices: {}",
            ribbon.vertices.len(),
            ribbon.indices.len()
        );
    }
}
