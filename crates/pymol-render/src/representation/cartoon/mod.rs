//! Cartoon representation for protein secondary structure visualization
//!
//! This module provides cartoon rendering of protein backbones, showing:
//! - Alpha helices as ribbons or cylinders
//! - Beta sheets as flat arrows
//! - Loops and coils as tubes
//!
//! # Architecture
//!
//! The cartoon representation is built in several stages:
//! 1. **Backbone extraction**: Extract CA and O atoms from protein residues
//! 2. **PyMOL-style smoothing**: Apply helix centering, sheet flattening, loop smoothing
//! 3. **Spline interpolation**: Smooth the backbone with Catmull-Rom splines
//! 4. **Frame generation**: Calculate reference frames along the spline
//! 5. **Geometry extrusion**: Extrude cross-section profiles along frames
//!
//! # PyMOL Compatibility
//!
//! This implementation includes PyMOL-compatible smoothing algorithms:
//! - Helix axis centering (ExtrudeShiftToAxis)
//! - Sheet flattening (RepCartoonFlattenSheets)
//! - Loop smoothing (RepCartoonSmoothLoops)
//! - Normal refinement (RepCartoonRefineNormals)
//!
//! # Example
//!
//! ```ignore
//! use pymol_render::{CartoonRep, RenderContext, ColorResolver};
//! use pymol_mol::ObjectMolecule;
//!
//! let mut cartoon = CartoonRep::new();
//! cartoon.build(&molecule, &coord_set, &color_resolver, &settings);
//! cartoon.upload(&device, &queue);
//! cartoon.render(&mut render_pass);
//! ```

pub mod backbone;
pub mod frame;
pub mod geometry;
pub mod spline;

use crate::buffer::GrowableBuffer;
use crate::color_resolver::ColorResolver;
use crate::representation::Representation;
use crate::vertex::MeshVertex;

use pymol_mol::{CoordSet, ObjectMolecule};
use pymol_settings::SettingResolver;

use self::backbone::{
    apply_pymol_smoothing, extract_backbone_segments, smooth_orientations, CartoonSmoothSettings,
};
use self::frame::generate_frames;
use self::geometry::{find_sheet_termini, generate_cartoon_mesh, CartoonGeometrySettings};
use self::spline::CatmullRomSpline;

/// Cartoon representation for protein secondary structure
///
/// Renders protein backbones as ribbons (helices), arrows (sheets), and tubes (loops).
/// Uses mesh-based rendering with the standard mesh shader for proper lighting.
pub struct CartoonRep {
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

impl CartoonRep {
    /// Create a new cartoon representation
    pub fn new() -> Self {
        Self {
            vertices: Vec::new(),
            indices: Vec::new(),
            vertex_buffer: GrowableBuffer::new("Cartoon Vertices", wgpu::BufferUsages::VERTEX),
            index_buffer: GrowableBuffer::new("Cartoon Indices", wgpu::BufferUsages::INDEX),
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

impl Default for CartoonRep {
    fn default() -> Self {
        Self::new()
    }
}

impl Representation for CartoonRep {
    fn build(
        &mut self,
        molecule: &ObjectMolecule,
        coord_set: &CoordSet,
        colors: &ColorResolver,
        settings: &SettingResolver,
    ) {
        self.vertices.clear();
        self.indices.clear();

        // Get settings - PyMOL setting IDs (from pymol-settings/src/definitions.rs)
        const CARTOON_GAP_CUTOFF: u16 = 750;
        const CARTOON_SAMPLING: u16 = 91;
        const CARTOON_POWER: u16 = 94;
        const CARTOON_POWER_B: u16 = 95;
        const CARTOON_THROW: u16 = 122;
        const CARTOON_SMOOTH_CYCLES: u16 = 259;
        const CARTOON_FLAT_CYCLES: u16 = 260;
        const CARTOON_SMOOTH_FIRST: u16 = 257;
        const CARTOON_SMOOTH_LAST: u16 = 258;
        const CARTOON_REFINE: u16 = 123;
        const CARTOON_REFINE_NORMALS: u16 = 112;

        let gap_cutoff = settings.get_int_if_defined(CARTOON_GAP_CUTOFF).unwrap_or(10);
        let sampling = settings.get_int_if_defined(CARTOON_SAMPLING).unwrap_or(7);
        let power = settings.get_float_if_defined(CARTOON_POWER).unwrap_or(2.0);

        // PyMOL-style smoothing settings (defaults match PyMOL exactly)
        let power_b = settings.get_float_if_defined(CARTOON_POWER_B).unwrap_or(0.52);
        let throw = settings.get_float_if_defined(CARTOON_THROW).unwrap_or(1.35);
        
        let smooth_settings = CartoonSmoothSettings {
            smooth_cycles: settings
                .get_int_if_defined(CARTOON_SMOOTH_CYCLES)
                .unwrap_or(2) as u32, // PyMOL default: 2
            flat_cycles: settings
                .get_int_if_defined(CARTOON_FLAT_CYCLES)
                .unwrap_or(4) as u32, // PyMOL default: 4
            smooth_first: settings
                .get_int_if_defined(CARTOON_SMOOTH_FIRST)
                .unwrap_or(1) as u32, // PyMOL default: 1
            smooth_last: settings
                .get_int_if_defined(CARTOON_SMOOTH_LAST)
                .unwrap_or(1) as u32, // PyMOL default: 1
            refine_cycles: settings
                .get_int_if_defined(CARTOON_REFINE)
                .unwrap_or(5) as u32, // PyMOL default: 5
            refine_normals: settings
                .get_bool_if_defined(CARTOON_REFINE_NORMALS)
                .unwrap_or(true),
            power_a: power,
            power_b,
            throw_factor: throw,
        };

        // Calculate subdivisions from sampling
        // Higher subdivisions = smoother curves, especially for helix turns
        let subdivisions = if sampling < 0 {
            // Auto: use higher value for smoother helix turns
            14u32
        } else {
            // Ensure minimum of 14 for smooth helix turns
            (sampling as u32).max(14)
        };

        // Create spline with low tension for smooth curves
        // Lower tension = smoother, more rounded curves
        // Standard Catmull-Rom is 0.5, but 0.0 gives smoother results for ribbons
        let tension = 0.0;
        let spline = CatmullRomSpline::with_tension(tension);

        // Get geometry settings with arrow length scaled by subdivisions
        let geom_settings =
            CartoonGeometrySettings::from_resolver(settings).with_subdivisions(subdivisions);

        // Extract backbone segments (using CARTOON visibility mask)
        let mut segments = extract_backbone_segments(molecule, coord_set, colors, gap_cutoff, pymol_mol::RepMask::CARTOON);

        // Process each segment
        for segment in &mut segments {
            if segment.len() < 2 {
                continue;
            }

            // Apply PyMOL-style smoothing operations:
            // 1. Helix axis centering
            // 2. Sheet flattening
            // 3. Loop smoothing
            // 4. Normal refinement
            apply_pymol_smoothing(segment, &smooth_settings);

            // Additional orientation smoothing (as before)
            smooth_orientations(segment, smooth_settings.smooth_cycles);

            // Generate reference frames with increased smoothing for better visual quality
            // Higher cycles = smoother ribbon appearance
            let frame_smooth_cycles = smooth_settings.smooth_cycles.max(8);
            let frames = generate_frames(
                &segment.guide_points,
                &spline,
                subdivisions,
                frame_smooth_cycles,
            );

            if frames.is_empty() {
                continue;
            }

            // Find sheet termini for arrow heads (only if fancy_sheets is enabled)
            let sheet_termini = if geom_settings.fancy_sheets {
                find_sheet_termini(&frames)
            } else {
                Vec::new()
            };

            // Generate mesh geometry
            let (mut seg_vertices, seg_indices) =
                generate_cartoon_mesh(&frames, &geom_settings, &sheet_termini);

            // Offset indices for this segment
            let base_index = self.vertices.len() as u32;
            let offset_indices: Vec<u32> = seg_indices.iter().map(|i| i + base_index).collect();

            // Append to global buffers
            self.vertices.append(&mut seg_vertices);
            self.indices.extend(offset_indices);
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
        self.index_count as usize / 3 // triangles
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::backbone::{extract_backbone_segments, GuidePoint};
    use super::frame::generate_frames;
    use super::spline::CatmullRomSpline;
    use lin_alg::f32::Vec3;
    use pymol_mol::{Atom, CoordSet, Element, ObjectMolecule};
    use pymol_settings::GlobalSettings;

    #[test]
    fn test_cartoon_rep_new() {
        let rep = CartoonRep::new();
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
            atom.resn = resn.to_string();
            atom.resv = *resv;
            atom.chain = "A".to_string();
            atom.visible_reps.set_visible(pymol_mol::RepMask::CARTOON);
            mol.add_atom(atom);
            coords.push(Vec3::new(coord.0, coord.1, coord.2));
        }

        mol.add_coord_set(CoordSet::from_vec3(&coords));
        mol.classify_atoms();
        mol
    }

    #[test]
    fn test_backbone_extraction() {
        use pymol_color::{ChainColors, ElementColors, NamedColors};
        
        let mol = create_test_peptide();
        let coord_set = mol.get_coord_set(0).unwrap();
        
        // Create color tables
        let named_colors = NamedColors::new();
        let element_colors = ElementColors::new();
        let chain_colors = ChainColors;
        let color_resolver = ColorResolver::new(&named_colors, &element_colors, &chain_colors);

        let segments = extract_backbone_segments(&mol, coord_set, &color_resolver, 10, pymol_mol::RepMask::CARTOON);
        
        assert!(!segments.is_empty(), "Should have at least one segment");
        assert_eq!(segments[0].len(), 3, "Should have 3 guide points (3 residues)");
    }

    #[test]
    fn test_frame_generation() {
        let guide_points = vec![
            GuidePoint::new(
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
                [1.0, 0.0, 0.0, 1.0],
                pymol_mol::SecondaryStructure::Loop,
                pymol_mol::AtomIndex(0),
                1,
                20.0,
            ),
            GuidePoint::new(
                Vec3::new(3.8, 0.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
                [0.0, 1.0, 0.0, 1.0],
                pymol_mol::SecondaryStructure::Loop,
                pymol_mol::AtomIndex(1),
                2,
                20.0,
            ),
            GuidePoint::new(
                Vec3::new(7.6, 0.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
                [0.0, 0.0, 1.0, 1.0],
                pymol_mol::SecondaryStructure::Loop,
                pymol_mol::AtomIndex(2),
                3,
                20.0,
            ),
        ];
        
        let spline = CatmullRomSpline::new();
        let frames = generate_frames(&guide_points, &spline, 7, 2);
        
        assert!(!frames.is_empty(), "Should have frames");
        assert_eq!(frames.len(), 17, "Should have 17 interpolated frames");
    }

    #[test]
    fn test_cartoon_build_simple_peptide() {
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

        let mut cartoon = CartoonRep::new();
        cartoon.build(&mol, coord_set, &color_resolver, &settings_resolver);

        assert!(!cartoon.vertices.is_empty(), "Cartoon should have vertices");
        assert!(!cartoon.indices.is_empty(), "Cartoon should have indices");
        assert!(!cartoon.is_empty(), "Cartoon should not be empty");
    }

    #[test]
    fn test_cartoon_with_secondary_structure() {
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
            n.resn = resn.to_string();
            n.resv = *resv;
            n.chain = "A".to_string();
            n.ss_type = *ss;
            n.visible_reps.set_visible(pymol_mol::RepMask::CARTOON);
            mol.add_atom(n);
            coords.push(Vec3::new(x - 1.0, 0.0, 0.0));
            
            // CA atom
            let mut ca = Atom::new("CA", Element::Carbon);
            ca.resn = resn.to_string();
            ca.resv = *resv;
            ca.chain = "A".to_string();
            ca.ss_type = *ss;
            ca.visible_reps.set_visible(pymol_mol::RepMask::CARTOON);
            mol.add_atom(ca);
            coords.push(Vec3::new(x, 0.0, 0.0));
            
            // C atom
            let mut c = Atom::new("C", Element::Carbon);
            c.resn = resn.to_string();
            c.resv = *resv;
            c.chain = "A".to_string();
            c.ss_type = *ss;
            c.visible_reps.set_visible(pymol_mol::RepMask::CARTOON);
            mol.add_atom(c);
            coords.push(Vec3::new(x + 0.5, 0.5, 0.0));
            
            // O atom
            let mut o = Atom::new("O", Element::Oxygen);
            o.resn = resn.to_string();
            o.resv = *resv;
            o.chain = "A".to_string();
            o.ss_type = *ss;
            o.visible_reps.set_visible(pymol_mol::RepMask::CARTOON);
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

        // Extract backbone segments
        let segments = extract_backbone_segments(&mol, coord_set, &color_resolver, 10, pymol_mol::RepMask::CARTOON);
        assert!(!segments.is_empty(), "Should have backbone segments");
        assert_eq!(segments[0].len(), 6, "Should have 6 guide points");

        // Build cartoon
        let mut cartoon = CartoonRep::new();
        cartoon.build(&mol, coord_set, &color_resolver, &settings_resolver);

        assert!(!cartoon.vertices.is_empty(), "Cartoon should have vertices");
        assert!(!cartoon.indices.is_empty(), "Cartoon should have indices");
        
        // Print some diagnostic info
        eprintln!("Cartoon vertices: {}, indices: {}", cartoon.vertices.len(), cartoon.indices.len());
    }

    #[test]
    fn test_cartoon_from_pdb_file() {
        use pymol_color::{ChainColors, ElementColors, NamedColors};
        use std::path::Path;
        
        // Try to load a real PDB file
        let path = Path::new("../../pymol-open-source/test/dat/1tii.pdb");
        if !path.exists() {
            eprintln!("Skipping test - PDB file not found at {:?}", path);
            return;
        }
        
        let mut mol = match pymol_io::read_file(path) {
            Ok(m) => m,
            Err(e) => {
                eprintln!("Failed to load PDB: {}", e);
                return;
            }
        };
        
        eprintln!("Loaded {} atoms from {:?}", mol.atom_count(), path);
        
        // Count atoms by secondary structure
        let mut helix_count = 0;
        let mut sheet_count = 0;
        let mut loop_count = 0;
        for atom in mol.atoms() {
            match atom.ss_type {
                pymol_mol::SecondaryStructure::Helix |
                pymol_mol::SecondaryStructure::Helix310 |
                pymol_mol::SecondaryStructure::HelixPi => helix_count += 1,
                pymol_mol::SecondaryStructure::Sheet => sheet_count += 1,
                _ => loop_count += 1,
            }
        }
        eprintln!("Secondary structure: Helix={}, Sheet={}, Loop={}", helix_count, sheet_count, loop_count);
        
        // Count protein residues
        let protein_count = mol.residues().filter(|r| r.is_protein()).count();
        eprintln!("Protein residues: {}", protein_count);
        
        // Check if atoms have CARTOON visibility - by default they don't
        let cartoon_visible_before = mol.atoms().filter(|a| a.visible_reps.is_visible(pymol_mol::RepMask::CARTOON)).count();
        eprintln!("Atoms with CARTOON visible (before): {}", cartoon_visible_before);
        
        // Set CARTOON visibility on all atoms (simulating what toggle() does)
        for atom in mol.atoms_mut() {
            atom.visible_reps.set_visible(pymol_mol::RepMask::CARTOON);
        }
        
        let cartoon_visible_after = mol.atoms().filter(|a| a.visible_reps.is_visible(pymol_mol::RepMask::CARTOON)).count();
        eprintln!("Atoms with CARTOON visible (after): {}", cartoon_visible_after);
        
        // Get coordinate set
        let coord_set = match mol.get_coord_set(0) {
            Some(cs) => cs,
            None => {
                eprintln!("No coordinate set found");
                return;
            }
        };
        
        let settings = GlobalSettings::new();
        let settings_resolver = pymol_settings::SettingResolver::global(&settings);
        
        let named_colors = NamedColors::new();
        let element_colors = ElementColors::new();
        let chain_colors = ChainColors;
        let color_resolver = ColorResolver::new(&named_colors, &element_colors, &chain_colors);

        // Extract backbone segments
        let segments = extract_backbone_segments(&mol, coord_set, &color_resolver, 10, pymol_mol::RepMask::CARTOON);
        eprintln!("Backbone segments: {}", segments.len());
        for (i, seg) in segments.iter().enumerate() {
            eprintln!("  Segment {}: {} guide points, chain={}", i, seg.len(), seg.chain_id);
        }
        
        // Build cartoon
        let mut cartoon = CartoonRep::new();
        cartoon.build(&mol, coord_set, &color_resolver, &settings_resolver);
        
        eprintln!("Cartoon vertices: {}, indices: {}", cartoon.vertices.len(), cartoon.indices.len());
        
        // The cartoon should have vertices for a protein file
        assert!(!cartoon.vertices.is_empty(), "Cartoon should have vertices for a protein PDB file");
    }
}
