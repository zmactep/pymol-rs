//! Surface representation for molecular surfaces
//!
//! Uses a marching cubes algorithm for generating molecular surfaces:
//!
//! - Computes distance field on a 3D grid
//! - Extracts isosurface at distance = 0
//! - O(n³) complexity where n is grid size

pub mod coloring;
pub mod distance_field;
pub mod grid;
pub mod marching_cubes;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

use crate::buffer::GrowableBuffer;
use crate::color_resolver::ColorResolver;
use crate::representation::Representation;
use crate::vertex::MeshVertex;

use pymol_mol::{CoordSet, ObjectMolecule, RepMask};
use pymol_settings::ResolvedSettings;

pub use coloring::AtomColor;
pub use distance_field::{SurfaceAtom, SurfaceType};
pub use grid::Grid3D;

/// Parameters used to compute the cached SDF grid.
/// If these match, the cached grid can be reused.
#[derive(PartialEq, Clone)]
struct SdfParams {
    probe_radius: f32,
    quality: i32,
    surface_type: SurfaceType,
    atom_count: usize,
}

/// Surface representation for molecular surfaces
///
/// Generates a triangulated molecular surface from atomic positions
/// using a marching cubes (volumetric) algorithm.
pub struct SurfaceRep {
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
    /// Surface type
    surface_type: SurfaceType,
    /// Solvent probe radius (Angstroms)
    probe_radius: f32,
    /// Surface quality (-4 to 4)
    quality: i32,
    /// Whether any vertices have transparency (alpha < 1.0)
    has_transparency: bool,
    /// Cached SDF grid (persists across rebuilds to avoid recomputation)
    cached_sdf: Option<Grid3D>,
    /// Parameters that were used when cached_sdf was computed
    cached_sdf_params: Option<SdfParams>,
}

impl SurfaceRep {
    /// Create a new surface representation
    pub fn new() -> Self {
        Self {
            vertices: Vec::new(),
            indices: Vec::new(),
            vertex_buffer: GrowableBuffer::new("Surface Vertices", wgpu::BufferUsages::VERTEX),
            index_buffer: GrowableBuffer::new("Surface Indices", wgpu::BufferUsages::INDEX),
            pipeline: None,
            bind_group: None,
            dirty: true,
            index_count: 0,
            surface_type: SurfaceType::default(),
            probe_radius: 1.4,
            quality: 0,
            has_transparency: false,
            cached_sdf: None,
            cached_sdf_params: None,
        }
    }

    /// Check if any vertices have transparency (alpha < 1.0)
    ///
    /// When true, the surface should be rendered with a transparent blend mode.
    pub fn has_transparency(&self) -> bool {
        self.has_transparency
    }

    /// Set the render pipeline and bind group
    pub fn set_pipeline(&mut self, pipeline: wgpu::RenderPipeline, bind_group: wgpu::BindGroup) {
        self.pipeline = Some(pipeline);
        self.bind_group = Some(bind_group);
    }

    /// Get the vertices (for raytracing)
    pub fn vertices(&self) -> &[MeshVertex] {
        &self.vertices
    }

    /// Get the indices (for raytracing)
    pub fn indices(&self) -> &[u32] {
        &self.indices
    }

    /// Set surface type
    pub fn set_surface_type(&mut self, surface_type: SurfaceType) {
        if self.surface_type != surface_type {
            self.surface_type = surface_type;
            self.dirty = true;
        }
    }

    /// Set probe radius
    pub fn set_probe_radius(&mut self, radius: f32) {
        if (self.probe_radius - radius).abs() > 0.001 {
            self.probe_radius = radius;
            self.dirty = true;
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

    /// Get grid spacing from quality setting
    fn quality_to_spacing(quality: i32) -> f32 {
        // Quality -4 to 4 maps to spacing 2.0 to 0.125
        // Default (0) = 0.5 Angstroms
        match quality {
            -4 => 2.0,
            -3 => 1.5,
            -2 => 1.0,
            -1 => 0.75,
            0 => 0.5,
            1 => 0.35,
            2 => 0.25,
            3 => 0.18,
            _ => 0.125, // 4 or higher
        }
    }

    /// Get adaptive grid spacing based on molecule size
    /// 
    /// Uses the quality setting directly without size-based offsets.
    /// The max_grid_points limit and min_spacing protect against excessive computation.
    /// 
    /// Grid point estimates for 70Å molecule:
    /// - 2.0Å spacing: ~35³ = 43K points (fast)
    /// - 1.5Å spacing: ~47³ = 100K points (acceptable)
    /// - 1.0Å spacing: ~70³ = 343K points (slow)
    /// - 0.5Å spacing: ~140³ = 2.7M points (very slow)
    fn adaptive_spacing(_n_atoms: usize, quality: i32) -> f32 {
        // Use quality directly - the max_grid_points limit protects large molecules
        Self::quality_to_spacing(quality)
    }

    /// Clear the surface mesh
    pub fn clear(&mut self) {
        self.vertices.clear();
        self.indices.clear();
        self.index_count = 0;
        self.has_transparency = false;
    }

    /// Build separate surfaces per chain, then merge into one mesh.
    fn build_per_chain(
        &mut self,
        atoms: &[SurfaceAtom],
        atom_colors: &[AtomColor],
        chain_ids: &[String],
    ) {
        use std::collections::BTreeMap;

        self.clear();

        if atoms.is_empty() {
            self.dirty = false;
            return;
        }

        // Group atoms by chain ID (BTreeMap for deterministic order)
        let mut groups: BTreeMap<&str, (Vec<SurfaceAtom>, Vec<AtomColor>)> = BTreeMap::new();
        for i in 0..atoms.len() {
            let entry = groups.entry(chain_ids[i].as_str()).or_default();
            entry.0.push(atoms[i].clone());
            entry.1.push(atom_colors[i].clone());
        }

        for (chain_atoms, chain_colors) in groups.values() {
            let index_offset = self.vertices.len() as u32;

            // Build this chain's surface in a temporary rep
            let mut temp = SurfaceRep::new();
            temp.surface_type = self.surface_type;
            temp.probe_radius = self.probe_radius;
            temp.quality = self.quality;
            temp.build_from_atoms(chain_atoms, chain_colors);

            // Merge into our buffers with index offset
            self.vertices.extend_from_slice(&temp.vertices);
            self.indices.extend(temp.indices.iter().map(|&i| i + index_offset));
        }

        self.index_count = self.indices.len() as u32;
        self.dirty = false;
    }

    /// Build surface from pre-computed atoms
    pub fn build_from_atoms(
        &mut self,
        atoms: &[SurfaceAtom],
        atom_colors: &[AtomColor],
    ) {
        self.clear();

        if atoms.is_empty() {
            self.dirty = false;
            return;
        }

        self.build_marching_cubes(atoms, atom_colors);
    }

    /// Recolor existing vertices without rebuilding geometry.
    /// Returns false if vertices are empty (full rebuild needed).
    pub fn recolor(&mut self, atom_colors: &[AtomColor], transparency: f32) -> bool {
        if self.vertices.is_empty() {
            return false;
        }
        let alpha = 1.0 - transparency.clamp(0.0, 1.0);
        let positions: Vec<[f32; 3]> = self.vertices.iter().map(|v| v.position).collect();
        let colors = coloring::color_vertices(&positions, atom_colors);
        #[cfg(feature = "parallel")]
        let iter = self.vertices.par_iter_mut().zip(colors.par_iter());
        #[cfg(not(feature = "parallel"))]
        let iter = self.vertices.iter_mut().zip(colors.iter());
        iter.for_each(|(v, c)| {
            v.color = [c[0], c[1], c[2], alpha];
        });
        self.has_transparency = transparency > 0.0;
        true
    }

    /// Build surface using marching cubes algorithm (volumetric)
    fn build_marching_cubes(
        &mut self,
        atoms: &[SurfaceAtom],
        atom_colors: &[AtomColor],
    ) {
        // Compute bounding box with padding for probe
        let padding = self.probe_radius + 2.0; // Extra padding for smooth surface
        let (min, max) = distance_field::compute_bounds_with_radii(atoms, padding);

        // Calculate bounding box dimensions
        let box_size = [
            max[0] - min[0],
            max[1] - min[1],
            max[2] - min[2],
        ];
        
        // Create grid with adaptive spacing based on molecule size
        // Balance between quality and performance
        let mut spacing = Self::adaptive_spacing(atoms.len(), self.quality);
        
        // Limit grid size based on quality setting
        // Higher quality settings allow more grid points (slower but finer)
        let max_grid_points = match self.quality {
            4 => 1_000_000.0_f32,   // Very fine: 1M points (100³)
            3 => 700_000.0_f32,     // Fine: 700K points (~89³)
            2 => 500_000.0_f32,     // High: 500K points (~79³)
            1 => 350_000.0_f32,     // Medium-high: 350K points (~70³)
            0 => 250_000.0_f32,     // Default: 250K points (~63³)
            -1 => 175_000.0_f32,    // Medium-low: 175K points (~56³)
            -2 => 125_000.0_f32,    // Low: 125K points (~50³)
            -3 => 80_000.0_f32,     // Very low: 80K points (~43³)
            _ => 50_000.0_f32,      // Coarse (-4): 50K points (~37³)
        };
        
        let max_dim = box_size[0].max(box_size[1]).max(box_size[2]);
        let min_spacing_for_limit = max_dim / max_grid_points.cbrt();
        
        if spacing < min_spacing_for_limit {
            spacing = min_spacing_for_limit;
        }
        
        // Apply a soft minimum spacing based on molecule size
        // The max_grid_points limit above already protects against excessive computation
        // So we use a modest minimum to prevent pathological cases
        let min_spacing = 0.5;  // Minimum 0.5Å regardless of quality
        spacing = spacing.max(min_spacing);

        // Check if cached SDF grid can be reused
        let current_params = SdfParams {
            probe_radius: self.probe_radius,
            quality: self.quality,
            surface_type: self.surface_type,
            atom_count: atoms.len(),
        };

        let sdf_grid = if self.cached_sdf_params.as_ref() == Some(&current_params) {
            self.cached_sdf.as_ref().unwrap()
        } else {
            let mut grid = Grid3D::from_bounds(min, max, spacing, 0.0);
            distance_field::compute_distance_field(
                &mut grid,
                atoms,
                self.surface_type,
                self.probe_radius,
            );
            self.cached_sdf = Some(grid);
            self.cached_sdf_params = Some(current_params);
            self.cached_sdf.as_ref().unwrap()
        };

        // Extract isosurface using marching cubes
        let result = marching_cubes::extract_isosurface_smooth(sdf_grid, 0.0);

        if result.positions.is_empty() {
            self.dirty = false;
            return;
        }

        // Weld duplicate vertices and average normals for smooth shading
        // Use spacing/2 as tolerance to catch vertices that should be the same
        let result = marching_cubes::weld_vertices(&result, spacing * 0.5);

        // Color vertices based on nearest atoms
        let colors = coloring::color_vertices(&result.positions, atom_colors);

        // Build mesh vertices
        self.vertices.reserve(result.positions.len());
        for ((position, normal), color) in result.positions.iter().zip(&result.normals).zip(&colors) {
            self.vertices.push(MeshVertex {
                position: *position,
                normal: *normal,
                color: *color,
            });
        }

        // Copy indices
        self.indices = result.indices;
        self.index_count = self.indices.len() as u32;
        self.dirty = false;
    }
}

impl Default for SurfaceRep {
    fn default() -> Self {
        Self::new()
    }
}

impl Representation for SurfaceRep {
    fn build(
        &mut self,
        molecule: &ObjectMolecule,
        coord_set: &CoordSet,
        colors: &ColorResolver,
        settings: &ResolvedSettings,
    ) {
        self.clear();

        // Get transparency and convert to alpha
        // PyMOL uses transparency (0=opaque, 1=invisible), we need alpha (1=opaque, 0=invisible)
        let transparency = settings.surface.transparency.clamp(0.0, 1.0);
        let alpha = 1.0 - transparency;

        let surface_type_setting = settings.surface.surface_type;
        let surface_solvent = settings.surface.solvent;
        self.surface_type = if surface_solvent {
            SurfaceType::SolventAccessible
        } else {
            SurfaceType::from_setting(surface_type_setting)
        };

        let individual_chains = settings.surface.individual_chains;
        let surface_color = settings.surface.color;

        // Collect atoms with SURFACE visibility
        let mut atoms: Vec<SurfaceAtom> = Vec::new();
        let mut atom_colors: Vec<AtomColor> = Vec::new();
        let mut chain_ids: Vec<String> = Vec::new();

        for (atom_idx, coord) in coord_set.iter_with_atoms() {
            let atom = match molecule.get_atom(atom_idx) {
                Some(a) => a,
                None => continue,
            };

            // Check if surface representation is visible for this atom
            if !atom.repr.visible_reps.is_visible(RepMask::SURFACE) {
                continue;
            }

            let position = [coord.x, coord.y, coord.z];
            let radius = atom.effective_vdw();
            let color = colors.resolve_rep_color(atom, atom.repr.colors.surface, surface_color);

            // Apply transparency to the color's alpha channel
            let color_with_alpha = [color[0], color[1], color[2], alpha];

            atoms.push(SurfaceAtom {
                position,
                radius,
                atom_index: atom_idx.into(),
            });

            atom_colors.push(AtomColor { position, color: color_with_alpha });

            if individual_chains {
                chain_ids.push(atom.residue.chain.clone());
            }
        }

        if individual_chains {
            self.build_per_chain(&atoms, &atom_colors, &chain_ids);
        } else {
            // Build from collected atoms (this calls clear() internally)
            self.build_from_atoms(&atoms, &atom_colors);
        }

        // Set has_transparency AFTER building since build methods call clear() which resets the flag
        self.has_transparency = transparency > 0.0;
        
        log::debug!("Surface built with transparency={}, has_transparency={}", transparency, self.has_transparency);
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

    fn clear_cpu_data(&mut self) {
        self.vertices = Vec::new();
        self.indices = Vec::new();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_surface_rep_new() {
        let rep = SurfaceRep::new();
        assert!(rep.is_dirty());
        assert!(rep.is_empty());
    }

    #[test]
    fn test_quality_to_spacing() {
        assert!((SurfaceRep::quality_to_spacing(-4) - 2.0).abs() < 0.001);
        assert!((SurfaceRep::quality_to_spacing(0) - 0.5).abs() < 0.001);
        assert!((SurfaceRep::quality_to_spacing(4) - 0.125).abs() < 0.001);
    }

    #[test]
    fn test_build_single_atom() {
        let mut rep = SurfaceRep::new();
        
        let atoms = vec![SurfaceAtom {
            position: [0.0, 0.0, 0.0],
            radius: 1.5,
            atom_index: 0,
        }];
        
        let colors = vec![AtomColor {
            position: [0.0, 0.0, 0.0],
            color: [1.0, 0.0, 0.0, 1.0],
        }];

        rep.set_quality(-2); // Lower quality for faster test
        rep.build_from_atoms(&atoms, &colors);

        assert!(!rep.is_empty(), "Surface should have triangles");
        assert!(!rep.vertices.is_empty());
        assert!(!rep.indices.is_empty());
        assert_eq!(rep.indices.len() % 3, 0, "Indices should be triangle triplets");
    }

    #[test]
    fn test_build_two_atoms() {
        let mut rep = SurfaceRep::new();
        
        let atoms = vec![
            SurfaceAtom {
                position: [0.0, 0.0, 0.0],
                radius: 1.5,
                atom_index: 0,
            },
            SurfaceAtom {
                position: [2.5, 0.0, 0.0],
                radius: 1.5,
                atom_index: 1,
            },
        ];
        
        let colors = vec![
            AtomColor {
                position: [0.0, 0.0, 0.0],
                color: [1.0, 0.0, 0.0, 1.0],
            },
            AtomColor {
                position: [2.5, 0.0, 0.0],
                color: [0.0, 0.0, 1.0, 1.0],
            },
        ];

        rep.set_quality(-2);
        rep.build_from_atoms(&atoms, &colors);

        assert!(!rep.is_empty(), "Surface should have triangles for two overlapping atoms");
    }

    #[test]
    fn test_surface_type_setting() {
        let mut rep = SurfaceRep::new();
        
        rep.set_surface_type(SurfaceType::SolventAccessible);
        assert_eq!(rep.surface_type, SurfaceType::SolventAccessible);
        
        rep.set_surface_type(SurfaceType::VanDerWaals);
        assert_eq!(rep.surface_type, SurfaceType::VanDerWaals);
    }

    #[test]
    fn test_marching_cubes_produces_output() {
        let atoms = vec![
            SurfaceAtom {
                position: [0.0, 0.0, 0.0],
                radius: 1.5,
                atom_index: 0,
            },
            SurfaceAtom {
                position: [2.5, 0.0, 0.0],
                radius: 1.5,
                atom_index: 1,
            },
        ];

        let colors = vec![
            AtomColor { position: [0.0, 0.0, 0.0], color: [1.0, 0.0, 0.0, 1.0] },
            AtomColor { position: [2.5, 0.0, 0.0], color: [0.0, 0.0, 1.0, 1.0] },
        ];

        let mut rep = SurfaceRep::new();
        rep.set_quality(-2);
        rep.build_from_atoms(&atoms, &colors);

        assert!(!rep.is_empty(), "Marching cubes should produce output");
        assert_eq!(rep.indices.len() % 3, 0);
    }
}
