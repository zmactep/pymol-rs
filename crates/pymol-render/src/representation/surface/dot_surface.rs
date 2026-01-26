//! Dot-based molecular surface generation.
//!
//! This module implements PyMOL's analytical surface algorithm:
//! 1. Generate dots on sphere surfaces (from tessellation)
//! 2. Filter dots to keep only exposed ones
//! 3. Output dots for triangulation
//!
//! This approach has O(n) complexity vs O(n³) for volumetric methods.

use super::distance_field::SurfaceAtom;
use super::grid::SpatialHash;
use super::sphere::{get_sphere, quality_to_level};
use rayon::prelude::*;

/// A surface dot with position, normal, and source atom
#[derive(Debug, Clone, Copy)]
pub struct SurfaceDot {
    /// World-space position of the dot
    pub position: [f32; 3],
    /// Surface normal (outward-pointing, same as unit sphere direction)
    pub normal: [f32; 3],
    /// Index of the atom this dot belongs to (for coloring, future use)
    #[allow(dead_code)]
    pub atom_index: usize,
}

/// Surface type for dot generation
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum DotSurfaceType {
    /// Van der Waals surface (union of atomic spheres)
    VanDerWaals,
    /// Solvent Accessible Surface (VdW + probe radius)
    SolventAccessible,
    /// Solvent Excluded Surface (Connolly surface)
    SolventExcluded,
}

/// Configuration for dot surface generation
#[derive(Debug, Clone)]
pub struct DotSurfaceConfig {
    /// Surface type to generate
    pub surface_type: DotSurfaceType,
    /// Probe radius for solvent surfaces (typically 1.4 Å for water)
    pub probe_radius: f32,
    /// Quality level (-4 to 4), affects sphere tessellation level
    pub quality: i32,
    /// Tolerance for exposure check (small positive value for robustness)
    pub tolerance: f32,
}

impl Default for DotSurfaceConfig {
    fn default() -> Self {
        Self {
            surface_type: DotSurfaceType::SolventAccessible,
            probe_radius: 1.4,
            quality: 0,
            tolerance: 0.01,
        }
    }
}

/// Generate exposed surface dots from atom positions
///
/// This is the main entry point for dot-based surface generation.
pub fn generate_surface_dots(
    atoms: &[SurfaceAtom],
    config: &DotSurfaceConfig,
) -> Vec<SurfaceDot> {
    if atoms.is_empty() {
        return Vec::new();
    }

    // Calculate effective radii based on surface type
    let effective_radii: Vec<f32> = atoms
        .iter()
        .map(|a| match config.surface_type {
            DotSurfaceType::VanDerWaals => a.radius,
            DotSurfaceType::SolventAccessible => a.radius + config.probe_radius,
            DotSurfaceType::SolventExcluded => a.radius + config.probe_radius,
        })
        .collect();

    // Find max radius for spatial hash cell sizing
    let max_radius = effective_radii.iter().cloned().fold(0.0_f32, f32::max);

    // Compute bounds with padding for spatial hash
    let padding = max_radius + 1.0;
    let (min, max) = compute_bounds(atoms, padding);

    // Build spatial hash for neighbor queries
    let positions: Vec<[f32; 3]> = atoms.iter().map(|a| a.position).collect();
    let spatial_hash = SpatialHash::new(&positions, &effective_radii, min, max);

    // Get sphere tessellation for the quality level
    let sphere_level = quality_to_level(config.quality);
    let sphere = get_sphere(sphere_level);

    // Generate dots in parallel
    let dots: Vec<SurfaceDot> = atoms
        .par_iter()
        .enumerate()
        .flat_map(|(atom_idx, atom)| {
            let mut atom_dots = Vec::with_capacity(sphere.num_dots() / 4); // Estimate ~25% exposed
            let mut nearby_buffer = Vec::with_capacity(64);
            
            let eff_radius = effective_radii[atom_idx];

            for unit_dot in &sphere.dots {
                // Transform dot to world space
                let world_pos = [
                    atom.position[0] + unit_dot[0] * eff_radius,
                    atom.position[1] + unit_dot[1] * eff_radius,
                    atom.position[2] + unit_dot[2] * eff_radius,
                ];

                // Check if dot is exposed (not buried by other atoms)
                let is_exposed = is_dot_exposed(
                    &world_pos,
                    atom_idx,
                    atoms,
                    &effective_radii,
                    &spatial_hash,
                    config.tolerance,
                    &mut nearby_buffer,
                );

                if is_exposed {
                    atom_dots.push(SurfaceDot {
                        position: world_pos,
                        normal: *unit_dot, // Unit sphere direction is the surface normal
                        atom_index: atom_idx,
                    });
                }
            }

            atom_dots
        })
        .collect();

    dots
}

/// Check if a dot is exposed (not buried inside another atom's sphere)
#[inline]
fn is_dot_exposed(
    dot_pos: &[f32; 3],
    source_atom_idx: usize,
    atoms: &[SurfaceAtom],
    effective_radii: &[f32],
    spatial_hash: &SpatialHash,
    tolerance: f32,
    nearby_buffer: &mut Vec<usize>,
) -> bool {
    // Query neighboring atoms
    spatial_hash.query_neighborhood(*dot_pos, nearby_buffer);

    for &neighbor_idx in nearby_buffer.iter() {
        // Skip the source atom
        if neighbor_idx == source_atom_idx {
            continue;
        }

        let neighbor = &atoms[neighbor_idx];
        let neighbor_radius = effective_radii[neighbor_idx];

        // Calculate distance from dot to neighbor center
        let dx = dot_pos[0] - neighbor.position[0];
        let dy = dot_pos[1] - neighbor.position[1];
        let dz = dot_pos[2] - neighbor.position[2];
        let dist_sq = dx * dx + dy * dy + dz * dz;

        // Dot is buried if it's inside neighbor's sphere (with tolerance)
        let threshold = neighbor_radius - tolerance;
        if dist_sq < threshold * threshold {
            return false;
        }
    }

    true
}

/// Compute bounding box of atoms with padding
fn compute_bounds(atoms: &[SurfaceAtom], padding: f32) -> ([f32; 3], [f32; 3]) {
    let mut min = [f32::MAX; 3];
    let mut max = [f32::MIN; 3];

    for atom in atoms {
        for i in 0..3 {
            min[i] = min[i].min(atom.position[i]);
            max[i] = max[i].max(atom.position[i]);
        }
    }

    // Add padding
    for i in 0..3 {
        min[i] -= padding;
        max[i] += padding;
    }

    (min, max)
}

/// Compute the maximum edge length cutoff for triangulation.
///
/// This is based on the sphere tessellation density and effective radii.
/// Edges longer than this should not form triangles.
pub fn compute_edge_cutoff(atoms: &[SurfaceAtom], config: &DotSurfaceConfig) -> f32 {
    if atoms.is_empty() {
        return 1.0;
    }

    // Calculate average effective radius
    let avg_radius: f32 = atoms
        .iter()
        .map(|a| match config.surface_type {
            DotSurfaceType::VanDerWaals => a.radius,
            DotSurfaceType::SolventAccessible | DotSurfaceType::SolventExcluded => {
                a.radius + config.probe_radius
            }
        })
        .sum::<f32>()
        / atoms.len() as f32;

    // Get sphere tessellation level
    let sphere_level = quality_to_level(config.quality);
    let sphere = get_sphere(sphere_level);

    // Estimate average edge length on a tessellated sphere
    // For an icosahedron subdivision, edge length ≈ 4πr² / (n_dots * avg_edges_per_dot / 2)
    // Simplified: edge_length ≈ 2 * r * sqrt(4π / n_dots)
    let n_dots = sphere.num_dots() as f32;
    let avg_edge_on_unit_sphere = 2.0 * (4.0 * std::f32::consts::PI / n_dots).sqrt();

    // Scale by average radius and add safety margin
    let base_cutoff = avg_edge_on_unit_sphere * avg_radius;

    // Add some margin (1.5x) to allow for variation in atom sizes
    base_cutoff * 1.5
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_atom(x: f32, y: f32, z: f32, radius: f32) -> SurfaceAtom {
        SurfaceAtom {
            position: [x, y, z],
            radius,
            atom_index: 0,
        }
    }

    #[test]
    fn test_single_atom_all_exposed() {
        let atoms = vec![make_atom(0.0, 0.0, 0.0, 1.5)];
        let config = DotSurfaceConfig {
            surface_type: DotSurfaceType::VanDerWaals,
            probe_radius: 0.0,
            quality: 0,
            tolerance: 0.01,
        };

        let dots = generate_surface_dots(&atoms, &config);

        // All dots should be exposed for a single atom
        let sphere = get_sphere(quality_to_level(config.quality));
        assert_eq!(dots.len(), sphere.num_dots());

        // All dots should be at the correct distance from center
        for dot in &dots {
            let dist_sq = dot.position[0].powi(2)
                + dot.position[1].powi(2)
                + dot.position[2].powi(2);
            let dist = dist_sq.sqrt();
            assert!((dist - 1.5).abs() < 0.01, "Dot at wrong distance: {}", dist);
        }
    }

    #[test]
    fn test_two_overlapping_atoms() {
        // Two atoms that overlap - some dots should be buried
        let atoms = vec![
            make_atom(0.0, 0.0, 0.0, 1.5),
            make_atom(2.0, 0.0, 0.0, 1.5), // 2.0 Å apart, radii 1.5 each = overlap
        ];
        let config = DotSurfaceConfig {
            surface_type: DotSurfaceType::VanDerWaals,
            probe_radius: 0.0,
            quality: 0,
            tolerance: 0.01,
        };

        let dots = generate_surface_dots(&atoms, &config);

        let sphere = get_sphere(quality_to_level(config.quality));
        let max_dots = sphere.num_dots() * 2;

        // Should have fewer dots than 2 full spheres due to overlap
        assert!(dots.len() < max_dots, "Expected fewer dots due to overlap");
        assert!(dots.len() > max_dots / 2, "Expected some dots to remain");
    }

    #[test]
    fn test_dot_normals_are_unit() {
        let atoms = vec![make_atom(5.0, 3.0, 2.0, 2.0)];
        let config = DotSurfaceConfig::default();

        let dots = generate_surface_dots(&atoms, &config);

        for dot in &dots {
            let len = (dot.normal[0].powi(2) + dot.normal[1].powi(2) + dot.normal[2].powi(2)).sqrt();
            assert!((len - 1.0).abs() < 1e-5, "Normal not unit length: {}", len);
        }
    }

    #[test]
    fn test_empty_atoms() {
        let atoms: Vec<SurfaceAtom> = vec![];
        let config = DotSurfaceConfig::default();

        let dots = generate_surface_dots(&atoms, &config);
        assert!(dots.is_empty());
    }

    #[test]
    fn test_edge_cutoff_calculation() {
        let atoms = vec![
            make_atom(0.0, 0.0, 0.0, 1.5),
            make_atom(3.0, 0.0, 0.0, 1.5),
        ];
        let config = DotSurfaceConfig {
            quality: 0,
            ..Default::default()
        };

        let cutoff = compute_edge_cutoff(&atoms, &config);

        // Cutoff should be reasonable (not too small, not too large)
        assert!(cutoff > 0.5, "Cutoff too small: {}", cutoff);
        assert!(cutoff < 5.0, "Cutoff too large: {}", cutoff);
    }
}
