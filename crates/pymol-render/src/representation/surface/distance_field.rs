//! Distance field computation for molecular surfaces
//!
//! Computes signed distance fields from atomic positions and radii
//! for various surface types (VdW, SAS, SES).

use super::grid::{Grid3D, SpatialHash};
use rayon::prelude::*;

/// Surface type for distance field computation
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SurfaceType {
    /// Van der Waals surface (union of atomic spheres)
    #[default]
    VanDerWaals,
    /// Solvent Accessible Surface (VdW + probe radius)
    SolventAccessible,
    /// Solvent Excluded Surface (Connolly surface)
    /// This is the surface traced by the probe sphere center
    SolventExcluded,
}

impl SurfaceType {
    /// Create from PyMOL surface_type setting value
    pub fn from_setting(value: i32) -> Self {
        match value {
            1 => SurfaceType::SolventAccessible,
            2 => SurfaceType::VanDerWaals,
            _ => SurfaceType::SolventExcluded, // Default (0)
        }
    }
}

/// Atom data for surface generation
#[derive(Debug, Clone)]
pub struct SurfaceAtom {
    /// Position
    pub position: [f32; 3],
    /// Van der Waals radius
    pub radius: f32,
    /// Original atom index (for coloring)
    pub atom_index: usize,
}

/// Compute the distance field on a grid from atoms using sparse computation
///
/// This implementation only computes exact distances for grid cells near atoms,
/// significantly reducing computation time for large molecules with internal cavities.
///
/// # Arguments
/// * `grid` - The grid to fill with distance values
/// * `atoms` - List of atoms with positions and radii
/// * `surface_type` - Type of surface to compute
/// * `probe_radius` - Solvent probe radius (typically 1.4 Å)
pub fn compute_distance_field(
    grid: &mut Grid3D,
    atoms: &[SurfaceAtom],
    surface_type: SurfaceType,
    probe_radius: f32,
) {
    if atoms.is_empty() {
        return;
    }

    // Extract radii for spatial hash and distance computation
    let radii: Vec<f32> = atoms.iter().map(|a| {
        match surface_type {
            SurfaceType::VanDerWaals => a.radius,
            SurfaceType::SolventAccessible => a.radius + probe_radius,
            SurfaceType::SolventExcluded => a.radius + probe_radius,
        }
    }).collect();

    let vd = grid.vertex_dims();
    let total_cells = vd[0] * vd[1] * vd[2];
    
    // Build active cell mask - cells near any atom that need exact distance computation
    // The narrowband extends slightly beyond the effective radius to ensure marching cubes
    // gets accurate values where the surface might pass
    // Minimum of 1.5Å ensures reliability with 3.0Å spatial hash cells
    let narrowband = (grid.spacing * 2.0).max(1.5);
    let active_mask = build_active_mask(grid, atoms, &radii, narrowband);
    
    // Count active cells for logging/debugging
    let active_count = active_mask.iter().filter(|&&x| x).count();
    let skip_ratio = 1.0 - (active_count as f64 / total_cells as f64);
    
    // Only use sparse computation if it saves significant work (>20% cells skipped)
    if skip_ratio > 0.2 {
        compute_distance_field_sparse(grid, atoms, &radii, surface_type, probe_radius, &active_mask);
    } else {
        compute_distance_field_dense(grid, atoms, &radii, surface_type, probe_radius);
    }
}

/// Build a mask of grid cells that need exact distance computation
fn build_active_mask(
    grid: &Grid3D,
    atoms: &[SurfaceAtom],
    effective_radii: &[f32],
    narrowband: f32,
) -> Vec<bool> {
    let vd = grid.vertex_dims();
    let total_cells = vd[0] * vd[1] * vd[2];
    let mut active = vec![false; total_cells];
    
    let spacing = grid.spacing;
    let origin = grid.origin;
    
    // For each atom, mark grid cells within (radius + narrowband) as active
    for (atom, &radius) in atoms.iter().zip(effective_radii.iter()) {
        let influence_radius = radius + narrowband;
        
        // Calculate bounding box of affected cells
        let min_x = ((atom.position[0] - influence_radius - origin[0]) / spacing).floor() as i32;
        let max_x = ((atom.position[0] + influence_radius - origin[0]) / spacing).ceil() as i32;
        let min_y = ((atom.position[1] - influence_radius - origin[1]) / spacing).floor() as i32;
        let max_y = ((atom.position[1] + influence_radius - origin[1]) / spacing).ceil() as i32;
        let min_z = ((atom.position[2] - influence_radius - origin[2]) / spacing).floor() as i32;
        let max_z = ((atom.position[2] + influence_radius - origin[2]) / spacing).ceil() as i32;
        
        // Clamp to grid bounds
        let min_x = min_x.max(0) as usize;
        let max_x = (max_x as usize).min(vd[0] - 1);
        let min_y = min_y.max(0) as usize;
        let max_y = (max_y as usize).min(vd[1] - 1);
        let min_z = min_z.max(0) as usize;
        let max_z = (max_z as usize).min(vd[2] - 1);
        
        // Mark all cells in the bounding box
        // (Could be optimized to check actual sphere intersection, but box is fast enough)
        for z in min_z..=max_z {
            for y in min_y..=max_y {
                let row_start = z * vd[1] * vd[0] + y * vd[0];
                for x in min_x..=max_x {
                    active[row_start + x] = true;
                }
            }
        }
    }
    
    active
}

/// Sparse distance field computation - only computes for active cells
fn compute_distance_field_sparse(
    grid: &mut Grid3D,
    atoms: &[SurfaceAtom],
    radii: &[f32],
    surface_type: SurfaceType,
    probe_radius: f32,
    active_mask: &[bool],
) {
    let vd = grid.vertex_dims();
    let grid_min = grid.origin;
    let grid_max = [
        grid.origin[0] + (vd[0] - 1) as f32 * grid.spacing,
        grid.origin[1] + (vd[1] - 1) as f32 * grid.spacing,
        grid.origin[2] + (vd[2] - 1) as f32 * grid.spacing,
    ];
    
    // Extend radii by narrowband for spatial hash insertion
    // This ensures atoms are found when querying any active cell
    // Minimum of 1.5Å ensures single-cell queries work reliably
    // (half the spatial hash cell size of 3.0Å)
    let narrowband = (grid.spacing * 2.0).max(1.5);
    let extended_radii: Vec<f32> = radii.iter().map(|r| r + narrowband).collect();
    
    let spatial_hash = SpatialHash::new(
        &atoms.iter().map(|a| a.position).collect::<Vec<_>>(),
        &extended_radii,
        grid_min,
        grid_max,
    );

    // Process in parallel by z-slices
    let values: Vec<f32> = (0..vd[2])
        .into_par_iter()
        .flat_map(|z| {
            let mut slice_values = Vec::with_capacity(vd[0] * vd[1]);
            let z_offset = z * vd[1] * vd[0];
            
            for y in 0..vd[1] {
                let row_offset = z_offset + y * vd[0];
                for x in 0..vd[0] {
                    let idx = row_offset + x;
                    
                    if active_mask[idx] {
                        // Compute exact distance for active cells using fast single-cell query
                        let pos = grid.vertex_position(x, y, z);
                        let dist = compute_point_distance_fast(
                            &pos, atoms, radii, &spatial_hash, surface_type, probe_radius
                        );
                        slice_values.push(dist);
                    } else {
                        // Inactive cells are far from surface - mark as outside
                        slice_values.push(f32::MAX);
                    }
                }
            }
            slice_values
        })
        .collect();

    grid.values = values;
}

/// Dense distance field computation - computes for all cells
fn compute_distance_field_dense(
    grid: &mut Grid3D,
    atoms: &[SurfaceAtom],
    radii: &[f32],
    surface_type: SurfaceType,
    probe_radius: f32,
) {
    let vd = grid.vertex_dims();
    let grid_min = grid.origin;
    let grid_max = [
        grid.origin[0] + (vd[0] - 1) as f32 * grid.spacing,
        grid.origin[1] + (vd[1] - 1) as f32 * grid.spacing,
        grid.origin[2] + (vd[2] - 1) as f32 * grid.spacing,
    ];
    
    // Extend radii to ensure atoms are found for all grid cells
    // Minimum of 1.5Å ensures single-cell queries work reliably
    // (half the spatial hash cell size of 3.0Å)
    let margin = (grid.spacing * 2.0).max(1.5);
    let extended_radii: Vec<f32> = radii.iter().map(|r| r + margin).collect();
    
    let spatial_hash = SpatialHash::new(
        &atoms.iter().map(|a| a.position).collect::<Vec<_>>(),
        &extended_radii,
        grid_min,
        grid_max,
    );

    let values: Vec<f32> = (0..vd[2])
        .into_par_iter()
        .flat_map(|z| {
            let mut slice_values = Vec::with_capacity(vd[0] * vd[1]);
            for y in 0..vd[1] {
                for x in 0..vd[0] {
                    let pos = grid.vertex_position(x, y, z);
                    let dist = compute_point_distance_fast(
                        &pos, atoms, radii, &spatial_hash, surface_type, probe_radius
                    );
                    slice_values.push(dist);
                }
            }
            slice_values
        })
        .collect();

    grid.values = values;
}

/// Fast distance computation using single-cell query
/// Since atoms are inserted into all cells they influence, single-cell query is sufficient
#[inline]
fn compute_point_distance_fast(
    point: &[f32; 3],
    atoms: &[SurfaceAtom],
    effective_radii: &[f32],
    spatial_hash: &SpatialHash,
    surface_type: SurfaceType,
    probe_radius: f32,
) -> f32 {
    // Use fast single-cell query - atoms are already inserted into all cells they influence
    let nearby = spatial_hash.query(*point);
    
    if nearby.is_empty() {
        return f32::MAX;
    }

    // Find minimum distance to any atom surface
    let mut min_dist = f32::MAX;
    
    for &atom_idx in nearby {
        let atom = &atoms[atom_idx];
        let r = effective_radii[atom_idx];
        
        // Squared distance from point to atom center
        let dx = point[0] - atom.position[0];
        let dy = point[1] - atom.position[1];
        let dz = point[2] - atom.position[2];
        let dist_to_center = (dx * dx + dy * dy + dz * dz).sqrt();
        
        // Distance to atom surface
        let dist_to_surface = dist_to_center - r;
        
        if dist_to_surface < min_dist {
            min_dist = dist_to_surface;
        }
    }

    // For SES, subtract probe radius at the end
    match surface_type {
        SurfaceType::SolventExcluded => min_dist - probe_radius,
        _ => min_dist,
    }
}

/// Compute distance from a point to the surface using a reusable buffer for neighbor queries
/// 
/// Optimized to avoid sqrt for atoms that are clearly not the closest.
/// Note: Currently unused in favor of compute_point_distance_fast, kept for potential future use.
#[allow(dead_code)]
fn compute_point_distance_with_buffer(
    point: &[f32; 3],
    atoms: &[SurfaceAtom],
    effective_radii: &[f32],
    spatial_hash: &SpatialHash,
    surface_type: SurfaceType,
    probe_radius: f32,
    nearby_buffer: &mut Vec<usize>,
) -> f32 {
    // Query 3x3x3 neighborhood of cells for robustness at cell boundaries
    spatial_hash.query_neighborhood(*point, nearby_buffer);
    
    if nearby_buffer.is_empty() {
        // No nearby atoms, return large positive distance
        return f32::MAX;
    }

    // Find minimum distance to any atom surface
    // Use squared distances where possible to avoid sqrt
    let mut min_dist = f32::MAX;
    let mut min_dist_plus_r_sq = f32::MAX; // (min_dist + r_max)^2 for early rejection
    
    for &atom_idx in nearby_buffer.iter() {
        let atom = &atoms[atom_idx];
        let r = effective_radii[atom_idx];
        
        // Squared distance from point to atom center
        let dx = point[0] - atom.position[0];
        let dy = point[1] - atom.position[1];
        let dz = point[2] - atom.position[2];
        let dist_sq = dx * dx + dy * dy + dz * dz;
        
        // Early rejection: if dist > min_dist + r, this atom can't be closer
        // Equivalently: dist_sq > (min_dist + r)^2
        if dist_sq > min_dist_plus_r_sq {
            continue;
        }
        
        // Must compute actual distance
        let dist_to_center = dist_sq.sqrt();
        let dist_to_surface = dist_to_center - r;
        
        if dist_to_surface < min_dist {
            min_dist = dist_to_surface;
            // Update threshold for early rejection of future atoms
            // We can reject atoms where dist > min_dist + max_radius_in_set
            // For simplicity, use a conservative estimate
            let threshold = min_dist + r + 5.0; // 5 Å margin
            min_dist_plus_r_sq = threshold * threshold;
        }
    }

    // For SES, we need to subtract probe radius at the end
    // (the effective radii already include probe for SAS)
    match surface_type {
        SurfaceType::SolventExcluded => min_dist - probe_radius,
        _ => min_dist,
    }
}

/// Compute distance from a point to the surface (legacy single-cell version)
#[allow(dead_code)]
fn compute_point_distance(
    point: &[f32; 3],
    atoms: &[SurfaceAtom],
    effective_radii: &[f32],
    spatial_hash: &SpatialHash,
    surface_type: SurfaceType,
    probe_radius: f32,
) -> f32 {
    let nearby = spatial_hash.query(*point);
    
    if nearby.is_empty() {
        // No nearby atoms, return large positive distance
        return f32::MAX;
    }

    // Find minimum distance to any atom surface
    let mut min_dist = f32::MAX;
    
    for &atom_idx in nearby {
        let atom = &atoms[atom_idx];
        let r = effective_radii[atom_idx];
        
        // Distance from point to atom center
        let dx = point[0] - atom.position[0];
        let dy = point[1] - atom.position[1];
        let dz = point[2] - atom.position[2];
        let dist_to_center = (dx * dx + dy * dy + dz * dz).sqrt();
        
        // Distance to atom surface
        let dist_to_surface = dist_to_center - r;
        
        if dist_to_surface < min_dist {
            min_dist = dist_to_surface;
        }
    }

    // For SES, we need to subtract probe radius at the end
    // (the effective radii already include probe for SAS)
    match surface_type {
        SurfaceType::SolventExcluded => min_dist - probe_radius,
        _ => min_dist,
    }
}

/// Find the nearest atom to a point
#[allow(dead_code)]
pub fn find_nearest_atom(
    point: &[f32; 3],
    atoms: &[SurfaceAtom],
    spatial_hash: &SpatialHash,
) -> Option<usize> {
    let nearby = spatial_hash.query(*point);
    
    if nearby.is_empty() {
        // Fall back to brute force for points outside hash coverage
        return atoms
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| {
                let da = distance_squared(point, &a.position);
                let db = distance_squared(point, &b.position);
                da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(i, _)| i);
    }

    nearby
        .iter()
        .min_by(|&&a, &&b| {
            let da = distance_squared(point, &atoms[a].position);
            let db = distance_squared(point, &atoms[b].position);
            da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
        })
        .copied()
}

/// Create a spatial hash for atom lookups
pub fn create_spatial_hash(atoms: &[SurfaceAtom]) -> SpatialHash {
    let positions: Vec<[f32; 3]> = atoms.iter().map(|a| a.position).collect();
    let radii: Vec<f32> = atoms.iter().map(|a| a.radius).collect();
    let (min, max) = compute_bounds(&positions);
    SpatialHash::new(&positions, &radii, min, max)
}

#[inline]
fn distance_squared(a: &[f32; 3], b: &[f32; 3]) -> f32 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    dx * dx + dy * dy + dz * dz
}

/// Compute bounding box of positions
pub fn compute_bounds(positions: &[[f32; 3]]) -> ([f32; 3], [f32; 3]) {
    if positions.is_empty() {
        return ([0.0; 3], [0.0; 3]);
    }

    let mut min = positions[0];
    let mut max = positions[0];

    for pos in positions.iter().skip(1) {
        min[0] = min[0].min(pos[0]);
        min[1] = min[1].min(pos[1]);
        min[2] = min[2].min(pos[2]);
        max[0] = max[0].max(pos[0]);
        max[1] = max[1].max(pos[1]);
        max[2] = max[2].max(pos[2]);
    }

    (min, max)
}

/// Compute bounding box with radius expansion
pub fn compute_bounds_with_radii(atoms: &[SurfaceAtom], padding: f32) -> ([f32; 3], [f32; 3]) {
    if atoms.is_empty() {
        return ([0.0; 3], [0.0; 3]);
    }

    let mut min = [f32::MAX; 3];
    let mut max = [f32::MIN; 3];

    for atom in atoms {
        let r = atom.radius + padding;
        min[0] = min[0].min(atom.position[0] - r);
        min[1] = min[1].min(atom.position[1] - r);
        min[2] = min[2].min(atom.position[2] - r);
        max[0] = max[0].max(atom.position[0] + r);
        max[1] = max[1].max(atom.position[1] + r);
        max[2] = max[2].max(atom.position[2] + r);
    }

    (min, max)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compute_bounds() {
        let positions = vec![
            [0.0, 0.0, 0.0],
            [10.0, 5.0, 3.0],
            [-2.0, 8.0, 1.0],
        ];
        
        let (min, max) = compute_bounds(&positions);
        assert_eq!(min, [-2.0, 0.0, 0.0]);
        assert_eq!(max, [10.0, 8.0, 3.0]);
    }

    #[test]
    fn test_single_atom_distance_field() {
        let atoms = vec![SurfaceAtom {
            position: [0.0, 0.0, 0.0],
            radius: 1.5,
            atom_index: 0,
        }];

        let mut grid = Grid3D::from_bounds(
            [-3.0, -3.0, -3.0],
            [3.0, 3.0, 3.0],
            0.5,
            0.0,
        );

        compute_distance_field(&mut grid, &atoms, SurfaceType::VanDerWaals, 1.4);

        // Center should be negative (inside)
        let vd = grid.vertex_dims();
        let center_x = vd[0] / 2;
        let center_y = vd[1] / 2;
        let center_z = vd[2] / 2;
        let center_val = grid.get(center_x, center_y, center_z);
        assert!(center_val < 0.0, "Center should be inside the surface");

        // Corner should be positive (outside)
        let corner_val = grid.get(0, 0, 0);
        assert!(corner_val > 0.0, "Corner should be outside the surface");
    }

    #[test]
    fn test_surface_type_from_setting() {
        assert_eq!(SurfaceType::from_setting(0), SurfaceType::SolventExcluded);
        assert_eq!(SurfaceType::from_setting(1), SurfaceType::SolventAccessible);
        assert_eq!(SurfaceType::from_setting(2), SurfaceType::VanDerWaals);
    }
}
