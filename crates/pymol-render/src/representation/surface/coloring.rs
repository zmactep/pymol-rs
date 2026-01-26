//! Vertex coloring for molecular surfaces
//!
//! Assigns colors to surface vertices based on the nearest atoms.

use super::distance_field::{create_spatial_hash, SurfaceAtom};
use super::grid::SpatialHash;

/// Atom color information for surface coloring
#[derive(Debug, Clone)]
pub struct AtomColor {
    /// Position
    pub position: [f32; 3],
    /// RGBA color
    pub color: [f32; 4],
}

/// Assign colors to surface vertices based on nearest atoms
///
/// # Arguments
/// * `positions` - Vertex positions
/// * `atoms` - Atoms with colors
///
/// # Returns
/// A vector of RGBA colors, one per vertex
pub fn color_vertices(
    positions: &[[f32; 3]],
    atoms: &[AtomColor],
) -> Vec<[f32; 4]> {
    if atoms.is_empty() || positions.is_empty() {
        return vec![[0.5, 0.5, 0.5, 1.0]; positions.len()];
    }

    // Create spatial hash for fast lookups
    let surface_atoms: Vec<SurfaceAtom> = atoms
        .iter()
        .enumerate()
        .map(|(i, a)| SurfaceAtom {
            position: a.position,
            radius: 2.0, // Use a reasonable search radius
            atom_index: i,
        })
        .collect();
    
    let spatial_hash = create_spatial_hash(&surface_atoms);

    positions
        .iter()
        .map(|pos| {
            find_nearest_color(pos, atoms, &spatial_hash)
        })
        .collect()
}

/// Find the color of the nearest atom to a point
fn find_nearest_color(
    point: &[f32; 3],
    atoms: &[AtomColor],
    spatial_hash: &SpatialHash,
) -> [f32; 4] {
    let nearby = spatial_hash.query(*point);
    
    if nearby.is_empty() {
        // Fall back to brute force
        return atoms
            .iter()
            .min_by(|a, b| {
                let da = distance_squared(point, &a.position);
                let db = distance_squared(point, &b.position);
                da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|a| a.color)
            .unwrap_or([0.5, 0.5, 0.5, 1.0]);
    }

    nearby
        .iter()
        .min_by(|&&a, &&b| {
            let da = distance_squared(point, &atoms[a].position);
            let db = distance_squared(point, &atoms[b].position);
            da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|&i| atoms[i].color)
        .unwrap_or([0.5, 0.5, 0.5, 1.0])
}

/// Assign a uniform color to all vertices
#[allow(dead_code)]
pub fn color_uniform(positions: &[[f32; 3]], color: [f32; 4]) -> Vec<[f32; 4]> {
    vec![color; positions.len()]
}

/// Blend colors smoothly based on distance to multiple nearby atoms
///
/// This produces a smoother color transition on the surface.
#[allow(dead_code)]
pub fn color_vertices_smooth(
    positions: &[[f32; 3]],
    atoms: &[AtomColor],
    blend_radius: f32,
) -> Vec<[f32; 4]> {
    if atoms.is_empty() || positions.is_empty() {
        return vec![[0.5, 0.5, 0.5, 1.0]; positions.len()];
    }

    // Create spatial hash
    let surface_atoms: Vec<SurfaceAtom> = atoms
        .iter()
        .enumerate()
        .map(|(i, a)| SurfaceAtom {
            position: a.position,
            radius: blend_radius,
            atom_index: i,
        })
        .collect();
    
    let spatial_hash = create_spatial_hash(&surface_atoms);

    positions
        .iter()
        .map(|pos| {
            blend_nearby_colors(pos, atoms, &spatial_hash, blend_radius)
        })
        .collect()
}

/// Blend colors from nearby atoms using inverse distance weighting
fn blend_nearby_colors(
    point: &[f32; 3],
    atoms: &[AtomColor],
    spatial_hash: &SpatialHash,
    blend_radius: f32,
) -> [f32; 4] {
    let nearby = spatial_hash.query(*point);
    
    if nearby.is_empty() {
        // Fall back to nearest atom color
        return atoms
            .iter()
            .min_by(|a, b| {
                let da = distance_squared(point, &a.position);
                let db = distance_squared(point, &b.position);
                da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|a| a.color)
            .unwrap_or([0.5, 0.5, 0.5, 1.0]);
    }

    let blend_radius_sq = blend_radius * blend_radius;
    let mut total_weight = 0.0_f32;
    let mut blended = [0.0_f32; 4];
    let mut closest_idx = 0;
    let mut closest_dist = f32::MAX;

    for &atom_idx in nearby {
        let atom = &atoms[atom_idx];
        let dist_sq = distance_squared(point, &atom.position);
        
        if dist_sq < closest_dist {
            closest_dist = dist_sq;
            closest_idx = atom_idx;
        }

        if dist_sq < blend_radius_sq && dist_sq > 1e-10 {
            let dist = dist_sq.sqrt();
            let weight = 1.0 / (dist + 0.01); // Add small epsilon to avoid division by zero
            
            blended[0] += atom.color[0] * weight;
            blended[1] += atom.color[1] * weight;
            blended[2] += atom.color[2] * weight;
            blended[3] += atom.color[3] * weight;
            total_weight += weight;
        }
    }

    if total_weight > 0.0 {
        [
            blended[0] / total_weight,
            blended[1] / total_weight,
            blended[2] / total_weight,
            blended[3] / total_weight,
        ]
    } else {
        // Use nearest atom color
        atoms[closest_idx].color
    }
}

#[inline]
fn distance_squared(a: &[f32; 3], b: &[f32; 3]) -> f32 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    dx * dx + dy * dy + dz * dz
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_color_vertices_single_atom() {
        let atoms = vec![AtomColor {
            position: [0.0, 0.0, 0.0],
            color: [1.0, 0.0, 0.0, 1.0], // Red
        }];

        let positions = vec![
            [0.5, 0.0, 0.0],
            [0.0, 0.5, 0.0],
            [0.0, 0.0, 0.5],
        ];

        let colors = color_vertices(&positions, &atoms);
        
        assert_eq!(colors.len(), 3);
        for color in colors {
            assert_eq!(color, [1.0, 0.0, 0.0, 1.0]);
        }
    }

    #[test]
    fn test_color_vertices_two_atoms() {
        let atoms = vec![
            AtomColor {
                position: [0.0, 0.0, 0.0],
                color: [1.0, 0.0, 0.0, 1.0], // Red
            },
            AtomColor {
                position: [10.0, 0.0, 0.0],
                color: [0.0, 0.0, 1.0, 1.0], // Blue
            },
        ];

        let positions = vec![
            [1.0, 0.0, 0.0],  // Closer to first atom
            [9.0, 0.0, 0.0],  // Closer to second atom
        ];

        let colors = color_vertices(&positions, &atoms);
        
        assert_eq!(colors[0], [1.0, 0.0, 0.0, 1.0]); // Red
        assert_eq!(colors[1], [0.0, 0.0, 1.0, 1.0]); // Blue
    }

    #[test]
    fn test_color_uniform() {
        let positions = vec![[0.0, 0.0, 0.0]; 10];
        let color = [0.5, 0.5, 0.5, 1.0];
        
        let colors = color_uniform(&positions, color);
        
        assert_eq!(colors.len(), 10);
        for c in colors {
            assert_eq!(c, color);
        }
    }
}
