//! Shared utility functions for bond visualization
//!
//! Used by both stick and line representations for valence bond display,
//! perpendicular calculations, and neighbor atom lookup.

use pymol_mol::{AtomIndex, BondOrder, CoordSet, ObjectMolecule};

/// Calculate perpendicular vector for bond offset using neighbor atom
///
/// Returns a normalized vector perpendicular to the bond direction that lies
/// in the plane defined by the bond and a neighboring atom. This ensures
/// double/triple bond cylinders stay within the molecular plane (e.g., in a ring).
///
/// If no neighbor is available, falls back to using a fixed reference axis.
pub fn calculate_perpendicular_with_neighbor(
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
pub fn calculate_perpendicular_fallback(bond_dir: [f32; 3]) -> [f32; 3] {
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
pub fn cross(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Normalize a 3D vector
#[inline]
pub fn normalize(v: [f32; 3]) -> [f32; 3] {
    let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if len > 0.0001 {
        [v[0] / len, v[1] / len, v[2] / len]
    } else {
        [0.0, 1.0, 0.0] // Fallback
    }
}

/// Get the offset factors for multiple bond display based on bond order
///
/// Returns factors that will be multiplied by valence_size to get actual offsets.
pub fn get_bond_offsets(order: BondOrder) -> &'static [f32] {
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
pub fn find_neighbor_position(
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
