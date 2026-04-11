//! Coordinate set data structure
//!
//! Provides the `CoordSet` struct for storing atomic coordinates in a molecular state.
//! Multi-state molecules have multiple coordinate sets, one per state/frame.

use serde::{Deserialize, Serialize};

use crate::index::{AtomIndex, CoordIndex, INVALID_INDEX};
use lin_alg::f32::{Mat4, Vec3};
use pymol_settings::GlobalSettings;

/// Mapping between atom indices and coordinate indices
///
/// For memory efficiency, uses identity mapping when all atoms have coordinates
/// (the common case), and sparse mapping only when some atoms are missing.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
enum CoordMapping {
    /// Identity mapping: coord_idx == atom_idx for all indices
    /// No arrays needed - just need to know the count
    #[default]
    Dense,
    /// Sparse mapping: only some atoms have coordinates
    Sparse {
        /// Mapping from coordinate index to atom index
        idx_to_atm: Vec<AtomIndex>,
        /// Mapping from atom index to coordinate index (INVALID_INDEX if missing)
        atm_to_idx: Vec<u32>,
    },
}


/// Iterator over (AtomIndex, Vec3) pairs from a CoordSet
pub struct CoordAtomIter<'a> {
    coords: &'a [f32],
    mapping: &'a CoordMapping,
    n_index: usize,
    current: usize,
}

impl<'a> Iterator for CoordAtomIter<'a> {
    type Item = (AtomIndex, Vec3);

    fn next(&mut self) -> Option<Self::Item> {
        if self.current >= self.n_index {
            return None;
        }

        let coord_idx = self.current;
        let base = coord_idx * 3;
        let coord = Vec3::new(
            self.coords[base],
            self.coords[base + 1],
            self.coords[base + 2],
        );

        let atom_idx = match self.mapping {
            CoordMapping::Dense => AtomIndex(coord_idx as u32),
            CoordMapping::Sparse { idx_to_atm, .. } => idx_to_atm[coord_idx],
        };

        self.current += 1;
        Some((atom_idx, coord))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = self.n_index - self.current;
        (remaining, Some(remaining))
    }
}

impl<'a> ExactSizeIterator for CoordAtomIter<'a> {}

/// Symmetry information for a coordinate set
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Symmetry {
    /// Space group name (e.g., "P 21 21 21")
    pub space_group: String,
    /// Unit cell parameters: a, b, c (in Angstroms)
    pub cell_lengths: [f32; 3],
    /// Unit cell angles: alpha, beta, gamma (in degrees)
    pub cell_angles: [f32; 3],
    /// Z value (number of molecules per unit cell)
    pub z_value: i32,
}

impl Symmetry {
    /// Create new symmetry with given space group and cell parameters
    pub fn new(
        space_group: impl Into<String>,
        lengths: [f32; 3],
        angles: [f32; 3],
    ) -> Self {
        Symmetry {
            space_group: space_group.into(),
            cell_lengths: lengths,
            cell_angles: angles,
            z_value: 1,
        }
    }

    /// Create P1 symmetry (no symmetry)
    pub fn p1(lengths: [f32; 3], angles: [f32; 3]) -> Self {
        Self::new("P 1", lengths, angles)
    }
}

/// Coordinate set containing atomic positions for one state/frame
///
/// A coordinate set stores the 3D coordinates of atoms in a particular
/// conformational state. Multi-state molecules (e.g., NMR ensembles, MD
/// trajectories) have multiple coordinate sets.
///
/// Not all atoms in the molecule necessarily have coordinates in every
/// coordinate set. The `mapping` field handles the atom-to-coordinate mapping,
/// using an efficient identity mapping when all atoms have coordinates.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CoordSet {
    /// Flat array of coordinates (3 floats per atom: x, y, z)
    coords: Vec<f32>,

    /// Number of atoms with coordinates in this set
    n_index: usize,

    /// Mapping between atom indices and coordinate indices
    /// Dense when all atoms have coordinates (identity mapping)
    /// Sparse when some atoms are missing coordinates
    mapping: CoordMapping,

    /// Name of this coordinate set (e.g., for trajectory frames)
    pub name: String,

    /// Symmetry information (unit cell, space group)
    pub symmetry: Option<Symmetry>,

    /// Per-state settings
    pub settings: Option<GlobalSettings>,
}

impl Default for CoordSet {
    fn default() -> Self {
        CoordSet {
            coords: Vec::new(),
            n_index: 0,
            mapping: CoordMapping::Dense,
            name: String::new(),
            symmetry: None,
            settings: None,
        }
    }
}

impl CoordSet {
    /// Create a new empty coordinate set
    pub fn new() -> Self {
        CoordSet::default()
    }

    /// Create a coordinate set with pre-allocated capacity
    ///
    /// Initially uses dense mapping. Call `add_coord` to switch to sparse
    /// mapping if atoms are added out of order.
    pub fn with_capacity(n_atoms: usize) -> Self {
        CoordSet {
            coords: Vec::with_capacity(n_atoms * 3),
            n_index: 0,
            mapping: CoordMapping::Dense,
            name: String::new(),
            symmetry: None,
            settings: None,
        }
    }

    /// Create a coordinate set from a flat array of coordinates
    ///
    /// Assumes all atoms have coordinates (1:1 mapping).
    /// Uses dense mapping for memory efficiency.
    pub fn from_coords(coords: Vec<f32>) -> Self {
        let n_atoms = coords.len() / 3;

        CoordSet {
            coords,
            n_index: n_atoms,
            mapping: CoordMapping::Dense,
            name: String::new(),
            symmetry: None,
            settings: None,
        }
    }

    /// Create a coordinate set from a vector of Vec3 positions
    ///
    /// Uses dense mapping for memory efficiency.
    pub fn from_vec3(positions: &[Vec3]) -> Self {
        let n_atoms = positions.len();
        let mut coords = Vec::with_capacity(n_atoms * 3);
        for pos in positions {
            coords.push(pos.x);
            coords.push(pos.y);
            coords.push(pos.z);
        }

        CoordSet {
            coords,
            n_index: n_atoms,
            mapping: CoordMapping::Dense,
            name: String::new(),
            symmetry: None,
            settings: None,
        }
    }

    /// Get the number of atoms with coordinates
    #[inline]
    pub fn len(&self) -> usize {
        self.n_index
    }

    /// Check if the coordinate set is empty
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.n_index == 0
    }



    /// Get coordinates by coordinate index
    #[inline]
    pub fn get_coord(&self, idx: CoordIndex) -> Option<Vec3> {
        let i = idx.as_usize();
        if i >= self.n_index {
            return None;
        }
        let base = i * 3;
        Some(Vec3::new(
            self.coords[base],
            self.coords[base + 1],
            self.coords[base + 2],
        ))
    }

    /// Set coordinates by coordinate index
    #[inline]
    pub fn set_coord(&mut self, idx: CoordIndex, coord: Vec3) {
        let i = idx.as_usize();
        if i < self.n_index {
            let base = i * 3;
            self.coords[base] = coord.x;
            self.coords[base + 1] = coord.y;
            self.coords[base + 2] = coord.z;
        }
    }

    /// Get coordinates by atom index
    pub fn get_atom_coord(&self, atom: AtomIndex) -> Option<Vec3> {
        let atom_idx = atom.as_usize();
        match &self.mapping {
            CoordMapping::Dense => {
                if atom_idx >= self.n_index {
                    return None;
                }
                self.get_coord(CoordIndex(atom_idx as u32))
            }
            CoordMapping::Sparse { atm_to_idx, .. } => {
                if atom_idx >= atm_to_idx.len() {
                    return None;
                }
                let coord_idx = atm_to_idx[atom_idx];
                if coord_idx == INVALID_INDEX {
                    return None;
                }
                self.get_coord(CoordIndex(coord_idx))
            }
        }
    }

    /// Set coordinates by atom index
    pub fn set_atom_coord(&mut self, atom: AtomIndex, coord: Vec3) -> bool {
        let atom_idx = atom.as_usize();
        match &self.mapping {
            CoordMapping::Dense => {
                if atom_idx >= self.n_index {
                    return false;
                }
                self.set_coord(CoordIndex(atom_idx as u32), coord);
                true
            }
            CoordMapping::Sparse { atm_to_idx, .. } => {
                if atom_idx >= atm_to_idx.len() {
                    return false;
                }
                let coord_idx = atm_to_idx[atom_idx];
                if coord_idx == INVALID_INDEX {
                    return false;
                }
                self.set_coord(CoordIndex(coord_idx), coord);
                true
            }
        }
    }

    /// Check if an atom has coordinates in this set
    #[inline]
    pub fn has_atom(&self, atom: AtomIndex) -> bool {
        let atom_idx = atom.as_usize();
        match &self.mapping {
            CoordMapping::Dense => atom_idx < self.n_index,
            CoordMapping::Sparse { atm_to_idx, .. } => {
                atom_idx < atm_to_idx.len() && atm_to_idx[atom_idx] != INVALID_INDEX
            }
        }
    }

    /// Get the coordinate index for an atom index
    #[inline]
    pub fn atom_to_coord(&self, atom: AtomIndex) -> Option<CoordIndex> {
        let atom_idx = atom.as_usize();
        match &self.mapping {
            CoordMapping::Dense => {
                if atom_idx < self.n_index {
                    Some(CoordIndex(atom_idx as u32))
                } else {
                    None
                }
            }
            CoordMapping::Sparse { atm_to_idx, .. } => {
                if atom_idx >= atm_to_idx.len() {
                    return None;
                }
                let coord_idx = atm_to_idx[atom_idx];
                if coord_idx == INVALID_INDEX {
                    None
                } else {
                    Some(CoordIndex(coord_idx))
                }
            }
        }
    }

    /// Iterate over all coordinates
    pub fn iter(&self) -> impl Iterator<Item = Vec3> + '_ {
        (0..self.n_index).map(move |i| {
            let base = i * 3;
            Vec3::new(
                self.coords[base],
                self.coords[base + 1],
                self.coords[base + 2],
            )
        })
    }

    /// Iterate over (atom_index, coord) pairs
    pub fn iter_with_atoms(&self) -> CoordAtomIter<'_> {
        CoordAtomIter {
            coords: &self.coords,
            mapping: &self.mapping,
            n_index: self.n_index,
            current: 0,
        }
    }

    /// Add coordinates for an atom
    ///
    /// This extends the coordinate set to include the given atom.
    /// If atoms are added out of order (not sequential), converts to sparse mapping.
    pub fn add_coord(&mut self, atom: AtomIndex, coord: Vec3) {
        let atom_idx = atom.as_usize();
        let coord_idx = self.n_index;

        // Check if we need to convert from dense to sparse
        // Dense mapping requires: atom_idx == coord_idx (sequential 0, 1, 2, ...)
        if matches!(self.mapping, CoordMapping::Dense) && atom_idx != coord_idx {
            // Convert to sparse mapping
            let idx_to_atm: Vec<AtomIndex> = (0..self.n_index as u32).map(AtomIndex).collect();
            let mut atm_to_idx: Vec<u32> = (0..self.n_index as u32).collect();

            // Extend atm_to_idx to accommodate the new atom
            while atm_to_idx.len() <= atom_idx {
                atm_to_idx.push(INVALID_INDEX);
            }

            self.mapping = CoordMapping::Sparse {
                idx_to_atm,
                atm_to_idx,
            };
        }

        // Add the coordinate
        self.coords.push(coord.x);
        self.coords.push(coord.y);
        self.coords.push(coord.z);
        self.n_index += 1;

        // Update mapping
        match &mut self.mapping {
            CoordMapping::Dense => {
                // Dense mapping: atom_idx == coord_idx, no explicit tracking needed
            }
            CoordMapping::Sparse {
                idx_to_atm,
                atm_to_idx,
            } => {
                // Extend atm_to_idx if needed
                while atm_to_idx.len() <= atom_idx {
                    atm_to_idx.push(INVALID_INDEX);
                }
                idx_to_atm.push(atom);
                atm_to_idx[atom_idx] = coord_idx as u32;
            }
        }
    }

    /// Compute the center of mass (geometric center)
    pub fn center(&self) -> Vec3 {
        if self.n_index == 0 {
            return Vec3::new(0.0, 0.0, 0.0);
        }

        let mut sum = Vec3::new(0.0, 0.0, 0.0);
        for coord in self.iter() {
            sum += coord;
        }
        sum * (1.0 / self.n_index as f32)
    }

    /// Compute the bounding box (min, max corners)
    pub fn bounding_box(&self) -> Option<(Vec3, Vec3)> {
        if self.n_index == 0 {
            return None;
        }

        let mut min = Vec3::new(f32::MAX, f32::MAX, f32::MAX);
        let mut max = Vec3::new(f32::MIN, f32::MIN, f32::MIN);

        for coord in self.iter() {
            min.x = min.x.min(coord.x);
            min.y = min.y.min(coord.y);
            min.z = min.z.min(coord.z);
            max.x = max.x.max(coord.x);
            max.y = max.y.max(coord.y);
            max.z = max.z.max(coord.z);
        }

        Some((min, max))
    }

    /// Translate all coordinates by a vector
    pub fn translate(&mut self, delta: Vec3) {
        for i in 0..self.n_index {
            let base = i * 3;
            self.coords[base] += delta.x;
            self.coords[base + 1] += delta.y;
            self.coords[base + 2] += delta.z;
        }
    }

    /// Center the coordinates at the origin
    pub fn center_origin(&mut self) {
        let center = self.center();
        self.translate(-center);
    }

    /// Translate specific atoms by atom indices
    pub fn translate_atoms(&mut self, atoms: &[AtomIndex], delta: Vec3) {
        for &atom in atoms {
            if let Some(coord_idx) = self.atom_to_coord(atom) {
                let i = coord_idx.as_usize();
                let base = i * 3;
                self.coords[base] += delta.x;
                self.coords[base + 1] += delta.y;
                self.coords[base + 2] += delta.z;
            }
        }
    }

    /// Rotate all coordinates about an axis passing through origin
    ///
    /// # Arguments
    /// * `axis` - Rotation axis (will be normalized)
    /// * `angle` - Rotation angle in radians
    /// * `origin` - Center of rotation
    pub fn rotate(&mut self, axis: Vec3, angle: f32, origin: Vec3) {
        let matrix = rotation_matrix(axis, angle, origin);
        self.transform(&matrix);
    }

    /// Rotate specific atoms by atom indices
    pub fn rotate_atoms(&mut self, atoms: &[AtomIndex], axis: Vec3, angle: f32, origin: Vec3) {
        let matrix = rotation_matrix(axis, angle, origin);
        self.transform_atoms(atoms, &matrix);
    }

    /// Apply a 4x4 homogeneous transformation matrix to all coordinates
    ///
    /// The matrix is in column-major format (lin_alg convention):
    /// ```text
    /// data = [col0, col1, col2, col3] flattened
    /// ```
    /// Translation is at data[12], data[13], data[14].
    pub fn transform(&mut self, matrix: &Mat4) {
        for i in 0..self.n_index {
            apply_mat4(&mut self.coords, i * 3, matrix);
        }
    }

    /// Return a clone with all coordinates transformed by the given matrix
    pub fn transformed(&self, matrix: &Mat4) -> CoordSet {
        let mut cs = self.clone();
        cs.transform(matrix);
        cs
    }

    /// Apply a 4x4 homogeneous transformation matrix to specific atoms
    pub fn transform_atoms(&mut self, atoms: &[AtomIndex], matrix: &Mat4) {
        for &atom in atoms {
            if let Some(coord_idx) = self.atom_to_coord(atom) {
                apply_mat4(&mut self.coords, coord_idx.as_usize() * 3, matrix);
            }
        }
    }

    /// Apply a TTT (Translate-Transform-Translate) matrix
    ///
    /// TTT format is a 16-element array where:
    /// - [0-2, 4-6, 8-10]: 3x3 rotation matrix (row-major)
    /// - [3, 7, 11]: post-rotation translation
    /// - [12, 13, 14]: pre-rotation translation
    /// - [15]: always 1.0
    ///
    /// Transform: y = R * (x + pre_trans) + post_trans
    pub fn transform_ttt(&mut self, ttt: &[f32; 16]) {
        for i in 0..self.n_index {
            apply_ttt(&mut self.coords, i * 3, ttt);
        }
    }

    /// Apply a TTT matrix to specific atoms
    pub fn transform_ttt_atoms(&mut self, atoms: &[AtomIndex], ttt: &[f32; 16]) {
        for &atom in atoms {
            if let Some(coord_idx) = self.atom_to_coord(atom) {
                apply_ttt(&mut self.coords, coord_idx.as_usize() * 3, ttt);
            }
        }
    }
}

// =============================================================================
// Transformation Helpers (private)
// =============================================================================

/// Apply a 4x4 transformation matrix (column-major) to a coordinate at the given base index.
/// This is an inline helper to eliminate code duplication in transform methods.
#[inline]
fn apply_mat4(coords: &mut [f32], base: usize, matrix: &Mat4) {
    let x = coords[base];
    let y = coords[base + 1];
    let z = coords[base + 2];

    // Column-major: data[col*4 + row], so M[row][col] = data[col*4 + row]
    coords[base] =
        matrix.data[0] * x + matrix.data[4] * y + matrix.data[8] * z + matrix.data[12];
    coords[base + 1] =
        matrix.data[1] * x + matrix.data[5] * y + matrix.data[9] * z + matrix.data[13];
    coords[base + 2] =
        matrix.data[2] * x + matrix.data[6] * y + matrix.data[10] * z + matrix.data[14];
}

/// Apply a TTT transformation to a coordinate at the given base index.
/// TTT format: y = R * (x + pre_trans) + post_trans
#[inline]
fn apply_ttt(coords: &mut [f32], base: usize, ttt: &[f32; 16]) {
    // Add pre-translation
    let x = coords[base] + ttt[12];
    let y = coords[base + 1] + ttt[13];
    let z = coords[base + 2] + ttt[14];

    // Apply rotation and post-translation
    coords[base] = ttt[0] * x + ttt[1] * y + ttt[2] * z + ttt[3];
    coords[base + 1] = ttt[4] * x + ttt[5] * y + ttt[6] * z + ttt[7];
    coords[base + 2] = ttt[8] * x + ttt[9] * y + ttt[10] * z + ttt[11];
}

/// 3x3 rotation matrix elements computed using Rodrigues' rotation formula.
/// Used internally by `rotation_matrix` and `rotation_ttt`.
struct RotationElements {
    r00: f32,
    r01: f32,
    r02: f32,
    r10: f32,
    r11: f32,
    r12: f32,
    r20: f32,
    r21: f32,
    r22: f32,
}

/// Compute rotation matrix elements using Rodrigues' rotation formula.
/// Returns None if the axis is degenerate (zero length).
fn rodrigues_rotation(axis: Vec3, angle: f32) -> Option<RotationElements> {
    let len = (axis.x * axis.x + axis.y * axis.y + axis.z * axis.z).sqrt();
    if len < 1e-10 {
        return None;
    }

    let x = axis.x / len;
    let y = axis.y / len;
    let z = axis.z / len;

    let c = angle.cos();
    let s = angle.sin();
    let t = 1.0 - c;

    let xx = x * x;
    let yy = y * y;
    let zz = z * z;
    let xy = x * y;
    let yz = y * z;
    let zx = z * x;
    let xs = x * s;
    let ys = y * s;
    let zs = z * s;

    Some(RotationElements {
        r00: t * xx + c,
        r01: t * xy - zs,
        r02: t * zx + ys,
        r10: t * xy + zs,
        r11: t * yy + c,
        r12: t * yz - xs,
        r20: t * zx - ys,
        r21: t * yz + xs,
        r22: t * zz + c,
    })
}

// =============================================================================
// Transformation Utilities (public)
// =============================================================================

/// Create a rotation matrix for rotation about an axis through an origin
///
/// Uses Rodrigues' rotation formula to construct the rotation matrix,
/// then combines it with translations to/from the origin.
///
/// # Arguments
/// * `axis` - Rotation axis (will be normalized)
/// * `angle` - Rotation angle in radians
/// * `origin` - Center of rotation
///
/// # Returns
/// A 4x4 homogeneous transformation matrix (column-major: data[col*4 + row])
pub fn rotation_matrix(axis: Vec3, angle: f32, origin: Vec3) -> Mat4 {
    let r = match rodrigues_rotation(axis, angle) {
        Some(r) => r,
        None => return Mat4::new_identity(),
    };

    // Combine with translation: T(origin) * R * T(-origin)
    // Result is: R * (v - origin) + origin = R*v + (origin - R*origin)
    let tx = origin.x - (r.r00 * origin.x + r.r01 * origin.y + r.r02 * origin.z);
    let ty = origin.y - (r.r10 * origin.x + r.r11 * origin.y + r.r12 * origin.z);
    let tz = origin.z - (r.r20 * origin.x + r.r21 * origin.y + r.r22 * origin.z);

    // Build matrix in column-major format: data[col*4 + row]
    let mut result = Mat4::new_identity();
    // Column 0
    result.data[0] = r.r00;
    result.data[1] = r.r10;
    result.data[2] = r.r20;
    // Column 1
    result.data[4] = r.r01;
    result.data[5] = r.r11;
    result.data[6] = r.r21;
    // Column 2
    result.data[8] = r.r02;
    result.data[9] = r.r12;
    result.data[10] = r.r22;
    // Column 3 (translation)
    result.data[12] = tx;
    result.data[13] = ty;
    result.data[14] = tz;

    result
}

/// Create a translation matrix
///
/// Creates a column-major 4x4 translation matrix where translation
/// is in column 3: data[12], data[13], data[14].
pub fn translation_matrix(delta: Vec3) -> Mat4 {
    let mut result = Mat4::new_identity();
    result.data[12] = delta.x;
    result.data[13] = delta.y;
    result.data[14] = delta.z;
    result
}

/// Convert a TTT matrix to a standard 4x4 homogeneous matrix (column-major)
///
/// TTT format: pre-translate, rotate, post-translate
/// Standard: just rotation + translation
pub fn ttt_to_mat4(ttt: &[f32; 16]) -> Mat4 {
    // TTT: y = R * (x + pre) + post
    // Standard: y = R * x + t where t = R * pre + post
    let r00 = ttt[0];
    let r01 = ttt[1];
    let r02 = ttt[2];
    let r10 = ttt[4];
    let r11 = ttt[5];
    let r12 = ttt[6];
    let r20 = ttt[8];
    let r21 = ttt[9];
    let r22 = ttt[10];

    let pre_x = ttt[12];
    let pre_y = ttt[13];
    let pre_z = ttt[14];

    let post_x = ttt[3];
    let post_y = ttt[7];
    let post_z = ttt[11];

    // t = R * pre + post
    let tx = r00 * pre_x + r01 * pre_y + r02 * pre_z + post_x;
    let ty = r10 * pre_x + r11 * pre_y + r12 * pre_z + post_y;
    let tz = r20 * pre_x + r21 * pre_y + r22 * pre_z + post_z;

    // Build column-major Mat4: data[col*4 + row]
    Mat4::new([
        r00, r10, r20, 0.0, // col 0
        r01, r11, r21, 0.0, // col 1
        r02, r12, r22, 0.0, // col 2
        tx,  ty,  tz,  1.0, // col 3
    ])
}

/// Convert a standard 4x4 column-major matrix to TTT format (with zero pre-translation)
pub fn mat4_to_ttt(matrix: &Mat4) -> [f32; 16] {
    let d = &matrix.data;
    // Column-major: M[row][col] = d[col*4 + row]
    // TTT row-major layout: [r00, r01, r02, post_x, r10, r11, r12, post_y, ...]
    [
        d[0],  d[4],  d[8],  d[12], // row 0: M[0,0], M[0,1], M[0,2], post_x
        d[1],  d[5],  d[9],  d[13], // row 1: M[1,0], M[1,1], M[1,2], post_y
        d[2],  d[6],  d[10], d[14], // row 2: M[2,0], M[2,1], M[2,2], post_z
        0.0,   0.0,   0.0,   1.0,   // pre-translation (zero)
    ]
}

/// Build a TTT matrix for rotation about an axis through an origin
///
/// TTT format layout:
/// - The rotation matrix is in positions [0-2, 4-6, 8-10]
/// - The origin goes in positions [3, 7, 11] (post-translation)
/// - The negative origin goes in positions [12, 13, 14] (pre-translation)
pub fn rotation_ttt(axis: Vec3, angle: f32, origin: Vec3) -> [f32; 16] {
    let r = match rodrigues_rotation(axis, angle) {
        Some(r) => r,
        None => {
            return [
                1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
            ];
        }
    };

    [
        r.r00, r.r01, r.r02, origin.x, r.r10, r.r11, r.r12, origin.y, r.r20, r.r21, r.r22,
        origin.z, -origin.x, -origin.y, -origin.z, 1.0,
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_coordset_from_coords() {
        let coords = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let cs = CoordSet::from_coords(coords);

        assert_eq!(cs.len(), 2);
        assert_eq!(
            cs.get_coord(CoordIndex(0)),
            Some(Vec3::new(1.0, 2.0, 3.0))
        );
        assert_eq!(
            cs.get_coord(CoordIndex(1)),
            Some(Vec3::new(4.0, 5.0, 6.0))
        );
    }

    #[test]
    fn test_coordset_from_vec3() {
        let positions = vec![Vec3::new(1.0, 2.0, 3.0), Vec3::new(4.0, 5.0, 6.0)];
        let cs = CoordSet::from_vec3(&positions);

        assert_eq!(cs.len(), 2);
        assert_eq!(cs.get_atom_coord(AtomIndex(0)), Some(Vec3::new(1.0, 2.0, 3.0)));
    }

    #[test]
    fn test_coordset_set_coord() {
        let mut cs = CoordSet::from_vec3(&[Vec3::new(0.0, 0.0, 0.0)]);
        cs.set_coord(CoordIndex(0), Vec3::new(1.0, 2.0, 3.0));
        assert_eq!(cs.get_coord(CoordIndex(0)), Some(Vec3::new(1.0, 2.0, 3.0)));
    }

    #[test]
    fn test_coordset_add_coord() {
        let mut cs = CoordSet::new();
        cs.add_coord(AtomIndex(0), Vec3::new(1.0, 2.0, 3.0));
        cs.add_coord(AtomIndex(2), Vec3::new(4.0, 5.0, 6.0));

        assert_eq!(cs.len(), 2);
        assert!(cs.has_atom(AtomIndex(0)));
        assert!(!cs.has_atom(AtomIndex(1)));
        assert!(cs.has_atom(AtomIndex(2)));
    }

    #[test]
    fn test_coordset_center() {
        let positions = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(2.0, 0.0, 0.0),
            Vec3::new(0.0, 2.0, 0.0),
            Vec3::new(0.0, 0.0, 2.0),
        ];
        let cs = CoordSet::from_vec3(&positions);

        let center = cs.center();
        assert!((center.x - 0.5).abs() < 0.001);
        assert!((center.y - 0.5).abs() < 0.001);
        assert!((center.z - 0.5).abs() < 0.001);
    }

    #[test]
    fn test_coordset_bounding_box() {
        let positions = vec![
            Vec3::new(1.0, 2.0, 3.0),
            Vec3::new(4.0, 5.0, 6.0),
        ];
        let cs = CoordSet::from_vec3(&positions);

        let (min, max) = cs.bounding_box().unwrap();
        assert_eq!(min, Vec3::new(1.0, 2.0, 3.0));
        assert_eq!(max, Vec3::new(4.0, 5.0, 6.0));
    }

    #[test]
    fn test_symmetry() {
        let sym = Symmetry::new("P 21 21 21", [10.0, 20.0, 30.0], [90.0, 90.0, 90.0]);
        assert_eq!(sym.space_group, "P 21 21 21");
        assert_eq!(sym.cell_lengths, [10.0, 20.0, 30.0]);
    }

    #[test]
    fn test_rotate_90_degrees_z() {
        let mut cs = CoordSet::from_vec3(&[Vec3::new(1.0, 0.0, 0.0)]);

        // Rotate 90 degrees around Z axis through origin
        let angle = std::f32::consts::FRAC_PI_2; // 90 degrees
        cs.rotate(Vec3::new(0.0, 0.0, 1.0), angle, Vec3::new(0.0, 0.0, 0.0));

        let result = cs.get_coord(CoordIndex(0)).unwrap();
        assert!((result.x - 0.0).abs() < 0.001);
        assert!((result.y - 1.0).abs() < 0.001);
        assert!((result.z - 0.0).abs() < 0.001);
    }

    #[test]
    fn test_rotate_about_origin() {
        let mut cs = CoordSet::from_vec3(&[Vec3::new(2.0, 1.0, 0.0)]);

        // Rotate 90 degrees around Z axis through point (1, 1, 0)
        let angle = std::f32::consts::FRAC_PI_2;
        cs.rotate(Vec3::new(0.0, 0.0, 1.0), angle, Vec3::new(1.0, 1.0, 0.0));

        let result = cs.get_coord(CoordIndex(0)).unwrap();
        // Point (2,1) relative to origin (1,1) is (1,0)
        // After 90 deg rotation: (0, 1)
        // Back to absolute: (1, 2)
        assert!((result.x - 1.0).abs() < 0.001);
        assert!((result.y - 2.0).abs() < 0.001);
    }

    #[test]
    fn test_transform_translation() {
        let mut cs = CoordSet::from_vec3(&[Vec3::new(1.0, 2.0, 3.0)]);

        let matrix = super::translation_matrix(Vec3::new(10.0, 20.0, 30.0));
        cs.transform(&matrix);

        let result = cs.get_coord(CoordIndex(0)).unwrap();
        assert!((result.x - 11.0).abs() < 0.001);
        assert!((result.y - 22.0).abs() < 0.001);
        assert!((result.z - 33.0).abs() < 0.001);
    }

    #[test]
    fn test_transform_ttt() {
        let mut cs = CoordSet::from_vec3(&[Vec3::new(1.0, 0.0, 0.0)]);

        // TTT for 90 degree rotation around Z through origin
        let angle = std::f32::consts::FRAC_PI_2;
        let c = angle.cos();
        let s = angle.sin();
        let ttt: [f32; 16] = [
            c, -s, 0.0, 0.0, // row 0: rotation + post_x
            s, c, 0.0, 0.0, // row 1: rotation + post_y
            0.0, 0.0, 1.0, 0.0, // row 2: rotation + post_z
            0.0, 0.0, 0.0, 1.0, // row 3: pre-translation + 1.0
        ];

        cs.transform_ttt(&ttt);

        let result = cs.get_coord(CoordIndex(0)).unwrap();
        assert!((result.x - 0.0).abs() < 0.001);
        assert!((result.y - 1.0).abs() < 0.001);
    }

    #[test]
    fn test_ttt_to_mat4_identity() {
        let ttt: [f32; 16] = [
            1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
        ];
        let mat = super::ttt_to_mat4(&ttt);

        // Should be identity
        assert!((mat.data[0] - 1.0).abs() < 0.001);
        assert!((mat.data[5] - 1.0).abs() < 0.001);
        assert!((mat.data[10] - 1.0).abs() < 0.001);
    }
}
