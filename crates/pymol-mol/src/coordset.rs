//! Coordinate set data structure
//!
//! Provides the `CoordSet` struct for storing atomic coordinates in a molecular state.
//! Multi-state molecules have multiple coordinate sets, one per state/frame.

use crate::index::{AtomIndex, CoordIndex, INVALID_INDEX};
use lin_alg::f32::Vec3;
use pymol_settings::GlobalSettings;

/// Symmetry information for a coordinate set
#[derive(Debug, Clone)]
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
/// coordinate set. The `idx_to_atm` and `atm_to_idx` mappings handle this.
#[derive(Debug, Clone)]
pub struct CoordSet {
    /// Flat array of coordinates (3 floats per atom: x, y, z)
    coords: Vec<f32>,

    /// Number of atoms with coordinates in this set
    n_index: usize,

    /// Mapping from coordinate index to atom index
    /// idx_to_atm[coord_idx] = atom_idx
    idx_to_atm: Vec<AtomIndex>,

    /// Mapping from atom index to coordinate index
    /// atm_to_idx[atom_idx] = coord_idx (INVALID_INDEX if not present)
    atm_to_idx: Vec<u32>,

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
            idx_to_atm: Vec::new(),
            atm_to_idx: Vec::new(),
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
    pub fn with_capacity(n_atoms: usize) -> Self {
        CoordSet {
            coords: Vec::with_capacity(n_atoms * 3),
            n_index: 0,
            idx_to_atm: Vec::with_capacity(n_atoms),
            atm_to_idx: Vec::new(),
            name: String::new(),
            symmetry: None,
            settings: None,
        }
    }

    /// Create a coordinate set from a flat array of coordinates
    ///
    /// Assumes all atoms have coordinates (1:1 mapping).
    pub fn from_coords(coords: Vec<f32>) -> Self {
        let n_atoms = coords.len() / 3;
        let idx_to_atm: Vec<AtomIndex> = (0..n_atoms as u32).map(AtomIndex).collect();
        let atm_to_idx: Vec<u32> = (0..n_atoms as u32).collect();

        CoordSet {
            coords,
            n_index: n_atoms,
            idx_to_atm,
            atm_to_idx,
            name: String::new(),
            symmetry: None,
            settings: None,
        }
    }

    /// Create a coordinate set from a vector of Vec3 positions
    pub fn from_vec3(positions: &[Vec3]) -> Self {
        let n_atoms = positions.len();
        let mut coords = Vec::with_capacity(n_atoms * 3);
        for pos in positions {
            coords.push(pos.x);
            coords.push(pos.y);
            coords.push(pos.z);
        }

        let idx_to_atm: Vec<AtomIndex> = (0..n_atoms as u32).map(AtomIndex).collect();
        let atm_to_idx: Vec<u32> = (0..n_atoms as u32).collect();

        CoordSet {
            coords,
            n_index: n_atoms,
            idx_to_atm,
            atm_to_idx,
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

    /// Get the number of atoms in the parent molecule
    /// (size of the atm_to_idx mapping)
    #[inline]
    pub fn atom_capacity(&self) -> usize {
        self.atm_to_idx.len()
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
        if atom_idx >= self.atm_to_idx.len() {
            return None;
        }
        let coord_idx = self.atm_to_idx[atom_idx];
        if coord_idx == INVALID_INDEX {
            return None;
        }
        self.get_coord(CoordIndex(coord_idx))
    }

    /// Set coordinates by atom index
    pub fn set_atom_coord(&mut self, atom: AtomIndex, coord: Vec3) -> bool {
        let atom_idx = atom.as_usize();
        if atom_idx >= self.atm_to_idx.len() {
            return false;
        }
        let coord_idx = self.atm_to_idx[atom_idx];
        if coord_idx == INVALID_INDEX {
            return false;
        }
        self.set_coord(CoordIndex(coord_idx), coord);
        true
    }

    /// Check if an atom has coordinates in this set
    #[inline]
    pub fn has_atom(&self, atom: AtomIndex) -> bool {
        let atom_idx = atom.as_usize();
        atom_idx < self.atm_to_idx.len() && self.atm_to_idx[atom_idx] != INVALID_INDEX
    }

    /// Get the atom index for a coordinate index
    #[inline]
    pub fn coord_to_atom(&self, coord: CoordIndex) -> Option<AtomIndex> {
        self.idx_to_atm.get(coord.as_usize()).copied()
    }

    /// Get the coordinate index for an atom index
    #[inline]
    pub fn atom_to_coord(&self, atom: AtomIndex) -> Option<CoordIndex> {
        let atom_idx = atom.as_usize();
        if atom_idx >= self.atm_to_idx.len() {
            return None;
        }
        let coord_idx = self.atm_to_idx[atom_idx];
        if coord_idx == INVALID_INDEX {
            None
        } else {
            Some(CoordIndex(coord_idx))
        }
    }

    /// Get a pointer to the coordinate data for a coordinate index
    #[inline]
    pub fn coord_ptr(&self, idx: CoordIndex) -> Option<&[f32; 3]> {
        let i = idx.as_usize();
        if i >= self.n_index {
            return None;
        }
        let base = i * 3;
        // SAFETY: We checked bounds above, and coords is always a multiple of 3
        Some(unsafe {
            &*(self.coords.as_ptr().add(base) as *const [f32; 3])
        })
    }

    /// Get a mutable pointer to the coordinate data
    #[inline]
    pub fn coord_ptr_mut(&mut self, idx: CoordIndex) -> Option<&mut [f32; 3]> {
        let i = idx.as_usize();
        if i >= self.n_index {
            return None;
        }
        let base = i * 3;
        // SAFETY: We checked bounds above, and coords is always a multiple of 3
        Some(unsafe {
            &mut *(self.coords.as_mut_ptr().add(base) as *mut [f32; 3])
        })
    }

    /// Get the raw coordinate slice
    #[inline]
    pub fn coords_raw(&self) -> &[f32] {
        &self.coords
    }

    /// Get the raw coordinate slice mutably
    #[inline]
    pub fn coords_raw_mut(&mut self) -> &mut [f32] {
        &mut self.coords
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
    pub fn iter_with_atoms(&self) -> impl Iterator<Item = (AtomIndex, Vec3)> + '_ {
        self.idx_to_atm.iter().enumerate().map(move |(i, &atom)| {
            let base = i * 3;
            let coord = Vec3::new(
                self.coords[base],
                self.coords[base + 1],
                self.coords[base + 2],
            );
            (atom, coord)
        })
    }

    /// Add coordinates for an atom
    ///
    /// This extends the coordinate set to include the given atom.
    /// Call `rebuild_mappings` after adding all atoms.
    pub fn add_coord(&mut self, atom: AtomIndex, coord: Vec3) {
        // Extend atm_to_idx if needed
        let atom_idx = atom.as_usize();
        while self.atm_to_idx.len() <= atom_idx {
            self.atm_to_idx.push(INVALID_INDEX);
        }

        // Add the coordinate
        let coord_idx = self.n_index;
        self.coords.push(coord.x);
        self.coords.push(coord.y);
        self.coords.push(coord.z);
        self.idx_to_atm.push(atom);
        self.atm_to_idx[atom_idx] = coord_idx as u32;
        self.n_index += 1;
    }

    /// Compute the center of mass (geometric center)
    pub fn center(&self) -> Vec3 {
        if self.n_index == 0 {
            return Vec3::new(0.0, 0.0, 0.0);
        }

        let mut sum = Vec3::new(0.0, 0.0, 0.0);
        for coord in self.iter() {
            sum = sum + coord;
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
}
