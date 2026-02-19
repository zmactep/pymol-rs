//! Molecular object container
//!
//! Provides the `ObjectMolecule` struct which is the main container for molecular data.
//! It holds atoms, bonds, coordinate sets, and settings.

use serde::{Deserialize, Serialize};
use std::sync::Arc;

use lin_alg::f32::{Mat4, Vec3};
use pymol_settings::{GlobalSettings, UniqueSettings};
use smallvec::SmallVec;

use crate::atom::{Atom, AtomResidue};
use crate::bond::{Bond, BondOrder};
use crate::coordset::{rotation_matrix, CoordSet, Symmetry};
use crate::error::{MolError, MolResult};
use crate::index::{AtomIndex, BondIndex, StateIndex};
use crate::residue::{ChainIterator, ResidueIterator};

// ============================================================================
// Protein Double Bond Lookup Table
// ============================================================================

/// Double bonds in standard amino acid residues.
/// Format: (residue_pattern, atom1, atom2)
/// - "*" matches all residues (for backbone carbonyl)
/// - Multiple residue names separated by "|" for variants
static PROTEIN_DOUBLE_BONDS: &[(&[&str], &str, &str)] = &[
    // Backbone carbonyl C=O (all amino acids)
    (&["*"], "C", "O"),
    // ARG: guanidinium group CZ=NH1
    (&["ARG", "ARGP"], "CZ", "NH1"),
    // ASP: carboxylate CG=OD1
    (&["ASP", "ASPM"], "CG", "OD1"),
    // ASN: amide CG=OD1
    (&["ASN"], "CG", "OD1"),
    // GLU: carboxylate CD=OE1
    (&["GLU", "GLUM"], "CD", "OE1"),
    // GLN: amide CD=OE1
    (&["GLN"], "CD", "OE1"),
    // HIS and variants: imidazole ring (Kekulé alternating pattern)
    (&["HIS", "HID", "HIE", "HIP"], "ND1", "CE1"),
    (&["HIS", "HID", "HIE", "HIP"], "NE2", "CD2"),
    // PHE: benzene ring (Kekulé alternating pattern)
    (&["PHE"], "CG", "CD1"),
    (&["PHE"], "CE1", "CZ"),
    (&["PHE"], "CE2", "CD2"),
    // TYR: benzene ring (same Kekulé pattern as PHE)
    (&["TYR"], "CG", "CD1"),
    (&["TYR"], "CE1", "CZ"),
    (&["TYR"], "CE2", "CD2"),
    // TRP: indole ring system (Kekulé alternating pattern)
    (&["TRP"], "CG", "CD1"),
    (&["TRP"], "CE2", "CZ2"),
    (&["TRP"], "CH2", "CZ3"),
    (&["TRP"], "CE3", "CD2"),
];

/// Get bond order for a bond between two atoms in a known protein residue
///
/// Based on PyMOL's assign_pdb_known_residue() function.
/// Uses atom names and residue names to determine bond orders.
fn get_protein_bond_order(resn: &str, name1: &str, name2: &str) -> BondOrder {
    // Helper to check if a bond connects two specific atoms (in either order)
    let is_bond = |a: &str, b: &str| -> bool {
        (name1 == a && name2 == b) || (name1 == b && name2 == a)
    };

    // Look up in the double bonds table
    for &(residues, atom1, atom2) in PROTEIN_DOUBLE_BONDS {
        if is_bond(atom1, atom2) {
            // Check if residue matches
            let matches = residues.iter().any(|&r| r == "*" || r == resn);
            if matches {
                return BondOrder::Double;
            }
        }
    }

    BondOrder::Single
}

/// A molecular object containing atoms, bonds, and coordinate sets
///
/// This is the Rust equivalent of PyMOL's `ObjectMolecule`.
/// It provides a container for all molecular data including:
/// - Atom properties (name, element, residue, etc.)
/// - Bond connectivity
/// - Coordinates for one or more states (conformations)
/// - Per-object and per-atom/bond settings
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ObjectMolecule {
    // =========================================================================
    // Core Data
    // =========================================================================
    /// All atoms in the molecule
    atoms: Vec<Atom>,

    /// All bonds in the molecule
    bonds: Vec<Bond>,

    /// Coordinate sets (one per state/frame)
    coord_sets: Vec<CoordSet>,

    // =========================================================================
    // Bond Lookup
    // =========================================================================
    /// Map from atom index to list of bond indices involving that atom
    /// Uses SmallVec to avoid heap allocation for atoms with 4 or fewer bonds (common case)
    /// Skipped during serialization — rebuilt from bonds via `rebuild_atom_bonds()`.
    #[serde(skip)]
    atom_bonds: Vec<SmallVec<[BondIndex; 4]>>,

    // =========================================================================
    // Metadata
    // =========================================================================
    /// Object name
    pub name: String,

    /// Title or description
    pub title: String,

    /// Current state index (0-based)
    pub current_state: usize,

    // =========================================================================
    // Discrete Mode
    // =========================================================================
    /// Whether this is a discrete object
    /// In discrete mode, each state has independent atoms (not shared)
    pub discrete: bool,

    // =========================================================================
    // Settings
    // =========================================================================
    /// Per-object settings
    pub settings: Option<GlobalSettings>,

    /// Per-atom/bond unique settings
    pub unique_settings: UniqueSettings,

    // =========================================================================
    // Symmetry
    // =========================================================================
    /// Crystallographic symmetry (applies to all states unless overridden)
    pub symmetry: Option<Symmetry>,
}

impl Default for ObjectMolecule {
    fn default() -> Self {
        ObjectMolecule {
            atoms: Vec::new(),
            bonds: Vec::new(),
            coord_sets: Vec::new(),
            atom_bonds: Vec::new(),
            name: String::new(),
            title: String::new(),
            current_state: 0,
            discrete: false,
            settings: None,
            unique_settings: UniqueSettings::new(),
            symmetry: None,
        }
    }
}

impl ObjectMolecule {
    /// Create a new empty molecule with the given name
    pub fn new(name: impl Into<String>) -> Self {
        ObjectMolecule {
            name: name.into(),
            ..Default::default()
        }
    }

    /// Rebuild the `atom_bonds` lookup table from the current bond list.
    ///
    /// This must be called after deserialization since `atom_bonds` is skipped
    /// during serialization (it is derived data).
    pub fn rebuild_atom_bonds(&mut self) {
        self.atom_bonds = vec![SmallVec::new(); self.atoms.len()];
        for (bi, bond) in self.bonds.iter().enumerate() {
            let bond_idx = BondIndex(bi as u32);
            if bond.atom1.as_usize() < self.atoms.len() {
                self.atom_bonds[bond.atom1.as_usize()].push(bond_idx);
            }
            if bond.atom2.as_usize() < self.atoms.len() {
                self.atom_bonds[bond.atom2.as_usize()].push(bond_idx);
            }
        }
    }

    /// Create a molecule with pre-allocated capacity
    pub fn with_capacity(name: impl Into<String>, atoms: usize, bonds: usize) -> Self {
        ObjectMolecule {
            atoms: Vec::with_capacity(atoms),
            bonds: Vec::with_capacity(bonds),
            coord_sets: Vec::new(),
            atom_bonds: Vec::with_capacity(atoms),
            name: name.into(),
            title: String::new(),
            current_state: 0,
            discrete: false,
            settings: None,
            unique_settings: UniqueSettings::new(),
            symmetry: None,
        }
    }

    // =========================================================================
    // Atom Operations
    // =========================================================================

    /// Add an atom to the molecule
    ///
    /// Returns the index of the newly added atom.
    pub fn add_atom(&mut self, atom: Atom) -> AtomIndex {
        let index = AtomIndex(self.atoms.len() as u32);
        self.atoms.push(atom);
        self.atom_bonds.push(SmallVec::new());
        index
    }

    /// Get an atom by index
    #[inline]
    pub fn get_atom(&self, index: AtomIndex) -> Option<&Atom> {
        self.atoms.get(index.as_usize())
    }

    /// Get a mutable reference to an atom by index
    #[inline]
    pub fn get_atom_mut(&mut self, index: AtomIndex) -> Option<&mut Atom> {
        self.atoms.get_mut(index.as_usize())
    }

    /// Get the number of atoms
    #[inline]
    pub fn atom_count(&self) -> usize {
        self.atoms.len()
    }

    /// Iterate over all atoms
    pub fn atoms(&self) -> impl Iterator<Item = &Atom> {
        self.atoms.iter()
    }

    /// Iterate over all atoms with indices
    pub fn atoms_indexed(&self) -> impl Iterator<Item = (AtomIndex, &Atom)> {
        self.atoms
            .iter()
            .enumerate()
            .map(|(i, a)| (AtomIndex(i as u32), a))
    }

    /// Get the atoms slice
    pub fn atoms_slice(&self) -> &[Atom] {
        &self.atoms
    }

    /// Get mutable atoms slice
    pub fn atoms_slice_mut(&mut self) -> &mut [Atom] {
        &mut self.atoms
    }

    /// Iterate over all atoms mutably
    pub fn atoms_mut(&mut self) -> impl Iterator<Item = &mut Atom> {
        self.atoms.iter_mut()
    }

    // =========================================================================
    // Bond Operations
    // =========================================================================

    /// Add a bond between two atoms
    ///
    /// Returns the index of the newly added bond.
    /// Returns an error if the atoms don't exist or a bond already exists.
    pub fn add_bond(
        &mut self,
        atom1: AtomIndex,
        atom2: AtomIndex,
        order: BondOrder,
    ) -> MolResult<BondIndex> {
        // Validate atoms exist
        let n_atoms = self.atoms.len();
        if atom1.as_usize() >= n_atoms {
            return Err(MolError::atom_out_of_bounds(atom1.0, n_atoms));
        }
        if atom2.as_usize() >= n_atoms {
            return Err(MolError::atom_out_of_bounds(atom2.0, n_atoms));
        }

        // Don't allow self-loops
        if atom1 == atom2 {
            return Err(MolError::InvalidBond(atom1.0, atom2.0));
        }

        // Check for duplicate bond
        if self.find_bond(atom1, atom2).is_some() {
            return Err(MolError::DuplicateBond(atom1.0, atom2.0));
        }

        // Create the bond
        let bond = Bond::new(atom1, atom2, order);
        let bond_index = BondIndex(self.bonds.len() as u32);
        self.bonds.push(bond);

        // Update atom-to-bond mapping
        self.atom_bonds[atom1.as_usize()].push(bond_index);
        self.atom_bonds[atom2.as_usize()].push(bond_index);

        // Mark atoms as bonded
        if let Some(a1) = self.atoms.get_mut(atom1.as_usize()) {
            a1.state.bonded = true;
        }
        if let Some(a2) = self.atoms.get_mut(atom2.as_usize()) {
            a2.state.bonded = true;
        }

        Ok(bond_index)
    }

    /// Add a bond (simplified version that ignores errors for duplicate bonds)
    pub fn add_bond_unchecked(&mut self, atom1: AtomIndex, atom2: AtomIndex, order: BondOrder) -> Option<BondIndex> {
        self.add_bond(atom1, atom2, order).ok()
    }

    /// Get a bond by index
    #[inline]
    pub fn get_bond(&self, index: BondIndex) -> Option<&Bond> {
        self.bonds.get(index.as_usize())
    }

    /// Get a mutable reference to a bond by index
    #[inline]
    pub fn get_bond_mut(&mut self, index: BondIndex) -> Option<&mut Bond> {
        self.bonds.get_mut(index.as_usize())
    }

    /// Find a bond between two atoms
    pub fn find_bond(&self, atom1: AtomIndex, atom2: AtomIndex) -> Option<BondIndex> {
        // Use the smaller atom's bond list for efficiency
        let smaller = if self.atom_bonds.get(atom1.as_usize())?.len()
            < self.atom_bonds.get(atom2.as_usize())?.len()
        {
            atom1
        } else {
            atom2
        };

        for &bond_idx in &self.atom_bonds[smaller.as_usize()] {
            if let Some(bond) = self.bonds.get(bond_idx.as_usize()) {
                if bond.connects(atom1, atom2) {
                    return Some(bond_idx);
                }
            }
        }
        None
    }

    /// Get all atoms bonded to a given atom
    pub fn bonded_atoms(&self, atom: AtomIndex) -> Vec<AtomIndex> {
        let mut result = Vec::new();
        if let Some(bond_indices) = self.atom_bonds.get(atom.as_usize()) {
            for &bond_idx in bond_indices {
                if let Some(bond) = self.bonds.get(bond_idx.as_usize()) {
                    if let Some(other) = bond.other(atom) {
                        result.push(other);
                    }
                }
            }
        }
        result
    }

    /// Get the bond indices for a given atom
    pub fn atom_bond_indices(&self, atom: AtomIndex) -> &[BondIndex] {
        self.atom_bonds
            .get(atom.as_usize())
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }

    /// Get the number of bonds
    #[inline]
    pub fn bond_count(&self) -> usize {
        self.bonds.len()
    }

    /// Iterate over all bonds
    pub fn bonds(&self) -> impl Iterator<Item = &Bond> {
        self.bonds.iter()
    }

    /// Iterate over all bonds with indices
    pub fn bonds_indexed(&self) -> impl Iterator<Item = (BondIndex, &Bond)> {
        self.bonds
            .iter()
            .enumerate()
            .map(|(i, b)| (BondIndex(i as u32), b))
    }

    /// Generate bonds based on inter-atomic distances
    ///
    /// Creates bonds between atoms that are within bonding distance based on
    /// their VdW radii. A bond is created if:
    ///   distance < (vdw1 + vdw2) * tolerance
    ///
    /// Typical tolerance values:
    /// - 0.6 for standard covalent bonds
    /// - 0.45 for very strict bonding (peptide bonds, etc.)
    ///
    /// This method operates on the first coordinate set (state 0).
    pub fn generate_bonds(&mut self, tolerance: f32) {
        self.generate_bonds_for_state(tolerance, 0);
    }

    /// Generate bonds for a specific coordinate state
    ///
    /// Creates bonds between atoms based on distance criteria.
    /// Uses a spatial hash grid for O(n) performance instead of O(n²) all-pairs.
    pub fn generate_bonds_for_state(&mut self, tolerance: f32, state: usize) {
        use crate::spatial::SpatialGrid;

        let coord_set = match self.coord_sets.get(state) {
            Some(cs) => cs,
            None => return,
        };

        let n_atoms = self.atoms.len();
        if n_atoms < 2 {
            return;
        }

        // Collect atom data for bond determination
        struct AtomData {
            coord: Vec3,
            vdw: f32,
            is_hydrogen: bool,
        }

        let atom_data: Vec<AtomData> = (0..n_atoms)
            .map(|i| {
                let coord = coord_set
                    .get_atom_coord(AtomIndex(i as u32))
                    .unwrap_or(Vec3::new(0.0, 0.0, 0.0));
                let atom = &self.atoms[i];
                AtomData {
                    coord,
                    vdw: atom.effective_vdw(),
                    is_hydrogen: atom.is_hydrogen(),
                }
            })
            .collect();

        // Compute the maximum possible bonding cutoff to use as cell size.
        // This ensures all potential bonds are found within the 3×3×3 neighborhood.
        let max_vdw = atom_data
            .iter()
            .map(|a| a.vdw)
            .fold(0.0f32, f32::max);
        let max_cutoff = 2.0 * max_vdw * tolerance;

        // Build spatial hash grid — O(n)
        let mut grid = SpatialGrid::with_capacity(max_cutoff, n_atoms);
        for (i, ad) in atom_data.iter().enumerate() {
            grid.insert(ad.coord, i);
        }

        // Find bonded pairs using spatial grid — O(n·k) where k is avg neighbors per cell
        let mut new_bonds: Vec<(AtomIndex, AtomIndex)> = Vec::new();
        let mut neighbors = Vec::new();

        for i in 0..n_atoms {
            let a1 = &atom_data[i];
            grid.query_neighbors(a1.coord, &mut neighbors);

            for &j in &neighbors {
                // Only check j > i to avoid duplicate pairs
                if j <= i {
                    continue;
                }

                let a2 = &atom_data[j];

                // Calculate distance squared
                let dx = a1.coord.x - a2.coord.x;
                let dy = a1.coord.y - a2.coord.y;
                let dz = a1.coord.z - a2.coord.z;
                let dist_sq = dx * dx + dy * dy + dz * dz;

                // Calculate bonding threshold
                // Use smaller tolerance for H-H bonds (rare)
                let effective_tolerance = if a1.is_hydrogen && a2.is_hydrogen {
                    tolerance * 0.5
                } else {
                    tolerance
                };

                let threshold = (a1.vdw + a2.vdw) * effective_tolerance;
                let threshold_sq = threshold * threshold;

                if dist_sq < threshold_sq {
                    new_bonds.push((AtomIndex(i as u32), AtomIndex(j as u32)));
                }
            }
        }

        // Bulk-insert bonds without duplicate checking (j > i guarantees uniqueness)
        for (a1, a2) in new_bonds {
            let bond = Bond::new(a1, a2, BondOrder::Single);
            let bond_index = BondIndex(self.bonds.len() as u32);
            self.bonds.push(bond);
            self.atom_bonds[a1.as_usize()].push(bond_index);
            self.atom_bonds[a2.as_usize()].push(bond_index);
            self.atoms[a1.as_usize()].state.bonded = true;
            self.atoms[a2.as_usize()].state.bonded = true;
        }
    }

    /// Clear all existing bonds
    pub fn clear_bonds(&mut self) {
        self.bonds.clear();
        for atom_bonds in &mut self.atom_bonds {
            atom_bonds.clear();
        }
        for atom in &mut self.atoms {
            atom.state.bonded = false;
        }
    }

    /// Assign bond orders for known protein residues based on atom names
    ///
    /// This function uses PDB atom name conventions to assign bond orders for
    /// standard amino acids. Unlike distance-based inference, this method is
    /// chemically accurate because it uses the known structure of amino acids.
    ///
    /// Bond orders assigned:
    /// - Backbone C=O (carbonyl): all amino acids
    /// - ARG: CZ=NH1 (guanidinium)
    /// - ASP: CG=OD1 (carboxylate)
    /// - ASN: CG=OD1 (amide)
    /// - GLU: CD=OE1 (carboxylate)
    /// - GLN: CD=OE1 (amide)
    pub fn assign_known_residue_bond_orders(&mut self) {
        use crate::residue::{atoms_same_residue, is_amino_acid};

        for bond in &mut self.bonds {
            let atom1 = match self.atoms.get(bond.atom1.as_usize()) {
                Some(a) => a,
                None => continue,
            };
            let atom2 = match self.atoms.get(bond.atom2.as_usize()) {
                Some(a) => a,
                None => continue,
            };

            // Only assign bond orders within the same residue
            if !atoms_same_residue(atom1, atom2) {
                continue;
            }

            // Only for amino acid residues
            if !is_amino_acid(&atom1.residue.resn) {
                continue;
            }

            let name1 = &*atom1.name;
            let name2 = &*atom2.name;
            let resn = atom1.residue.resn.as_str();

            // Check for double bonds based on atom names
            let order = get_protein_bond_order(resn, name1, name2);
            if order != BondOrder::Single {
                bond.order = order;
            }
        }
    }

    // =========================================================================
    // Coordinate Set Operations
    // =========================================================================

    /// Add a coordinate set
    ///
    /// Returns the state index of the newly added coordinate set.
    pub fn add_coord_set(&mut self, coord_set: CoordSet) -> StateIndex {
        let index = StateIndex(self.coord_sets.len() as u32);
        self.coord_sets.push(coord_set);
        index
    }

    /// Get a coordinate set by state index
    #[inline]
    pub fn get_coord_set(&self, state: usize) -> Option<&CoordSet> {
        self.coord_sets.get(state)
    }

    /// Get a mutable reference to a coordinate set
    #[inline]
    pub fn get_coord_set_mut(&mut self, state: usize) -> Option<&mut CoordSet> {
        self.coord_sets.get_mut(state)
    }

    /// Get the current coordinate set
    pub fn current_coord_set(&self) -> Option<&CoordSet> {
        self.coord_sets.get(self.current_state)
    }

    /// Get the number of states (coordinate sets)
    #[inline]
    pub fn state_count(&self) -> usize {
        self.coord_sets.len()
    }

    /// Check if the molecule has any coordinates
    #[inline]
    pub fn has_coords(&self) -> bool {
        !self.coord_sets.is_empty()
    }

    // =========================================================================
    // Coordinate Access (Convenience)
    // =========================================================================

    /// Get coordinates for an atom in a specific state
    pub fn get_coord(&self, atom: AtomIndex, state: usize) -> Option<Vec3> {
        self.coord_sets.get(state)?.get_atom_coord(atom)
    }

    /// Set coordinates for an atom in a specific state
    pub fn set_coord(&mut self, atom: AtomIndex, state: usize, coord: Vec3) -> bool {
        if let Some(cs) = self.coord_sets.get_mut(state) {
            cs.set_atom_coord(atom, coord)
        } else {
            false
        }
    }

    /// Get coordinates for an atom in the current state
    pub fn get_current_coord(&self, atom: AtomIndex) -> Option<Vec3> {
        self.get_coord(atom, self.current_state)
    }

    // =========================================================================
    // Residue/Chain Iteration
    // =========================================================================

    /// Iterate over residues in the molecule
    pub fn residues(&self) -> ResidueIterator<'_> {
        ResidueIterator::new(&self.atoms, 0)
    }

    /// Iterate over chains in the molecule
    pub fn chains(&self) -> ChainIterator<'_> {
        ChainIterator::new(&self.atoms, 0)
    }

    /// Count the number of residues
    pub fn residue_count(&self) -> usize {
        self.residues().count()
    }

    /// Count the number of chains
    pub fn chain_count(&self) -> usize {
        self.chains().count()
    }

    /// Get unique chain IDs
    pub fn chain_ids(&self) -> Vec<String> {
        let mut ids: Vec<String> = self.atoms.iter().map(|a| a.residue.chain.clone()).collect();
        ids.sort();
        ids.dedup();
        ids
    }

    // =========================================================================
    // Analysis
    // =========================================================================

    /// Calculate the center of mass (geometric center) for a state
    pub fn center(&self, state: usize) -> Option<Vec3> {
        self.coord_sets.get(state).map(|cs| cs.center())
    }

    /// Calculate the bounding box for a state
    pub fn bounding_box(&self, state: usize) -> Option<(Vec3, Vec3)> {
        self.coord_sets.get(state)?.bounding_box()
    }

    /// Get the maximum VdW radius of any atom
    pub fn max_vdw(&self) -> f32 {
        self.atoms
            .iter()
            .map(|a| a.effective_vdw())
            .fold(0.0f32, f32::max)
    }

    // =========================================================================
    // Atom Classification
    // =========================================================================

    /// Classify atoms as protein, nucleic acid, solvent, organic, or inorganic
    ///
    /// This method sets `AtomFlags` on atoms based on their residue names and
    /// other properties. It follows PyMOL's `SelectorClassifyAtoms` logic:
    ///
    /// - Known amino acid residues (non-HETATM) → `PROTEIN | POLYMER`
    /// - Known nucleic acid residues (non-HETATM) → `NUCLEIC | POLYMER`
    /// - Water residues → `SOLVENT`
    /// - CA atoms in protein residues → `GUIDE` flag (for cartoon)
    ///
    /// This method should be called after loading atoms from a file or
    /// after programmatic molecule creation to enable cartoon representation.
    pub fn classify_atoms(&mut self) {
        use crate::flags::AtomFlags;
        use crate::residue::{is_amino_acid, is_nucleotide, is_water};

        // We need to collect residue information first since we can't borrow
        // self.atoms mutably while iterating through residues
        let residue_info: Vec<(std::ops::Range<usize>, AtomFlags, Option<usize>)> = self
            .residues()
            .map(|residue| {
                let resn = residue.resn();
                let is_hetatm = residue
                    .atoms
                    .first()
                    .map(|a| a.state.hetatm)
                    .unwrap_or(true);

                // Determine classification flags
                // Note: nucleotides are recognised even with HETATM (common for
                // modified bases), but proteins keep the !is_hetatm guard because
                // small peptides can legitimately be HETATM ligands.
                let mask = if !is_hetatm && is_amino_acid(resn) {
                    AtomFlags::POLYMER | AtomFlags::PROTEIN
                } else if is_nucleotide(resn) {
                    AtomFlags::POLYMER | AtomFlags::NUCLEIC
                } else if is_water(resn) {
                    AtomFlags::SOLVENT
                } else {
                    // Non-polymer atoms: ligands, ions, etc.
                    // Classify based on content (carbon → organic, otherwise inorganic)
                    if residue.atoms.iter().any(|a| a.element.is_carbon()) {
                        AtomFlags::ORGANIC
                    } else {
                        AtomFlags::INORGANIC
                    }
                };

                // Find guide atom (CA for proteins, C4' for nucleic acids)
                let guide_atom_local_idx = if mask.contains(AtomFlags::PROTEIN) {
                    // Find CA atom
                    residue.atoms.iter().position(|a| &*a.name == "CA")
                } else if mask.contains(AtomFlags::NUCLEIC) {
                    // Find C4' or C4* atom
                    residue
                        .atoms
                        .iter()
                        .position(|a| &*a.name == "C4'" || &*a.name == "C4*")
                } else {
                    None
                };

                // Convert local index to global index
                let guide_atom_global_idx =
                    guide_atom_local_idx.map(|local| residue.atom_range.start + local);

                (residue.atom_range, mask, guide_atom_global_idx)
            })
            .collect();

        // Apply flags to atoms
        for (range, mask, guide_idx) in residue_info {
            let is_non_polymer = !mask.contains(AtomFlags::POLYMER);
            for idx in range {
                if let Some(atom) = self.atoms.get_mut(idx) {
                    atom.state.flags |= mask;
                    // Mark non-polymer atoms as hetatm for proper representation
                    // (mirrors PyMOL's need_hetatm_classification post-processing)
                    if is_non_polymer {
                        atom.state.hetatm = true;
                    }
                }
            }
            // Set guide flag on the guide atom
            if let Some(guide) = guide_idx {
                if let Some(atom) = self.atoms.get_mut(guide) {
                    atom.state.flags |= AtomFlags::GUIDE;
                }
            }
        }
    }

    // =========================================================================
    // Chain Assignment
    // =========================================================================

    /// Assign chain IDs to atoms that have no chain information
    ///
    /// This is designed for formats like GRO that lack chain identifiers.
    /// The algorithm classifies residues by type (protein, nucleic, solvent,
    /// ion, lipid, other) and assigns sequential chain letters.
    ///
    /// Chain breaks occur when:
    /// - Residue category changes
    /// - Residue number decreases (GROMACS wraps at 99999)
    /// - No backbone bond between consecutive polymer residues
    ///
    /// All solvent residues share one chain. All ions share one chain.
    ///
    /// Must be called after `classify_atoms()` and `generate_bonds()`.
    pub fn assign_chains(&mut self) {
        use crate::residue::{classify_residue, ResidueCategory};

        // Only act if all atoms have empty chain IDs
        if self.atoms.iter().any(|a| !a.residue.chain.is_empty()) {
            return;
        }

        // Phase 1: Collect residue info
        let residues: Vec<(std::ops::Range<usize>, ResidueCategory, i32)> = self
            .residues()
            .map(|r| {
                let cat = classify_residue(r.resn());
                let resv = r.resv();
                (r.atom_range, cat, resv)
            })
            .collect();

        if residues.is_empty() {
            return;
        }

        // Phase 2: Assign chain labels
        let mut chain_labels: Vec<String> = vec![String::new(); residues.len()];
        let mut chain_gen = ChainIdGenerator::new();

        let mut solvent_indices: Vec<usize> = Vec::new();
        let mut ion_indices: Vec<usize> = Vec::new();

        let mut i = 0;
        while i < residues.len() {
            let cat = residues[i].1;

            match cat {
                ResidueCategory::Protein | ResidueCategory::Nucleic => {
                    let chain_id = chain_gen.next();
                    chain_labels[i] = chain_id.clone();

                    let mut j = i + 1;
                    while j < residues.len() {
                        let next_cat = residues[j].1;
                        let next_resv = residues[j].2;
                        let prev_resv = residues[j - 1].2;

                        // Break on category change
                        if next_cat != cat {
                            break;
                        }
                        // Break on resnum decrease (GROMACS wrap)
                        if next_resv < prev_resv {
                            break;
                        }
                        // Break if no backbone bond
                        if !self.has_backbone_bond(&residues[j - 1].0, &residues[j].0, cat) {
                            break;
                        }

                        chain_labels[j] = chain_id.clone();
                        j += 1;
                    }
                    i = j;
                }
                ResidueCategory::Solvent => {
                    solvent_indices.push(i);
                    i += 1;
                }
                ResidueCategory::Ion => {
                    ion_indices.push(i);
                    i += 1;
                }
                ResidueCategory::Lipid => {
                    let chain_id = chain_gen.next();
                    chain_labels[i] = chain_id.clone();

                    let mut j = i + 1;
                    while j < residues.len() {
                        if residues[j].1 != ResidueCategory::Lipid {
                            break;
                        }
                        if residues[j].2 < residues[j - 1].2 {
                            break;
                        }
                        chain_labels[j] = chain_id.clone();
                        j += 1;
                    }
                    i = j;
                }
                ResidueCategory::Other => {
                    let chain_id = chain_gen.next();
                    chain_labels[i] = chain_id.clone();

                    let mut j = i + 1;
                    while j < residues.len() {
                        if residues[j].1 != ResidueCategory::Other {
                            break;
                        }
                        if residues[j].2 < residues[j - 1].2 {
                            break;
                        }
                        chain_labels[j] = chain_id.clone();
                        j += 1;
                    }
                    i = j;
                }
            }
        }

        // Assign shared chain for solvent and ions
        if !solvent_indices.is_empty() {
            let solvent_chain = chain_gen.next();
            for &idx in &solvent_indices {
                chain_labels[idx] = solvent_chain.clone();
            }
        }
        if !ion_indices.is_empty() {
            let ion_chain = chain_gen.next();
            for &idx in &ion_indices {
                chain_labels[idx] = ion_chain.clone();
            }
        }

        // Phase 3: Apply chain labels to atoms
        let residue_ranges: Vec<std::ops::Range<usize>> =
            self.residues().map(|r| r.atom_range).collect();

        for (res_idx, range) in residue_ranges.iter().enumerate() {
            let chain = &chain_labels[res_idx];
            if chain.is_empty() {
                continue;
            }

            let first_atom = &self.atoms[range.start];
            let new_residue = Arc::new(AtomResidue::from_parts(
                chain.as_str(),
                &first_atom.residue.resn,
                first_atom.residue.resv,
                first_atom.residue.inscode,
                &first_atom.residue.segi,
            ));

            for atom_idx in range.clone() {
                self.atoms[atom_idx].residue = new_residue.clone();
            }
        }
    }

    /// Check if there's a backbone bond between two consecutive residues
    fn has_backbone_bond(
        &self,
        prev_range: &std::ops::Range<usize>,
        curr_range: &std::ops::Range<usize>,
        cat: crate::residue::ResidueCategory,
    ) -> bool {
        use crate::residue::ResidueCategory;

        let (prev_name, curr_name) = match cat {
            ResidueCategory::Protein => ("C", "N"),
            ResidueCategory::Nucleic => ("O3'", "P"),
            _ => return false,
        };

        let prev_atom_idx = self.atoms[prev_range.clone()]
            .iter()
            .enumerate()
            .find(|(_, a)| &*a.name == prev_name)
            .map(|(local, _)| AtomIndex((prev_range.start + local) as u32));

        let curr_atom_idx = self.atoms[curr_range.clone()]
            .iter()
            .enumerate()
            .find(|(_, a)| &*a.name == curr_name)
            .map(|(local, _)| AtomIndex((curr_range.start + local) as u32));

        match (prev_atom_idx, curr_atom_idx) {
            (Some(p), Some(c)) => self.find_bond(p, c).is_some(),
            _ => false,
        }
    }

    // =========================================================================
    // Transformation
    // =========================================================================

    /// Translate all coordinates in a state
    pub fn translate(&mut self, state: usize, delta: Vec3) {
        if let Some(cs) = self.coord_sets.get_mut(state) {
            cs.translate(delta);
        }
    }

    /// Center coordinates at the origin for a state
    pub fn center_origin(&mut self, state: usize) {
        if let Some(cs) = self.coord_sets.get_mut(state) {
            cs.center_origin();
        }
    }

    /// Translate specific atoms in a state
    pub fn translate_atoms(&mut self, state: usize, atoms: &[AtomIndex], delta: Vec3) {
        if let Some(cs) = self.coord_sets.get_mut(state) {
            cs.translate_atoms(atoms, delta);
        }
    }

    /// Translate all atoms in all states
    pub fn translate_all_states(&mut self, delta: Vec3) {
        for cs in &mut self.coord_sets {
            cs.translate(delta);
        }
    }

    /// Rotate all coordinates in a state
    ///
    /// # Arguments
    /// * `state` - State index (0-based)
    /// * `axis` - Rotation axis (will be normalized)
    /// * `angle` - Rotation angle in radians
    /// * `origin` - Center of rotation
    pub fn rotate(&mut self, state: usize, axis: Vec3, angle: f32, origin: Vec3) {
        if let Some(cs) = self.coord_sets.get_mut(state) {
            cs.rotate(axis, angle, origin);
        }
    }

    /// Rotate specific atoms in a state
    pub fn rotate_atoms(
        &mut self,
        state: usize,
        atoms: &[AtomIndex],
        axis: Vec3,
        angle: f32,
        origin: Vec3,
    ) {
        if let Some(cs) = self.coord_sets.get_mut(state) {
            cs.rotate_atoms(atoms, axis, angle, origin);
        }
    }

    /// Rotate all atoms in all states
    pub fn rotate_all_states(&mut self, axis: Vec3, angle: f32, origin: Vec3) {
        let matrix = rotation_matrix(axis, angle, origin);
        for cs in &mut self.coord_sets {
            cs.transform(&matrix);
        }
    }

    /// Apply a 4x4 transformation matrix to all coordinates in a state
    pub fn transform(&mut self, state: usize, matrix: &Mat4) {
        if let Some(cs) = self.coord_sets.get_mut(state) {
            cs.transform(matrix);
        }
    }

    /// Apply a 4x4 transformation matrix to specific atoms in a state
    pub fn transform_atoms(&mut self, state: usize, atoms: &[AtomIndex], matrix: &Mat4) {
        if let Some(cs) = self.coord_sets.get_mut(state) {
            cs.transform_atoms(atoms, matrix);
        }
    }

    /// Apply a 4x4 transformation matrix to all atoms in all states
    pub fn transform_all_states(&mut self, matrix: &Mat4) {
        for cs in &mut self.coord_sets {
            cs.transform(matrix);
        }
    }

    /// Apply a PyMOL-style TTT matrix to all coordinates in a state
    pub fn transform_ttt(&mut self, state: usize, ttt: &[f32; 16]) {
        if let Some(cs) = self.coord_sets.get_mut(state) {
            cs.transform_ttt(ttt);
        }
    }

    /// Apply a PyMOL-style TTT matrix to specific atoms in a state
    pub fn transform_ttt_atoms(&mut self, state: usize, atoms: &[AtomIndex], ttt: &[f32; 16]) {
        if let Some(cs) = self.coord_sets.get_mut(state) {
            cs.transform_ttt_atoms(atoms, ttt);
        }
    }

    /// Get the geometric center of the molecule for a state
    #[deprecated(since = "0.2.0", note = "Use center() instead")]
    pub fn get_center(&self, state: usize) -> Option<Vec3> {
        self.center(state)
    }

    /// Get the geometric center of specific atoms for a state
    pub fn get_atoms_center(&self, state: usize, atoms: &[AtomIndex]) -> Option<Vec3> {
        let cs = self.coord_sets.get(state)?;
        if atoms.is_empty() {
            return Some(Vec3::new(0.0, 0.0, 0.0));
        }

        let mut sum = Vec3::new(0.0, 0.0, 0.0);
        let mut count = 0;
        for &atom in atoms {
            if let Some(coord) = cs.get_atom_coord(atom) {
                sum += coord;
                count += 1;
            }
        }

        if count > 0 {
            Some(sum * (1.0 / count as f32))
        } else {
            None
        }
    }

    // =========================================================================
    // Atom Removal
    // =========================================================================

    /// Remove atoms by indices, reindexing remaining atoms and bonds.
    ///
    /// This removes the specified atoms, any bonds involving those atoms,
    /// and updates all coordinate sets. Remaining atom/bond indices are
    /// compacted so they remain contiguous.
    pub fn remove_atoms(&mut self, indices: &[AtomIndex]) {
        if indices.is_empty() {
            return;
        }

        let n_atoms = self.atoms.len();

        // Build a set of atoms to remove for O(1) lookup
        let mut remove_set = vec![false; n_atoms];
        for &idx in indices {
            if idx.as_usize() < n_atoms {
                remove_set[idx.as_usize()] = true;
            }
        }

        // Build old→new index mapping (-1 means removed)
        let mut old_to_new: Vec<Option<u32>> = vec![None; n_atoms];
        let mut new_idx = 0u32;
        for old_idx in 0..n_atoms {
            if !remove_set[old_idx] {
                old_to_new[old_idx] = Some(new_idx);
                new_idx += 1;
            }
        }

        // Filter atoms (keep only non-removed)
        let mut new_atoms = Vec::with_capacity(new_idx as usize);
        for (i, atom) in self.atoms.drain(..).enumerate() {
            if !remove_set[i] {
                new_atoms.push(atom);
            }
        }
        self.atoms = new_atoms;

        // Filter bonds: keep only bonds where both atoms survived, and remap indices
        let mut new_bonds = Vec::new();
        for bond in &self.bonds {
            if let (Some(new_a1), Some(new_a2)) = (
                old_to_new[bond.atom1.as_usize()],
                old_to_new[bond.atom2.as_usize()],
            ) {
                let mut new_bond = bond.clone();
                new_bond.atom1 = AtomIndex(new_a1);
                new_bond.atom2 = AtomIndex(new_a2);
                new_bonds.push(new_bond);
            }
        }
        self.bonds = new_bonds;

        // Rebuild atom_bonds mapping
        self.atom_bonds = vec![SmallVec::new(); self.atoms.len()];
        for (bi, bond) in self.bonds.iter().enumerate() {
            let bond_idx = BondIndex(bi as u32);
            self.atom_bonds[bond.atom1.as_usize()].push(bond_idx);
            self.atom_bonds[bond.atom2.as_usize()].push(bond_idx);
        }

        // Update coordinate sets: extract coords for surviving atoms only
        for cs in &mut self.coord_sets {
            let mut new_coords = Vec::with_capacity(self.atoms.len() * 3);
            for old_idx in 0..n_atoms {
                if remove_set[old_idx] {
                    continue;
                }
                if let Some(v) = cs.get_atom_coord(AtomIndex(old_idx as u32)) {
                    new_coords.push(v.x);
                    new_coords.push(v.y);
                    new_coords.push(v.z);
                } else {
                    new_coords.extend_from_slice(&[0.0, 0.0, 0.0]);
                }
            }
            *cs = CoordSet::from_coords(new_coords);
        }
    }

    // =========================================================================
    // Selection Support
    // =========================================================================

    /// Get atoms matching a predicate
    pub fn select<F>(&self, predicate: F) -> Vec<AtomIndex>
    where
        F: Fn(&Atom) -> bool,
    {
        self.atoms
            .iter()
            .enumerate()
            .filter(|(_, a)| predicate(a))
            .map(|(i, _)| AtomIndex(i as u32))
            .collect()
    }

    /// Get heavy atoms (non-hydrogen)
    pub fn heavy_atoms(&self) -> Vec<AtomIndex> {
        self.select(|a| a.is_heavy())
    }

    /// Get C-alpha atoms
    pub fn ca_atoms(&self) -> Vec<AtomIndex> {
        self.select(|a| a.is_ca())
    }

    // =========================================================================
    // Validation
    // =========================================================================

    /// Validate the molecule structure
    pub fn validate(&self) -> Vec<String> {
        let mut issues = Vec::new();

        // Check bond atom indices
        for (i, bond) in self.bonds.iter().enumerate() {
            if bond.atom1.as_usize() >= self.atoms.len() {
                issues.push(format!(
                    "Bond {} references invalid atom1: {}",
                    i, bond.atom1
                ));
            }
            if bond.atom2.as_usize() >= self.atoms.len() {
                issues.push(format!(
                    "Bond {} references invalid atom2: {}",
                    i, bond.atom2
                ));
            }
        }

        // Check atom_bonds consistency
        if self.atom_bonds.len() != self.atoms.len() {
            issues.push(format!(
                "atom_bonds length ({}) doesn't match atom count ({})",
                self.atom_bonds.len(),
                self.atoms.len()
            ));
        }

        issues
    }
}

impl std::fmt::Display for ObjectMolecule {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "ObjectMolecule({}: {} atoms, {} bonds, {} states)",
            self.name,
            self.atom_count(),
            self.bond_count(),
            self.state_count()
        )
    }
}

/// Builder for creating molecules with a fluent interface
#[derive(Debug, Default)]
pub struct MoleculeBuilder {
    mol: ObjectMolecule,
    pending_coords: Vec<Vec3>,
}

impl MoleculeBuilder {
    /// Create a new molecule builder
    pub fn new(name: impl Into<String>) -> Self {
        MoleculeBuilder {
            mol: ObjectMolecule::new(name),
            pending_coords: Vec::new(),
        }
    }

    /// Add an atom with coordinates
    pub fn add_atom(mut self, atom: Atom, coord: Vec3) -> Self {
        self.mol.add_atom(atom);
        self.pending_coords.push(coord);
        self
    }

    /// Add a bond
    pub fn add_bond(mut self, atom1: AtomIndex, atom2: AtomIndex, order: BondOrder) -> Self {
        let _ = self.mol.add_bond(atom1, atom2, order);
        self
    }

    /// Set the title
    pub fn title(mut self, title: impl Into<String>) -> Self {
        self.mol.title = title.into();
        self
    }

    /// Build the molecule
    pub fn build(mut self) -> ObjectMolecule {
        // Create coordinate set from pending coords
        if !self.pending_coords.is_empty() {
            let cs = CoordSet::from_vec3(&self.pending_coords);
            self.mol.add_coord_set(cs);
        }
        self.mol
    }
}

// ============================================================================
// Chain ID Generator
// ============================================================================

/// Generates chain IDs: A, B, C, ..., Z, AA, AB, ..., AZ, BA, ...
struct ChainIdGenerator {
    counter: usize,
}

impl ChainIdGenerator {
    fn new() -> Self {
        Self { counter: 0 }
    }

    fn next(&mut self) -> String {
        let id = self.counter;
        self.counter += 1;

        if id < 26 {
            String::from((b'A' + id as u8) as char)
        } else {
            // Bijective base-26: 26=AA, 27=AB, ..., 51=AZ, 52=BA, ...
            let mut result = String::new();
            let mut n = id;
            loop {
                result.insert(0, (b'A' + (n % 26) as u8) as char);
                n /= 26;
                if n == 0 {
                    break;
                }
                n -= 1;
            }
            result
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::element::Element;

    fn create_water() -> ObjectMolecule {
        let mut mol = ObjectMolecule::new("water");

        // Add atoms
        let o = mol.add_atom(Atom::new("O", Element::Oxygen));
        let h1 = mol.add_atom(Atom::new("H1", Element::Hydrogen));
        let h2 = mol.add_atom(Atom::new("H2", Element::Hydrogen));

        // Add bonds
        mol.add_bond(o, h1, BondOrder::Single).unwrap();
        mol.add_bond(o, h2, BondOrder::Single).unwrap();

        // Add coordinates
        let cs = CoordSet::from_vec3(&[
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.96, 0.0, 0.0),
            Vec3::new(-0.24, 0.93, 0.0),
        ]);
        mol.add_coord_set(cs);

        mol
    }

    #[test]
    fn test_molecule_creation() {
        let mol = create_water();
        assert_eq!(mol.name, "water");
        assert_eq!(mol.atom_count(), 3);
        assert_eq!(mol.bond_count(), 2);
        assert_eq!(mol.state_count(), 1);
    }

    #[test]
    fn test_atom_access() {
        let mol = create_water();
        let atom = mol.get_atom(AtomIndex(0)).unwrap();
        assert_eq!(&*atom.name, "O");
        assert_eq!(atom.element, Element::Oxygen);
    }

    #[test]
    fn test_bond_lookup() {
        let mol = create_water();

        // Find bond between O and H1
        let bond_idx = mol.find_bond(AtomIndex(0), AtomIndex(1)).unwrap();
        let bond = mol.get_bond(bond_idx).unwrap();
        assert_eq!(bond.order, BondOrder::Single);

        // No bond between H1 and H2
        assert!(mol.find_bond(AtomIndex(1), AtomIndex(2)).is_none());
    }

    #[test]
    fn test_bonded_atoms() {
        let mol = create_water();
        let bonded = mol.bonded_atoms(AtomIndex(0));
        assert_eq!(bonded.len(), 2);
        assert!(bonded.contains(&AtomIndex(1)));
        assert!(bonded.contains(&AtomIndex(2)));
    }

    #[test]
    fn test_coordinate_access() {
        let mol = create_water();
        let coord = mol.get_coord(AtomIndex(0), 0).unwrap();
        assert_eq!(coord, Vec3::new(0.0, 0.0, 0.0));
    }

    #[test]
    fn test_duplicate_bond() {
        let mut mol = ObjectMolecule::new("test");
        let a1 = mol.add_atom(Atom::new("A", Element::Carbon));
        let a2 = mol.add_atom(Atom::new("B", Element::Carbon));

        mol.add_bond(a1, a2, BondOrder::Single).unwrap();
        let result = mol.add_bond(a1, a2, BondOrder::Single);
        assert!(matches!(result, Err(MolError::DuplicateBond(_, _))));
    }

    #[test]
    fn test_builder() {
        let mol = MoleculeBuilder::new("test")
            .add_atom(Atom::new("C", Element::Carbon), Vec3::new(0.0, 0.0, 0.0))
            .add_atom(Atom::new("C", Element::Carbon), Vec3::new(1.5, 0.0, 0.0))
            .add_bond(AtomIndex(0), AtomIndex(1), BondOrder::Single)
            .title("Test molecule")
            .build();

        assert_eq!(mol.atom_count(), 2);
        assert_eq!(mol.bond_count(), 1);
        assert_eq!(mol.state_count(), 1);
        assert_eq!(mol.title, "Test molecule");
    }

    #[test]
    fn test_molecule_display() {
        let mol = create_water();
        let display = format!("{}", mol);
        assert!(display.contains("water"));
        assert!(display.contains("3 atoms"));
        assert!(display.contains("2 bonds"));
    }

    #[test]
    fn test_validation() {
        let mol = create_water();
        let issues = mol.validate();
        assert!(issues.is_empty());
    }

    #[test]
    fn test_classify_atoms_protein() {
        use crate::atom::AtomResidue;
        use crate::flags::AtomFlags;
        use std::sync::Arc;

        // Create a simple peptide with amino acid residues
        let mut mol = ObjectMolecule::new("peptide");

        // Create shared residue
        let ala_res = Arc::new(AtomResidue::from_parts("A", "ALA", 1, ' ', ""));

        // Add ALA residue atoms
        let mut n = Atom::new("N", Element::Nitrogen);
        n.residue = ala_res.clone();
        mol.add_atom(n);

        let mut ca = Atom::new("CA", Element::Carbon);
        ca.residue = ala_res.clone();
        mol.add_atom(ca);

        let mut c = Atom::new("C", Element::Carbon);
        c.residue = ala_res.clone();
        mol.add_atom(c);

        let mut o = Atom::new("O", Element::Oxygen);
        o.residue = ala_res.clone();
        mol.add_atom(o);

        // Add coordinate set
        let cs = CoordSet::from_vec3(&[
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.45, 0.0, 0.0),
            Vec3::new(2.0, 1.4, 0.0),
            Vec3::new(1.4, 2.4, 0.0),
        ]);
        mol.add_coord_set(cs);

        // Before classification, atoms should NOT have PROTEIN flag
        assert!(!mol.get_atom(AtomIndex(0)).unwrap().state.flags.contains(AtomFlags::PROTEIN));

        // Classify atoms
        mol.classify_atoms();

        // After classification, all atoms should have PROTEIN and POLYMER flags
        for atom in mol.atoms() {
            assert!(atom.state.flags.contains(AtomFlags::PROTEIN), "Atom {} should have PROTEIN flag", atom.name);
            assert!(atom.state.flags.contains(AtomFlags::POLYMER), "Atom {} should have POLYMER flag", atom.name);
        }

        // CA atom should also have GUIDE flag
        let ca_atom = mol.get_atom(AtomIndex(1)).unwrap();
        assert!(ca_atom.state.flags.contains(AtomFlags::GUIDE), "CA atom should have GUIDE flag");
    }

    #[test]
    fn test_classify_atoms_water() {
        use crate::atom::AtomResidue;
        use crate::flags::AtomFlags;
        use std::sync::Arc;

        // Create water molecule
        let mut mol = ObjectMolecule::new("water");

        // Create shared residue for water
        let hoh_res = Arc::new(AtomResidue::from_parts("", "HOH", 1, ' ', ""));

        let mut o = Atom::new("O", Element::Oxygen);
        o.residue = hoh_res.clone();
        mol.add_atom(o);

        let mut h1 = Atom::new("H1", Element::Hydrogen);
        h1.residue = hoh_res.clone();
        mol.add_atom(h1);

        let cs = CoordSet::from_vec3(&[
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.96, 0.0, 0.0),
        ]);
        mol.add_coord_set(cs);

        mol.classify_atoms();

        // Water should have SOLVENT flag, not PROTEIN
        for atom in mol.atoms() {
            assert!(atom.state.flags.contains(AtomFlags::SOLVENT), "Atom {} should have SOLVENT flag", atom.name);
            assert!(!atom.state.flags.contains(AtomFlags::PROTEIN), "Atom {} should NOT have PROTEIN flag", atom.name);
        }
    }

    #[test]
    fn test_residue_is_protein_after_classification() {
        use crate::atom::AtomResidue;
        use std::sync::Arc;

        // Create a peptide and verify that is_protein() returns true after classification
        let mut mol = ObjectMolecule::new("peptide");

        let gly_res = Arc::new(AtomResidue::from_parts("A", "GLY", 1, ' ', ""));

        let mut ca = Atom::new("CA", Element::Carbon);
        ca.residue = gly_res.clone();
        mol.add_atom(ca);

        let cs = CoordSet::from_vec3(&[Vec3::new(0.0, 0.0, 0.0)]);
        mol.add_coord_set(cs);

        // Before classification
        let residue = mol.residues().next().unwrap();
        assert!(!residue.is_protein(), "Residue should NOT be protein before classification");

        // After classification
        mol.classify_atoms();
        let residue = mol.residues().next().unwrap();
        assert!(residue.is_protein(), "Residue should be protein after classification");
    }

    #[test]
    fn test_chain_id_generator() {
        let mut gen = ChainIdGenerator::new();
        assert_eq!(gen.next(), "A");
        assert_eq!(gen.next(), "B");
        // Skip to Z
        for _ in 2..26 {
            gen.next();
        }
        assert_eq!(gen.next(), "AA");
        assert_eq!(gen.next(), "AB");
    }

    #[test]
    fn test_assign_chains_protein_solvent_ions() {
        use crate::atom::AtomResidue;
        use std::sync::Arc;

        let mut mol = ObjectMolecule::new("test");

        // Protein: ALA(1) - GLY(2)
        let ala_res = Arc::new(AtomResidue::from_parts("", "ALA", 1, ' ', ""));
        for name in &["N", "CA", "C", "O"] {
            let mut atom = Atom::new(*name, if *name == "N" { Element::Nitrogen } else if *name == "O" { Element::Oxygen } else { Element::Carbon });
            atom.residue = ala_res.clone();
            mol.add_atom(atom);
        }
        let gly_res = Arc::new(AtomResidue::from_parts("", "GLY", 2, ' ', ""));
        for name in &["N", "CA", "C", "O"] {
            let mut atom = Atom::new(*name, if *name == "N" { Element::Nitrogen } else if *name == "O" { Element::Oxygen } else { Element::Carbon });
            atom.residue = gly_res.clone();
            mol.add_atom(atom);
        }

        // Ion: NA(3)
        let na_res = Arc::new(AtomResidue::from_parts("", "NA", 3, ' ', ""));
        let mut na = Atom::new("NA", Element::Sodium);
        na.residue = na_res.clone();
        mol.add_atom(na);

        // Solvent: SOL(4), SOL(5)
        for resv in [4, 5] {
            let sol_res = Arc::new(AtomResidue::from_parts("", "SOL", resv, ' ', ""));
            for name in &["OW", "HW1", "HW2"] {
                let mut atom = Atom::new(*name, if *name == "OW" { Element::Oxygen } else { Element::Hydrogen });
                atom.residue = sol_res.clone();
                mol.add_atom(atom);
            }
        }

        // Add coordinates (needed for bond generation)
        // ALA: N(0,0,0) CA(1.47,0,0) C(2.5,1.2,0) O(2.5,2.4,0)
        // GLY: N(3.8,0.5,0) CA(5.2,0.5,0) C(6.3,1.7,0) O(6.3,2.9,0)
        // Place C of ALA and N of GLY close enough for a bond (~1.33 Å)
        let coords = vec![
            // ALA
            Vec3::new(0.0, 0.0, 0.0),   // N
            Vec3::new(1.47, 0.0, 0.0),   // CA
            Vec3::new(2.5, 1.2, 0.0),    // C
            Vec3::new(2.5, 2.4, 0.0),    // O
            // GLY — N close to ALA's C for backbone bond
            Vec3::new(3.8, 1.2, 0.0),    // N (1.3 Å from C)
            Vec3::new(5.2, 1.2, 0.0),    // CA
            Vec3::new(6.3, 2.4, 0.0),    // C
            Vec3::new(6.3, 3.6, 0.0),    // O
            // NA ion far away
            Vec3::new(20.0, 20.0, 20.0), // NA
            // SOL 1
            Vec3::new(30.0, 30.0, 30.0), // OW
            Vec3::new(30.5, 30.5, 30.0), // HW1
            Vec3::new(30.5, 29.5, 30.0), // HW2
            // SOL 2
            Vec3::new(35.0, 35.0, 35.0), // OW
            Vec3::new(35.5, 35.5, 35.0), // HW1
            Vec3::new(35.5, 34.5, 35.0), // HW2
        ];
        let cs = CoordSet::from_vec3(&coords);
        mol.add_coord_set(cs);

        mol.classify_atoms();
        mol.generate_bonds(0.6);
        mol.assign_chains();

        // Protein should be chain A
        assert_eq!(mol.atoms_slice()[0].residue.chain, "A"); // ALA
        assert_eq!(mol.atoms_slice()[4].residue.chain, "A"); // GLY

        // Ion should share a chain
        let ion_chain = &mol.atoms_slice()[8].residue.chain;
        assert!(!ion_chain.is_empty());

        // Solvent should share a chain
        let sol_chain = &mol.atoms_slice()[9].residue.chain;
        assert!(!sol_chain.is_empty());

        // Ion and solvent chains should be different from protein and each other
        assert_ne!(ion_chain, "A");
        assert_ne!(sol_chain, "A");
        assert_ne!(ion_chain, sol_chain);
    }

    #[test]
    fn test_assign_chains_noop_when_chains_exist() {
        use crate::atom::AtomResidue;
        use std::sync::Arc;

        let mut mol = ObjectMolecule::new("test");

        let res = Arc::new(AtomResidue::from_parts("X", "ALA", 1, ' ', ""));
        let mut atom = Atom::new("CA", Element::Carbon);
        atom.residue = res;
        mol.add_atom(atom);

        let cs = CoordSet::from_vec3(&[Vec3::new(0.0, 0.0, 0.0)]);
        mol.add_coord_set(cs);

        mol.assign_chains();

        // Chain should remain "X" (not overwritten)
        assert_eq!(mol.atoms_slice()[0].residue.chain, "X");
    }

    #[test]
    fn test_assign_chains_resnum_wrap() {
        use crate::atom::AtomResidue;
        use std::sync::Arc;

        let mut mol = ObjectMolecule::new("test");
        let mut coords = Vec::new();

        // First protein block: ALA(99998), GLY(99999) — bonded backbone
        let ala1 = Arc::new(AtomResidue::from_parts("", "ALA", 99998, ' ', ""));
        let ala1_names = [("N", Element::Nitrogen), ("CA", Element::Carbon), ("C", Element::Carbon), ("O", Element::Oxygen)];
        let ala1_coords = [
            Vec3::new(0.0, 0.0, 0.0),   // N
            Vec3::new(1.47, 0.0, 0.0),   // CA
            Vec3::new(2.52, 1.25, 0.0),  // C
            Vec3::new(2.52, 2.40, 0.0),  // O
        ];
        for ((name, elem), coord) in ala1_names.iter().zip(&ala1_coords) {
            let mut atom = Atom::new(*name, *elem);
            atom.residue = ala1.clone();
            mol.add_atom(atom);
            coords.push(*coord);
        }

        let gly1 = Arc::new(AtomResidue::from_parts("", "GLY", 99999, ' ', ""));
        let gly1_names = [("N", Element::Nitrogen), ("CA", Element::Carbon), ("C", Element::Carbon), ("O", Element::Oxygen)];
        let gly1_coords = [
            Vec3::new(3.80, 1.25, 0.0),  // N — 1.28 Å from ALA's C (bonded)
            Vec3::new(5.24, 1.25, 0.0),  // CA
            Vec3::new(6.30, 2.45, 0.0),  // C
            Vec3::new(6.30, 3.60, 0.0),  // O
        ];
        for ((name, elem), coord) in gly1_names.iter().zip(&gly1_coords) {
            let mut atom = Atom::new(*name, *elem);
            atom.residue = gly1.clone();
            mol.add_atom(atom);
            coords.push(*coord);
        }

        // Second protein block (resnum wraps): ALA(1) — far from first block
        let ala2 = Arc::new(AtomResidue::from_parts("", "ALA", 1, ' ', ""));
        let ala2_names = [("N", Element::Nitrogen), ("CA", Element::Carbon), ("C", Element::Carbon), ("O", Element::Oxygen)];
        let ala2_coords = [
            Vec3::new(20.0, 0.0, 0.0),
            Vec3::new(21.47, 0.0, 0.0),
            Vec3::new(22.52, 1.25, 0.0),
            Vec3::new(22.52, 2.40, 0.0),
        ];
        for ((name, elem), coord) in ala2_names.iter().zip(&ala2_coords) {
            let mut atom = Atom::new(*name, *elem);
            atom.residue = ala2.clone();
            mol.add_atom(atom);
            coords.push(*coord);
        }

        let cs = CoordSet::from_vec3(&coords);
        mol.add_coord_set(cs);

        mol.classify_atoms();
        mol.generate_bonds(0.6);
        mol.assign_chains();

        // First block = chain A, second block = chain B (due to resnum decrease)
        assert_eq!(mol.atoms_slice()[0].residue.chain, "A");
        assert_eq!(mol.atoms_slice()[4].residue.chain, "A");
        assert_eq!(mol.atoms_slice()[8].residue.chain, "B");
    }
}
