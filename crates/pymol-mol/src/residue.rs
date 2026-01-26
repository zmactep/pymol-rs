//! Residue and chain utilities
//!
//! Provides types and functions for working with residues and chains.
//! Since PyMOL stores residue/chain information inline in atoms rather than
//! as separate structures, this module provides view types for convenient access.

use crate::atom::Atom;
use crate::index::AtomIndex;
use std::ops::Range;

/// Key for uniquely identifying a residue within a molecule
///
/// A residue is uniquely identified by its chain, name, number, and insertion code.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct ResidueKey {
    /// Chain identifier
    pub chain: String,
    /// Residue name (e.g., "ALA", "GLY")
    pub resn: String,
    /// Residue sequence number
    pub resv: i32,
    /// Insertion code
    pub inscode: char,
}

impl ResidueKey {
    /// Create a new residue key
    pub fn new(
        chain: impl Into<String>,
        resn: impl Into<String>,
        resv: i32,
        inscode: char,
    ) -> Self {
        ResidueKey {
            chain: chain.into(),
            resn: resn.into(),
            resv,
            inscode,
        }
    }

    /// Create a residue key from an atom
    pub fn from_atom(atom: &Atom) -> Self {
        ResidueKey {
            chain: atom.chain.clone(),
            resn: atom.resn.clone(),
            resv: atom.resv,
            inscode: atom.inscode,
        }
    }

    /// Get a display string for the residue (e.g., "A/ALA`1")
    pub fn display(&self) -> String {
        if self.inscode != ' ' && self.inscode != '\0' {
            format!("{}/{}'{}{}", self.chain, self.resn, self.resv, self.inscode)
        } else {
            format!("{}/{}'{}",self.chain, self.resn, self.resv)
        }
    }
}

impl std::fmt::Display for ResidueKey {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.display())
    }
}

/// A view into the atoms of a single residue
///
/// This is a borrowed view that does not own the atom data.
/// It provides convenient access to residue-level operations.
#[derive(Debug)]
pub struct ResidueView<'a> {
    /// The residue key
    pub key: ResidueKey,
    /// Slice of atoms in this residue
    pub atoms: &'a [Atom],
    /// Range of atom indices in the parent molecule
    pub atom_range: Range<usize>,
}

impl<'a> ResidueView<'a> {
    /// Create a new residue view
    pub fn new(key: ResidueKey, atoms: &'a [Atom], atom_range: Range<usize>) -> Self {
        ResidueView {
            key,
            atoms,
            atom_range,
        }
    }

    /// Get the chain ID
    #[inline]
    pub fn chain(&self) -> &str {
        &self.key.chain
    }

    /// Get the residue name
    #[inline]
    pub fn resn(&self) -> &str {
        &self.key.resn
    }

    /// Get the residue number
    #[inline]
    pub fn resv(&self) -> i32 {
        self.key.resv
    }

    /// Get the number of atoms in this residue
    #[inline]
    pub fn len(&self) -> usize {
        self.atoms.len()
    }

    /// Check if the residue is empty
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.atoms.is_empty()
    }

    /// Get an atom by index within the residue
    pub fn get(&self, index: usize) -> Option<&'a Atom> {
        self.atoms.get(index)
    }

    /// Iterate over atoms in the residue
    pub fn iter(&self) -> impl Iterator<Item = &'a Atom> {
        self.atoms.iter()
    }

    /// Iterate over atom indices and atoms
    pub fn iter_indexed(&self) -> impl Iterator<Item = (AtomIndex, &'a Atom)> {
        self.atom_range
            .clone()
            .zip(self.atoms.iter())
            .map(|(idx, atom)| (AtomIndex(idx as u32), atom))
    }

    /// Find an atom by name
    pub fn find_by_name(&self, name: &str) -> Option<(AtomIndex, &'a Atom)> {
        for (idx, atom) in self.iter_indexed() {
            if atom.name == name {
                return Some((idx, atom));
            }
        }
        None
    }

    /// Get the C-alpha atom (for proteins)
    pub fn ca(&self) -> Option<(AtomIndex, &'a Atom)> {
        self.find_by_name("CA")
    }

    /// Check if this is a protein residue
    pub fn is_protein(&self) -> bool {
        self.atoms
            .first()
            .map(|a| a.flags.contains(crate::flags::AtomFlags::PROTEIN))
            .unwrap_or(false)
    }

    /// Check if this is a nucleic acid residue
    pub fn is_nucleic(&self) -> bool {
        self.atoms
            .first()
            .map(|a| a.flags.contains(crate::flags::AtomFlags::NUCLEIC))
            .unwrap_or(false)
    }
}

/// A view into the atoms of a single chain
#[derive(Debug)]
pub struct ChainView<'a> {
    /// Chain identifier
    pub chain_id: String,
    /// Slice of atoms in this chain
    pub atoms: &'a [Atom],
    /// Range of atom indices in the parent molecule
    pub atom_range: Range<usize>,
}

impl<'a> ChainView<'a> {
    /// Create a new chain view
    pub fn new(chain_id: impl Into<String>, atoms: &'a [Atom], atom_range: Range<usize>) -> Self {
        ChainView {
            chain_id: chain_id.into(),
            atoms,
            atom_range,
        }
    }

    /// Get the chain ID
    #[inline]
    pub fn id(&self) -> &str {
        &self.chain_id
    }

    /// Get the number of atoms in this chain
    #[inline]
    pub fn len(&self) -> usize {
        self.atoms.len()
    }

    /// Check if the chain is empty
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.atoms.is_empty()
    }

    /// Iterate over atoms in the chain
    pub fn iter(&self) -> impl Iterator<Item = &'a Atom> {
        self.atoms.iter()
    }

    /// Iterate over atom indices and atoms
    pub fn iter_indexed(&self) -> impl Iterator<Item = (AtomIndex, &'a Atom)> {
        self.atom_range
            .clone()
            .zip(self.atoms.iter())
            .map(|(idx, atom)| (AtomIndex(idx as u32), atom))
    }

    /// Iterate over residues in this chain
    pub fn residues(&self) -> ResidueIterator<'a> {
        ResidueIterator::new(self.atoms, self.atom_range.start)
    }
}

/// Check if two atoms are in the same residue
#[inline]
pub fn atoms_same_residue(a: &Atom, b: &Atom) -> bool {
    a.chain == b.chain
        && a.resv == b.resv
        && a.inscode == b.inscode
        && a.resn == b.resn
}

/// Check if two atoms are in the same chain
#[inline]
pub fn atoms_same_chain(a: &Atom, b: &Atom) -> bool {
    a.chain == b.chain
}

/// Check if two atoms are in the same segment
#[inline]
pub fn atoms_same_segment(a: &Atom, b: &Atom) -> bool {
    a.segi == b.segi
}

/// Iterator over residues in a slice of atoms
pub struct ResidueIterator<'a> {
    atoms: &'a [Atom],
    base_index: usize,
    current: usize,
}

impl<'a> ResidueIterator<'a> {
    /// Create a new residue iterator
    pub fn new(atoms: &'a [Atom], base_index: usize) -> Self {
        ResidueIterator {
            atoms,
            base_index,
            current: 0,
        }
    }
}

impl<'a> Iterator for ResidueIterator<'a> {
    type Item = ResidueView<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.current >= self.atoms.len() {
            return None;
        }

        let start = self.current;
        let first_atom = &self.atoms[start];
        let key = ResidueKey::from_atom(first_atom);

        // Find the end of this residue
        let mut end = start + 1;
        while end < self.atoms.len() {
            if !atoms_same_residue(&self.atoms[end], first_atom) {
                break;
            }
            end += 1;
        }

        self.current = end;

        Some(ResidueView {
            key,
            atoms: &self.atoms[start..end],
            atom_range: (self.base_index + start)..(self.base_index + end),
        })
    }
}

/// Iterator over chains in a slice of atoms
pub struct ChainIterator<'a> {
    atoms: &'a [Atom],
    base_index: usize,
    current: usize,
}

impl<'a> ChainIterator<'a> {
    /// Create a new chain iterator
    pub fn new(atoms: &'a [Atom], base_index: usize) -> Self {
        ChainIterator {
            atoms,
            base_index,
            current: 0,
        }
    }
}

impl<'a> Iterator for ChainIterator<'a> {
    type Item = ChainView<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.current >= self.atoms.len() {
            return None;
        }

        let start = self.current;
        let first_atom = &self.atoms[start];
        let chain_id = first_atom.chain.clone();

        // Find the end of this chain
        let mut end = start + 1;
        while end < self.atoms.len() {
            if !atoms_same_chain(&self.atoms[end], first_atom) {
                break;
            }
            end += 1;
        }

        self.current = end;

        Some(ChainView {
            chain_id,
            atoms: &self.atoms[start..end],
            atom_range: (self.base_index + start)..(self.base_index + end),
        })
    }
}

/// Standard amino acid residue names (3-letter codes)
pub const AMINO_ACIDS: &[&str] = &[
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    // Non-standard but common
    "MSE", "SEC", "PYL",
];

/// Standard nucleotide residue names
pub const NUCLEOTIDES: &[&str] = &[
    // DNA
    "DA", "DC", "DG", "DT",
    // RNA
    "A", "C", "G", "U",
    // Alternative names
    "ADE", "CYT", "GUA", "THY", "URA",
];

/// Check if a residue name is a standard amino acid
pub fn is_amino_acid(resn: &str) -> bool {
    AMINO_ACIDS.contains(&resn)
}

/// Check if a residue name is a nucleotide
pub fn is_nucleotide(resn: &str) -> bool {
    NUCLEOTIDES.contains(&resn)
}

/// Common water residue names
pub const WATER_NAMES: &[&str] = &["HOH", "WAT", "H2O", "DOD", "TIP", "TIP3", "SPC"];

/// Check if a residue name is water
pub fn is_water(resn: &str) -> bool {
    WATER_NAMES.contains(&resn)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::element::Element;

    fn make_test_atoms() -> Vec<Atom> {
        let mut atoms = Vec::new();

        // Chain A, residue 1 (ALA)
        for name in &["N", "CA", "C", "O", "CB"] {
            let mut atom = Atom::new(*name, Element::Carbon);
            atom.chain = "A".to_string();
            atom.resn = "ALA".to_string();
            atom.resv = 1;
            atoms.push(atom);
        }

        // Chain A, residue 2 (GLY)
        for name in &["N", "CA", "C", "O"] {
            let mut atom = Atom::new(*name, Element::Carbon);
            atom.chain = "A".to_string();
            atom.resn = "GLY".to_string();
            atom.resv = 2;
            atoms.push(atom);
        }

        // Chain B, residue 1 (SER)
        for name in &["N", "CA", "C", "O", "CB", "OG"] {
            let mut atom = Atom::new(*name, Element::Carbon);
            atom.chain = "B".to_string();
            atom.resn = "SER".to_string();
            atom.resv = 1;
            atoms.push(atom);
        }

        atoms
    }

    #[test]
    fn test_residue_key() {
        let key = ResidueKey::new("A", "ALA", 1, ' ');
        assert_eq!(key.chain, "A");
        assert_eq!(key.resn, "ALA");
        assert_eq!(key.resv, 1);
    }

    #[test]
    fn test_residue_iterator() {
        let atoms = make_test_atoms();
        let residues: Vec<_> = ResidueIterator::new(&atoms, 0).collect();

        assert_eq!(residues.len(), 3);
        assert_eq!(residues[0].resn(), "ALA");
        assert_eq!(residues[0].len(), 5);
        assert_eq!(residues[1].resn(), "GLY");
        assert_eq!(residues[1].len(), 4);
        assert_eq!(residues[2].resn(), "SER");
        assert_eq!(residues[2].len(), 6);
    }

    #[test]
    fn test_chain_iterator() {
        let atoms = make_test_atoms();
        let chains: Vec<_> = ChainIterator::new(&atoms, 0).collect();

        assert_eq!(chains.len(), 2);
        assert_eq!(chains[0].id(), "A");
        assert_eq!(chains[0].len(), 9);
        assert_eq!(chains[1].id(), "B");
        assert_eq!(chains[1].len(), 6);
    }

    #[test]
    fn test_residue_view_find() {
        let atoms = make_test_atoms();
        let residue = ResidueIterator::new(&atoms, 0).next().unwrap();

        let (idx, ca) = residue.ca().unwrap();
        assert_eq!(ca.name, "CA");
        assert_eq!(idx.as_usize(), 1);
    }

    #[test]
    fn test_atoms_same_residue() {
        let atoms = make_test_atoms();
        assert!(atoms_same_residue(&atoms[0], &atoms[1]));
        assert!(!atoms_same_residue(&atoms[0], &atoms[5]));
    }

    #[test]
    fn test_is_amino_acid() {
        assert!(is_amino_acid("ALA"));
        assert!(is_amino_acid("GLY"));
        assert!(!is_amino_acid("HOH"));
    }

    #[test]
    fn test_is_water() {
        assert!(is_water("HOH"));
        assert!(is_water("WAT"));
        assert!(!is_water("ALA"));
    }
}
