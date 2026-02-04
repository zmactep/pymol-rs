//! Iterators for molecular data structures
//!
//! Provides iterators for traversing atoms grouped by residue or chain.
//! Uses a generic `GroupingIterator` foundation to eliminate code duplication.

use crate::atom::Atom;
use crate::residue::{atoms_same_chain, atoms_same_residue, ChainView, ResidueKey, ResidueView};

// =============================================================================
// Generic Grouping Iterator
// =============================================================================

/// A generic iterator that groups consecutive atoms based on a comparison function.
///
/// This is the foundation for both `ResidueIterator` and `ChainIterator`, eliminating
/// code duplication between them. The iterator advances through a slice of atoms,
/// grouping consecutive atoms that satisfy the `same_group` predicate.
///
/// # Type Parameters
///
/// * `'a` - Lifetime of the atom slice
/// * `F` - Comparison function `Fn(&Atom, &Atom) -> bool` that determines if two atoms
///   belong to the same group
pub(crate) struct GroupingIterator<'a, F> {
    atoms: &'a [Atom],
    base_index: usize,
    current: usize,
    same_group: F,
}

impl<'a, F> GroupingIterator<'a, F>
where
    F: Fn(&Atom, &Atom) -> bool,
{
    /// Create a new grouping iterator
    ///
    /// # Arguments
    ///
    /// * `atoms` - Slice of atoms to iterate over
    /// * `base_index` - Starting atom index in the parent molecule (for correct range calculation)
    /// * `same_group` - Predicate that returns `true` if two atoms belong to the same group
    pub(crate) fn new(atoms: &'a [Atom], base_index: usize, same_group: F) -> Self {
        GroupingIterator {
            atoms,
            base_index,
            current: 0,
            same_group,
        }
    }

    /// Advance to the next group and return (start, end, first_atom)
    ///
    /// Returns `None` when all atoms have been consumed.
    pub(crate) fn next_group(&mut self) -> Option<(usize, usize, &'a Atom)> {
        if self.current >= self.atoms.len() {
            return None;
        }

        let start = self.current;
        let first_atom = &self.atoms[start];

        // Find the end of this group
        let mut end = start + 1;
        while end < self.atoms.len() {
            if !(self.same_group)(&self.atoms[end], first_atom) {
                break;
            }
            end += 1;
        }

        self.current = end;
        Some((start, end, first_atom))
    }

    /// Get the base index (starting atom index in parent molecule)
    #[inline]
    pub(crate) fn base_index(&self) -> usize {
        self.base_index
    }

    /// Get the atoms slice
    #[inline]
    pub(crate) fn atoms(&self) -> &'a [Atom] {
        self.atoms
    }
}

// =============================================================================
// Residue Iterator
// =============================================================================

/// Iterator over residues in a slice of atoms
///
/// Groups consecutive atoms by residue (chain, residue name, residue number, insertion code).
/// Each iteration yields a `ResidueView` containing all atoms in that residue.
///
/// # Example
///
/// ```ignore
/// for residue in molecule.residues() {
///     println!("Residue {} has {} atoms", residue.key, residue.len());
/// }
/// ```
pub struct ResidueIterator<'a> {
    inner: GroupingIterator<'a, fn(&Atom, &Atom) -> bool>,
}

impl<'a> ResidueIterator<'a> {
    /// Create a new residue iterator
    ///
    /// # Arguments
    ///
    /// * `atoms` - Slice of atoms to iterate over
    /// * `base_index` - Starting atom index in the parent molecule
    pub fn new(atoms: &'a [Atom], base_index: usize) -> Self {
        ResidueIterator {
            inner: GroupingIterator::new(atoms, base_index, atoms_same_residue),
        }
    }
}

impl<'a> Iterator for ResidueIterator<'a> {
    type Item = ResidueView<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        let (start, end, first_atom) = self.inner.next_group()?;
        let base = self.inner.base_index();
        let atoms = self.inner.atoms();

        Some(ResidueView {
            key: ResidueKey::from_atom(first_atom),
            atoms: &atoms[start..end],
            atom_range: (base + start)..(base + end),
        })
    }
}

// =============================================================================
// Chain Iterator
// =============================================================================

/// Iterator over chains in a slice of atoms
///
/// Groups consecutive atoms by chain identifier. Each iteration yields a `ChainView`
/// containing all atoms in that chain.
///
/// # Example
///
/// ```ignore
/// for chain in molecule.chains() {
///     println!("Chain {} has {} atoms", chain.id(), chain.len());
/// }
/// ```
pub struct ChainIterator<'a> {
    inner: GroupingIterator<'a, fn(&Atom, &Atom) -> bool>,
}

impl<'a> ChainIterator<'a> {
    /// Create a new chain iterator
    ///
    /// # Arguments
    ///
    /// * `atoms` - Slice of atoms to iterate over
    /// * `base_index` - Starting atom index in the parent molecule
    pub fn new(atoms: &'a [Atom], base_index: usize) -> Self {
        ChainIterator {
            inner: GroupingIterator::new(atoms, base_index, atoms_same_chain),
        }
    }
}

impl<'a> Iterator for ChainIterator<'a> {
    type Item = ChainView<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        let (start, end, first_atom) = self.inner.next_group()?;
        let base = self.inner.base_index();
        let atoms = self.inner.atoms();

        Some(ChainView {
            chain_id: first_atom.residue.chain.clone(),
            atoms: &atoms[start..end],
            atom_range: (base + start)..(base + end),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::AtomResidue;
    use crate::element::Element;
    use std::sync::Arc;

    fn make_test_atoms() -> Vec<Atom> {
        let mut atoms = Vec::new();

        // Chain A, residue 1 (ALA)
        let res_ala = Arc::new(AtomResidue::from_parts("A", "ALA", 1, ' ', ""));
        for name in &["N", "CA", "C", "O", "CB"] {
            let mut atom = Atom::new(*name, Element::Carbon);
            atom.residue = res_ala.clone();
            atoms.push(atom);
        }

        // Chain A, residue 2 (GLY)
        let res_gly = Arc::new(AtomResidue::from_parts("A", "GLY", 2, ' ', ""));
        for name in &["N", "CA", "C", "O"] {
            let mut atom = Atom::new(*name, Element::Carbon);
            atom.residue = res_gly.clone();
            atoms.push(atom);
        }

        // Chain B, residue 1 (SER)
        let res_ser = Arc::new(AtomResidue::from_parts("B", "SER", 1, ' ', ""));
        for name in &["N", "CA", "C", "O", "CB", "OG"] {
            let mut atom = Atom::new(*name, Element::Carbon);
            atom.residue = res_ser.clone();
            atoms.push(atom);
        }

        atoms
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
    fn test_residue_iterator_with_base_index() {
        let atoms = make_test_atoms();
        let residues: Vec<_> = ResidueIterator::new(&atoms, 100).collect();

        // Check that atom ranges are correctly offset
        assert_eq!(residues[0].atom_range, 100..105);
        assert_eq!(residues[1].atom_range, 105..109);
        assert_eq!(residues[2].atom_range, 109..115);
    }

    #[test]
    fn test_chain_iterator_with_base_index() {
        let atoms = make_test_atoms();
        let chains: Vec<_> = ChainIterator::new(&atoms, 50).collect();

        assert_eq!(chains[0].atom_range, 50..59);
        assert_eq!(chains[1].atom_range, 59..65);
    }

    #[test]
    fn test_empty_iterator() {
        let atoms: Vec<Atom> = Vec::new();
        let residues: Vec<_> = ResidueIterator::new(&atoms, 0).collect();
        assert!(residues.is_empty());

        let chains: Vec<_> = ChainIterator::new(&atoms, 0).collect();
        assert!(chains.is_empty());
    }
}
