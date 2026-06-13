//! Iterators for molecular data structures
//!
//! Provides iterators for traversing atoms grouped by residue, chain, or
//! chain subchain. All iterators are built on a shared [`GroupingIterator`]
//! foundation that walks consecutive atoms grouped by a comparison predicate.
//!
//! Hierarchy:
//! - [`ChainIterator`] — one [`ChainView`] per unique chain identifier.
//!   Chain "A" with mixed polymer + HEM + waters yields a single chain view.
//! - [`SubchainIterator`] — one [`SubchainView`] per homogeneous run within a
//!   chain. Chain "A" with polymer + HEM splits into two subchains.
//! - [`ResidueIterator`] — one [`ResidueView`] per residue.
//!
//! Use `ChainView::subchains()` (or [`SubchainIterator`]) when you need the
//! polymer / HET / solvent split (Objects panel, popover commands). Use
//! `ChainView::residues()` directly when chain-scoped iteration is enough
//! (sequence viewer). Use `ChainView::polymer_subchains()` for representations
//! that only operate on biopolymer runs (cartoon).

use crate::atom::Atom;
use crate::residue::{atoms_same_residue, ChainView, ResidueKey, ResidueView};
use crate::subchain::{atoms_same_subchain, SubchainKind, SubchainLabel, SubchainView};

// =============================================================================
// Generic Grouping Iterator
// =============================================================================

/// A generic iterator that groups consecutive atoms based on a comparison function.
///
/// This is the foundation for residue, chain, and subchain iterators, eliminating
/// code duplication. The iterator advances through a slice of atoms, grouping
/// consecutive atoms that satisfy the `same_group` predicate.
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

/// Iterator over chains in a slice of atoms.
///
/// Groups consecutive atoms by chain identifier — and **only** by chain
/// identifier. A chain "A" containing both biopolymer and HET atoms yields
/// a single [`ChainView`] covering all of them. Use [`ChainView::subchains`]
/// to descend into polymer / HET / solvent sub-groups.
pub struct ChainIterator<'a> {
    inner: GroupingIterator<'a, fn(&Atom, &Atom) -> bool>,
}

impl<'a> ChainIterator<'a> {
    /// Create a new chain iterator.
    pub fn new(atoms: &'a [Atom], base_index: usize) -> Self {
        ChainIterator {
            inner: GroupingIterator::new(atoms, base_index, atoms_same_chain_id),
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

/// Predicate: do two atoms belong to the same chain identifier?
///
/// This is the predicate behind [`ChainIterator`]. The HET / resn split is
/// handled one level down by [`SubchainIterator`].
#[inline]
pub fn atoms_same_chain_id(a: &Atom, b: &Atom) -> bool {
    a.residue.chain == b.residue.chain
}

// =============================================================================
// Subchain Iterator
// =============================================================================

/// Iterator over chain subchains — homogeneous runs of atoms within a chain.
///
/// Groups consecutive atoms by the [`atoms_same_subchain`] predicate:
/// same `chain_id`, same `hetatm` flag, and (for non-polymer atoms) same
/// residue name. Each iteration yields a [`SubchainView`].
pub struct SubchainIterator<'a> {
    inner: GroupingIterator<'a, fn(&Atom, &Atom) -> bool>,
}

impl<'a> SubchainIterator<'a> {
    /// Create a new subchain iterator.
    pub fn new(atoms: &'a [Atom], base_index: usize) -> Self {
        SubchainIterator {
            inner: GroupingIterator::new(atoms, base_index, atoms_same_subchain),
        }
    }
}

impl<'a> Iterator for SubchainIterator<'a> {
    type Item = SubchainView<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        let (start, end, first_atom) = self.inner.next_group()?;
        let base = self.inner.base_index();
        let atoms = self.inner.atoms();

        let kind = SubchainKind::from_atom(first_atom);
        let label = if first_atom.state.hetatm {
            SubchainLabel::Single(first_atom.residue.resn.to_string())
        } else {
            SubchainLabel::Polymer
        };

        Some(SubchainView::from_range(
            first_atom.residue.chain.clone(),
            kind,
            label,
            atoms,
            base as u32,
            start as u32,
            end as u32,
        ))
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

    fn make_mixed_chain_atoms() -> Vec<Atom> {
        use crate::flags::AtomFlags;
        let mut atoms = Vec::new();

        // Chain A biopolymer (5 atoms)
        let res_ala = Arc::new(AtomResidue::from_parts("A", "ALA", 1, ' ', ""));
        for name in &["N", "CA", "C", "O", "CB"] {
            let mut atom = Atom::new(*name, Element::Carbon);
            atom.residue = res_ala.clone();
            atom.state.flags |= AtomFlags::PROTEIN;
            atoms.push(atom);
        }

        // Chain A PO4 hetatm (1 atom)
        let res_po4 = Arc::new(AtomResidue::from_parts("A", "PO4", 300, ' ', ""));
        let mut po4 = Atom::new("P", Element::Phosphorus);
        po4.residue = res_po4;
        po4.state.hetatm = true;
        po4.state.flags |= AtomFlags::ORGANIC;
        atoms.push(po4);

        // Chain A HEM hetatm (4 atoms)
        let res_hem = Arc::new(AtomResidue::from_parts("A", "HEM", 301, ' ', ""));
        for name in &["FE", "NA", "NB", "NC"] {
            let mut atom = Atom::new(*name, Element::Iron);
            atom.residue = res_hem.clone();
            atom.state.hetatm = true;
            atom.state.flags |= AtomFlags::ORGANIC;
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
    fn test_chain_iterator_groups_by_id() {
        let atoms = make_test_atoms();
        let chains: Vec<_> = ChainIterator::new(&atoms, 0).collect();

        assert_eq!(chains.len(), 2);
        assert_eq!(chains[0].id(), "A");
        assert_eq!(chains[0].len(), 9);
        assert_eq!(chains[1].id(), "B");
        assert_eq!(chains[1].len(), 6);
    }

    #[test]
    fn test_chain_iterator_does_not_split_hetatm() {
        let atoms = make_mixed_chain_atoms();
        let chains: Vec<_> = ChainIterator::new(&atoms, 0).collect();

        // All atoms share chain "A" — exactly one ChainView regardless of HET split.
        assert_eq!(chains.len(), 1);
        assert_eq!(chains[0].id(), "A");
        assert_eq!(chains[0].len(), 10);
    }

    #[test]
    fn test_subchain_iterator_splits_hetatm() {
        let atoms = make_mixed_chain_atoms();
        let subchains: Vec<_> = SubchainIterator::new(&atoms, 0).collect();

        assert_eq!(subchains.len(), 3);

        assert_eq!(subchains[0].chain_id(), "A");
        assert_eq!(subchains[0].kind, SubchainKind::Biopolymer);
        assert_eq!(subchains[0].label, SubchainLabel::Polymer);
        assert_eq!(subchains[0].len(), 5);

        assert_eq!(subchains[1].chain_id(), "A");
        assert_eq!(subchains[1].kind, SubchainKind::Organic);
        assert_eq!(subchains[1].label, SubchainLabel::Single("PO4".to_string()));
        assert_eq!(subchains[1].len(), 1);

        assert_eq!(subchains[2].chain_id(), "A");
        assert_eq!(subchains[2].kind, SubchainKind::Organic);
        assert_eq!(subchains[2].label, SubchainLabel::Single("HEM".to_string()));
        assert_eq!(subchains[2].len(), 4);
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
    fn test_subchain_iterator_with_base_index() {
        let atoms = make_mixed_chain_atoms();
        let subchains: Vec<_> = SubchainIterator::new(&atoms, 100).collect();

        let collect_global =
            |s: &SubchainView| -> Vec<u32> { s.iter_indexed().map(|(i, _)| i.0).collect() };
        assert_eq!(
            collect_global(&subchains[0]),
            (100u32..105).collect::<Vec<_>>()
        );
        assert_eq!(collect_global(&subchains[1]), vec![105u32]);
        assert_eq!(
            collect_global(&subchains[2]),
            (106u32..110).collect::<Vec<_>>()
        );
    }

    #[test]
    fn test_chain_view_subchains() {
        let atoms = make_mixed_chain_atoms();
        let chains: Vec<_> = ChainIterator::new(&atoms, 0).collect();
        let subchains: Vec<_> = chains[0].subchains().collect();

        assert_eq!(subchains.len(), 3);
        assert!(subchains[0].kind.is_biopolymer());
        assert_eq!(subchains[1].label, SubchainLabel::Single("PO4".to_string()));
        assert_eq!(subchains[2].label, SubchainLabel::Single("HEM".to_string()));
    }

    #[test]
    fn test_chain_view_polymer_subchains() {
        let atoms = make_mixed_chain_atoms();
        let chains: Vec<_> = ChainIterator::new(&atoms, 0).collect();
        let polymer: Vec<_> = chains[0].polymer_subchains().collect();

        assert_eq!(polymer.len(), 1);
        assert!(polymer[0].kind.is_biopolymer());
        assert_eq!(polymer[0].len(), 5);
    }

    #[test]
    fn test_empty_iterator() {
        let atoms: Vec<Atom> = Vec::new();
        let residues: Vec<_> = ResidueIterator::new(&atoms, 0).collect();
        assert!(residues.is_empty());

        let chains: Vec<_> = ChainIterator::new(&atoms, 0).collect();
        assert!(chains.is_empty());

        let subchains: Vec<_> = SubchainIterator::new(&atoms, 0).collect();
        assert!(subchains.is_empty());
    }
}
