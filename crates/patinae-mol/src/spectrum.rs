//! Per-residue rank tables for spectrum / rainbow coloring.
//!
//! [`ResidueRankTable`] maps each polymer residue to `(rank, total)` where
//! `rank` is the 0-based polymer-residue index within the chain and `total`
//! is the count of polymer residues in that chain. Non-polymer residues
//! (waters, ions, ligands, capping groups) surface as `None` from
//! [`ResidueRankTable::get`], letting the renderer and the sequence panel
//! fall back to element / residue-type coloring without staining HET atoms
//! with the polymer rainbow.

use std::collections::HashMap;

use crate::atom::Atom;
use crate::molecule::ObjectMolecule;
use crate::residue::ResidueKey;

/// Per-residue rank lookup, shared between the renderer's `ColorResolver`
/// and the sequence panel's color path.
#[derive(Debug, Default, Clone)]
pub struct ResidueRankTable {
    by_residue: HashMap<ResidueKey, (u32, u32)>,
}

impl ResidueRankTable {
    /// Look up `(rank, total)` for an atom's residue, or `None` if the
    /// residue is not part of a polymer chain.
    #[inline]
    pub fn get(&self, atom: &Atom) -> Option<(u32, u32)> {
        self.by_residue.get(&atom.residue.key).copied()
    }

    /// Sample position in `[0, 1]` for a polymer atom, or `None` for
    /// non-polymer. Single-residue chains map to `0.5` (mid-spectrum) to
    /// avoid division by zero.
    #[inline]
    pub fn t_for(&self, atom: &Atom) -> Option<f32> {
        let (rank, total) = self.get(atom)?;
        Some(if total <= 1 {
            0.5
        } else {
            rank as f32 / (total - 1) as f32
        })
    }
}

/// Build a [`ResidueRankTable`] for `mol`.
///
/// A residue is "polymer" if its first atom carries [`crate::flags::AtomFlags::POLYMER`]
/// (set during [`ObjectMolecule::classify_atoms`]). Ranks are assigned in
/// residue iteration order, per chain identifier.
pub fn polymer_residue_ranks(mol: &ObjectMolecule) -> ResidueRankTable {
    let mut by_residue: HashMap<ResidueKey, (u32, u32)> = HashMap::new();
    let mut chain_counts: HashMap<String, u32> = HashMap::new();

    // Pass 1: assign rank per polymer residue (per chain).
    for residue in mol.residues() {
        if !residue.is_polymer() {
            continue;
        }
        let counter = chain_counts.entry(residue.chain().to_string()).or_insert(0);
        let rank = *counter;
        *counter += 1;
        by_residue.insert(residue.key.clone(), (rank, 0));
    }

    // Pass 2: fill in chain totals.
    for (key, entry) in by_residue.iter_mut() {
        if let Some(total) = chain_counts.get(&key.chain) {
            entry.1 = *total;
        }
    }

    ResidueRankTable { by_residue }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::AtomResidue;
    use crate::element::Element;
    use crate::flags::AtomFlags;
    use crate::index::AtomIndex;
    use std::sync::Arc;

    fn polymer_atom(resn: &str, resv: i32, chain: &str) -> Atom {
        let mut atom = Atom::new("CA", Element::Carbon);
        atom.residue = Arc::new(AtomResidue::from_parts(chain, resn, resv, ' ', ""));
        atom.state.flags = AtomFlags::PROTEIN | AtomFlags::POLYMER;
        atom
    }

    fn het_atom(resn: &str, resv: i32, chain: &str) -> Atom {
        let mut atom = Atom::new("FE", Element::Iron);
        atom.residue = Arc::new(AtomResidue::from_parts(chain, resn, resv, ' ', ""));
        atom.state.flags = AtomFlags::ORGANIC;
        atom.state.hetatm = true;
        atom
    }

    #[test]
    fn ranks_polymer_residues_per_chain() {
        let mut mol = ObjectMolecule::new("t");
        mol.add_atom(polymer_atom("ALA", 1, "A"));
        mol.add_atom(polymer_atom("GLY", 2, "A"));
        mol.add_atom(polymer_atom("SER", 3, "A"));
        mol.add_atom(het_atom("HEM", 100, "A"));
        mol.add_atom(polymer_atom("ALA", 1, "B"));
        mol.add_atom(polymer_atom("GLY", 2, "B"));

        let table = polymer_residue_ranks(&mol);

        let a0 = mol.get_atom(AtomIndex(0)).unwrap();
        let a1 = mol.get_atom(AtomIndex(1)).unwrap();
        let a2 = mol.get_atom(AtomIndex(2)).unwrap();
        let a3 = mol.get_atom(AtomIndex(3)).unwrap();
        let a4 = mol.get_atom(AtomIndex(4)).unwrap();
        let a5 = mol.get_atom(AtomIndex(5)).unwrap();

        assert_eq!(table.get(a0), Some((0, 3)));
        assert_eq!(table.get(a1), Some((1, 3)));
        assert_eq!(table.get(a2), Some((2, 3)));
        assert_eq!(table.get(a3), None, "HET ligand is not polymer");
        assert_eq!(table.get(a4), Some((0, 2)));
        assert_eq!(table.get(a5), Some((1, 2)));

        assert_eq!(table.t_for(a0), Some(0.0));
        assert_eq!(table.t_for(a1), Some(0.5));
        assert_eq!(table.t_for(a2), Some(1.0));
        assert_eq!(table.t_for(a3), None);
        assert_eq!(table.t_for(a4), Some(0.0));
        assert_eq!(table.t_for(a5), Some(1.0));
    }

    #[test]
    fn single_residue_chain_maps_to_midpoint() {
        let mut mol = ObjectMolecule::new("t");
        mol.add_atom(polymer_atom("ALA", 1, "A"));
        let table = polymer_residue_ranks(&mol);
        let atom = mol.get_atom(AtomIndex(0)).unwrap();
        assert_eq!(table.t_for(atom), Some(0.5));
    }
}
