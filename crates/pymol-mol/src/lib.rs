//! PyMOL-RS Molecular Data Structures
//!
//! This crate provides the core molecular data structures for PyMOL-RS, including:
//!
//! - [`Atom`] - Atom properties (name, element, residue info, charges, etc.)
//! - [`Bond`] - Bond connectivity and properties
//! - [`CoordSet`] - Coordinate sets for multi-state molecules
//! - [`ObjectMolecule`] - Main molecular container
//!
//! # Architecture
//!
//! The crate follows PyMOL's original data model where:
//! - Atoms are stored in a flat array with residue/chain info inline
//! - Bonds reference atoms by index
//! - Multiple coordinate sets (states) can exist for trajectory/ensemble data
//! - Settings can be applied at per-atom/bond level via unique IDs
//!
//! # Example
//!
//! ```rust
//! use pymol_mol::{ObjectMolecule, Atom, BondOrder, Element, AtomIndex};
//! use lin_alg::f32::Vec3;
//!
//! // Create a simple water molecule
//! let mut mol = ObjectMolecule::new("water");
//!
//! // Add atoms
//! let o = mol.add_atom(Atom::new("O", Element::Oxygen));
//! let h1 = mol.add_atom(Atom::new("H1", Element::Hydrogen));
//! let h2 = mol.add_atom(Atom::new("H2", Element::Hydrogen));
//!
//! // Add bonds
//! mol.add_bond(o, h1, BondOrder::Single).unwrap();
//! mol.add_bond(o, h2, BondOrder::Single).unwrap();
//!
//! assert_eq!(mol.atom_count(), 3);
//! assert_eq!(mol.bond_count(), 2);
//! ```

// Module declarations
mod atom;
mod bond;
mod coordset;
pub mod dss;
mod element;
mod error;
mod flags;
mod index;
mod iterator;
mod molecule;
mod residue;
mod spatial;
mod secondary;

// Re-export main types
pub use atom::{Atom, AtomBuilder, AtomColors, AtomRepresentation, AtomResidue, AtomState, RepMask, COLOR_UNSET};
pub use bond::{Bond, BondOrder, BondStereo, SymOp};
pub use coordset::{
    mat4_to_ttt, rotation_matrix, rotation_ttt, translation_matrix, ttt_to_mat4, CoordSet,
    Symmetry,
};
pub use element::{Element, DEFAULT_VDW_RADIUS, ELEMENT_COUNT};
pub use error::{MolError, MolResult};
pub use flags::{AtomFlags, AtomGeometry, Chirality, Stereo};
pub use index::{AtomIndex, BondIndex, CoordIndex, StateIndex, INVALID_INDEX};
pub use molecule::{MoleculeBuilder, ObjectMolecule};
pub use residue::{
    atoms_same_chain, atoms_same_residue, atoms_same_segment, is_amino_acid, is_ion, is_lipid,
    is_nucleotide, is_water, nucleotide_to_char, residue_to_char, three_to_one, ChainIterator,
    ChainView, ResidueIterator, ResidueKey, ResidueView, AMINO_ACIDS, ION_NAMES, LIPID_NAMES,
    NUCLEOTIDES, WATER_NAMES,
};
pub use secondary::SecondaryStructure;

// Re-export DSS types
pub use dss::{assign_secondary_structure, DssSettings, PhiPsi};

/// Re-export commonly used types for convenience
pub mod prelude {
    pub use crate::atom::{Atom, AtomBuilder, AtomColors, AtomRepresentation, AtomResidue, AtomState};
    pub use crate::bond::{Bond, BondOrder};
    pub use crate::coordset::CoordSet;
    pub use crate::element::Element;
    pub use crate::error::{MolError, MolResult};
    pub use crate::flags::{AtomFlags, AtomGeometry};
    pub use crate::index::{AtomIndex, BondIndex, StateIndex};
    pub use crate::molecule::{MoleculeBuilder, ObjectMolecule};
    pub use crate::residue::{ChainView, ResidueKey, ResidueView};
    pub use crate::secondary::SecondaryStructure;
}

#[cfg(test)]
mod tests {
    use super::*;
    use lin_alg::f32::Vec3;

    #[test]
    fn test_create_molecule() {
        let mut mol = ObjectMolecule::new("test");

        let c1 = mol.add_atom(Atom::new("C1", Element::Carbon));
        let c2 = mol.add_atom(Atom::new("C2", Element::Carbon));
        mol.add_bond(c1, c2, BondOrder::Single).unwrap();

        assert_eq!(mol.atom_count(), 2);
        assert_eq!(mol.bond_count(), 1);
    }

    #[test]
    fn test_builder_pattern() {
        let mol = MoleculeBuilder::new("ethane")
            .add_atom(Atom::new("C1", Element::Carbon), Vec3::new(0.0, 0.0, 0.0))
            .add_atom(Atom::new("C2", Element::Carbon), Vec3::new(1.54, 0.0, 0.0))
            .add_bond(AtomIndex(0), AtomIndex(1), BondOrder::Single)
            .build();

        assert_eq!(mol.atom_count(), 2);
        assert!(mol.has_coords());
    }

    #[test]
    fn test_element_lookup() {
        assert_eq!(Element::from_symbol("C"), Some(Element::Carbon));
        assert_eq!(Element::Carbon.symbol(), "C");
        assert!((Element::Carbon.vdw_radius() - 1.70).abs() < 0.01);
    }

    #[test]
    fn test_residue_iteration() {
        use std::sync::Arc;

        let mut mol = ObjectMolecule::new("peptide");

        // Create a simple dipeptide
        for (resn, resv) in &[("ALA", 1), ("GLY", 2)] {
            let residue = Arc::new(AtomResidue::from_parts("A", *resn, *resv, ' ', ""));
            for name in &["N", "CA", "C", "O"] {
                let mut atom = Atom::new(*name, Element::Carbon);
                atom.residue = residue.clone();
                mol.add_atom(atom);
            }
        }

        let residue_count = mol.residues().count();
        assert_eq!(residue_count, 2);
    }
}
