//! MOL2 file writer
//!
//! Writes molecular structures in TRIPOS MOL2 format.

use std::io::Write;

use pymol_mol::{BondOrder, ObjectMolecule};

use crate::error::IoResult;
use crate::traits::MoleculeWriter;

/// MOL2 file writer
pub struct Mol2Writer<W> {
    writer: W,
    state: Option<usize>,
}

impl<W: Write> Mol2Writer<W> {
    /// Create a new MOL2 writer
    pub fn new(writer: W) -> Self {
        Mol2Writer {
            writer,
            state: None,
        }
    }

    /// Create a MOL2 writer that writes a specific state
    pub fn with_state(writer: W, state: usize) -> Self {
        Mol2Writer {
            writer,
            state: Some(state),
        }
    }

    /// Get SYBYL atom type for an atom
    fn get_sybyl_type(atom: &pymol_mol::Atom) -> String {
        // Simplified SYBYL type assignment
        // A full implementation would consider hybridization, bonds, etc.
        let base = atom.element.symbol();

        // Common types
        match atom.element {
            pymol_mol::Element::Carbon => "C.3".to_string(),
            pymol_mol::Element::Nitrogen => "N.3".to_string(),
            pymol_mol::Element::Oxygen => "O.3".to_string(),
            pymol_mol::Element::Sulfur => "S.3".to_string(),
            pymol_mol::Element::Phosphorus => "P.3".to_string(),
            pymol_mol::Element::Hydrogen => "H".to_string(),
            _ => base.to_string(),
        }
    }

    /// Write a single molecule
    fn write_molecule(&mut self, mol: &ObjectMolecule) -> IoResult<()> {
        let state = self.state.unwrap_or(0);

        // @<TRIPOS>MOLECULE
        writeln!(self.writer, "@<TRIPOS>MOLECULE")?;
        writeln!(self.writer, "{}", if mol.name.is_empty() { "molecule" } else { &mol.name })?;

        // Counts line
        let n_atoms = mol.atom_count();
        let n_bonds = mol.bond_count();
        let n_subst = mol.residue_count().max(1);
        writeln!(self.writer, " {} {} {} 0 0", n_atoms, n_bonds, n_subst)?;

        // Molecule type
        writeln!(self.writer, "SMALL")?;

        // Charge type
        let has_charges = mol.atoms().any(|a| a.partial_charge != 0.0);
        if has_charges {
            writeln!(self.writer, "USER_CHARGES")?;
        } else {
            writeln!(self.writer, "NO_CHARGES")?;
        }
        writeln!(self.writer)?;

        // @<TRIPOS>ATOM
        writeln!(self.writer, "@<TRIPOS>ATOM")?;

        for (idx, atom) in mol.atoms_indexed() {
            let coord = mol.get_coord(idx, state).unwrap_or_default();
            let sybyl_type = Self::get_sybyl_type(atom);

            // subst_name: residue name + residue number
            let subst_name = if atom.resn.is_empty() {
                format!("UNK{}", atom.resv.max(1))
            } else {
                format!("{}{}", atom.resn, atom.resv.max(1))
            };

            writeln!(
                self.writer,
                "{:7} {:<8} {:10.4} {:10.4} {:10.4} {:<8} {:3} {:<8} {:8.4}",
                idx.0 + 1,  // 1-indexed
                atom.name,
                coord.x,
                coord.y,
                coord.z,
                sybyl_type,
                atom.resv.max(1),
                subst_name,
                atom.partial_charge
            )?;
        }

        // @<TRIPOS>BOND
        writeln!(self.writer, "@<TRIPOS>BOND")?;

        for (bond_idx, bond) in mol.bonds_indexed() {
            let bond_type = match bond.order {
                BondOrder::Single => "1",
                BondOrder::Double => "2",
                BondOrder::Triple => "3",
                BondOrder::Aromatic => "ar",
                BondOrder::Unknown => "un",
            };

            writeln!(
                self.writer,
                "{:6} {:5} {:5} {}",
                bond_idx.0 + 1,  // 1-indexed
                bond.atom1.0 + 1,
                bond.atom2.0 + 1,
                bond_type
            )?;
        }

        Ok(())
    }
}

impl<W: Write> MoleculeWriter for Mol2Writer<W> {
    fn write(&mut self, mol: &ObjectMolecule) -> IoResult<()> {
        self.write_molecule(mol)
    }

    fn write_all(&mut self, molecules: &[ObjectMolecule]) -> IoResult<()> {
        for mol in molecules {
            self.write_molecule(mol)?;
        }
        Ok(())
    }

    fn flush(&mut self) -> IoResult<()> {
        self.writer.flush()?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use lin_alg::f32::Vec3;
    use pymol_mol::{Atom, CoordSet, Element};
    use crate::traits::MoleculeReader;

    fn create_water() -> ObjectMolecule {
        let mut mol = ObjectMolecule::new("water");

        let o = mol.add_atom(Atom::new("O", Element::Oxygen));
        let h1 = mol.add_atom(Atom::new("H1", Element::Hydrogen));
        let h2 = mol.add_atom(Atom::new("H2", Element::Hydrogen));

        mol.add_bond(o, h1, BondOrder::Single).unwrap();
        mol.add_bond(o, h2, BondOrder::Single).unwrap();

        let coords = CoordSet::from_vec3(&[
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.96, 0.0, 0.0),
            Vec3::new(-0.24, 0.93, 0.0),
        ]);
        mol.add_coord_set(coords);

        mol
    }

    #[test]
    fn test_write_mol2() {
        let mol = create_water();
        let mut output = Vec::new();

        {
            let mut writer = Mol2Writer::new(&mut output);
            writer.write(&mol).unwrap();
        }

        let mol2_string = String::from_utf8(output).unwrap();

        assert!(mol2_string.contains("@<TRIPOS>MOLECULE"));
        assert!(mol2_string.contains("water"));
        assert!(mol2_string.contains("@<TRIPOS>ATOM"));
        assert!(mol2_string.contains("@<TRIPOS>BOND"));
    }

    #[test]
    fn test_roundtrip() {
        let mol = create_water();
        let mut output = Vec::new();

        {
            let mut writer = Mol2Writer::new(&mut output);
            writer.write(&mol).unwrap();
        }

        // Parse it back
        let mut reader = crate::mol2::Mol2Reader::new(output.as_slice());
        let parsed = reader.read().unwrap();

        assert_eq!(parsed.atom_count(), mol.atom_count());
        assert_eq!(parsed.bond_count(), mol.bond_count());
    }
}
