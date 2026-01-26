//! SDF/MOL file writer
//!
//! Writes molecular structures in MDL SDF/MOL V2000 format.

use std::io::Write;

use pymol_mol::{BondOrder, ObjectMolecule};

use crate::error::{IoError, IoResult};
use crate::traits::MoleculeWriter;

/// SDF/MOL file writer
pub struct SdfWriter<W> {
    writer: W,
    state: Option<usize>,
}

impl<W: Write> SdfWriter<W> {
    /// Create a new SDF writer
    pub fn new(writer: W) -> Self {
        SdfWriter {
            writer,
            state: None,
        }
    }

    /// Create an SDF writer that writes a specific state
    pub fn with_state(writer: W, state: usize) -> Self {
        SdfWriter {
            writer,
            state: Some(state),
        }
    }

    /// Write a single molecule
    fn write_molecule(&mut self, mol: &ObjectMolecule) -> IoResult<()> {
        let state = self.state.unwrap_or(0);

        // Line 1: Molecule name
        writeln!(self.writer, "{}", mol.name)?;

        // Line 2: Program/timestamp (blank or informational)
        writeln!(self.writer, "  pymol-io          3D")?;

        // Line 3: Comment
        if mol.title.is_empty() {
            writeln!(self.writer)?;
        } else {
            writeln!(self.writer, "{}", &mol.title[..mol.title.len().min(80)])?;
        }

        // Line 4: Counts line
        let n_atoms = mol.atom_count();
        let n_bonds = mol.bond_count();

        if n_atoms > 999 || n_bonds > 999 {
            return Err(IoError::unsupported(
                "Molecule too large for V2000 format (max 999 atoms/bonds)",
            ));
        }

        writeln!(
            self.writer,
            "{:3}{:3}  0  0  0  0  0  0  0  0999 V2000",
            n_atoms, n_bonds
        )?;

        // Atom block
        for (idx, atom) in mol.atoms_indexed() {
            let coord = mol.get_coord(idx, state).unwrap_or_default();

            writeln!(
                self.writer,
                "{:10.4}{:10.4}{:10.4} {:<3} 0  0  0  0  0  0  0  0  0  0  0  0",
                coord.x,
                coord.y,
                coord.z,
                atom.element.symbol()
            )?;
        }

        // Bond block
        for (_, bond) in mol.bonds_indexed() {
            let bond_type = match bond.order {
                BondOrder::Single => 1,
                BondOrder::Double => 2,
                BondOrder::Triple => 3,
                BondOrder::Aromatic => 4,
                BondOrder::Unknown => 1,
            };

            writeln!(
                self.writer,
                "{:3}{:3}{:3}  0  0  0  0",
                bond.atom1.0 + 1, // Convert to 1-indexed
                bond.atom2.0 + 1,
                bond_type
            )?;
        }

        // Collect atoms with non-zero formal charges
        let charged_atoms: Vec<_> = mol
            .atoms_indexed()
            .filter(|(_, atom)| atom.formal_charge != 0)
            .collect();

        // Write M  CHG lines (max 8 per line)
        for chunk in charged_atoms.chunks(8) {
            write!(self.writer, "M  CHG{:3}", chunk.len())?;
            for (idx, atom) in chunk {
                write!(self.writer, " {:3} {:3}", idx.0 + 1, atom.formal_charge)?;
            }
            writeln!(self.writer)?;
        }

        // End of molecule
        writeln!(self.writer, "M  END")?;

        Ok(())
    }
}

impl<W: Write> MoleculeWriter for SdfWriter<W> {
    fn write(&mut self, mol: &ObjectMolecule) -> IoResult<()> {
        self.write_molecule(mol)?;
        // Write SDF separator
        writeln!(self.writer, "$$$$")?;
        Ok(())
    }

    fn write_all(&mut self, molecules: &[ObjectMolecule]) -> IoResult<()> {
        for mol in molecules {
            self.write_molecule(mol)?;
            writeln!(self.writer, "$$$$")?;
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
        let h1 = mol.add_atom(Atom::new("H", Element::Hydrogen));
        let h2 = mol.add_atom(Atom::new("H", Element::Hydrogen));

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
    fn test_write_sdf() {
        let mol = create_water();
        let mut output = Vec::new();

        {
            let mut writer = SdfWriter::new(&mut output);
            writer.write(&mol).unwrap();
        }

        let sdf_string = String::from_utf8(output).unwrap();

        assert!(sdf_string.contains("water"));
        assert!(sdf_string.contains("V2000"));
        assert!(sdf_string.contains("M  END"));
        assert!(sdf_string.contains("$$$$"));
        assert!(sdf_string.contains("  3  2")); // 3 atoms, 2 bonds
    }

    #[test]
    fn test_roundtrip() {
        let mol = create_water();
        let mut output = Vec::new();

        {
            let mut writer = SdfWriter::new(&mut output);
            writer.write(&mol).unwrap();
        }

        // Parse it back
        let mut reader = crate::sdf::SdfReader::new(output.as_slice());
        let parsed = reader.read().unwrap();

        assert_eq!(parsed.atom_count(), mol.atom_count());
        assert_eq!(parsed.bond_count(), mol.bond_count());
    }
}
