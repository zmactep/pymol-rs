//! XYZ file writer
//!
//! Writes molecular structures in XYZ coordinate format.

use std::io::Write;

use pymol_mol::ObjectMolecule;

use crate::error::IoResult;
use crate::traits::MoleculeWriter;

/// XYZ file writer
pub struct XyzWriter<W> {
    writer: W,
    state: Option<usize>,
}

impl<W: Write> XyzWriter<W> {
    /// Create a new XYZ writer
    pub fn new(writer: W) -> Self {
        XyzWriter {
            writer,
            state: None,
        }
    }

    /// Create an XYZ writer that writes a specific state
    pub fn with_state(writer: W, state: usize) -> Self {
        XyzWriter {
            writer,
            state: Some(state),
        }
    }

    /// Write a single frame
    fn write_frame(&mut self, mol: &ObjectMolecule, state: usize) -> IoResult<()> {
        let n_atoms = mol.atom_count();

        // Line 1: Number of atoms
        writeln!(self.writer, "{}", n_atoms)?;

        // Line 2: Comment line (use title or name)
        let comment = if !mol.title.is_empty() {
            &mol.title
        } else if !mol.name.is_empty() {
            &mol.name
        } else {
            "molecule"
        };
        writeln!(self.writer, "{}", comment)?;

        // Atom lines
        for (idx, atom) in mol.atoms_indexed() {
            let coord = mol.get_coord(idx, state).unwrap_or_default();
            writeln!(
                self.writer,
                "{:2}  {:14.8}  {:14.8}  {:14.8}",
                atom.element.symbol(),
                coord.x,
                coord.y,
                coord.z
            )?;
        }

        Ok(())
    }

    /// Write a molecule (all states as trajectory or single state)
    fn write_molecule(&mut self, mol: &ObjectMolecule) -> IoResult<()> {
        let num_states = mol.state_count();

        if let Some(state) = self.state {
            // Write specific state
            let state = state.min(num_states.saturating_sub(1));
            self.write_frame(mol, state)?;
        } else {
            // Write all states (trajectory)
            for state in 0..num_states.max(1) {
                self.write_frame(mol, state)?;
            }
        }

        Ok(())
    }
}

impl<W: Write> MoleculeWriter for XyzWriter<W> {
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

        mol.add_atom(Atom::new("O", Element::Oxygen));
        mol.add_atom(Atom::new("H", Element::Hydrogen));
        mol.add_atom(Atom::new("H", Element::Hydrogen));

        let coords = CoordSet::from_vec3(&[
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.96, 0.0, 0.0),
            Vec3::new(-0.24, 0.93, 0.0),
        ]);
        mol.add_coord_set(coords);

        mol
    }

    #[test]
    fn test_write_xyz() {
        let mol = create_water();
        let mut output = Vec::new();

        {
            let mut writer = XyzWriter::new(&mut output);
            writer.write(&mol).unwrap();
        }

        let xyz_string = String::from_utf8(output).unwrap();

        // Check structure
        let lines: Vec<&str> = xyz_string.lines().collect();
        assert_eq!(lines[0], "3");
        assert_eq!(lines[1], "water");
        assert!(lines[2].starts_with("O ") || lines[2].starts_with(" O"));
        assert!(lines[3].starts_with("H ") || lines[3].starts_with(" H"));
        assert!(lines[4].starts_with("H ") || lines[4].starts_with(" H"));
    }

    #[test]
    fn test_roundtrip() {
        let mol = create_water();
        let mut output = Vec::new();

        {
            let mut writer = XyzWriter::new(&mut output);
            writer.write(&mol).unwrap();
        }

        // Parse it back
        let mut reader = crate::xyz::XyzReader::new(output.as_slice());
        let parsed = reader.read().unwrap();

        assert_eq!(parsed.atom_count(), mol.atom_count());
        assert_eq!(parsed.state_count(), mol.state_count());
    }
}
