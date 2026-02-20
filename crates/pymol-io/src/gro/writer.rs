//! GRO format writer (GROMACS structure file)
//!
//! Writes molecular structures in GROMACS GRO coordinate format.
//! Coordinates are converted from Angstroms to nanometers (÷10).

use std::io::Write;

use pymol_mol::ObjectMolecule;

use crate::error::IoResult;
use crate::traits::MoleculeWriter;

/// Angstroms to nanometers conversion factor
const ANGSTROM_TO_NM: f32 = 0.1;

/// GRO file writer
pub struct GroWriter<W> {
    writer: W,
    state: Option<usize>,
}

impl<W: Write> GroWriter<W> {
    /// Create a new GRO writer
    pub fn new(writer: W) -> Self {
        GroWriter {
            writer,
            state: None,
        }
    }

    /// Create a GRO writer that writes a specific state
    pub fn with_state(writer: W, state: usize) -> Self {
        GroWriter {
            writer,
            state: Some(state),
        }
    }

    /// Write a single frame
    fn write_frame(&mut self, mol: &ObjectMolecule, state: usize) -> IoResult<()> {
        // Title line
        let title = if !mol.title.is_empty() {
            &mol.title
        } else if !mol.name.is_empty() {
            &mol.name
        } else {
            "molecule"
        };
        writeln!(self.writer, "{}", title)?;

        // Atom count
        let n_atoms = mol.atom_count();
        writeln!(self.writer, "{:5}", n_atoms)?;

        // Atom lines: %5d%-5s%5s%5d%8.3f%8.3f%8.3f
        let mut atom_serial: i32 = 1;
        for (idx, atom) in mol.atoms_indexed() {
            let coord = mol.get_coord(idx, state).unwrap_or_default();

            // Residue number wraps at 99999
            let resid = ((atom.residue.resv - 1).rem_euclid(99999)) + 1;

            // Atom serial wraps at 99999
            let serial = ((atom_serial - 1) % 99999) + 1;

            writeln!(
                self.writer,
                "{:5}{:<5}{:>5}{:5}{:8.3}{:8.3}{:8.3}",
                resid,
                truncate(&atom.residue.resn, 5),
                truncate(&atom.name, 5),
                serial,
                coord.x * ANGSTROM_TO_NM,
                coord.y * ANGSTROM_TO_NM,
                coord.z * ANGSTROM_TO_NM,
            )?;

            atom_serial += 1;
        }

        // Box vectors (last line)
        if let Some(sym) = &mol.symmetry {
            let [a, b, c] = sym.cell_lengths;
            writeln!(
                self.writer,
                "{:10.5}{:10.5}{:10.5}",
                a * ANGSTROM_TO_NM,
                b * ANGSTROM_TO_NM,
                c * ANGSTROM_TO_NM,
            )?;
        } else {
            writeln!(self.writer, "   0.00000   0.00000   0.00000")?;
        }

        Ok(())
    }

    /// Write a molecule (all states as trajectory or single state)
    fn write_molecule(&mut self, mol: &ObjectMolecule) -> IoResult<()> {
        let num_states = mol.state_count();

        if let Some(state) = self.state {
            let state = state.min(num_states.saturating_sub(1));
            self.write_frame(mol, state)?;
        } else {
            for state in 0..num_states.max(1) {
                self.write_frame(mol, state)?;
            }
        }

        Ok(())
    }
}

impl<W: Write> MoleculeWriter for GroWriter<W> {
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

/// Truncate a string to at most `max` bytes, respecting UTF-8 char boundaries.
fn truncate(s: &str, max: usize) -> &str {
    if s.len() <= max {
        return s;
    }
    // Find the last char boundary at or before `max`
    let mut end = max;
    while end > 0 && !s.is_char_boundary(end) {
        end -= 1;
    }
    &s[..end]
}

#[cfg(test)]
mod tests {
    use super::*;
    use lin_alg::f32::Vec3;
    use pymol_mol::{Atom, AtomIndex, CoordSet, Element, Symmetry};

    fn create_water() -> ObjectMolecule {
        let mut mol = ObjectMolecule::new("water");

        let mut o = Atom::new("OW", Element::Oxygen);
        o.set_residue("SOL", 1, "");
        mol.add_atom(o);

        let mut h1 = Atom::new("HW1", Element::Hydrogen);
        h1.set_residue("SOL", 1, "");
        mol.add_atom(h1);

        let mut h2 = Atom::new("HW2", Element::Hydrogen);
        h2.set_residue("SOL", 1, "");
        mol.add_atom(h2);

        // Coordinates in Angstroms
        let coords = CoordSet::from_vec3(&[
            Vec3::new(2.30, 6.28, 1.13),
            Vec3::new(1.37, 6.28, 1.50),
            Vec3::new(2.31, 5.89, 0.21),
        ]);
        mol.add_coord_set(coords);

        mol.symmetry = Some(Symmetry::new("P 1", [10.0, 10.0, 10.0], [90.0, 90.0, 90.0]));

        mol
    }

    #[test]
    fn test_write_gro() {
        let mol = create_water();
        let mut output = Vec::new();

        {
            let mut writer = GroWriter::new(&mut output);
            writer.write(&mol).unwrap();
        }

        let gro_string = String::from_utf8(output).unwrap();
        let lines: Vec<&str> = gro_string.lines().collect();

        // Title
        assert_eq!(lines[0], "water");
        // Atom count (right-justified in 5 chars)
        assert_eq!(lines[1].trim(), "3");
        // First atom line should contain SOL, OW, and nm coordinates
        assert!(lines[2].contains("SOL"));
        assert!(lines[2].contains("OW"));
        // Box vectors
        assert_eq!(lines.len(), 6); // title + count + 3 atoms + box
        assert!(lines[5].contains("1.00000"));
    }

    #[test]
    fn test_coordinate_conversion() {
        let mol = create_water();
        let mut output = Vec::new();

        {
            let mut writer = GroWriter::new(&mut output);
            writer.write(&mol).unwrap();
        }

        let gro_string = String::from_utf8(output).unwrap();
        let lines: Vec<&str> = gro_string.lines().collect();

        // First atom: OW at (2.30, 6.28, 1.13) Å → (0.230, 0.628, 0.113) nm
        let atom_line = lines[2];
        // Coordinates start at column 20
        let x: f32 = atom_line[20..28].trim().parse().unwrap();
        let y: f32 = atom_line[28..36].trim().parse().unwrap();
        let z: f32 = atom_line[36..44].trim().parse().unwrap();
        assert!((x - 0.230).abs() < 0.001);
        assert!((y - 0.628).abs() < 0.001);
        assert!((z - 0.113).abs() < 0.001);
    }

    #[test]
    fn test_roundtrip() {
        let mol = create_water();
        let mut output = Vec::new();

        {
            let mut writer = GroWriter::new(&mut output);
            writer.write(&mol).unwrap();
        }

        // Parse it back
        let mut reader = crate::gro::GroReader::new(output.as_slice());
        use crate::traits::MoleculeReader;
        let parsed = reader.read().unwrap();

        assert_eq!(parsed.atom_count(), mol.atom_count());
        assert_eq!(parsed.state_count(), mol.state_count());

        // Check coordinates survive roundtrip (within GRO precision: 3 decimal places in nm)
        for i in 0..mol.atom_count() {
            let idx = AtomIndex::from(i);
            let orig = mol.get_coord(idx, 0).unwrap();
            let rt = parsed.get_coord(idx, 0).unwrap();
            assert!((orig.x - rt.x).abs() < 0.1, "x mismatch for atom {}", i);
            assert!((orig.y - rt.y).abs() < 0.1, "y mismatch for atom {}", i);
            assert!((orig.z - rt.z).abs() < 0.1, "z mismatch for atom {}", i);
        }
    }

    #[test]
    fn test_truncate() {
        assert_eq!(truncate("hello", 5), "hello");
        assert_eq!(truncate("hello world", 5), "hello");
        assert_eq!(truncate("hi", 5), "hi");
        assert_eq!(truncate("", 5), "");
    }
}
