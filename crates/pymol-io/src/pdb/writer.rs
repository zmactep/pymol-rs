//! PDB file writer
//!
//! Writes molecular structures in PDB format.

use std::io::Write;

use pymol_mol::{AtomIndex, ObjectMolecule};

use crate::error::{IoError, IoResult};
use crate::traits::MoleculeWriter;

/// PDB file writer
pub struct PdbWriter<W> {
    writer: W,
    state: Option<usize>,
}

impl<W: Write> PdbWriter<W> {
    /// Create a new PDB writer
    pub fn new(writer: W) -> Self {
        PdbWriter {
            writer,
            state: None,
        }
    }

    /// Create a PDB writer that writes a specific state
    pub fn with_state(writer: W, state: usize) -> Self {
        PdbWriter {
            writer,
            state: Some(state),
        }
    }

    /// Write title record
    fn write_title(&mut self, title: &str) -> IoResult<()> {
        if !title.is_empty() {
            // Split title into 70-character chunks
            for (i, chunk) in title.as_bytes().chunks(70).enumerate() {
                let cont = if i > 0 {
                    format!("{:>2} ", i + 1)
                } else {
                    "   ".to_string()
                };
                writeln!(
                    self.writer,
                    "TITLE  {}{}",
                    cont,
                    String::from_utf8_lossy(chunk)
                )?;
            }
        }
        Ok(())
    }

    /// Write CRYST1 record
    fn write_cryst1(&mut self, mol: &ObjectMolecule) -> IoResult<()> {
        if let Some(ref sym) = mol.symmetry {
            writeln!(
                self.writer,
                "CRYST1{:9.3}{:9.3}{:9.3}{:7.2}{:7.2}{:7.2} {:11}{:4}",
                sym.cell_lengths[0],
                sym.cell_lengths[1],
                sym.cell_lengths[2],
                sym.cell_angles[0],
                sym.cell_angles[1],
                sym.cell_angles[2],
                sym.space_group,
                sym.z_value
            )?;
        }
        Ok(())
    }

    /// Write ATOM/HETATM record
    fn write_atom(
        &mut self,
        serial: i32,
        atom: &pymol_mol::Atom,
        x: f32,
        y: f32,
        z: f32,
    ) -> IoResult<()> {
        let record_type = if atom.state.hetatm { "HETATM" } else { "ATOM  " };

        // Format atom name with proper justification
        let name = format_atom_name(&atom.name, atom.element.symbol());

        // Format chain (single character)
        let chain = if atom.residue.chain.is_empty() {
            " "
        } else {
            &atom.residue.chain[..1.min(atom.residue.chain.len())]
        };

        // Format insertion code
        let icode = if atom.residue.inscode == ' ' {
            ' '
        } else {
            atom.residue.inscode
        };

        // Format element symbol (right-justified in 2 characters)
        let element = format!("{:>2}", atom.element.symbol());

        // Format charge
        let charge = format_charge(atom.formal_charge);

        writeln!(
            self.writer,
            "{}{:5} {:4}{}{:3} {}{:4}{}   {:8.3}{:8.3}{:8.3}{:6.2}{:6.2}          {}{}",
            record_type,
            serial % 100000,
            name,
            if atom.alt == ' ' { ' ' } else { atom.alt },
            if atom.residue.resn.len() > 3 {
                &atom.residue.resn[..3]
            } else {
                &atom.residue.resn
            },
            chain,
            atom.residue.resv,
            icode,
            x,
            y,
            z,
            atom.occupancy,
            atom.b_factor,
            element,
            charge
        )?;

        Ok(())
    }

    /// Write TER record
    fn write_ter(&mut self, serial: i32, resn: &str, chain: &str, resv: i32) -> IoResult<()> {
        writeln!(
            self.writer,
            "TER   {:5}      {:3} {:1}{:4}",
            serial % 100000,
            if resn.len() > 3 { &resn[..3] } else { resn },
            if chain.is_empty() {
                " "
            } else {
                &chain[..1.min(chain.len())]
            },
            resv
        )?;
        Ok(())
    }

    /// Write CONECT records
    fn write_conect(&mut self, mol: &ObjectMolecule) -> IoResult<()> {
        // Group bonds by atom
        let mut bonds_by_atom: std::collections::HashMap<i32, Vec<i32>> =
            std::collections::HashMap::new();

        for (_, bond) in mol.bonds_indexed() {
            if let (Some(atom1), Some(atom2)) =
                (mol.get_atom(bond.atom1), mol.get_atom(bond.atom2))
            {
                bonds_by_atom
                    .entry(atom1.id)
                    .or_default()
                    .push(atom2.id);
                bonds_by_atom
                    .entry(atom2.id)
                    .or_default()
                    .push(atom1.id);
            }
        }

        // Write CONECT records
        let mut serials: Vec<i32> = bonds_by_atom.keys().copied().collect();
        serials.sort();

        for serial in serials {
            if let Some(bonded) = bonds_by_atom.get(&serial) {
                // Write in chunks of 4 bonded atoms per CONECT record
                for chunk in bonded.chunks(4) {
                    write!(self.writer, "CONECT{:5}", serial)?;
                    for &bonded_serial in chunk {
                        write!(self.writer, "{:5}", bonded_serial)?;
                    }
                    writeln!(self.writer)?;
                }
            }
        }

        Ok(())
    }

    /// Write END record
    fn write_end(&mut self) -> IoResult<()> {
        writeln!(self.writer, "END")?;
        Ok(())
    }

    /// Write a molecule
    fn write_molecule(&mut self, mol: &ObjectMolecule) -> IoResult<()> {
        // Write header records
        self.write_title(&mol.title)?;
        self.write_cryst1(mol)?;

        let state = self.state.unwrap_or(0);
        let num_states = mol.state_count();

        if num_states == 0 {
            return Err(IoError::EmptyFile);
        }

        // Write multi-model file if more than one state and no specific state requested
        let write_models = num_states > 1 && self.state.is_none();

        let states_to_write: Vec<usize> = if write_models {
            (0..num_states).collect()
        } else {
            vec![state.min(num_states.saturating_sub(1))]
        };

        for (model_num, state_idx) in states_to_write.iter().enumerate() {
            if write_models {
                writeln!(self.writer, "MODEL     {:4}", model_num + 1)?;
            }

            let mut serial = 1;
            let mut last_chain = String::new();

            for (atom_idx, atom) in mol.atoms_indexed() {
                // Get coordinates for this state
                let coord = mol.get_coord(atom_idx, *state_idx).unwrap_or_default();

                // Write TER between chains
                if !last_chain.is_empty() && atom.residue.chain != last_chain {
                    // Get the previous atom's info for TER record
                    if let Some(prev_atom) = mol.get_atom(AtomIndex(atom_idx.0.saturating_sub(1))) {
                        self.write_ter(serial, &prev_atom.residue.resn, &prev_atom.residue.chain, prev_atom.residue.resv)?;
                        serial += 1;
                    }
                }
                last_chain = atom.residue.chain.clone();

                self.write_atom(serial, atom, coord.x, coord.y, coord.z)?;
                serial += 1;
            }

            // Write final TER
            if let Some(last_atom) = mol
                .atom_count()
                .checked_sub(1)
                .and_then(|i| mol.get_atom(AtomIndex(i as u32)))
            {
                self.write_ter(serial, &last_atom.residue.resn, &last_atom.residue.chain, last_atom.residue.resv)?;
            }

            if write_models {
                writeln!(self.writer, "ENDMDL")?;
            }
        }

        // Write CONECT records (only for single model)
        if !write_models && mol.bond_count() > 0 {
            self.write_conect(mol)?;
        }

        self.write_end()?;

        Ok(())
    }
}

impl<W: Write> MoleculeWriter for PdbWriter<W> {
    fn write(&mut self, mol: &ObjectMolecule) -> IoResult<()> {
        self.write_molecule(mol)
    }

    fn flush(&mut self) -> IoResult<()> {
        self.writer.flush()?;
        Ok(())
    }
}

/// Format atom name according to PDB conventions
fn format_atom_name(name: &str, element: &str) -> String {
    let name = name.trim();
    let element_len = element.len();

    // PDB convention: 1-letter elements start in column 14 (index 1 of 4-char field)
    // 2-letter elements start in column 13 (index 0 of 4-char field)
    if name.len() >= 4 {
        name[..4].to_string()
    } else if element_len == 1 && !name.starts_with(char::is_numeric) {
        format!(" {:<3}", name)
    } else {
        format!("{:<4}", name)
    }
}

/// Format charge for PDB format (e.g., "2+" or "1-")
fn format_charge(charge: i8) -> String {
    match charge {
        0 => "  ".to_string(),
        c if c > 0 => format!("{}+", c.abs()),
        c => format!("{}-", c.abs()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use pymol_mol::{Atom, CoordSet, Element};
    use lin_alg::f32::Vec3;

    fn create_test_molecule() -> ObjectMolecule {
        let mut mol = ObjectMolecule::new("test");

        let mut atom1 = Atom::new("N", Element::Nitrogen);
        atom1.set_residue("ALA", 1, "A");
        atom1.id = 1;
        mol.add_atom(atom1);

        let mut atom2 = Atom::new("CA", Element::Carbon);
        atom2.set_residue("ALA", 1, "A");
        atom2.id = 2;
        mol.add_atom(atom2);

        let coords = CoordSet::from_vec3(&[Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.5, 0.0, 0.0)]);
        mol.add_coord_set(coords);

        mol
    }

    #[test]
    fn test_format_atom_name() {
        assert_eq!(format_atom_name("N", "N"), " N  ");
        assert_eq!(format_atom_name("CA", "C"), " CA ");
        assert_eq!(format_atom_name("FE", "Fe"), "FE  ");
    }

    #[test]
    fn test_format_charge() {
        assert_eq!(format_charge(0), "  ");
        assert_eq!(format_charge(1), "1+");
        assert_eq!(format_charge(2), "2+");
        assert_eq!(format_charge(-1), "1-");
        assert_eq!(format_charge(-2), "2-");
    }

    #[test]
    fn test_write_pdb() {
        let mol = create_test_molecule();
        let mut output = Vec::new();

        {
            let mut writer = PdbWriter::new(&mut output);
            writer.write(&mol).unwrap();
        }

        let pdb_string = String::from_utf8(output).unwrap();
        assert!(pdb_string.contains("ATOM"));
        assert!(pdb_string.contains("ALA"));
        assert!(pdb_string.contains("END"));
    }
}
