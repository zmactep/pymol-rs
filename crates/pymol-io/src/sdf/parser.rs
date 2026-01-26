//! SDF/MOL file parser
//!
//! Parses MDL SDF/MOL V2000 format files.

use std::io::{BufRead, BufReader, Read};

use lin_alg::f32::Vec3;
use pymol_mol::{Atom, BondOrder, CoordSet, Element, ObjectMolecule};

use crate::error::{IoError, IoResult};
use crate::traits::MoleculeReader;

/// SDF/MOL file reader
pub struct SdfReader<R> {
    reader: BufReader<R>,
    line_number: usize,
}

impl<R: Read> SdfReader<R> {
    /// Create a new SDF reader
    pub fn new(reader: R) -> Self {
        SdfReader {
            reader: BufReader::new(reader),
            line_number: 0,
        }
    }

    /// Read a single line from the file
    fn read_line(&mut self) -> IoResult<Option<String>> {
        let mut line = String::new();
        match self.reader.read_line(&mut line) {
            Ok(0) => Ok(None),
            Ok(_) => {
                self.line_number += 1;
                // Remove trailing newline
                if line.ends_with('\n') {
                    line.pop();
                    if line.ends_with('\r') {
                        line.pop();
                    }
                }
                Ok(Some(line))
            }
            Err(e) => Err(IoError::Io(e)),
        }
    }

    /// Parse a single molecule from the file
    fn parse_molecule(&mut self) -> IoResult<Option<ObjectMolecule>> {
        // Line 1: Molecule name
        let name = match self.read_line()? {
            Some(line) => line.trim().to_string(),
            None => return Ok(None),
        };

        // Check for end of SDF file
        if name == "$$$$" {
            return self.parse_molecule(); // Try next molecule
        }

        // Line 2: Program/timestamp line (ignored)
        self.read_line()?;

        // Line 3: Comment line (ignored)
        self.read_line()?;

        // Line 4: Counts line
        let counts_line = match self.read_line()? {
            Some(line) => line,
            None => return Err(IoError::parse(self.line_number, "Unexpected end of file")),
        };

        let (n_atoms, n_bonds, version) = parse_counts_line(&counts_line)?;

        // Check for V3000 format
        if version == "V3000" {
            return Err(IoError::unsupported("V3000 format not yet supported"));
        }

        let mut mol = ObjectMolecule::with_capacity(&name, n_atoms as usize, n_bonds as usize);
        mol.name = name;

        let mut coords = Vec::with_capacity(n_atoms as usize);

        // Parse atom block
        for i in 0..n_atoms {
            let line = match self.read_line()? {
                Some(line) => line,
                None => {
                    return Err(IoError::parse(
                        self.line_number,
                        format!("Expected atom {}, got end of file", i + 1),
                    ))
                }
            };

            let (atom, coord) = parse_atom_line(&line, self.line_number)?;
            mol.add_atom(atom);
            coords.push(coord);
        }

        // Parse bond block
        for i in 0..n_bonds {
            let line = match self.read_line()? {
                Some(line) => line,
                None => {
                    return Err(IoError::parse(
                        self.line_number,
                        format!("Expected bond {}, got end of file", i + 1),
                    ))
                }
            };

            let (atom1, atom2, order) = parse_bond_line(&line, self.line_number)?;

            // Convert from 1-indexed to 0-indexed
            let idx1 = pymol_mol::AtomIndex((atom1 - 1) as u32);
            let idx2 = pymol_mol::AtomIndex((atom2 - 1) as u32);

            let _ = mol.add_bond_unchecked(idx1, idx2, order);
        }

        // Parse properties block (M lines)
        let mut charges: Vec<(u32, i8)> = Vec::new();

        loop {
            let line = match self.read_line()? {
                Some(line) => line,
                None => break,
            };

            if line.starts_with("M  END") {
                break;
            }

            if line.starts_with("M  CHG") {
                // Parse charge line: M  CHG  n  aaa vvv ...
                if let Some(charge_data) = parse_charge_line(&line) {
                    charges.extend(charge_data);
                }
            }
            // Other M lines (ISO, RAD, etc.) could be parsed here
        }

        // Apply charges
        for (atom_idx, charge) in charges {
            if let Some(atom) = mol.get_atom_mut(pymol_mol::AtomIndex(atom_idx)) {
                atom.formal_charge = charge;
            }
        }

        // Add coordinate set
        let coord_set = CoordSet::from_vec3(&coords);
        mol.add_coord_set(coord_set);

        // Skip to end of molecule (for SDF files)
        loop {
            let line = match self.read_line()? {
                Some(line) => line,
                None => break,
            };

            if line.starts_with("$$$$") {
                break;
            }
            // Skip data fields in SDF format
        }

        Ok(Some(mol))
    }
}

impl<R: Read> MoleculeReader for SdfReader<R> {
    fn read(&mut self) -> IoResult<ObjectMolecule> {
        match self.parse_molecule()? {
            Some(mol) => Ok(mol),
            None => Err(IoError::EmptyFile),
        }
    }

    fn read_all(&mut self) -> IoResult<Vec<ObjectMolecule>> {
        let mut molecules = Vec::new();
        while let Some(mol) = self.parse_molecule()? {
            molecules.push(mol);
        }
        if molecules.is_empty() {
            Err(IoError::EmptyFile)
        } else {
            Ok(molecules)
        }
    }
}

/// Parse the counts line (line 4 of MOL file)
fn parse_counts_line(line: &str) -> IoResult<(u32, u32, String)> {
    // Format: aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv
    // aaa = number of atoms (3 chars)
    // bbb = number of bonds (3 chars)
    // ... other fields ...
    // vvvvvv = version (V2000 or V3000)

    if line.len() < 6 {
        return Err(IoError::parse_msg("Counts line too short"));
    }

    let n_atoms: u32 = line
        .get(0..3)
        .and_then(|s| s.trim().parse().ok())
        .ok_or_else(|| IoError::parse_msg("Invalid atom count"))?;

    let n_bonds: u32 = line
        .get(3..6)
        .and_then(|s| s.trim().parse().ok())
        .ok_or_else(|| IoError::parse_msg("Invalid bond count"))?;

    // Check for version string
    let version = if line.contains("V3000") {
        "V3000".to_string()
    } else {
        "V2000".to_string()
    };

    Ok((n_atoms, n_bonds, version))
}

/// Parse an atom line
fn parse_atom_line(line: &str, line_number: usize) -> IoResult<(Atom, Vec3)> {
    // Format: xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee
    // x, y, z = coordinates (10.4 format each)
    // aaa = atom symbol (3 chars)
    // dd = mass difference
    // ccc = charge (0-7, mapped: 0=0, 1=+3, 2=+2, 3=+1, 4=doublet, 5=-1, 6=-2, 7=-3)
    // ... other fields ...

    if line.len() < 34 {
        return Err(IoError::parse(line_number, "Atom line too short"));
    }

    let x: f32 = line
        .get(0..10)
        .and_then(|s| s.trim().parse().ok())
        .ok_or_else(|| IoError::parse(line_number, "Invalid x coordinate"))?;

    let y: f32 = line
        .get(10..20)
        .and_then(|s| s.trim().parse().ok())
        .ok_or_else(|| IoError::parse(line_number, "Invalid y coordinate"))?;

    let z: f32 = line
        .get(20..30)
        .and_then(|s| s.trim().parse().ok())
        .ok_or_else(|| IoError::parse(line_number, "Invalid z coordinate"))?;

    let symbol = line
        .get(31..34)
        .map(|s| s.trim())
        .ok_or_else(|| IoError::parse(line_number, "Invalid atom symbol"))?;

    let element = Element::from_symbol(symbol).unwrap_or(Element::Unknown);
    let atom = Atom::new(symbol, element);

    // Parse charge field if present
    // Note: Actual charge is typically specified in M  CHG line, not here

    Ok((atom, Vec3::new(x, y, z)))
}

/// Parse a bond line
fn parse_bond_line(line: &str, line_number: usize) -> IoResult<(u32, u32, BondOrder)> {
    // Format: 111222tttsssxxxrrrccc
    // 111 = first atom number (3 chars)
    // 222 = second atom number (3 chars)
    // ttt = bond type (1=single, 2=double, 3=triple, 4=aromatic, etc.)
    // sss = bond stereo (not used here)

    if line.len() < 9 {
        return Err(IoError::parse(line_number, "Bond line too short"));
    }

    let atom1: u32 = line
        .get(0..3)
        .and_then(|s| s.trim().parse().ok())
        .ok_or_else(|| IoError::parse(line_number, "Invalid first atom in bond"))?;

    let atom2: u32 = line
        .get(3..6)
        .and_then(|s| s.trim().parse().ok())
        .ok_or_else(|| IoError::parse(line_number, "Invalid second atom in bond"))?;

    let bond_type: u32 = line
        .get(6..9)
        .and_then(|s| s.trim().parse().ok())
        .unwrap_or(1);

    let order = match bond_type {
        1 => BondOrder::Single,
        2 => BondOrder::Double,
        3 => BondOrder::Triple,
        4 => BondOrder::Aromatic,
        _ => BondOrder::Single,
    };

    Ok((atom1, atom2, order))
}

/// Parse an M  CHG charge line
fn parse_charge_line(line: &str) -> Option<Vec<(u32, i8)>> {
    // Format: M  CHG  n aaa vvv aaa vvv ...
    // n = count of charge specifications
    // aaa = atom number
    // vvv = charge value

    if !line.starts_with("M  CHG") {
        return None;
    }

    let parts: Vec<&str> = line.split_whitespace().collect();
    if parts.len() < 4 {
        return None;
    }

    let count: usize = parts.get(2)?.parse().ok()?;
    let mut charges = Vec::with_capacity(count);

    for i in 0..count {
        let atom_idx = 3 + i * 2;
        let value_idx = 4 + i * 2;

        if value_idx >= parts.len() {
            break;
        }

        let atom: u32 = parts.get(atom_idx)?.parse().ok()?;
        let value: i8 = parts.get(value_idx)?.parse().ok()?;

        // Convert to 0-indexed
        charges.push((atom - 1, value));
    }

    Some(charges)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_counts_line() {
        let line = "  3  2  0  0  0  0  0  0  0  0999 V2000";
        let (atoms, bonds, version) = parse_counts_line(line).unwrap();
        assert_eq!(atoms, 3);
        assert_eq!(bonds, 2);
        assert_eq!(version, "V2000");
    }

    #[test]
    fn test_parse_atom_line() {
        let line = "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0";
        let (atom, coord) = parse_atom_line(line, 1).unwrap();
        assert_eq!(atom.element, Element::Carbon);
        assert!((coord.x).abs() < 0.0001);
    }

    #[test]
    fn test_parse_bond_line() {
        let line = "  1  2  1  0  0  0  0";
        let (a1, a2, order) = parse_bond_line(line, 1).unwrap();
        assert_eq!(a1, 1);
        assert_eq!(a2, 2);
        assert_eq!(order, BondOrder::Single);
    }

    #[test]
    fn test_read_simple_mol() {
        let mol_data = r#"Water
  test
Comment
  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.9572    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2400    0.9266    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
M  END
"#;

        let mut reader = SdfReader::new(mol_data.as_bytes());
        let mol = reader.read().unwrap();

        assert_eq!(mol.name, "Water");
        assert_eq!(mol.atom_count(), 3);
        assert_eq!(mol.bond_count(), 2);
    }

    #[test]
    fn test_read_sdf_multiple() {
        let sdf_data = r#"Mol1
  test
Comment
  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
Mol2
  test
Comment
  1  0  0  0  0  0  0  0  0  0999 V2000
    1.0000    1.0000    1.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
"#;

        let mut reader = SdfReader::new(sdf_data.as_bytes());
        let molecules = reader.read_all().unwrap();

        assert_eq!(molecules.len(), 2);
        assert_eq!(molecules[0].name, "Mol1");
        assert_eq!(molecules[1].name, "Mol2");
    }
}
