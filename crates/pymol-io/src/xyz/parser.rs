//! XYZ file parser
//!
//! Parses simple XYZ coordinate format files.

use std::io::{BufRead, BufReader, Read};

use lin_alg::f32::Vec3;
use pymol_mol::{Atom, CoordSet, Element, ObjectMolecule};

use crate::error::{IoError, IoResult};
use crate::traits::MoleculeReader;

/// XYZ file reader
pub struct XyzReader<R> {
    reader: BufReader<R>,
    line_number: usize,
}

impl<R: Read> XyzReader<R> {
    /// Create a new XYZ reader
    pub fn new(reader: R) -> Self {
        XyzReader {
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

    /// Parse a single frame/molecule from the file
    fn parse_frame(&mut self) -> IoResult<Option<ObjectMolecule>> {
        // Line 1: Number of atoms
        let n_atoms_line = match self.read_line()? {
            Some(line) => line,
            None => return Ok(None),
        };

        let n_atoms: usize = n_atoms_line
            .trim()
            .parse()
            .map_err(|_| IoError::parse(self.line_number, "Invalid atom count"))?;

        if n_atoms == 0 {
            return Err(IoError::parse(self.line_number, "Zero atoms in XYZ file"));
        }

        // Line 2: Comment/title line
        let title = match self.read_line()? {
            Some(line) => line.trim().to_string(),
            None => {
                return Err(IoError::parse(
                    self.line_number,
                    "Expected comment line after atom count",
                ))
            }
        };

        let mut mol = ObjectMolecule::with_capacity("", n_atoms, 0);
        mol.title = title.clone();
        mol.name = if title.is_empty() {
            "molecule".to_string()
        } else {
            // Use first word of title as name
            title.split_whitespace().next().unwrap_or("molecule").to_string()
        };

        let mut coords = Vec::with_capacity(n_atoms);

        // Parse atom lines
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

        // Add coordinate set
        let coord_set = CoordSet::from_vec3(&coords);
        mol.add_coord_set(coord_set);

        Ok(Some(mol))
    }
}

impl<R: Read> MoleculeReader for XyzReader<R> {
    fn read(&mut self) -> IoResult<ObjectMolecule> {
        // Try to read all frames and merge them as coordinate sets
        let mut mol = match self.parse_frame()? {
            Some(mol) => mol,
            None => return Err(IoError::EmptyFile),
        };

        // Read additional frames as coordinate sets (trajectory)
        while let Some(frame) = self.parse_frame()? {
            if frame.atom_count() != mol.atom_count() {
                // Different number of atoms - treat as separate molecule
                // For now, just stop reading
                break;
            }
            // Add coordinates as new state
            if let Some(cs) = frame.get_coord_set(0) {
                mol.add_coord_set(cs.clone());
            }
        }

        Ok(mol)
    }

    fn read_all(&mut self) -> IoResult<Vec<ObjectMolecule>> {
        // For XYZ, we typically want trajectory mode where frames become states
        // read() already handles this
        let mol = self.read()?;
        Ok(vec![mol])
    }
}

/// Parse an XYZ atom line
fn parse_atom_line(line: &str, line_number: usize) -> IoResult<(Atom, Vec3)> {
    // Format: element x y z [optional additional columns]
    let parts: Vec<&str> = line.split_whitespace().collect();

    if parts.len() < 4 {
        return Err(IoError::parse(
            line_number,
            format!("Atom line too short: expected 'element x y z', got '{}'", line),
        ));
    }

    let symbol = parts[0];
    let x: f32 = parts[1]
        .parse()
        .map_err(|_| IoError::parse(line_number, "Invalid x coordinate"))?;
    let y: f32 = parts[2]
        .parse()
        .map_err(|_| IoError::parse(line_number, "Invalid y coordinate"))?;
    let z: f32 = parts[3]
        .parse()
        .map_err(|_| IoError::parse(line_number, "Invalid z coordinate"))?;

    let element = Element::from_symbol(symbol).unwrap_or(Element::Unknown);
    let atom = Atom::new(symbol, element);

    Ok((atom, Vec3::new(x, y, z)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_atom_line() {
        let line = "O  0.0000  0.0000  0.0000";
        let (atom, coord) = parse_atom_line(line, 1).unwrap();

        assert_eq!(atom.element, Element::Oxygen);
        assert!((coord.x).abs() < 0.0001);
        assert!((coord.y).abs() < 0.0001);
        assert!((coord.z).abs() < 0.0001);
    }

    #[test]
    fn test_read_simple_xyz() {
        let xyz_data = r#"3
Water molecule
O     0.0000    0.0000    0.0000
H     0.9572    0.0000    0.0000
H    -0.2400    0.9266    0.0000
"#;

        let mut reader = XyzReader::new(xyz_data.as_bytes());
        let mol = reader.read().unwrap();

        assert_eq!(mol.atom_count(), 3);
        assert_eq!(mol.state_count(), 1);
        assert_eq!(mol.title, "Water molecule");
    }

    #[test]
    fn test_read_trajectory_xyz() {
        let xyz_data = r#"2
Frame 1
C  0.0  0.0  0.0
C  1.5  0.0  0.0
2
Frame 2
C  0.1  0.0  0.0
C  1.6  0.0  0.0
"#;

        let mut reader = XyzReader::new(xyz_data.as_bytes());
        let mol = reader.read().unwrap();

        assert_eq!(mol.atom_count(), 2);
        assert_eq!(mol.state_count(), 2);
    }
}
