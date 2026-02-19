//! GROMACS GRO file parser
//!
//! Parses fixed-width GRO coordinate format files.
//! Coordinates are converted from nanometers to Angstroms (×10).

use std::collections::HashMap;
use std::io::{BufRead, BufReader, Read};
use std::sync::Arc;

use lin_alg::f32::Vec3;
use pymol_mol::{Atom, AtomResidue, CoordSet, ObjectMolecule, Symmetry};

use crate::error::{IoError, IoResult};
use crate::pdb::infer_element_from_name;
use crate::traits::MoleculeReader;

/// Nanometers to Angstroms conversion factor
const NM_TO_ANGSTROM: f32 = 10.0;

/// GRO file reader
pub struct GroReader<R> {
    reader: BufReader<R>,
    line_number: usize,
}

impl<R: Read> GroReader<R> {
    /// Create a new GRO reader
    pub fn new(reader: R) -> Self {
        GroReader {
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

    /// Parse a single frame from the file
    fn parse_frame(&mut self) -> IoResult<Option<ObjectMolecule>> {
        // Line 1: Title
        let title = match self.read_line()? {
            Some(line) => line.trim().to_string(),
            None => return Ok(None),
        };

        // Line 2: Number of atoms
        let n_atoms_line = match self.read_line()? {
            Some(line) => line,
            None => {
                return Err(IoError::parse(
                    self.line_number,
                    "Expected atom count after title",
                ))
            }
        };

        let n_atoms: usize = n_atoms_line
            .trim()
            .parse()
            .map_err(|_| IoError::parse(self.line_number, "Invalid atom count"))?;

        if n_atoms == 0 {
            return Err(IoError::parse(self.line_number, "Zero atoms in GRO file"));
        }

        let mut mol = ObjectMolecule::with_capacity("", n_atoms, 0);
        mol.title = title.clone();
        mol.name = if title.is_empty() {
            "molecule".to_string()
        } else {
            title
                .split_whitespace()
                .next()
                .unwrap_or("molecule")
                .to_string()
        };

        let mut coords = Vec::with_capacity(n_atoms);
        let mut residue_cache: HashMap<AtomResidue, Arc<AtomResidue>> = HashMap::new();

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

            let (atom, coord) = parse_atom_line(&line, self.line_number, &mut residue_cache)?;
            mol.add_atom(atom);
            coords.push(coord);
        }

        // Box vectors line
        let box_line = match self.read_line()? {
            Some(line) => line,
            None => {
                return Err(IoError::parse(
                    self.line_number,
                    "Expected box vectors line",
                ))
            }
        };

        if let Some(symmetry) = parse_box_vectors(&box_line) {
            mol.symmetry = Some(symmetry);
        }

        // Add coordinate set
        let coord_set = CoordSet::from_vec3(&coords);
        mol.add_coord_set(coord_set);

        // Post-processing: classify atoms and generate bonds (GRO has no bond info)
        mol.classify_atoms();
        mol.generate_bonds(0.6);
        mol.assign_known_residue_bond_orders();
        mol.assign_chains();

        Ok(Some(mol))
    }
}

impl<R: Read> MoleculeReader for GroReader<R> {
    fn read(&mut self) -> IoResult<ObjectMolecule> {
        let mut mol = match self.parse_frame()? {
            Some(mol) => mol,
            None => return Err(IoError::EmptyFile),
        };

        // Read additional frames as coordinate sets (trajectory)
        while let Some(frame) = self.parse_frame()? {
            if frame.atom_count() != mol.atom_count() {
                break;
            }
            if let Some(cs) = frame.get_coord_set(0) {
                mol.add_coord_set(cs.clone());
            }
        }

        Ok(mol)
    }

    fn read_all(&mut self) -> IoResult<Vec<ObjectMolecule>> {
        let mol = self.read()?;
        Ok(vec![mol])
    }
}

/// Parse a GRO atom line (fixed-width columns)
///
/// Format: `%5d%-5s%5s%5d%8.3f%8.3f%8.3f[%8.4f%8.4f%8.4f]`
///   - columns  0.. 5: residue number (i32)
///   - columns  5..10: residue name
///   - columns 10..15: atom name
///   - columns 15..20: atom number (i32)
///   - columns 20..28: x (nm)
///   - columns 28..36: y (nm)
///   - columns 36..44: z (nm)
///   - columns 44..52: vx (optional, ignored)
///   - columns 52..60: vy (optional, ignored)
///   - columns 60..68: vz (optional, ignored)
fn parse_atom_line(
    line: &str,
    line_number: usize,
    residue_cache: &mut HashMap<AtomResidue, Arc<AtomResidue>>,
) -> IoResult<(Atom, Vec3)> {
    if line.len() < 44 {
        return Err(IoError::parse(
            line_number,
            format!(
                "GRO atom line too short ({}): expected ≥44 chars",
                line.len()
            ),
        ));
    }

    let resnum: i32 = line[0..5]
        .trim()
        .parse()
        .map_err(|_| IoError::parse(line_number, "Invalid residue number"))?;

    let resname = line[5..10].trim().to_string();
    let atom_name = line[10..15].trim().to_string();

    let atom_number: i32 = line[15..20]
        .trim()
        .parse()
        .map_err(|_| IoError::parse(line_number, "Invalid atom number"))?;

    let x: f32 = line[20..28]
        .trim()
        .parse::<f32>()
        .map_err(|_| IoError::parse(line_number, "Invalid x coordinate"))?
        * NM_TO_ANGSTROM;

    let y: f32 = line[28..36]
        .trim()
        .parse::<f32>()
        .map_err(|_| IoError::parse(line_number, "Invalid y coordinate"))?
        * NM_TO_ANGSTROM;

    let z: f32 = line[36..44]
        .trim()
        .parse::<f32>()
        .map_err(|_| IoError::parse(line_number, "Invalid z coordinate"))?
        * NM_TO_ANGSTROM;

    // Infer element from atom name (reuse PDB logic)
    let element = infer_element_from_name(&atom_name);

    let mut atom = Atom::new(&atom_name, element);
    atom.id = atom_number;

    // Set residue via cache (same pattern as PDB parser)
    let residue_data = AtomResidue::from_parts("", &resname, resnum, ' ', "");
    atom.residue = residue_cache
        .entry(residue_data.clone())
        .or_insert_with(|| Arc::new(residue_data))
        .clone();

    Ok((atom, Vec3::new(x, y, z)))
}

/// Parse box vectors line into Symmetry
///
/// GRO box line has 3 values (orthogonal) or 9 values (triclinic):
///   v1(x) v2(y) v3(z) [v1(y) v1(z) v2(x) v2(z) v3(x) v3(y)]
/// All values in nanometers.
fn parse_box_vectors(line: &str) -> Option<Symmetry> {
    let parts: Vec<f32> = line
        .split_whitespace()
        .filter_map(|s| s.parse::<f32>().ok())
        .collect();

    if parts.len() < 3 {
        return None;
    }

    // Convert box vectors from nm to Angstroms
    if parts.len() >= 9 {
        // Triclinic: full 3x3 matrix
        // v1 = (parts[0], parts[3], parts[4])
        // v2 = (parts[5], parts[1], parts[6])
        // v3 = (parts[7], parts[8], parts[2])
        let v1 = Vec3::new(parts[0], parts[3], parts[4]) * NM_TO_ANGSTROM;
        let v2 = Vec3::new(parts[5], parts[1], parts[6]) * NM_TO_ANGSTROM;
        let v3 = Vec3::new(parts[7], parts[8], parts[2]) * NM_TO_ANGSTROM;

        let a = v1.magnitude();
        let b = v2.magnitude();
        let c = v3.magnitude();

        let alpha = (v2.dot(v3) / (b * c)).clamp(-1.0, 1.0).acos().to_degrees();
        let beta = (v1.dot(v3) / (a * c)).clamp(-1.0, 1.0).acos().to_degrees();
        let gamma = (v1.dot(v2) / (a * b)).clamp(-1.0, 1.0).acos().to_degrees();

        Some(Symmetry::new("P 1", [a, b, c], [alpha, beta, gamma]))
    } else {
        // Orthogonal box: only diagonal elements
        let a = parts[0] * NM_TO_ANGSTROM;
        let b = parts[1] * NM_TO_ANGSTROM;
        let c = parts[2] * NM_TO_ANGSTROM;

        Some(Symmetry::new("P 1", [a, b, c], [90.0, 90.0, 90.0]))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use pymol_mol::{AtomIndex, Element};

    #[test]
    fn test_parse_single_frame() {
        let gro_data = "\
GROningen MAchine for Chemical Simulation
3
    1SOL     OW    1   0.230   0.628   0.113
    1SOL    HW1    2   0.137   0.628   0.150
    1SOL    HW2    3   0.231   0.589   0.021
   1.00000   1.00000   1.00000
";

        let mut reader = GroReader::new(gro_data.as_bytes());
        let mol = reader.read().unwrap();

        assert_eq!(mol.atom_count(), 3);
        assert_eq!(mol.state_count(), 1);

        // Check coordinates are in Angstroms (nm × 10)
        let coord = mol.get_coord(AtomIndex::from(0usize), 0).unwrap();
        assert!((coord.x - 2.30).abs() < 0.01);
        assert!((coord.y - 6.28).abs() < 0.01);
        assert!((coord.z - 1.13).abs() < 0.01);

        // Check residue info
        let atom = mol.get_atom(AtomIndex::from(0usize)).unwrap();
        assert_eq!(atom.residue.resn, "SOL");
        assert_eq!(atom.residue.resv, 1);
    }

    #[test]
    fn test_parse_with_velocities() {
        // Velocities in columns 44..68 should be ignored
        let gro_data = "\
Water with velocities
3
    1SOL     OW    1   0.230   0.628   0.113  0.1234  0.5678  0.9012
    1SOL    HW1    2   0.137   0.628   0.150 -0.1234 -0.5678  0.0000
    1SOL    HW2    3   0.231   0.589   0.021  0.0000  0.0000  0.0000
   1.00000   1.00000   1.00000
";

        let mut reader = GroReader::new(gro_data.as_bytes());
        let mol = reader.read().unwrap();

        assert_eq!(mol.atom_count(), 3);

        // Coordinates should be correct (velocities ignored)
        let coord = mol.get_coord(AtomIndex::from(0usize), 0).unwrap();
        assert!((coord.x - 2.30).abs() < 0.01);
    }

    #[test]
    fn test_parse_trajectory() {
        let gro_data = "\
Frame 1
2
    1ALA     CA    1   0.100   0.200   0.300
    1ALA      N    2   0.150   0.250   0.350
   1.00000   1.00000   1.00000
Frame 2
2
    1ALA     CA    1   0.110   0.210   0.310
    1ALA      N    2   0.160   0.260   0.360
   1.00000   1.00000   1.00000
";

        let mut reader = GroReader::new(gro_data.as_bytes());
        let mol = reader.read().unwrap();

        assert_eq!(mol.atom_count(), 2);
        assert_eq!(mol.state_count(), 2);

        // Check frame 2 coordinates
        let coord = mol.get_coord(AtomIndex::from(0usize), 1).unwrap();
        assert!((coord.x - 1.10).abs() < 0.01);
    }

    #[test]
    fn test_box_vectors_orthogonal() {
        let symmetry = parse_box_vectors("   1.00000   2.00000   3.00000").unwrap();
        assert!((symmetry.cell_lengths[0] - 10.0).abs() < 0.01);
        assert!((symmetry.cell_lengths[1] - 20.0).abs() < 0.01);
        assert!((symmetry.cell_lengths[2] - 30.0).abs() < 0.01);
        assert!((symmetry.cell_angles[0] - 90.0).abs() < 0.01);
        assert!((symmetry.cell_angles[1] - 90.0).abs() < 0.01);
        assert!((symmetry.cell_angles[2] - 90.0).abs() < 0.01);
    }

    #[test]
    fn test_box_vectors_triclinic() {
        // Orthogonal box expressed in 9-value format (off-diag = 0)
        let symmetry =
            parse_box_vectors("   1.00000   2.00000   3.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000")
                .unwrap();
        assert!((symmetry.cell_lengths[0] - 10.0).abs() < 0.01);
        assert!((symmetry.cell_lengths[1] - 20.0).abs() < 0.01);
        assert!((symmetry.cell_lengths[2] - 30.0).abs() < 0.01);
    }

    #[test]
    fn test_element_inference() {
        let gro_data = "\
Element test
5
    1ALA     CA    1   0.100   0.200   0.300
    1ALA      N    2   0.150   0.250   0.350
    1ALA      O    3   0.200   0.300   0.400
    1ALA      C    4   0.250   0.350   0.450
    1ALA      H    5   0.300   0.400   0.500
   1.00000   1.00000   1.00000
";

        let mut reader = GroReader::new(gro_data.as_bytes());
        let mol = reader.read().unwrap();

        let elements: Vec<Element> = mol.atoms().map(|a| a.element).collect();
        assert_eq!(elements[0], Element::Carbon); // CA
        assert_eq!(elements[1], Element::Nitrogen); // N
        assert_eq!(elements[2], Element::Oxygen); // O
        assert_eq!(elements[3], Element::Carbon); // C
        assert_eq!(elements[4], Element::Hydrogen); // H
    }

    #[test]
    fn test_format_detection() {
        assert_eq!(
            crate::FileFormat::from_extension("gro"),
            crate::FileFormat::Gro
        );
        assert_eq!(crate::FileFormat::Gro.extension(), "gro");
        assert_eq!(crate::FileFormat::Gro.name(), "GRO");
    }

    #[test]
    fn test_chain_assignment_protein_ions_solvent() {
        // A small system: 2 ALA residues + 1 NA ion + 2 SOL waters
        let gro_data = "\
Protein with ions and solvent
15
    1ALA      N    1   0.100   0.200   0.300
    1ALA     CA    2   0.247   0.200   0.300
    1ALA      C    3   0.380   0.200   0.300
    1ALA      O    4   0.380   0.320   0.300
    2ALA      N    5   0.513   0.200   0.300
    2ALA     CA    6   0.660   0.200   0.300
    2ALA      C    7   0.793   0.200   0.300
    2ALA      O    8   0.793   0.320   0.300
    3NA      NA    9   2.000   2.000   2.000
    4SOL     OW   10   3.000   3.000   3.000
    4SOL    HW1   11   3.050   3.050   3.000
    4SOL    HW2   12   3.050   2.950   3.000
    5SOL     OW   13   3.500   3.500   3.500
    5SOL    HW1   14   3.550   3.550   3.500
    5SOL    HW2   15   3.550   3.450   3.500
   5.00000   5.00000   5.00000
";

        let mut reader = GroReader::new(gro_data.as_bytes());
        let mol = reader.read().unwrap();

        assert_eq!(mol.atom_count(), 15);

        // Protein atoms should have chain "A"
        assert_eq!(mol.atoms_slice()[0].residue.chain, "A"); // ALA 1
        assert_eq!(mol.atoms_slice()[4].residue.chain, "A"); // ALA 2

        // Ion should have its own chain (not "A")
        let ion_chain = &mol.atoms_slice()[8].residue.chain;
        assert!(!ion_chain.is_empty());
        assert_ne!(ion_chain, "A");

        // Solvent should have its own shared chain
        let sol_chain = &mol.atoms_slice()[9].residue.chain;
        assert!(!sol_chain.is_empty());
        assert_eq!(mol.atoms_slice()[12].residue.chain, *sol_chain); // Both SOL same chain

        // All three groups have different chains
        assert_ne!(ion_chain, sol_chain);
    }
}
