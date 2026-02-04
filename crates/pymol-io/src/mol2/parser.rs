//! MOL2 file parser
//!
//! Parses TRIPOS MOL2 format files.

use std::collections::HashMap;
use std::io::{BufRead, BufReader, Read};
use std::sync::Arc;

use lin_alg::f32::Vec3;
use pymol_mol::{Atom, AtomIndex, AtomResidue, BondOrder, CoordSet, Element, ObjectMolecule};

use crate::error::{IoError, IoResult};
use crate::traits::MoleculeReader;

/// MOL2 file reader
pub struct Mol2Reader<R> {
    reader: BufReader<R>,
    line_number: usize,
    /// Buffer for peeking at lines
    peeked_line: Option<String>,
}

impl<R: Read> Mol2Reader<R> {
    /// Create a new MOL2 reader
    pub fn new(reader: R) -> Self {
        Mol2Reader {
            reader: BufReader::new(reader),
            line_number: 0,
            peeked_line: None,
        }
    }

    /// Read a single line from the file
    fn read_line(&mut self) -> IoResult<Option<String>> {
        // Check for peeked line first
        if let Some(line) = self.peeked_line.take() {
            return Ok(Some(line));
        }

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

    /// Peek at the next line without consuming it
    fn peek_line(&mut self) -> IoResult<Option<&str>> {
        if self.peeked_line.is_none() {
            self.peeked_line = self.read_line()?;
        }
        Ok(self.peeked_line.as_deref())
    }

    /// Skip to a specific section
    fn skip_to_section(&mut self, section: &str) -> IoResult<bool> {
        loop {
            match self.read_line()? {
                Some(line) => {
                    if line.starts_with("@<TRIPOS>") {
                        if line.contains(section) {
                            return Ok(true);
                        }
                    }
                }
                None => return Ok(false),
            }
        }
    }

    /// Read until end of section (next @<TRIPOS> or EOF)
    fn read_section_lines(&mut self) -> IoResult<Vec<String>> {
        let mut lines = Vec::new();
        loop {
            match self.peek_line()? {
                Some(line) => {
                    if line.starts_with("@<TRIPOS>") {
                        break;
                    }
                    if let Some(line) = self.read_line()? {
                        if !line.trim().is_empty() {
                            lines.push(line);
                        }
                    }
                }
                None => break,
            }
        }
        Ok(lines)
    }

    /// Parse a single molecule from the file
    fn parse_molecule(&mut self) -> IoResult<Option<ObjectMolecule>> {
        // Find @<TRIPOS>MOLECULE section
        if !self.skip_to_section("MOLECULE")? {
            return Ok(None);
        }

        // Line 1: Molecule name
        let name = match self.read_line()? {
            Some(line) => line.trim().to_string(),
            None => return Err(IoError::parse(self.line_number, "Expected molecule name")),
        };

        // Line 2: Counts (num_atoms num_bonds num_subst num_feat num_sets)
        let counts_line = match self.read_line()? {
            Some(line) => line,
            None => return Err(IoError::parse(self.line_number, "Expected counts line")),
        };

        let counts: Vec<&str> = counts_line.split_whitespace().collect();
        let n_atoms: usize = counts
            .first()
            .and_then(|s| s.parse().ok())
            .unwrap_or(0);
        let n_bonds: usize = counts
            .get(1)
            .and_then(|s| s.parse().ok())
            .unwrap_or(0);

        // Skip remaining molecule info lines (mol_type, charge_type, etc.)
        // until we hit another section
        loop {
            match self.peek_line()? {
                Some(line) => {
                    if line.starts_with("@<TRIPOS>") {
                        break;
                    }
                    self.read_line()?;
                }
                None => break,
            }
        }

        let mut mol = ObjectMolecule::with_capacity(&name, n_atoms, n_bonds);
        mol.name = name;

        // Parse @<TRIPOS>ATOM section
        if !self.skip_to_section("ATOM")? {
            return Err(IoError::parse(self.line_number, "Expected ATOM section"));
        }

        let atom_lines = self.read_section_lines()?;
        let mut coords = Vec::with_capacity(atom_lines.len());

        // Cache for sharing AtomResidue instances among atoms of the same residue
        let mut residue_cache: HashMap<AtomResidue, Arc<AtomResidue>> = HashMap::new();

        for line in &atom_lines {
            let (mut atom, coord) = parse_atom_line(line, self.line_number)?;

            // Use cache to share AtomResidue instances
            let residue_data = (*atom.residue).clone();
            atom.residue = residue_cache
                .entry(residue_data.clone())
                .or_insert_with(|| Arc::new(residue_data))
                .clone();

            mol.add_atom(atom);
            coords.push(coord);
        }

        // Add coordinate set
        let coord_set = CoordSet::from_vec3(&coords);
        mol.add_coord_set(coord_set);

        // Parse @<TRIPOS>BOND section
        if self.skip_to_section("BOND")? {
            let bond_lines = self.read_section_lines()?;

            for line in &bond_lines {
                if let Some((atom1, atom2, order)) = parse_bond_line(line) {
                    // Convert from 1-indexed to 0-indexed
                    let idx1 = AtomIndex((atom1 - 1) as u32);
                    let idx2 = AtomIndex((atom2 - 1) as u32);
                    let _ = mol.add_bond_unchecked(idx1, idx2, order);
                }
            }
        }

        // Skip other sections (SUBSTRUCTURE, etc.)
        // They will be handled by skip_to_section in next call

        Ok(Some(mol))
    }
}

impl<R: Read> MoleculeReader for Mol2Reader<R> {
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

/// Parse a MOL2 atom line
fn parse_atom_line(line: &str, line_number: usize) -> IoResult<(Atom, Vec3)> {
    // Format: atom_id atom_name x y z atom_type [subst_id [subst_name [charge [status_bit]]]]
    let parts: Vec<&str> = line.split_whitespace().collect();

    if parts.len() < 6 {
        return Err(IoError::parse(line_number, "Atom line too short"));
    }

    let atom_name = parts[1].to_string();
    let x: f32 = parts[2]
        .parse()
        .map_err(|_| IoError::parse(line_number, "Invalid x coordinate"))?;
    let y: f32 = parts[3]
        .parse()
        .map_err(|_| IoError::parse(line_number, "Invalid y coordinate"))?;
    let z: f32 = parts[4]
        .parse()
        .map_err(|_| IoError::parse(line_number, "Invalid z coordinate"))?;
    let atom_type = parts[5];

    // Parse element from SYBYL atom type (e.g., "C.3" -> "C", "N.ar" -> "N")
    let element_symbol = atom_type.split('.').next().unwrap_or(atom_type);
    let element = Element::from_symbol(element_symbol).unwrap_or(Element::Unknown);

    let mut atom = Atom::new(&atom_name, element);

    // Parse substructure info if present
    let mut resv = 1;
    let mut resn = String::new();

    if parts.len() > 6 {
        // subst_id
        if let Ok(parsed_resv) = parts[6].parse::<i32>() {
            resv = parsed_resv;
        }
    }

    if parts.len() > 7 {
        // subst_name (residue name)
        let resn_raw = parts[7];
        // Handle residue names like "ALA1" -> "ALA"
        let resn_clean: String = resn_raw.chars().take_while(|c| c.is_alphabetic()).collect();
        resn = if resn_clean.is_empty() {
            resn_raw.to_string()
        } else {
            resn_clean
        };
    }

    // Create AtomResidue
    atom.residue = Arc::new(AtomResidue::from_parts("", resn, resv, ' ', ""));

    if parts.len() > 8 {
        // partial charge
        if let Ok(charge) = parts[8].parse::<f32>() {
            atom.partial_charge = charge;
        }
    }

    Ok((atom, Vec3::new(x, y, z)))
}

/// Parse a MOL2 bond line
fn parse_bond_line(line: &str) -> Option<(usize, usize, BondOrder)> {
    // Format: bond_id origin_atom_id target_atom_id bond_type [status_bits]
    let parts: Vec<&str> = line.split_whitespace().collect();

    if parts.len() < 4 {
        return None;
    }

    let atom1: usize = parts[1].parse().ok()?;
    let atom2: usize = parts[2].parse().ok()?;
    let bond_type = parts[3];

    let order = match bond_type {
        "1" | "single" | "un" => BondOrder::Single,
        "2" | "double" => BondOrder::Double,
        "3" | "triple" => BondOrder::Triple,
        "ar" | "aromatic" => BondOrder::Aromatic,
        "am" => BondOrder::Single, // Amide bond treated as single
        "nc" => BondOrder::Unknown, // Not connected
        _ => BondOrder::Single,
    };

    Some((atom1, atom2, order))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_atom_line() {
        let line = "      1 N           1.0000    2.0000    3.0000 N.3       1 ALA1       -0.3000";
        let (atom, coord) = parse_atom_line(line, 1).unwrap();

        assert_eq!(&*atom.name, "N");
        assert_eq!(atom.element, Element::Nitrogen);
        assert!((coord.x - 1.0).abs() < 0.0001);
        assert!((coord.y - 2.0).abs() < 0.0001);
        assert!((coord.z - 3.0).abs() < 0.0001);
        assert_eq!(atom.residue.resn, "ALA");
        assert!((atom.partial_charge - (-0.3)).abs() < 0.0001);
    }

    #[test]
    fn test_parse_bond_line() {
        let line = "     1     1     2 1";
        let (a1, a2, order) = parse_bond_line(line).unwrap();
        assert_eq!(a1, 1);
        assert_eq!(a2, 2);
        assert_eq!(order, BondOrder::Single);

        let line2 = "     2     1     3 ar";
        let (_, _, order2) = parse_bond_line(line2).unwrap();
        assert_eq!(order2, BondOrder::Aromatic);
    }

    #[test]
    fn test_read_simple_mol2() {
        let mol2_data = r#"@<TRIPOS>MOLECULE
Water
 3 2 0 0 0
SMALL
NO_CHARGES


@<TRIPOS>ATOM
      1 O           0.0000    0.0000    0.0000 O.3       1 WAT1        0.0000
      2 H1          0.9572    0.0000    0.0000 H         1 WAT1        0.0000
      3 H2         -0.2400    0.9266    0.0000 H         1 WAT1        0.0000
@<TRIPOS>BOND
     1     1     2 1
     2     1     3 1
"#;

        let mut reader = Mol2Reader::new(mol2_data.as_bytes());
        let mol = reader.read().unwrap();

        assert_eq!(mol.name, "Water");
        assert_eq!(mol.atom_count(), 3);
        assert_eq!(mol.bond_count(), 2);
    }
}
