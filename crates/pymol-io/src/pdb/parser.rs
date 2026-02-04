//! PDB file parser
//!
//! Parses PDB format files using nom combinators.

use std::collections::HashMap;
use std::io::{BufRead, BufReader, Read};
use std::sync::Arc;

use lin_alg::f32::Vec3;
use nom::IResult;
use pymol_mol::{
    Atom, AtomIndex, AtomResidue, BondOrder, CoordSet, ObjectMolecule, SecondaryStructure,
};

use crate::error::{IoError, IoResult};
use crate::traits::MoleculeReader;

use super::records::{AtomRecord, ConectRecord, Cryst1Record, HelixRecord, SheetRecord};

/// PDB file reader
pub struct PdbReader<R> {
    reader: BufReader<R>,
    line_number: usize,
}

impl<R: Read> PdbReader<R> {
    /// Create a new PDB reader
    pub fn new(reader: R) -> Self {
        PdbReader {
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
                Ok(Some(line))
            }
            Err(e) => Err(IoError::Io(e)),
        }
    }

    /// Parse the PDB file
    fn parse(&mut self) -> IoResult<ObjectMolecule> {
        let mut atoms: Vec<AtomRecord> = Vec::new();
        let mut coords: Vec<Vec3> = Vec::new();
        let mut conects: Vec<ConectRecord> = Vec::new();
        let mut cryst1: Option<Cryst1Record> = None;
        let mut helices: Vec<HelixRecord> = Vec::new();
        let mut sheets: Vec<SheetRecord> = Vec::new();
        let mut title = String::new();
        let mut current_model = 0;
        let mut models: Vec<(Vec<AtomRecord>, Vec<Vec3>)> = Vec::new();
        let mut in_model = false;

        while let Some(line) = self.read_line()? {
            let line = line.trim_end();

            // Skip empty lines
            if line.is_empty() {
                continue;
            }

            // Get record type (first 6 characters)
            let record_type = if line.len() >= 6 {
                &line[0..6]
            } else {
                line
            };

            match record_type {
                "ATOM  " | "HETATM" => {
                    if let Ok((_, record)) = parse_atom_record(line) {
                        let coord = Vec3::new(record.x, record.y, record.z);
                        if in_model && current_model > 0 {
                            // Multi-model: store coordinates for this model
                            if let Some((_, model_coords)) = models.last_mut() {
                                model_coords.push(coord);
                            }
                        } else {
                            atoms.push(record);
                            coords.push(coord);
                        }
                    }
                }
                "CONECT" => {
                    if let Ok((_, record)) = parse_conect_record(line) {
                        conects.push(record);
                    }
                }
                "CRYST1" => {
                    if let Ok((_, record)) = parse_cryst1_record(line) {
                        cryst1 = Some(record);
                    }
                }
                "HELIX " => {
                    if let Ok((_, record)) = parse_helix_record(line) {
                        helices.push(record);
                    }
                }
                "SHEET " => {
                    if let Ok((_, record)) = parse_sheet_record(line) {
                        sheets.push(record);
                    }
                }
                "TITLE " => {
                    if line.len() > 10 {
                        if !title.is_empty() {
                            title.push(' ');
                        }
                        title.push_str(line[10..].trim());
                    }
                }
                "MODEL " => {
                    in_model = true;
                    current_model += 1;
                    if current_model > 1 {
                        // Start collecting coordinates for new model
                        models.push((Vec::new(), Vec::new()));
                    }
                }
                "ENDMDL" => {
                    // Model ended, coordinates for next model will be collected separately
                }
                "TER   " | "TER" => {
                    // Chain termination - we handle this implicitly via chain changes
                }
                "END   " | "END" => {
                    // End of file
                    break;
                }
                _ => {
                    // Ignore other record types for now
                }
            }
        }

        // Build the molecule
        self.build_molecule(atoms, coords, models, conects, cryst1, helices, sheets, title)
    }

    /// Build ObjectMolecule from parsed data
    fn build_molecule(
        &self,
        atom_records: Vec<AtomRecord>,
        coords: Vec<Vec3>,
        models: Vec<(Vec<AtomRecord>, Vec<Vec3>)>,
        conects: Vec<ConectRecord>,
        cryst1: Option<Cryst1Record>,
        helices: Vec<HelixRecord>,
        sheets: Vec<SheetRecord>,
        title: String,
    ) -> IoResult<ObjectMolecule> {
        if atom_records.is_empty() {
            return Err(IoError::EmptyFile);
        }

        let mut mol = ObjectMolecule::with_capacity("", atom_records.len(), 0);
        mol.title = title;

        // Map from PDB serial number to atom index
        let mut serial_to_index: HashMap<i32, AtomIndex> = HashMap::new();

        // Add atoms
        for record in &atom_records {
            let element = record.get_element();
            let mut atom = Atom::new(&record.name, element);

            // Create AtomResidue
            atom.residue = Arc::new(AtomResidue::from_parts(
                record.chain.clone(),
                record.resn.clone(),
                record.resv,
                record.icode,
                record.segi.clone(),
            ));
            atom.alt = record.alt_loc;
            atom.b_factor = record.b_factor;
            atom.occupancy = record.occupancy;
            atom.state.hetatm = record.hetatm;
            atom.formal_charge = record.get_formal_charge();
            atom.id = record.serial;

            let idx = mol.add_atom(atom);
            serial_to_index.insert(record.serial, idx);
        }

        // Add primary coordinate set
        let coord_set = CoordSet::from_vec3(&coords);
        mol.add_coord_set(coord_set);

        // Add additional models as coordinate sets
        for (_, model_coords) in models {
            if model_coords.len() == atom_records.len() {
                let coord_set = CoordSet::from_vec3(&model_coords);
                mol.add_coord_set(coord_set);
            }
        }

        // Add bonds from CONECT records
        for conect in &conects {
            if let Some(&atom1_idx) = serial_to_index.get(&conect.atom) {
                for &bonded_serial in &conect.bonded {
                    if let Some(&atom2_idx) = serial_to_index.get(&bonded_serial) {
                        // Avoid duplicate bonds (CONECT records are symmetric)
                        if atom1_idx.0 < atom2_idx.0 {
                            let _ = mol.add_bond_unchecked(atom1_idx, atom2_idx, BondOrder::Single);
                        }
                    }
                }
            }
        }

        // Apply secondary structure annotations
        apply_secondary_structure(&mut mol, &helices, &sheets);

        // Apply crystallographic symmetry
        if let Some(cryst) = cryst1 {
            use pymol_mol::Symmetry;
            mol.symmetry = Some(Symmetry::new(
                cryst.space_group,
                [cryst.a, cryst.b, cryst.c],
                [cryst.alpha, cryst.beta, cryst.gamma],
            ));
        }

        // Classify atoms (protein, nucleic, solvent, etc.)
        // This is needed for cartoon representation to work
        mol.classify_atoms();

        // Generate covalent bonds from distances
        // CONECT records only contain special bonds (disulfide, etc.), not standard covalent bonds
        mol.generate_bonds(0.6);

        // Assign bond orders for known protein residues using atom name templates
        mol.assign_known_residue_bond_orders();

        Ok(mol)
    }
}

impl<R: Read> MoleculeReader for PdbReader<R> {
    fn read(&mut self) -> IoResult<ObjectMolecule> {
        self.parse()
    }
}

/// Apply secondary structure annotations to atoms
fn apply_secondary_structure(
    mol: &mut ObjectMolecule,
    helices: &[HelixRecord],
    sheets: &[SheetRecord],
) {
    // This is a simplified implementation
    // A full implementation would iterate through atoms and check if they fall within
    // the residue ranges specified in HELIX/SHEET records

    // Collect helix atom indices first
    let helix_indices: Vec<AtomIndex> = helices
        .iter()
        .flat_map(|helix| {
            mol.atoms_indexed()
                .filter(|(_, atom)| {
                    atom.residue.chain == helix.init_chain
                        && atom.residue.resv >= helix.init_seq
                        && atom.residue.resv <= helix.end_seq
                })
                .map(|(idx, _)| idx)
                .collect::<Vec<_>>()
        })
        .collect();

    // Apply helix secondary structure
    for idx in helix_indices {
        if let Some(atom_mut) = mol.get_atom_mut(idx) {
            atom_mut.ss_type = SecondaryStructure::Helix;
        }
    }

    // Collect sheet atom indices
    let sheet_indices: Vec<AtomIndex> = sheets
        .iter()
        .flat_map(|sheet| {
            mol.atoms_indexed()
                .filter(|(_, atom)| {
                    atom.residue.chain == sheet.init_chain
                        && atom.residue.resv >= sheet.init_seq
                        && atom.residue.resv <= sheet.end_seq
                })
                .map(|(idx, _)| idx)
                .collect::<Vec<_>>()
        })
        .collect();

    // Apply sheet secondary structure
    for idx in sheet_indices {
        if let Some(atom_mut) = mol.get_atom_mut(idx) {
            atom_mut.ss_type = SecondaryStructure::Sheet;
        }
    }
}

// ============================================================================
// Nom parsers for PDB records
// ============================================================================

/// Parse an ATOM or HETATM record
fn parse_atom_record(input: &str) -> IResult<&str, AtomRecord> {
    // Ensure line is at least 54 characters (minimum for coordinates)
    if input.len() < 54 {
        return Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::TooLarge,
        )));
    }

    let hetatm = input.starts_with("HETATM");

    // Fixed column positions (0-indexed):
    // 6-10: serial (5)
    // 11: space
    // 12-15: name (4)
    // 16: altLoc (1)
    // 17-19: resName (3)
    // 20: space
    // 21: chainID (1)
    // 22-25: resSeq (4)
    // 26: iCode (1)
    // 27-29: spaces
    // 30-37: x (8)
    // 38-45: y (8)
    // 46-53: z (8)
    // 54-59: occupancy (6)
    // 60-65: tempFactor (6)
    // 66-71: spaces
    // 72-75: segment (4)
    // 76-77: element (2)
    // 78-79: charge (2)

    let serial: i32 = input
        .get(6..11)
        .and_then(|s| s.trim().parse().ok())
        .unwrap_or(0);
    let name = input.get(12..16).unwrap_or("    ").trim().to_string();
    let alt_loc = input.chars().nth(16).unwrap_or(' ');
    let resn = input.get(17..20).unwrap_or("   ").trim().to_string();
    let chain = input.get(21..22).unwrap_or(" ").trim().to_string();
    let resv: i32 = input
        .get(22..26)
        .and_then(|s| s.trim().parse().ok())
        .unwrap_or(0);
    let icode = input.chars().nth(26).unwrap_or(' ');

    let x: f32 = input
        .get(30..38)
        .and_then(|s| s.trim().parse().ok())
        .unwrap_or(0.0);
    let y: f32 = input
        .get(38..46)
        .and_then(|s| s.trim().parse().ok())
        .unwrap_or(0.0);
    let z: f32 = input
        .get(46..54)
        .and_then(|s| s.trim().parse().ok())
        .unwrap_or(0.0);

    let occupancy: f32 = input
        .get(54..60)
        .and_then(|s| s.trim().parse().ok())
        .unwrap_or(1.0);
    let b_factor: f32 = input
        .get(60..66)
        .and_then(|s| s.trim().parse().ok())
        .unwrap_or(0.0);

    let segi = input.get(72..76).unwrap_or("    ").trim().to_string();
    let element = input.get(76..78).unwrap_or("  ").trim().to_string();
    let charge = input.get(78..80).unwrap_or("  ").trim().to_string();

    let record = AtomRecord {
        hetatm,
        serial,
        name,
        alt_loc,
        resn,
        chain,
        resv,
        icode,
        x,
        y,
        z,
        occupancy,
        b_factor,
        segi,
        element,
        charge,
    };

    Ok(("", record))
}

/// Parse a CONECT record
fn parse_conect_record(input: &str) -> IResult<&str, ConectRecord> {
    // CONECT record format:
    // 6-10: atom serial number
    // 11-15, 16-20, 21-25, 26-30: bonded atom serial numbers

    let atom: i32 = input
        .get(6..11)
        .and_then(|s| s.trim().parse().ok())
        .unwrap_or(0);

    let mut bonded = Vec::new();

    for start in [11, 16, 21, 26] {
        if let Some(field) = input.get(start..start + 5) {
            if let Ok(serial) = field.trim().parse::<i32>() {
                if serial > 0 {
                    bonded.push(serial);
                }
            }
        }
    }

    Ok(("", ConectRecord { atom, bonded }))
}

/// Parse a CRYST1 record
fn parse_cryst1_record(input: &str) -> IResult<&str, Cryst1Record> {
    // CRYST1 record format:
    // 6-14: a (9)
    // 15-23: b (9)
    // 24-32: c (9)
    // 33-39: alpha (7)
    // 40-46: beta (7)
    // 47-53: gamma (7)
    // 55-65: space group (11)
    // 66-69: z (4)

    let a: f32 = input
        .get(6..15)
        .and_then(|s| s.trim().parse().ok())
        .unwrap_or(1.0);
    let b: f32 = input
        .get(15..24)
        .and_then(|s| s.trim().parse().ok())
        .unwrap_or(1.0);
    let c: f32 = input
        .get(24..33)
        .and_then(|s| s.trim().parse().ok())
        .unwrap_or(1.0);
    let alpha: f32 = input
        .get(33..40)
        .and_then(|s| s.trim().parse().ok())
        .unwrap_or(90.0);
    let beta: f32 = input
        .get(40..47)
        .and_then(|s| s.trim().parse().ok())
        .unwrap_or(90.0);
    let gamma: f32 = input
        .get(47..54)
        .and_then(|s| s.trim().parse().ok())
        .unwrap_or(90.0);
    let space_group = input.get(55..66).unwrap_or("P 1").trim().to_string();
    let z: i32 = input
        .get(66..70)
        .and_then(|s| s.trim().parse().ok())
        .unwrap_or(1);

    Ok((
        "",
        Cryst1Record {
            a,
            b,
            c,
            alpha,
            beta,
            gamma,
            space_group,
            z,
        },
    ))
}

/// Parse a HELIX record
fn parse_helix_record(input: &str) -> IResult<&str, HelixRecord> {
    // HELIX record format (simplified)
    let serial: i32 = input
        .get(7..10)
        .and_then(|s| s.trim().parse().ok())
        .unwrap_or(0);
    let helix_id = input.get(11..14).unwrap_or("").trim().to_string();
    let init_resn = input.get(15..18).unwrap_or("").trim().to_string();
    let init_chain = input.get(19..20).unwrap_or("").trim().to_string();
    let init_seq: i32 = input
        .get(21..25)
        .and_then(|s| s.trim().parse().ok())
        .unwrap_or(0);
    let init_icode = input.chars().nth(25).unwrap_or(' ');
    let end_resn = input.get(27..30).unwrap_or("").trim().to_string();
    let end_chain = input.get(31..32).unwrap_or("").trim().to_string();
    let end_seq: i32 = input
        .get(33..37)
        .and_then(|s| s.trim().parse().ok())
        .unwrap_or(0);
    let end_icode = input.chars().nth(37).unwrap_or(' ');
    let helix_class: i32 = input
        .get(38..40)
        .and_then(|s| s.trim().parse().ok())
        .unwrap_or(1);

    Ok((
        "",
        HelixRecord {
            serial,
            helix_id,
            init_resn,
            init_chain,
            init_seq,
            init_icode,
            end_resn,
            end_chain,
            end_seq,
            end_icode,
            helix_class,
        },
    ))
}

/// Parse a SHEET record
fn parse_sheet_record(input: &str) -> IResult<&str, SheetRecord> {
    // SHEET record format (simplified)
    let strand: i32 = input
        .get(7..10)
        .and_then(|s| s.trim().parse().ok())
        .unwrap_or(0);
    let sheet_id = input.get(11..14).unwrap_or("").trim().to_string();
    let num_strands: i32 = input
        .get(14..16)
        .and_then(|s| s.trim().parse().ok())
        .unwrap_or(0);
    let init_resn = input.get(17..20).unwrap_or("").trim().to_string();
    let init_chain = input.get(21..22).unwrap_or("").trim().to_string();
    let init_seq: i32 = input
        .get(22..26)
        .and_then(|s| s.trim().parse().ok())
        .unwrap_or(0);
    let init_icode = input.chars().nth(26).unwrap_or(' ');
    let end_resn = input.get(28..31).unwrap_or("").trim().to_string();
    let end_chain = input.get(32..33).unwrap_or("").trim().to_string();
    let end_seq: i32 = input
        .get(33..37)
        .and_then(|s| s.trim().parse().ok())
        .unwrap_or(0);
    let end_icode = input.chars().nth(37).unwrap_or(' ');
    let sense: i32 = input
        .get(38..40)
        .and_then(|s| s.trim().parse().ok())
        .unwrap_or(0);

    Ok((
        "",
        SheetRecord {
            strand,
            sheet_id,
            num_strands,
            init_resn,
            init_chain,
            init_seq,
            init_icode,
            end_resn,
            end_chain,
            end_seq,
            end_icode,
            sense,
        },
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_atom_record() {
        let line =
            "ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00 20.00           N  ";
        let (_, record) = parse_atom_record(line).unwrap();

        assert_eq!(record.serial, 1);
        assert_eq!(record.name.trim(), "N");
        assert_eq!(record.resn, "ALA");
        assert_eq!(record.chain, "A");
        assert_eq!(record.resv, 1);
        assert!((record.x - 1.0).abs() < 0.001);
        assert!((record.y - 2.0).abs() < 0.001);
        assert!((record.z - 3.0).abs() < 0.001);
        assert!((record.b_factor - 20.0).abs() < 0.001);
        assert_eq!(record.element, "N");
    }

    #[test]
    fn test_parse_hetatm_record() {
        let line =
            "HETATM    1  O   HOH A   1       1.000   2.000   3.000  1.00 20.00           O  ";
        let (_, record) = parse_atom_record(line).unwrap();

        assert!(record.hetatm);
        assert_eq!(record.resn, "HOH");
    }

    #[test]
    fn test_parse_conect_record() {
        let line = "CONECT    1    2    3    4    5";
        let (_, record) = parse_conect_record(line).unwrap();

        assert_eq!(record.atom, 1);
        assert_eq!(record.bonded, vec![2, 3, 4, 5]);
    }

    #[test]
    fn test_read_simple_pdb() {
        let pdb = r#"ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00 20.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00 20.00           C
ATOM      4  O   ALA A   1       1.251   2.390   0.000  1.00 20.00           O
END
"#;

        let mut reader = PdbReader::new(pdb.as_bytes());
        let mol = reader.read().unwrap();

        assert_eq!(mol.atom_count(), 4);
        assert_eq!(mol.state_count(), 1);
    }
}
