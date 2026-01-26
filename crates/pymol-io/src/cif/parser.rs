//! mmCIF file parser
//!
//! Parses macromolecular CIF format files.

use std::collections::HashMap;
use std::io::{BufReader, Read};

use lin_alg::f32::Vec3;
use pymol_mol::{Atom, AtomIndex, CoordSet, Element, ObjectMolecule, SecondaryStructure};

use crate::error::{IoError, IoResult};
use crate::traits::MoleculeReader;

use super::lexer::{tokenize, Token};

/// Secondary structure annotation from mmCIF
#[derive(Debug, Clone)]
struct SecondaryStructureRange {
    /// Type: helix or sheet
    ss_type: SecondaryStructure,
    /// Chain ID (label_asym_id)
    chain: String,
    /// Auth chain ID (auth_asym_id) - used for matching atoms
    auth_chain: Option<String>,
    /// Start residue sequence number
    start_seq: i32,
    /// End residue sequence number
    end_seq: i32,
}

/// mmCIF file reader
pub struct CifReader<R> {
    reader: BufReader<R>,
}

impl<R: Read> CifReader<R> {
    /// Create a new CIF reader
    pub fn new(reader: R) -> Self {
        CifReader {
            reader: BufReader::new(reader),
        }
    }

    /// Read and parse the CIF file
    fn parse(&mut self) -> IoResult<ObjectMolecule> {
        // Read entire content
        let mut content = String::new();
        self.reader
            .read_to_string(&mut content)
            .map_err(IoError::Io)?;

        // Tokenize
        let tokens = tokenize(&content);

        // Parse tokens
        parse_cif_tokens(&tokens)
    }
}

impl<R: Read> MoleculeReader for CifReader<R> {
    fn read(&mut self) -> IoResult<ObjectMolecule> {
        self.parse()
    }
}

/// Parse CIF tokens into a molecule
fn parse_cif_tokens(tokens: &[Token]) -> IoResult<ObjectMolecule> {
    let mut mol = ObjectMolecule::new("");
    let mut pos = 0;
    let mut ss_ranges: Vec<SecondaryStructureRange> = Vec::new();

    // Find data block name
    while pos < tokens.len() {
        match &tokens[pos] {
            Token::DataBlock(name) => {
                mol.name = name.clone();
                pos += 1;
                break;
            }
            Token::Eof => return Err(IoError::EmptyFile),
            _ => pos += 1,
        }
    }

    // Parse data items and loops
    while pos < tokens.len() {
        match &tokens[pos] {
            Token::Loop => {
                pos += 1;
                pos = parse_loop(tokens, pos, &mut mol, &mut ss_ranges)?;
            }
            Token::DataName(name) if name.starts_with("_cell.") => {
                pos = parse_cell(tokens, pos, &mut mol)?;
            }
            Token::DataName(name) if name.starts_with("_entry.") => {
                pos += 1;
                if let Some(Token::Value(val) | Token::SingleQuoted(val) | Token::DoubleQuoted(val)) =
                    tokens.get(pos)
                {
                    if mol.title.is_empty() {
                        mol.title = val.clone();
                    }
                    pos += 1;
                }
            }
            Token::DataBlock(_) => {
                // New data block - stop parsing this one
                break;
            }
            Token::Eof => break,
            _ => pos += 1,
        }
    }

    if mol.atom_count() == 0 {
        return Err(IoError::EmptyFile);
    }

    // Apply secondary structure annotations
    apply_secondary_structure(&mut mol, &ss_ranges);

    // Classify atoms (protein, nucleic, solvent, etc.)
    // This is needed for cartoon representation to work
    mol.classify_atoms();

    // Generate covalent bonds from distances
    // CIF files typically don't include explicit bond information
    mol.generate_bonds(0.6);

    Ok(mol)
}

/// Parse a loop_ construct
fn parse_loop(
    tokens: &[Token],
    mut pos: usize,
    mol: &mut ObjectMolecule,
    ss_ranges: &mut Vec<SecondaryStructureRange>,
) -> IoResult<usize> {
    // Collect column names
    let mut columns: Vec<String> = Vec::new();
    while pos < tokens.len() {
        match &tokens[pos] {
            Token::DataName(name) => {
                columns.push(name.clone());
                pos += 1;
            }
            _ => break,
        }
    }

    if columns.is_empty() {
        return Ok(pos);
    }

    // Check loop category
    let is_atom_site = columns.iter().any(|c| c.starts_with("_atom_site."));
    let is_struct_conf = columns.iter().any(|c| c.starts_with("_struct_conf."));
    let is_struct_sheet = columns.iter().any(|c| c.starts_with("_struct_sheet_range."));

    if is_atom_site {
        pos = parse_atom_site_loop(tokens, pos, &columns, mol)?;
    } else if is_struct_conf {
        pos = parse_struct_conf_loop(tokens, pos, &columns, ss_ranges)?;
    } else if is_struct_sheet {
        pos = parse_struct_sheet_range_loop(tokens, pos, &columns, ss_ranges)?;
    } else {
        // Skip other loops
        while pos < tokens.len() {
            match &tokens[pos] {
                Token::Loop | Token::DataName(_) | Token::DataBlock(_) | Token::Eof => break,
                _ => pos += 1,
            }
        }
    }

    Ok(pos)
}

/// Parse _atom_site loop
fn parse_atom_site_loop(
    tokens: &[Token],
    mut pos: usize,
    columns: &[String],
    mol: &mut ObjectMolecule,
) -> IoResult<usize> {
    // Map column names to indices
    let col_idx: HashMap<&str, usize> = columns
        .iter()
        .enumerate()
        .map(|(i, s)| (s.as_str(), i))
        .collect();

    let n_cols = columns.len();
    let mut coords: Vec<Vec3> = Vec::new();
    let mut models: HashMap<i32, Vec<Vec3>> = HashMap::new();

    // Parse rows
    loop {
        // Check if we've reached end of loop data
        if pos >= tokens.len() {
            break;
        }
        match &tokens[pos] {
            Token::Loop | Token::DataName(_) | Token::DataBlock(_) | Token::Eof => break,
            _ => {}
        }

        // Read one row
        let mut row: Vec<Option<String>> = Vec::with_capacity(n_cols);
        for _ in 0..n_cols {
            if pos >= tokens.len() {
                break;
            }
            match &tokens[pos] {
                Token::Value(v) | Token::SingleQuoted(v) | Token::DoubleQuoted(v) | Token::TextField(v) => {
                    row.push(Some(v.clone()));
                }
                Token::Missing | Token::Unknown => {
                    row.push(None);
                }
                Token::Loop | Token::DataName(_) | Token::DataBlock(_) | Token::Eof => {
                    break;
                }
            }
            pos += 1;
        }

        if row.len() < n_cols {
            // Incomplete row - end of loop
            break;
        }

        // Parse atom from row
        let get_col = |name: &str| -> Option<&str> {
            col_idx.get(name).and_then(|&i| row.get(i).and_then(|v| v.as_deref()))
        };

        // Get element (prefer type_symbol over label_atom_id)
        let element_str = get_col("_atom_site.type_symbol")
            .or_else(|| get_col("_atom_site.label_atom_id"))
            .unwrap_or("X");
        let element = Element::from_symbol(element_str).unwrap_or(Element::Unknown);

        // Get atom name (prefer auth_atom_id over label_atom_id)
        let atom_name = get_col("_atom_site.auth_atom_id")
            .or_else(|| get_col("_atom_site.label_atom_id"))
            .unwrap_or("X");

        let mut atom = Atom::new(atom_name, element);

        // Residue info (prefer auth over label)
        atom.resn = get_col("_atom_site.auth_comp_id")
            .or_else(|| get_col("_atom_site.label_comp_id"))
            .unwrap_or("")
            .to_string();

        atom.resv = get_col("_atom_site.auth_seq_id")
            .or_else(|| get_col("_atom_site.label_seq_id"))
            .and_then(|s| s.parse().ok())
            .unwrap_or(0);

        atom.chain = get_col("_atom_site.auth_asym_id")
            .or_else(|| get_col("_atom_site.label_asym_id"))
            .unwrap_or("")
            .to_string();

        // Insertion code
        atom.inscode = get_col("_atom_site.pdbx_PDB_ins_code")
            .and_then(|s| s.chars().next())
            .unwrap_or(' ');

        // Alt loc
        atom.alt = get_col("_atom_site.label_alt_id")
            .and_then(|s| s.chars().next())
            .unwrap_or(' ');

        // B-factor
        atom.b_factor = get_col("_atom_site.B_iso_or_equiv")
            .and_then(|s| s.parse().ok())
            .unwrap_or(0.0);

        // Occupancy
        atom.occupancy = get_col("_atom_site.occupancy")
            .and_then(|s| s.parse().ok())
            .unwrap_or(1.0);

        // HETATM flag
        atom.hetatm = get_col("_atom_site.group_PDB")
            .map(|s| s == "HETATM")
            .unwrap_or(false);

        // Formal charge
        atom.formal_charge = get_col("_atom_site.pdbx_formal_charge")
            .and_then(|s| s.parse().ok())
            .unwrap_or(0);

        // Serial number
        atom.id = get_col("_atom_site.id")
            .and_then(|s| s.parse().ok())
            .unwrap_or(0);

        // Coordinates
        let x: f32 = get_col("_atom_site.Cartn_x")
            .and_then(|s| s.parse().ok())
            .unwrap_or(0.0);
        let y: f32 = get_col("_atom_site.Cartn_y")
            .and_then(|s| s.parse().ok())
            .unwrap_or(0.0);
        let z: f32 = get_col("_atom_site.Cartn_z")
            .and_then(|s| s.parse().ok())
            .unwrap_or(0.0);

        // Model number
        let model_num: i32 = get_col("_atom_site.pdbx_PDB_model_num")
            .and_then(|s| s.parse().ok())
            .unwrap_or(1);

        // Add atom (only for first model, subsequent models add coordinates only)
        if model_num == 1 {
            mol.add_atom(atom);
            coords.push(Vec3::new(x, y, z));
        } else {
            models.entry(model_num).or_default().push(Vec3::new(x, y, z));
        }
    }

    // Add primary coordinate set
    if !coords.is_empty() {
        let coord_set = CoordSet::from_vec3(&coords);
        mol.add_coord_set(coord_set);
    }

    // Add additional models as coordinate sets
    let mut model_nums: Vec<i32> = models.keys().copied().collect();
    model_nums.sort();
    for model_num in model_nums {
        if let Some(model_coords) = models.get(&model_num) {
            if model_coords.len() == mol.atom_count() {
                let coord_set = CoordSet::from_vec3(model_coords);
                mol.add_coord_set(coord_set);
            }
        }
    }

    Ok(pos)
}

/// Parse _cell data items
fn parse_cell(tokens: &[Token], mut pos: usize, mol: &mut ObjectMolecule) -> IoResult<usize> {
    let mut a = 1.0f32;
    let mut b = 1.0f32;
    let mut c = 1.0f32;
    let mut alpha = 90.0f32;
    let mut beta = 90.0f32;
    let mut gamma = 90.0f32;

    while pos < tokens.len() {
        match &tokens[pos] {
            Token::DataName(name) if name.starts_with("_cell.") => {
                let field = &name[6..];
                pos += 1;
                if let Some(Token::Value(val)) = tokens.get(pos) {
                    match field {
                        "length_a" => a = val.parse().unwrap_or(1.0),
                        "length_b" => b = val.parse().unwrap_or(1.0),
                        "length_c" => c = val.parse().unwrap_or(1.0),
                        "angle_alpha" => alpha = val.parse().unwrap_or(90.0),
                        "angle_beta" => beta = val.parse().unwrap_or(90.0),
                        "angle_gamma" => gamma = val.parse().unwrap_or(90.0),
                        _ => {}
                    }
                    pos += 1;
                }
            }
            _ => break,
        }
    }

    // Only set symmetry if we have non-default values
    if (a - 1.0).abs() > 0.001 || (b - 1.0).abs() > 0.001 || (c - 1.0).abs() > 0.001 {
        use pymol_mol::Symmetry;
        mol.symmetry = Some(Symmetry::new("P 1", [a, b, c], [alpha, beta, gamma]));
    }

    Ok(pos)
}

/// Parse _struct_conf loop (helices)
///
/// This parses the mmCIF _struct_conf category which contains helix annotations.
fn parse_struct_conf_loop(
    tokens: &[Token],
    mut pos: usize,
    columns: &[String],
    ss_ranges: &mut Vec<SecondaryStructureRange>,
) -> IoResult<usize> {
    let col_idx: HashMap<&str, usize> = columns
        .iter()
        .enumerate()
        .map(|(i, s)| (s.as_str(), i))
        .collect();

    let n_cols = columns.len();

    loop {
        if pos >= tokens.len() {
            break;
        }
        match &tokens[pos] {
            Token::Loop | Token::DataName(_) | Token::DataBlock(_) | Token::Eof => break,
            _ => {}
        }

        // Read one row
        let mut row: Vec<Option<String>> = Vec::with_capacity(n_cols);
        for _ in 0..n_cols {
            if pos >= tokens.len() {
                break;
            }
            match &tokens[pos] {
                Token::Value(v) | Token::SingleQuoted(v) | Token::DoubleQuoted(v) | Token::TextField(v) => {
                    row.push(Some(v.clone()));
                }
                Token::Missing | Token::Unknown => {
                    row.push(None);
                }
                Token::Loop | Token::DataName(_) | Token::DataBlock(_) | Token::Eof => {
                    break;
                }
            }
            pos += 1;
        }

        if row.len() < n_cols {
            break;
        }

        let get_col = |name: &str| -> Option<&str> {
            col_idx.get(name).and_then(|&i| row.get(i).and_then(|v| v.as_deref()))
        };

        // Get helix type - check conf_type_id for HELX variants
        let conf_type = get_col("_struct_conf.conf_type_id").unwrap_or("");
        let ss_type = if conf_type.starts_with("HELX") {
            // HELX_P = alpha helix, HELX_OT_P = other helix types
            if conf_type.contains("310") || conf_type == "HELX_LH_3T_P" {
                SecondaryStructure::Helix310
            } else if conf_type.contains("PI") || conf_type == "HELX_LH_PI_P" {
                SecondaryStructure::HelixPi
            } else {
                SecondaryStructure::Helix
            }
        } else {
            continue; // Skip non-helix entries
        };

        // Get chain and residue range
        // Prefer auth_ columns if available (matches atom coordinates)
        let chain = get_col("_struct_conf.beg_auth_asym_id")
            .or_else(|| get_col("_struct_conf.beg_label_asym_id"))
            .unwrap_or("")
            .to_string();

        let auth_chain = get_col("_struct_conf.beg_auth_asym_id").map(|s| s.to_string());

        let start_seq = get_col("_struct_conf.beg_auth_seq_id")
            .or_else(|| get_col("_struct_conf.beg_label_seq_id"))
            .and_then(|s| s.parse().ok())
            .unwrap_or(0);

        let end_seq = get_col("_struct_conf.end_auth_seq_id")
            .or_else(|| get_col("_struct_conf.end_label_seq_id"))
            .and_then(|s| s.parse().ok())
            .unwrap_or(0);

        if !chain.is_empty() && start_seq > 0 && end_seq > 0 {
            ss_ranges.push(SecondaryStructureRange {
                ss_type,
                chain,
                auth_chain,
                start_seq,
                end_seq,
            });
        }
    }

    Ok(pos)
}

/// Parse _struct_sheet_range loop (beta sheets)
///
/// This parses the mmCIF _struct_sheet_range category which contains sheet annotations.
fn parse_struct_sheet_range_loop(
    tokens: &[Token],
    mut pos: usize,
    columns: &[String],
    ss_ranges: &mut Vec<SecondaryStructureRange>,
) -> IoResult<usize> {
    let col_idx: HashMap<&str, usize> = columns
        .iter()
        .enumerate()
        .map(|(i, s)| (s.as_str(), i))
        .collect();

    let n_cols = columns.len();

    loop {
        if pos >= tokens.len() {
            break;
        }
        match &tokens[pos] {
            Token::Loop | Token::DataName(_) | Token::DataBlock(_) | Token::Eof => break,
            _ => {}
        }

        // Read one row
        let mut row: Vec<Option<String>> = Vec::with_capacity(n_cols);
        for _ in 0..n_cols {
            if pos >= tokens.len() {
                break;
            }
            match &tokens[pos] {
                Token::Value(v) | Token::SingleQuoted(v) | Token::DoubleQuoted(v) | Token::TextField(v) => {
                    row.push(Some(v.clone()));
                }
                Token::Missing | Token::Unknown => {
                    row.push(None);
                }
                Token::Loop | Token::DataName(_) | Token::DataBlock(_) | Token::Eof => {
                    break;
                }
            }
            pos += 1;
        }

        if row.len() < n_cols {
            break;
        }

        let get_col = |name: &str| -> Option<&str> {
            col_idx.get(name).and_then(|&i| row.get(i).and_then(|v| v.as_deref()))
        };

        // Get chain and residue range
        // Prefer auth_ columns if available (matches atom coordinates)
        let chain = get_col("_struct_sheet_range.beg_auth_asym_id")
            .or_else(|| get_col("_struct_sheet_range.beg_label_asym_id"))
            .unwrap_or("")
            .to_string();

        let auth_chain = get_col("_struct_sheet_range.beg_auth_asym_id").map(|s| s.to_string());

        let start_seq = get_col("_struct_sheet_range.beg_auth_seq_id")
            .or_else(|| get_col("_struct_sheet_range.beg_label_seq_id"))
            .and_then(|s| s.parse().ok())
            .unwrap_or(0);

        let end_seq = get_col("_struct_sheet_range.end_auth_seq_id")
            .or_else(|| get_col("_struct_sheet_range.end_label_seq_id"))
            .and_then(|s| s.parse().ok())
            .unwrap_or(0);

        if !chain.is_empty() && start_seq > 0 && end_seq > 0 {
            ss_ranges.push(SecondaryStructureRange {
                ss_type: SecondaryStructure::Sheet,
                chain,
                auth_chain,
                start_seq,
                end_seq,
            });
        }
    }

    Ok(pos)
}

/// Apply secondary structure annotations to atoms in the molecule
fn apply_secondary_structure(mol: &mut ObjectMolecule, ss_ranges: &[SecondaryStructureRange]) {
    if ss_ranges.is_empty() {
        return;
    }

    // Collect atom indices that match each SS range
    let mut updates: Vec<(AtomIndex, SecondaryStructure)> = Vec::new();

    for range in ss_ranges {
        for (idx, atom) in mol.atoms_indexed() {
            // Match by chain and residue number
            // Try auth_chain first if available, otherwise use chain
            let chain_match = if let Some(ref auth_chain) = range.auth_chain {
                atom.chain == *auth_chain
            } else {
                atom.chain == range.chain
            };

            if chain_match && atom.resv >= range.start_seq && atom.resv <= range.end_seq {
                updates.push((idx, range.ss_type));
            }
        }
    }

    // Apply updates
    for (idx, ss_type) in updates {
        if let Some(atom) = mol.get_atom_mut(idx) {
            atom.ss_type = ss_type;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_simple_cif() {
        let cif_data = r#"data_TEST
_entry.id TEST
loop_
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
1 N N ALA A 1 0.000 0.000 0.000 1.00 20.00
2 C CA ALA A 1 1.458 0.000 0.000 1.00 20.00
3 C C ALA A 1 2.009 1.420 0.000 1.00 20.00
"#;

        let mut reader = CifReader::new(cif_data.as_bytes());
        let mol = reader.read().unwrap();

        assert_eq!(mol.name, "TEST");
        assert_eq!(mol.atom_count(), 3);
        assert_eq!(mol.state_count(), 1);
    }

    #[test]
    fn test_parse_with_cell() {
        let cif_data = r#"data_1ABC
_cell.length_a 50.0
_cell.length_b 60.0
_cell.length_c 70.0
_cell.angle_alpha 90.0
_cell.angle_beta 90.0
_cell.angle_gamma 90.0
loop_
_atom_site.id
_atom_site.type_symbol
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
1 C 0.0 0.0 0.0
"#;

        let mut reader = CifReader::new(cif_data.as_bytes());
        let mol = reader.read().unwrap();

        assert!(mol.symmetry.is_some());
        let sym = mol.symmetry.as_ref().unwrap();
        assert!((sym.cell_lengths[0] - 50.0).abs() < 0.001);
    }
}
