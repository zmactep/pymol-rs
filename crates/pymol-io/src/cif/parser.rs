//! mmCIF file parser
//!
//! Parses macromolecular CIF format files.

use std::collections::HashMap;
use std::io::{BufReader, Read};
use std::sync::Arc;

use lin_alg::f32::Vec3;
use pymol_mol::{Atom, AtomIndex, AtomResidue, CoordSet, Element, ObjectMolecule, SecondaryStructure};

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

/// Category of secondary structure annotation in mmCIF
#[derive(Debug, Clone, Copy)]
enum SsCategory {
    /// _struct_conf (helices)
    StructConf,
    /// _struct_sheet_range (beta sheets)
    StructSheetRange,
}

impl SsCategory {
    /// Returns the mmCIF column prefix for this category
    fn prefix(&self) -> &'static str {
        match self {
            SsCategory::StructConf => "_struct_conf.",
            SsCategory::StructSheetRange => "_struct_sheet_range.",
        }
    }
}

/// Parse a single secondary structure record from a field map.
///
/// This is the unified function that extracts chain, residue range, and SS type
/// from either loop rows or single-item data.
fn parse_ss_record(
    category: SsCategory,
    fields: &HashMap<String, String>,
) -> Option<SecondaryStructureRange> {
    let get_field = |name: &str| fields.get(name).map(|s| s.as_str());

    // Determine SS type based on category
    let ss_type = match category {
        SsCategory::StructConf => {
            let conf_type = get_field("conf_type_id").unwrap_or("");
            if !conf_type.starts_with("HELX") {
                return None; // Skip non-helix entries
            }
            // HELX_P = alpha helix, HELX_OT_P = other helix types
            if conf_type.contains("310") || conf_type == "HELX_LH_3T_P" {
                SecondaryStructure::Helix310
            } else if conf_type.contains("PI") || conf_type == "HELX_LH_PI_P" {
                SecondaryStructure::HelixPi
            } else {
                SecondaryStructure::Helix
            }
        }
        SsCategory::StructSheetRange => SecondaryStructure::Sheet,
    };

    // Extract chain and residue range (same logic for both categories)
    // Prefer auth_ columns if available (matches atom coordinates)
    let chain = get_field("beg_auth_asym_id")
        .or_else(|| get_field("beg_label_asym_id"))
        .filter(|s| !s.is_empty())?
        .to_string();

    let auth_chain = get_field("beg_auth_asym_id").map(|s| s.to_string());

    let start_seq = get_field("beg_auth_seq_id")
        .or_else(|| get_field("beg_label_seq_id"))
        .and_then(|s| s.parse().ok())
        .filter(|&n: &i32| n > 0)?;

    let end_seq = get_field("end_auth_seq_id")
        .or_else(|| get_field("end_label_seq_id"))
        .and_then(|s| s.parse().ok())
        .filter(|&n: &i32| n > 0)?;

    Some(SecondaryStructureRange {
        ss_type,
        chain,
        auth_chain,
        start_seq,
        end_seq,
    })
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
            // Handle single-item secondary structure format (when not in a loop)
            Token::DataName(name) if name.starts_with("_struct_conf.") => {
                pos = parse_ss_single(tokens, pos, SsCategory::StructConf, &mut ss_ranges);
            }
            Token::DataName(name) if name.starts_with("_struct_sheet_range.") => {
                pos = parse_ss_single(tokens, pos, SsCategory::StructSheetRange, &mut ss_ranges);
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

    // Assign bond orders for known protein residues using atom name templates
    mol.assign_known_residue_bond_orders();

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
        pos = parse_ss_loop(tokens, pos, &columns, SsCategory::StructConf, ss_ranges)?;
    } else if is_struct_sheet {
        pos = parse_ss_loop(tokens, pos, &columns, SsCategory::StructSheetRange, ss_ranges)?;
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

    // Cache for sharing AtomResidue instances among atoms of the same residue
    let mut residue_cache: HashMap<AtomResidue, Arc<AtomResidue>> = HashMap::new();

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
        let resn = get_col("_atom_site.auth_comp_id")
            .or_else(|| get_col("_atom_site.label_comp_id"))
            .unwrap_or("")
            .to_string();

        let resv = get_col("_atom_site.auth_seq_id")
            .or_else(|| get_col("_atom_site.label_seq_id"))
            .and_then(|s| s.parse().ok())
            .unwrap_or(0);

        let chain = get_col("_atom_site.auth_asym_id")
            .or_else(|| get_col("_atom_site.label_asym_id"))
            .unwrap_or("")
            .to_string();

        // Insertion code
        let inscode = get_col("_atom_site.pdbx_PDB_ins_code")
            .and_then(|s| s.chars().next())
            .unwrap_or(' ');

        // Create or reuse shared AtomResidue
        let residue_data = AtomResidue::from_parts(chain, resn, resv, inscode, "");
        atom.residue = residue_cache
            .entry(residue_data.clone())
            .or_insert_with(|| Arc::new(residue_data))
            .clone();

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
        atom.state.hetatm = get_col("_atom_site.group_PDB")
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

/// Parse secondary structure loop (_struct_conf or _struct_sheet_range)
///
/// This unified function handles both helix and sheet annotations in loop format.
fn parse_ss_loop(
    tokens: &[Token],
    mut pos: usize,
    columns: &[String],
    category: SsCategory,
    ss_ranges: &mut Vec<SecondaryStructureRange>,
) -> IoResult<usize> {
    let prefix = category.prefix();

    // Map column indices to short field names (without prefix)
    let col_to_field: Vec<Option<String>> = columns
        .iter()
        .map(|col| col.strip_prefix(prefix).map(|s| s.to_string()))
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

        // Build field map from row values
        let mut fields: HashMap<String, String> = HashMap::new();
        for (i, field_name) in col_to_field.iter().enumerate() {
            if let (Some(name), Some(Some(value))) = (field_name, row.get(i)) {
                fields.insert(name.clone(), value.clone());
            }
        }

        // Use unified record parser
        if let Some(range) = parse_ss_record(category, &fields) {
            ss_ranges.push(range);
        }
    }

    Ok(pos)
}

/// Parse secondary structure from single-item (non-loop) format
///
/// This handles cases where there's only one helix or sheet entry,
/// stored as individual data items rather than a loop.
fn parse_ss_single(
    tokens: &[Token],
    mut pos: usize,
    category: SsCategory,
    ss_ranges: &mut Vec<SecondaryStructureRange>,
) -> usize {
    let prefix = category.prefix();
    let mut fields: HashMap<String, String> = HashMap::new();

    // Collect consecutive data items with matching prefix
    while pos < tokens.len() {
        if let Token::DataName(name) = &tokens[pos] {
            if let Some(field) = name.strip_prefix(prefix) {
                pos += 1;
                // Read the value token
                if pos < tokens.len() {
                    match &tokens[pos] {
                        Token::Value(v)
                        | Token::SingleQuoted(v)
                        | Token::DoubleQuoted(v)
                        | Token::TextField(v) => {
                            fields.insert(field.to_string(), v.clone());
                            pos += 1;
                        }
                        Token::Missing | Token::Unknown => {
                            // Skip missing/unknown values
                            pos += 1;
                        }
                        _ => {}
                    }
                }
                continue;
            }
        }
        break;
    }

    // Use unified record parser
    if !fields.is_empty() {
        if let Some(range) = parse_ss_record(category, &fields) {
            ss_ranges.push(range);
        }
    }

    pos
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
                atom.residue.chain == *auth_chain
            } else {
                atom.residue.chain == range.chain
            };

            if chain_match && atom.residue.resv >= range.start_seq && atom.residue.resv <= range.end_seq {
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

    #[test]
    fn test_parse_single_item_secondary_structure() {
        // Test single-item (non-loop) secondary structure format
        let cif_data = r#"data_TEST
loop_
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
1 N N ALA A 1 0.0 0.0 0.0
2 C CA ALA A 1 1.5 0.0 0.0
3 N N GLU A 2 3.0 0.0 0.0
4 C CA GLU A 2 4.5 0.0 0.0
5 N N LEU A 3 6.0 0.0 0.0
6 C CA LEU A 3 7.5 0.0 0.0
#
_struct_conf.conf_type_id HELX_P
_struct_conf.id HELX_P1
_struct_conf.beg_auth_asym_id A
_struct_conf.beg_auth_seq_id 1
_struct_conf.end_auth_asym_id A
_struct_conf.end_auth_seq_id 3
"#;

        let mut reader = CifReader::new(cif_data.as_bytes());
        let mol = reader.read().unwrap();

        assert_eq!(mol.atom_count(), 6);

        // All atoms should have helix secondary structure
        for atom in mol.atoms() {
            assert_eq!(
                atom.ss_type,
                SecondaryStructure::Helix,
                "Atom {} resv={} should be Helix",
                atom.name,
                atom.residue.resv
            );
        }
    }

    #[test]
    fn test_parse_dna_nucleotides() {
        // Test that DNA nucleotides are parsed and classified correctly
        let cif_data = r#"data_TEST
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1  P P   DA I 1 DA 0.0 0.0 0.0
ATOM 2  C "C4'" DA I 1 DA 1.0 0.0 0.0
ATOM 3  P P   DT I 2 DT 2.0 0.0 0.0
ATOM 4  C "C4'" DT I 2 DT 3.0 0.0 0.0
ATOM 5  P P   DG I 3 DG 4.0 0.0 0.0
ATOM 6  C "C4'" DG I 3 DG 5.0 0.0 0.0
ATOM 7  P P   DC I 4 DC 6.0 0.0 0.0
ATOM 8  C "C4'" DC I 4 DC 7.0 0.0 0.0
"#;

        let mut reader = CifReader::new(cif_data.as_bytes());
        let mol = reader.read().unwrap();

        assert_eq!(mol.atom_count(), 8);

        // All atoms should be classified as nucleic
        use pymol_mol::AtomFlags;
        for atom in mol.atoms() {
            assert!(
                atom.state.flags.contains(AtomFlags::NUCLEIC),
                "Atom {} in residue {} should have NUCLEIC flag",
                atom.name,
                atom.residue.resn
            );
            assert!(
                atom.state.flags.contains(AtomFlags::POLYMER),
                "Atom {} in residue {} should have POLYMER flag",
                atom.name,
                atom.residue.resn
            );
        }

        // Check chain iteration: should have 1 chain (I) with 4 residues
        let chains: Vec<_> = mol.chains().collect();
        assert_eq!(chains.len(), 1, "Should have 1 chain");
        assert_eq!(chains[0].id(), "I");

        let residues: Vec<_> = chains[0].residues().collect();
        assert_eq!(residues.len(), 4, "Chain I should have 4 residues");
        assert_eq!(residues[0].resn(), "DA");
        assert_eq!(residues[1].resn(), "DT");
        assert_eq!(residues[2].resn(), "DG");
        assert_eq!(residues[3].resn(), "DC");

        // All should be classified as nucleic
        for r in &residues {
            assert!(r.is_nucleic(), "Residue {} should be nucleic", r.resn());
        }
    }

    #[test]
    fn test_parse_loop_secondary_structure() {
        // Test loop format secondary structure (existing behavior)
        let cif_data = r#"data_TEST
loop_
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
1 N N ALA A 1 0.0 0.0 0.0
2 C CA ALA A 1 1.5 0.0 0.0
3 N N GLU A 2 3.0 0.0 0.0
4 C CA GLU A 2 4.5 0.0 0.0
#
loop_
_struct_conf.conf_type_id
_struct_conf.id
_struct_conf.beg_auth_asym_id
_struct_conf.beg_auth_seq_id
_struct_conf.end_auth_asym_id
_struct_conf.end_auth_seq_id
HELX_P HELX_P1 A 1 A 2
"#;

        let mut reader = CifReader::new(cif_data.as_bytes());
        let mol = reader.read().unwrap();

        assert_eq!(mol.atom_count(), 4);

        // All atoms should have helix secondary structure
        for atom in mol.atoms() {
            assert_eq!(
                atom.ss_type,
                SecondaryStructure::Helix,
                "Atom {} resv={} should be Helix",
                atom.name,
                atom.residue.resv
            );
        }
    }
}
