//! mmCIF file parser
//!
//! Parses macromolecular CIF format files.

use std::collections::HashMap;
use std::io::{BufReader, Read};
use std::sync::Arc;

use lin_alg::f32::Vec3;
use pymol_mol::{Atom, AtomResidue, CoordSet, Element, ObjectMolecule};

use crate::error::{IoError, IoResult};
use crate::traits::MoleculeReader;

use super::common::{
    apply_secondary_structure, parse_ss_record, SecondaryStructureRange, SsCategory,
};
use super::lexer::{tokenize, Token};

/// Pre-resolved column indices for _atom_site loop — avoids per-atom HashMap lookups
#[derive(Default)]
struct AtomSiteColumns {
    type_symbol: Option<usize>,
    label_atom_id: Option<usize>,
    auth_atom_id: Option<usize>,
    label_comp_id: Option<usize>,
    auth_comp_id: Option<usize>,
    label_asym_id: Option<usize>,
    auth_asym_id: Option<usize>,
    label_seq_id: Option<usize>,
    auth_seq_id: Option<usize>,
    pdbx_pdb_ins_code: Option<usize>,
    label_alt_id: Option<usize>,
    b_iso_or_equiv: Option<usize>,
    occupancy: Option<usize>,
    group_pdb: Option<usize>,
    pdbx_formal_charge: Option<usize>,
    id: Option<usize>,
    cartn_x: Option<usize>,
    cartn_y: Option<usize>,
    cartn_z: Option<usize>,
    pdbx_pdb_model_num: Option<usize>,
}

impl AtomSiteColumns {
    /// Build column index struct from the loop column names
    fn from_columns(columns: &[String]) -> Self {
        let mut cols = AtomSiteColumns::default();
        for (i, name) in columns.iter().enumerate() {
            match name.as_str() {
                "_atom_site.type_symbol" => cols.type_symbol = Some(i),
                "_atom_site.label_atom_id" => cols.label_atom_id = Some(i),
                "_atom_site.auth_atom_id" => cols.auth_atom_id = Some(i),
                "_atom_site.label_comp_id" => cols.label_comp_id = Some(i),
                "_atom_site.auth_comp_id" => cols.auth_comp_id = Some(i),
                "_atom_site.label_asym_id" => cols.label_asym_id = Some(i),
                "_atom_site.auth_asym_id" => cols.auth_asym_id = Some(i),
                "_atom_site.label_seq_id" => cols.label_seq_id = Some(i),
                "_atom_site.auth_seq_id" => cols.auth_seq_id = Some(i),
                "_atom_site.pdbx_PDB_ins_code" => cols.pdbx_pdb_ins_code = Some(i),
                "_atom_site.label_alt_id" => cols.label_alt_id = Some(i),
                "_atom_site.B_iso_or_equiv" => cols.b_iso_or_equiv = Some(i),
                "_atom_site.occupancy" => cols.occupancy = Some(i),
                "_atom_site.group_PDB" => cols.group_pdb = Some(i),
                "_atom_site.pdbx_formal_charge" => cols.pdbx_formal_charge = Some(i),
                "_atom_site.id" => cols.id = Some(i),
                "_atom_site.Cartn_x" => cols.cartn_x = Some(i),
                "_atom_site.Cartn_y" => cols.cartn_y = Some(i),
                "_atom_site.Cartn_z" => cols.cartn_z = Some(i),
                "_atom_site.pdbx_PDB_model_num" => cols.pdbx_pdb_model_num = Some(i),
                _ => {}
            }
        }
        cols
    }
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
                if let Some(token) = tokens.get(pos) {
                    if let Some(val) = token.as_str() {
                        if mol.title.is_empty() && !matches!(token, Token::DataBlock(_) | Token::DataName(_)) {
                            mol.title = val.to_string();
                        }
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

/// Extract a `&str` from a value token (zero-copy)
fn token_value_str<'a>(token: &'a Token<'a>) -> Option<&'a str> {
    match token {
        Token::Value(v) | Token::SingleQuoted(v) | Token::DoubleQuoted(v) => Some(v),
        Token::TextField(v) => Some(v.as_str()),
        _ => None,
    }
}

/// Parse _atom_site loop
fn parse_atom_site_loop(
    tokens: &[Token],
    mut pos: usize,
    columns: &[String],
    mol: &mut ObjectMolecule,
) -> IoResult<usize> {
    // Pre-resolve column indices — O(1) lookups per field instead of HashMap
    let cols = AtomSiteColumns::from_columns(columns);

    let n_cols = columns.len();
    let mut coords: Vec<Vec3> = Vec::new();
    let mut models: HashMap<i32, Vec<Vec3>> = HashMap::new();

    // Cache for sharing AtomResidue instances — keyed by borrowed (&str, &str, i32, char)
    // to avoid allocating Strings for cache lookups (only allocate on cache miss)
    let mut residue_cache: HashMap<(&str, &str, i32, char), Arc<AtomResidue>> = HashMap::new();

    // Row buffer — allocated once and reused for every atom (zero-copy borrows)
    let mut row: Vec<Option<&str>> = vec![None; n_cols];

    // Helper: get column value by pre-resolved index
    macro_rules! col {
        ($idx:expr) => {
            $idx.and_then(|i| row[i])
        };
    }

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

        // Reset row
        for cell in &mut row {
            *cell = None;
        }

        // Read one row — zero-copy borrows from tokens
        let mut n_read = 0;
        for _ in 0..n_cols {
            if pos >= tokens.len() {
                break;
            }
            match &tokens[pos] {
                Token::Value(v) | Token::SingleQuoted(v) | Token::DoubleQuoted(v) => {
                    row[n_read] = Some(v);
                }
                Token::TextField(v) => {
                    row[n_read] = Some(v.as_str());
                }
                Token::Missing | Token::Unknown => {
                    row[n_read] = None;
                }
                Token::Loop | Token::DataName(_) | Token::DataBlock(_) | Token::Eof => {
                    break;
                }
            }
            n_read += 1;
            pos += 1;
        }

        if n_read < n_cols {
            // Incomplete row - end of loop
            break;
        }

        // Get element (prefer type_symbol over label_atom_id)
        let element_str = col!(cols.type_symbol)
            .or_else(|| col!(cols.label_atom_id))
            .unwrap_or("X");
        let element = Element::from_symbol(element_str).unwrap_or(Element::Unknown);

        // Get atom name (prefer auth_atom_id over label_atom_id)
        let atom_name = col!(cols.auth_atom_id)
            .or_else(|| col!(cols.label_atom_id))
            .unwrap_or("X");

        let mut atom = Atom::new(atom_name, element);

        // Residue info — zero-copy borrows; only allocate on cache miss
        let resn_str = col!(cols.auth_comp_id)
            .or_else(|| col!(cols.label_comp_id))
            .unwrap_or("");

        let resv: i32 = col!(cols.auth_seq_id)
            .or_else(|| col!(cols.label_seq_id))
            .and_then(|s| s.parse().ok())
            .unwrap_or(0);

        let chain_str = col!(cols.auth_asym_id)
            .or_else(|| col!(cols.label_asym_id))
            .unwrap_or("");

        let inscode = col!(cols.pdbx_pdb_ins_code)
            .and_then(|s| s.chars().next())
            .unwrap_or(' ');

        // Create or reuse shared AtomResidue — borrowed key avoids String allocation on cache hit
        let cache_key = (chain_str, resn_str, resv, inscode);
        atom.residue = residue_cache
            .entry(cache_key)
            .or_insert_with(|| {
                Arc::new(AtomResidue::from_parts(
                    chain_str.to_string(),
                    resn_str.to_string(),
                    resv,
                    inscode,
                    "",
                ))
            })
            .clone();

        // Alt loc
        atom.alt = col!(cols.label_alt_id)
            .and_then(|s| s.chars().next())
            .unwrap_or(' ');

        // B-factor
        atom.b_factor = col!(cols.b_iso_or_equiv)
            .and_then(|s| s.parse().ok())
            .unwrap_or(0.0);

        // Occupancy
        atom.occupancy = col!(cols.occupancy)
            .and_then(|s| s.parse().ok())
            .unwrap_or(1.0);

        // HETATM flag
        atom.state.hetatm = col!(cols.group_pdb)
            .map(|s| s == "HETATM")
            .unwrap_or(false);

        // Formal charge
        atom.formal_charge = col!(cols.pdbx_formal_charge)
            .and_then(|s| s.parse().ok())
            .unwrap_or(0);

        // Serial number
        atom.id = col!(cols.id)
            .and_then(|s| s.parse().ok())
            .unwrap_or(0);

        // Coordinates
        let x: f32 = col!(cols.cartn_x)
            .and_then(|s| s.parse().ok())
            .unwrap_or(0.0);
        let y: f32 = col!(cols.cartn_y)
            .and_then(|s| s.parse().ok())
            .unwrap_or(0.0);
        let z: f32 = col!(cols.cartn_z)
            .and_then(|s| s.parse().ok())
            .unwrap_or(0.0);

        // Model number
        let model_num: i32 = col!(cols.pdbx_pdb_model_num)
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
                if let Some(token) = tokens.get(pos) {
                    if let Some(val) = token_value_str(token) {
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
    let col_to_field: Vec<Option<&str>> = columns
        .iter()
        .map(|col| col.strip_prefix(prefix))
        .collect();

    let n_cols = columns.len();

    // Row buffer — allocated once and reused
    let mut row: Vec<Option<&str>> = vec![None; n_cols];

    loop {
        if pos >= tokens.len() {
            break;
        }
        match &tokens[pos] {
            Token::Loop | Token::DataName(_) | Token::DataBlock(_) | Token::Eof => break,
            _ => {}
        }

        // Reset row
        for cell in &mut row {
            *cell = None;
        }

        // Read one row — zero-copy
        let mut col = 0;
        for _ in 0..n_cols {
            if pos >= tokens.len() {
                break;
            }
            match &tokens[pos] {
                Token::Value(v) | Token::SingleQuoted(v) | Token::DoubleQuoted(v) => {
                    row[col] = Some(v);
                }
                Token::TextField(v) => {
                    row[col] = Some(v.as_str());
                }
                Token::Missing | Token::Unknown => {
                    row[col] = None;
                }
                Token::Loop | Token::DataName(_) | Token::DataBlock(_) | Token::Eof => {
                    break;
                }
            }
            col += 1;
            pos += 1;
        }

        if col < n_cols {
            break;
        }

        // Build field map from row values (zero-copy borrows)
        let mut fields: HashMap<&str, &str> = HashMap::new();
        for (i, field_name) in col_to_field.iter().enumerate() {
            if let (Some(name), Some(value)) = (field_name, row[i]) {
                fields.insert(name, value);
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
    let mut fields: HashMap<&str, &str> = HashMap::new();

    // Collect consecutive data items with matching prefix
    while pos < tokens.len() {
        if let Token::DataName(name) = &tokens[pos] {
            if let Some(field) = name.strip_prefix(prefix) {
                pos += 1;
                // Read the value token
                if pos < tokens.len() {
                    if let Some(val) = token_value_str(&tokens[pos]) {
                        fields.insert(field, val);
                        pos += 1;
                    } else if matches!(&tokens[pos], Token::Missing | Token::Unknown) {
                        // Skip missing/unknown values
                        pos += 1;
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

#[cfg(test)]
mod tests {
    use super::*;
    use pymol_mol::SecondaryStructure;

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
