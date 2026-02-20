//! BinaryCIF file parser
//!
//! Parses BinaryCIF format files into ObjectMolecule.

use std::collections::HashMap;
use std::io::Read;
use std::sync::Arc;

use lin_alg::f32::Vec3;
use pymol_mol::{Atom, AtomResidue, CoordSet, Element, ObjectMolecule, SecondaryStructure};

use crate::cif::common::{
    apply_secondary_structure, SecondaryStructureRange, SsCategory,
};
use crate::error::{IoError, IoResult};
use crate::traits::MoleculeReader;

use super::decode::{decode_column, decode_mask, ColumnMask, DecodedColumn};
use super::types::{BcifCategory, BcifDataBlock, BcifFile};

/// BinaryCIF file reader
pub struct BcifReader<R> {
    reader: R,
}

impl<R: Read> BcifReader<R> {
    pub fn new(reader: R) -> Self {
        BcifReader { reader }
    }

    fn parse(&mut self) -> IoResult<ObjectMolecule> {
        let mut bytes = Vec::new();
        self.reader
            .read_to_end(&mut bytes)
            .map_err(IoError::Io)?;

        let file: BcifFile = rmp_serde::from_slice(&bytes)
            .map_err(|e| IoError::parse_msg(format!("MessagePack error: {}", e)))?;

        let block = file
            .data_blocks
            .into_iter()
            .next()
            .ok_or(IoError::EmptyFile)?;

        parse_data_block(block)
    }
}

impl<R: Read> MoleculeReader for BcifReader<R> {
    fn read(&mut self) -> IoResult<ObjectMolecule> {
        self.parse()
    }
}

/// Decoded columns for a category, keyed by column name
struct CategoryColumns {
    columns: HashMap<String, (DecodedColumn, Option<ColumnMask>)>,
}

impl CategoryColumns {
    fn decode(category: &BcifCategory) -> IoResult<Self> {
        let mut columns = HashMap::new();
        for col in &category.columns {
            let decoded = decode_column(&col.data)?;
            let mask = decode_mask(&col.mask)?;
            columns.insert(col.name.clone(), (decoded, mask));
        }
        Ok(CategoryColumns { columns })
    }

    /// Decode only specific columns by name
    fn decode_selected(category: &BcifCategory, names: &[&str]) -> IoResult<Self> {
        let mut columns = HashMap::new();
        for col in &category.columns {
            if names.contains(&col.name.as_str()) {
                let decoded = decode_column(&col.data)?;
                let mask = decode_mask(&col.mask)?;
                columns.insert(col.name.clone(), (decoded, mask));
            }
        }
        Ok(CategoryColumns { columns })
    }

    fn str_at(&self, name: &str, i: usize) -> Option<&str> {
        let (col, mask) = self.columns.get(name)?;
        if let Some(m) = mask {
            if !m.is_present(i) {
                return None;
            }
        }
        col.str_at(i).filter(|s| !s.is_empty())
    }

    fn int_at(&self, name: &str, i: usize) -> Option<i32> {
        let (col, mask) = self.columns.get(name)?;
        if let Some(m) = mask {
            if !m.is_present(i) {
                return None;
            }
        }
        col.int_at(i)
    }

    fn float_at(&self, name: &str, i: usize) -> Option<f32> {
        let (col, mask) = self.columns.get(name)?;
        if let Some(m) = mask {
            if !m.is_present(i) {
                return None;
            }
        }
        col.float_at(i)
    }
}

fn parse_data_block(block: BcifDataBlock) -> IoResult<ObjectMolecule> {
    let mut mol = ObjectMolecule::new(&block.header);
    let mut ss_ranges: Vec<SecondaryStructureRange> = Vec::new();

    for category in &block.categories {
        match category.name.as_str() {
            "_atom_site" => parse_atom_site(category, &mut mol)?,
            "_cell" => parse_cell(category, &mut mol)?,
            "_struct_conf" => parse_ss(category, SsCategory::StructConf, &mut ss_ranges)?,
            "_struct_sheet_range" => {
                parse_ss(category, SsCategory::StructSheetRange, &mut ss_ranges)?;
            }
            "_entry" => parse_entry(category, &mut mol)?,
            _ => {}
        }
    }

    if mol.atom_count() == 0 {
        return Err(IoError::EmptyFile);
    }

    apply_secondary_structure(&mut mol, &ss_ranges);
    mol.classify_atoms();
    mol.generate_bonds(0.6);
    mol.assign_known_residue_bond_orders();

    Ok(mol)
}

fn parse_atom_site(category: &BcifCategory, mol: &mut ObjectMolecule) -> IoResult<()> {
    let row_count = category.row_count as usize;
    let cols = CategoryColumns::decode(category)?;

    let mut coords: Vec<Vec3> = Vec::with_capacity(row_count);
    let mut models: HashMap<i32, Vec<Vec3>> = HashMap::new();
    let mut residue_cache: HashMap<(String, String, i32, char), Arc<AtomResidue>> = HashMap::new();

    for i in 0..row_count {
        // Element
        let element_str = cols
            .str_at("type_symbol", i)
            .or_else(|| cols.str_at("label_atom_id", i))
            .unwrap_or("X");
        let element = Element::from_symbol(element_str).unwrap_or(Element::Unknown);

        // Atom name
        let atom_name = cols
            .str_at("auth_atom_id", i)
            .or_else(|| cols.str_at("label_atom_id", i))
            .unwrap_or("X");

        let mut atom = Atom::new(atom_name, element);

        // Residue info
        let resn = cols
            .str_at("auth_comp_id", i)
            .or_else(|| cols.str_at("label_comp_id", i))
            .unwrap_or("");

        let resv = cols
            .int_at("auth_seq_id", i)
            .or_else(|| cols.int_at("label_seq_id", i))
            .unwrap_or(0);

        let chain = cols
            .str_at("auth_asym_id", i)
            .or_else(|| cols.str_at("label_asym_id", i))
            .unwrap_or("");

        let inscode = cols
            .str_at("pdbx_PDB_ins_code", i)
            .and_then(|s| s.chars().next())
            .unwrap_or(' ');

        let cache_key = (chain.to_string(), resn.to_string(), resv, inscode);
        atom.residue = residue_cache
            .entry(cache_key)
            .or_insert_with(|| {
                Arc::new(AtomResidue::from_parts(
                    chain.to_string(),
                    resn.to_string(),
                    resv,
                    inscode,
                    "",
                ))
            })
            .clone();

        // Alt loc
        atom.alt = cols
            .str_at("label_alt_id", i)
            .and_then(|s| s.chars().next())
            .unwrap_or(' ');

        // B-factor
        atom.b_factor = cols.float_at("B_iso_or_equiv", i).unwrap_or(0.0);

        // Occupancy
        atom.occupancy = cols.float_at("occupancy", i).unwrap_or(1.0);

        // HETATM flag
        atom.state.hetatm = cols
            .str_at("group_PDB", i)
            .map(|s| s == "HETATM")
            .unwrap_or(false);

        // Formal charge
        atom.formal_charge = cols.int_at("pdbx_formal_charge", i).unwrap_or(0) as i8;

        // Serial number
        atom.id = cols.int_at("id", i).unwrap_or(0);

        // Coordinates
        let x = cols.float_at("Cartn_x", i).unwrap_or(0.0);
        let y = cols.float_at("Cartn_y", i).unwrap_or(0.0);
        let z = cols.float_at("Cartn_z", i).unwrap_or(0.0);

        // Model number
        let model_num = cols.int_at("pdbx_PDB_model_num", i).unwrap_or(1);

        if model_num == 1 {
            mol.add_atom(atom);
            coords.push(Vec3::new(x, y, z));
        } else {
            models
                .entry(model_num)
                .or_default()
                .push(Vec3::new(x, y, z));
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

    Ok(())
}

fn parse_cell(category: &BcifCategory, mol: &mut ObjectMolecule) -> IoResult<()> {
    let cols = CategoryColumns::decode_selected(
        category,
        &[
            "length_a",
            "length_b",
            "length_c",
            "angle_alpha",
            "angle_beta",
            "angle_gamma",
        ],
    )?;

    let a = cols.float_at("length_a", 0).unwrap_or(1.0);
    let b = cols.float_at("length_b", 0).unwrap_or(1.0);
    let c = cols.float_at("length_c", 0).unwrap_or(1.0);
    let alpha = cols.float_at("angle_alpha", 0).unwrap_or(90.0);
    let beta = cols.float_at("angle_beta", 0).unwrap_or(90.0);
    let gamma = cols.float_at("angle_gamma", 0).unwrap_or(90.0);

    if (a - 1.0).abs() > 0.001 || (b - 1.0).abs() > 0.001 || (c - 1.0).abs() > 0.001 {
        use pymol_mol::Symmetry;
        mol.symmetry = Some(Symmetry::new("P 1", [a, b, c], [alpha, beta, gamma]));
    }

    Ok(())
}

fn parse_ss(
    category: &BcifCategory,
    ss_category: SsCategory,
    ss_ranges: &mut Vec<SecondaryStructureRange>,
) -> IoResult<()> {
    let cols = CategoryColumns::decode(category)?;
    let row_count = category.row_count as usize;

    // Column names depend on category
    let (conf_type_col, beg_chain_col, beg_seq_col, _end_chain_col, end_seq_col) =
        match ss_category {
            SsCategory::StructConf => (
                Some("conf_type_id"),
                "beg_auth_asym_id",
                "beg_auth_seq_id",
                "end_auth_asym_id",
                "end_auth_seq_id",
            ),
            SsCategory::StructSheetRange => (
                None,
                "beg_auth_asym_id",
                "beg_auth_seq_id",
                "end_auth_asym_id",
                "end_auth_seq_id",
            ),
        };

    for i in 0..row_count {
        let ss_type = match ss_category {
            SsCategory::StructConf => {
                let conf_type = cols.str_at(conf_type_col.unwrap(), i).unwrap_or("");
                if !conf_type.starts_with("HELX") {
                    continue;
                }
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

        let chain = cols
            .str_at(beg_chain_col, i)
            .or_else(|| cols.str_at("beg_label_asym_id", i))
            .filter(|s| !s.is_empty());

        let chain = match chain {
            Some(c) => c.to_string(),
            None => continue,
        };

        let auth_chain = cols.str_at(beg_chain_col, i).map(|s| s.to_string());

        let start_seq = cols
            .int_at(beg_seq_col, i)
            .or_else(|| cols.int_at("beg_label_seq_id", i))
            .filter(|&n| n > 0);
        let end_seq = cols
            .int_at(end_seq_col, i)
            .or_else(|| cols.int_at("end_label_seq_id", i))
            .filter(|&n| n > 0);

        if let (Some(start), Some(end)) = (start_seq, end_seq) {
            ss_ranges.push(SecondaryStructureRange {
                ss_type,
                chain,
                auth_chain,
                start_seq: start,
                end_seq: end,
            });
        }
    }

    Ok(())
}

fn parse_entry(category: &BcifCategory, mol: &mut ObjectMolecule) -> IoResult<()> {
    let cols = CategoryColumns::decode_selected(category, &["id"])?;
    if let Some(id) = cols.str_at("id", 0) {
        if mol.title.is_empty() {
            mol.title = id.to_string();
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    /// Verify that CategoryColumns::str_at filters out empty strings,
    /// treating them as None (missing) — matching text CIF behavior
    /// where `.` and `?` tokens map to None.
    #[test]
    fn test_str_at_filters_empty_strings() {
        let mut columns = HashMap::new();
        columns.insert(
            "type_symbol".to_string(),
            (
                DecodedColumn::String(vec![
                    "".to_string(),    // index 0: empty (placeholder)
                    "O".to_string(),   // index 1: oxygen
                    "H".to_string(),   // index 2: hydrogen
                ]),
                None, // no mask
            ),
        );
        let cols = CategoryColumns { columns };

        // Empty string at index 0 should be filtered to None
        assert_eq!(cols.str_at("type_symbol", 0), None);
        // Valid strings should pass through
        assert_eq!(cols.str_at("type_symbol", 1), Some("O"));
        assert_eq!(cols.str_at("type_symbol", 2), Some("H"));
        // Out of bounds should be None
        assert_eq!(cols.str_at("type_symbol", 99), None);
        // Nonexistent column should be None
        assert_eq!(cols.str_at("nonexistent", 0), None);
    }

    /// Verify that masked entries return None regardless of string value
    #[test]
    fn test_str_at_respects_mask() {
        let mut columns = HashMap::new();
        columns.insert(
            "label_seq_id".to_string(),
            (
                DecodedColumn::String(vec![
                    "1".to_string(),
                    "2".to_string(),
                    "".to_string(), // masked entry
                ]),
                Some(ColumnMask(vec![0, 0, 1])), // 0=present, 1=missing
            ),
        );
        let cols = CategoryColumns { columns };

        assert_eq!(cols.str_at("label_seq_id", 0), Some("1"));
        assert_eq!(cols.str_at("label_seq_id", 1), Some("2"));
        // Masked entry returns None (mask takes priority)
        assert_eq!(cols.str_at("label_seq_id", 2), None);
    }
}
