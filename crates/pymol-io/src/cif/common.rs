//! Shared types and utilities for CIF-based parsers (text CIF and BinaryCIF)

use std::collections::HashMap;

use pymol_mol::{ObjectMolecule, SecondaryStructure};

/// Secondary structure annotation from mmCIF
#[derive(Debug, Clone)]
pub(crate) struct SecondaryStructureRange {
    /// Type: helix or sheet
    pub ss_type: SecondaryStructure,
    /// Chain ID (label_asym_id)
    pub chain: String,
    /// Auth chain ID (auth_asym_id) - used for matching atoms
    pub auth_chain: Option<String>,
    /// Start residue sequence number
    pub start_seq: i32,
    /// End residue sequence number
    pub end_seq: i32,
}

/// Category of secondary structure annotation in mmCIF
#[derive(Debug, Clone, Copy)]
pub(crate) enum SsCategory {
    /// _struct_conf (helices)
    StructConf,
    /// _struct_sheet_range (beta sheets)
    StructSheetRange,
}

impl SsCategory {
    /// Returns the mmCIF column prefix for this category
    pub fn prefix(&self) -> &'static str {
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
pub(crate) fn parse_ss_record(
    category: SsCategory,
    fields: &HashMap<&str, &str>,
) -> Option<SecondaryStructureRange> {
    let get_field = |name: &str| fields.get(name).copied();

    // Determine SS type based on category
    let ss_type = match category {
        SsCategory::StructConf => {
            let conf_type = get_field("conf_type_id").unwrap_or("");
            if !conf_type.starts_with("HELX") {
                return None; // Skip non-helix entries
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

/// Apply secondary structure annotations to atoms in the molecule
///
/// Uses a HashMap keyed by (chain, resv) for O(ranges x avg_span + atoms) instead of
/// the naive O(ranges x atoms) nested loop.
pub(crate) fn apply_secondary_structure(
    mol: &mut ObjectMolecule,
    ss_ranges: &[SecondaryStructureRange],
) {
    if ss_ranges.is_empty() {
        return;
    }

    // Build lookup: (chain, resv) -> SS type
    let mut ss_map: HashMap<(&str, i32), SecondaryStructure> = HashMap::new();
    for range in ss_ranges {
        let chain = range.auth_chain.as_deref().unwrap_or(&range.chain);
        for resv in range.start_seq..=range.end_seq {
            ss_map.insert((chain, resv), range.ss_type);
        }
    }

    // Single pass over atoms
    for atom in mol.atoms_mut() {
        if let Some(&ss) = ss_map.get(&(atom.residue.chain.as_str(), atom.residue.resv)) {
            atom.ss_type = ss;
        }
    }
}
