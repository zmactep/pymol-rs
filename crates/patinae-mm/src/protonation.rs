//! pH-dependent protonation-state assignment for titratable residues.
//!
//! Produces a map of [`residue_key`](crate::parametrize::residue_key) → forced
//! force-field template name, which [`crate::parametrize::parameterize_with`]
//! uses to select the right variant so the HDB step rebuilds the correct
//! hydrogens (e.g. `HIP` vs `HIE`, `ASH`, `GLH`, `LYN`).

use std::collections::HashMap;

use patinae_mol::AtomIndex;

use crate::parametrize::residue_key;
use crate::topology::SelectedMolecule;

// Standard side-chain pKa values (Nozaki–Tanford / common MD defaults).
const ASP_PKA: f64 = 3.65;
const GLU_PKA: f64 = 4.25;
const HIS_PKA: f64 = 6.0;
const LYS_PKA: f64 = 10.5;

/// Returns the force-field template a titratable residue should resolve to at
/// `ph`, or `None` to keep its default protonation.
///
/// - ASP/GLU: protonated (`ASH`/`GLH`) only below their pKa.
/// - HIS: always pinned — `HIP` (charged) below pKa, otherwise neutral `HIE`.
/// - LYS: neutral (`LYN`) only above its pKa.
///
/// CYS is intentionally left to the parametrizer (disulfide → CYX detection),
/// and ARG/TYR stay protonated across the physiological range.
pub fn variant_for(resn: &str, ph: f64) -> Option<&'static str> {
    match resn {
        "ASP" => (ph < ASP_PKA).then_some("ASH"),
        "GLU" => (ph < GLU_PKA).then_some("GLH"),
        "HIS" => Some(if ph < HIS_PKA { "HIP" } else { "HIE" }),
        "LYS" => (ph > LYS_PKA).then_some("LYN"),
        _ => None,
    }
}

/// Forced template overrides for every titratable residue in the selection at
/// the given `ph`. Keyed by [`residue_key`].
pub fn protonation_variants(molecules: &[SelectedMolecule], ph: f64) -> HashMap<String, String> {
    let mut variants = HashMap::new();
    let mut seen = std::collections::HashSet::new();
    for selected in molecules {
        for atom_idx in selected.selection.raw_indices() {
            let Some(atom) = selected.molecule.get_atom(AtomIndex(atom_idx as u32)) else {
                continue;
            };
            let key = residue_key(&selected.object, atom);
            if !seen.insert(key.clone()) {
                continue;
            }
            if let Some(variant) = variant_for(atom.residue.resn.as_str(), ph) {
                variants.insert(key, variant.to_string());
            }
        }
    }
    variants
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn his_is_always_pinned_by_ph() {
        assert_eq!(variant_for("HIS", 5.0), Some("HIP"));
        assert_eq!(variant_for("HIS", 7.4), Some("HIE"));
    }

    #[test]
    fn acids_protonate_only_below_pka() {
        assert_eq!(variant_for("ASP", 2.0), Some("ASH"));
        assert_eq!(variant_for("ASP", 7.4), None);
        assert_eq!(variant_for("GLU", 7.4), None);
    }

    #[test]
    fn lysine_neutralizes_only_above_pka() {
        assert_eq!(variant_for("LYS", 7.4), None);
        assert_eq!(variant_for("LYS", 11.5), Some("LYN"));
    }
}
