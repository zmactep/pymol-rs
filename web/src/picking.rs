//! Atom picking from screen coordinates.
//!
//! Provides the serializable `PickHitInfo` type returned from `pick_at_screen`
//! to the JavaScript caller, and a constructor that extracts the relevant
//! fields from a [`PickHit`] using the current `mouse_selection_mode`.

use serde::Serialize;

use pymol_mol::ObjectMolecule;
use pymol_scene::{pick_expression_for_hit, PickHit};

/// Pick result returned to JavaScript as JSON.
#[derive(Serialize)]
pub struct PickHitInfo {
    /// Name of the picked object.
    pub object_name: String,
    /// Zero-based atom index within the object, or `null` for non-atom hits.
    pub atom_index: Option<usize>,
    /// Chain identifier of the hit atom, or `null`.
    pub chain: Option<String>,
    /// Residue sequence number of the hit atom, or `null`.
    pub residue: Option<i32>,
    /// PyMOL selection expression (depends on `mouse_selection_mode`), or `null`.
    pub expression: Option<String>,
}

impl PickHitInfo {
    /// Build a `PickHitInfo` from a raw `PickHit`.
    ///
    /// `mode` is the value of the `mouse_selection_mode` setting (0–6).
    /// `mol` is the molecule that owns the hit atom.
    pub fn from_hit(hit: &PickHit, mode: i32, mol: &ObjectMolecule) -> Self {
        let expression = pick_expression_for_hit(hit, mode, mol);

        let (chain, residue) = hit
            .atom_index
            .and_then(|idx| mol.get_atom(idx))
            .map(|atom| {
                (
                    Some(atom.residue.key.chain.clone()),
                    Some(atom.residue.key.resv),
                )
            })
            .unwrap_or((None, None));

        PickHitInfo {
            object_name: hit.object_name.clone(),
            atom_index: hit.atom_index.map(|idx| idx.as_usize()),
            chain,
            residue,
            expression,
        }
    }
}
