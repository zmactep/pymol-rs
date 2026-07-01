//! Structure-building helpers derived from a parameterized system.
//!
//! Parametrization rebuilds missing hydrogens from the force field's HDB rules
//! as synthetic topology atoms. This module turns those synthetic atoms into a
//! concrete list of hydrogens to graft onto a loaded molecule, so a frontend can
//! produce a protonated copy of an object.

use lin_alg::f32::Vec3;

use crate::topology::{CoordinateSource, ParameterizedSystem};

/// One hydrogen to add to a loaded molecule.
#[derive(Debug, Clone)]
pub struct HydrogenAddition {
    /// Object the parent heavy atom belongs to.
    pub object: String,
    /// Index of the parent atom within that object's molecule.
    pub parent_atom_index: usize,
    /// Hydrogen atom name (e.g. `HA`, `HB1`).
    pub name: String,
    /// Residue name of the parent (for reference / display).
    pub residue_name: String,
    /// Rebuilt position.
    pub coord: Vec3,
}

/// Collects the hydrogens that parametrization rebuilt, paired with the loaded
/// parent atom they attach to. `coords` must be the per-topology-atom positions
/// from [`crate::energy::coordinates_for_frame`].
///
/// Synthetic hydrogens whose parent is itself synthetic (no loaded parent) are
/// skipped — they cannot be attached to an existing atom.
pub fn hydrogens_to_add(
    topology: &ParameterizedSystem,
    coords: &[Vec3],
) -> Vec<HydrogenAddition> {
    let mut additions = Vec::new();
    for (idx, atom) in topology.atoms.iter().enumerate() {
        if !matches!(atom.source, CoordinateSource::RebuiltHydrogen(_)) {
            continue;
        }
        let Some(parent_idx) = atom.synthetic_parent else {
            continue;
        };
        let Some(parent) = topology.atoms.get(parent_idx) else {
            continue;
        };
        let CoordinateSource::Loaded(parent_key) = &parent.source else {
            continue;
        };
        let Some(coord) = coords.get(idx).copied() else {
            continue;
        };
        additions.push(HydrogenAddition {
            object: parent_key.object.clone(),
            parent_atom_index: parent_key.atom_index,
            name: atom.atom_name.clone(),
            residue_name: atom.residue_name.clone(),
            coord,
        });
    }
    additions
}
