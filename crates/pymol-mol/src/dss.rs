//! Secondary structure assignment bridge
//!
//! Extracts backbone data from an [`ObjectMolecule`], delegates to the DSS
//! algorithm in `pymol-algos`, and applies the results back to the molecule.

use lin_alg::f32::Vec3;
use pymol_algos::dss::{BackboneResidue, SsType};

use crate::coordset::CoordSet;
use crate::index::AtomIndex;
use crate::molecule::ObjectMolecule;
use crate::residue::ResidueView;
use crate::secondary::SecondaryStructure;

// Re-export algorithm types for backward compatibility
pub use pymol_algos::dss::{AngleWindow, DssParams as DssSettings};

// ============================================================================
// Public API
// ============================================================================

/// Assign secondary structure to a molecule using the DSS algorithm
///
/// Extracts backbone residue data, runs the PyMOL DSS algorithm, and writes
/// the resulting secondary structure assignments back to each atom's `ss_type`.
///
/// Returns the number of atoms whose secondary structure was changed.
pub fn assign_secondary_structure(
    molecule: &mut ObjectMolecule,
    state: usize,
    settings: &DssSettings,
) -> usize {
    let backbone = extract_backbone(molecule, state);
    if backbone.is_empty() {
        return 0;
    }

    let assignments = pymol_algos::dss::dss(&backbone, settings);

    apply_ss(molecule, &backbone, &assignments)
}

// ============================================================================
// Backbone Extraction
// ============================================================================

/// Intermediate data from residue extraction, before backbone positions are resolved
struct RawResidue {
    ca: Vec3,
    n: Vec3,
    c: Vec3,
    o: Vec3,
    n_idx: AtomIndex,
    c_idx: AtomIndex,
    nh_direction: Option<Vec3>,
    chain: String,
    resv: i32,
}

fn extract_backbone(molecule: &ObjectMolecule, state: usize) -> Vec<BackboneResidue> {
    let coord_set = match molecule.get_coord_set(state) {
        Some(cs) => cs,
        None => return Vec::new(),
    };

    let raw: Vec<RawResidue> = molecule
        .residues()
        .filter(|r| crate::residue::is_amino_acid(&r.key.resn))
        .filter_map(|residue| {
            let ca_idx = find_atom_by_name(&residue, "CA")?;
            let n_idx = find_atom_by_name(&residue, "N")?;
            let c_idx = find_atom_by_name(&residue, "C")?;
            let o_idx = find_atom_by_name(&residue, "O")?;

            Some(RawResidue {
                ca: coord_set.get_atom_coord(ca_idx)?,
                n: coord_set.get_atom_coord(n_idx)?,
                c: coord_set.get_atom_coord(c_idx)?,
                o: coord_set.get_atom_coord(o_idx)?,
                n_idx,
                c_idx,
                nh_direction: compute_nh_direction(molecule, coord_set, n_idx),
                chain: residue.key.chain.clone(),
                resv: residue.key.resv,
            })
        })
        .collect();

    raw.iter()
        .enumerate()
        .map(|(i, r)| {
            let bonded_to_prev = i > 0
                && raw[i - 1].chain == r.chain
                && molecule.find_bond(raw[i - 1].c_idx, r.n_idx).is_some();

            BackboneResidue {
                ca: r.ca,
                n: r.n,
                c: r.c,
                o: r.o,
                chain: r.chain.clone(),
                resv: r.resv,
                nh_direction: r.nh_direction,
                bonded_to_prev,
            }
        })
        .collect()
}

fn find_atom_by_name(residue: &ResidueView, name: &str) -> Option<AtomIndex> {
    residue
        .iter_indexed()
        .find_map(|(idx, atom)| (&*atom.name == name).then_some(idx))
}

/// Compute the N-H bond direction for a nitrogen atom.
///
/// If an explicit hydrogen is bonded to N, uses the N→H direction.
/// Otherwise approximates from the average of heavy-atom bond directions.
fn compute_nh_direction(
    molecule: &ObjectMolecule,
    coord_set: &CoordSet,
    n_idx: AtomIndex,
) -> Option<Vec3> {
    let n_pos = coord_set.get_atom_coord(n_idx)?;
    let bonded: smallvec::SmallVec<[AtomIndex; 4]> = molecule.bonded_atoms_iter(n_idx).collect();
    if bonded.is_empty() {
        return None;
    }

    // Check if there's a bonded hydrogen — use its position directly
    for &b_idx in &bonded {
        if let Some(b_atom) = molecule.get_atom(b_idx) {
            if b_atom.is_hydrogen() {
                if let Some(h_pos) = coord_set.get_atom_coord(b_idx) {
                    let d = h_pos - n_pos;
                    let len = d.magnitude();
                    if len > 1e-6 {
                        return Some(d * (1.0 / len));
                    }
                }
            }
        }
    }

    // No hydrogen found — approximate from heavy atom geometry
    let mut avg_dir = Vec3::new(0.0, 0.0, 0.0);
    let mut count = 0;
    for &b_idx in &bonded {
        if let Some(b_pos) = coord_set.get_atom_coord(b_idx) {
            let d = b_pos - n_pos;
            let len = d.magnitude();
            if len > 1e-6 {
                avg_dir += d * (1.0 / len);
                count += 1;
            }
        }
    }
    if count == 0 {
        return None;
    }
    let h_dir = avg_dir * -1.0;
    let len = h_dir.magnitude();
    if len < 1e-6 {
        return None;
    }
    Some(h_dir * (1.0 / len))
}

// ============================================================================
// Apply Assignments
// ============================================================================

fn apply_ss(
    molecule: &mut ObjectMolecule,
    backbone: &[BackboneResidue],
    assignments: &[SsType],
) -> usize {
    let mut ss_map = ahash::AHashMap::new();
    for (bb, &ss) in backbone.iter().zip(assignments.iter()) {
        let ss = match ss {
            SsType::Helix => SecondaryStructure::Helix,
            SsType::Sheet => SecondaryStructure::Sheet,
            SsType::Loop => SecondaryStructure::Loop,
        };
        ss_map.insert((bb.chain.clone(), bb.resv), ss);
    }

    let mut updated_count = 0;
    for atom_idx in 0..molecule.atom_count() {
        let idx = AtomIndex(atom_idx as u32);
        if let Some(atom) = molecule.get_atom(idx) {
            let key = (atom.residue.chain.clone(), atom.residue.resv);
            if let Some(&ss) = ss_map.get(&key) {
                if let Some(atom_mut) = molecule.get_atom_mut(idx) {
                    if atom_mut.ss_type != ss {
                        atom_mut.ss_type = ss;
                        updated_count += 1;
                    }
                }
            }
        }
    }

    updated_count
}
