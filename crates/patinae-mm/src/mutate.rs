//! Residue substitution: replace a residue's side chain with a target type.
//!
//! Keeps the backbone (N/CA/C/O + amide/HA hydrogens), grafts the target
//! residue's ideal side chain onto the backbone frame, relabels the residue,
//! and regroups atoms so the new ones sit with their residue. Used by the
//! `mutate` tool and the `scan` tool (per candidate).

use std::collections::HashMap;
use std::sync::Arc;

use patinae_mol::{AtomIndex, AtomResidue, BondOrder, ObjectMolecule};

use crate::residue_geometry::ideal_side_chain_on;

const BACKBONE_KEEP: &[&str] = &[
    "N", "CA", "C", "O", "OXT", "H", "H1", "H2", "H3", "HA", "HA2", "HA3",
];

/// Builds a copy of `source` with the residue at (`chain`, `resv`, `inscode`)
/// mutated to `target_resn`. Returns an error if the residue lacks a backbone
/// or the target type has no ideal geometry.
pub fn build_mutant(
    source: &ObjectMolecule,
    chain: &str,
    resv: i32,
    inscode: char,
    target_resn: &str,
) -> Result<ObjectMolecule, String> {
    let mut mol = source.clone();

    let in_residue = |atom: &patinae_mol::Atom| {
        atom.residue.chain == chain && atom.residue.resv == resv && atom.residue.inscode == inscode
    };

    let remove: Vec<AtomIndex> = mol
        .atoms_indexed()
        .filter(|(_, atom)| in_residue(atom) && !BACKBONE_KEEP.contains(&atom.name.as_ref()))
        .map(|(idx, _)| idx)
        .collect();
    mol.remove_atoms(&remove);

    // Re-find the residue's backbone (indices shifted after removal).
    let mut n = None;
    let mut ca = None;
    let mut c = None;
    let mut survivors = Vec::new();
    let mut template = None;
    for (idx, atom) in mol.atoms_indexed() {
        if !in_residue(atom) {
            continue;
        }
        survivors.push(idx);
        if template.is_none() {
            template = Some(atom.clone());
        }
        if let Some(coord) = mol.get_coord(idx, 0) {
            match atom.name.as_ref() {
                "N" => n = Some(coord),
                "CA" => ca = Some(coord),
                "C" => c = Some(coord),
                _ => {}
            }
        }
    }
    let (Some(n), Some(ca), Some(c), Some(template)) = (n, ca, c, template) else {
        return Err("target residue is missing backbone N/CA/C".to_string());
    };

    let side_chain = ideal_side_chain_on(target_resn, n, ca, c)
        .ok_or_else(|| format!("residue {target_resn} is not supported"))?;

    let residue = Arc::new(AtomResidue::from_parts(
        chain,
        target_resn,
        resv,
        inscode,
        template.residue.segi.clone(),
    ));
    let mut name_to_idx: HashMap<String, AtomIndex> = HashMap::new();
    for idx in &survivors {
        if let Some(atom) = mol.get_atom_mut(*idx) {
            atom.residue = residue.clone();
            name_to_idx.insert(atom.name.to_string(), *idx);
        }
    }

    let mut next_id = mol.atoms().map(|a| a.id).max().unwrap_or(0) + 1;
    for (name, element, coord) in &side_chain.atoms {
        let mut atom = template.clone();
        atom.name = name.as_str().into();
        atom.element = *element;
        atom.residue = residue.clone();
        atom.partial_charge = 0.0;
        atom.id = next_id;
        next_id += 1;
        let idx = mol.push_atom_with_coord(atom, *coord);
        name_to_idx.insert(name.clone(), idx);
    }
    for (a, b) in &side_chain.bonds {
        if let (Some(&ia), Some(&ib)) = (name_to_idx.get(a), name_to_idx.get(b)) {
            let _ = mol.add_bond(ia, ib, BondOrder::Single);
        }
    }
    mol.regroup_by_residue();
    Ok(mol)
}
