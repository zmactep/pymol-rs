use std::collections::{HashMap, HashSet, VecDeque};
use std::path::PathBuf;

use patinae_mol::{AtomIndex, ObjectMolecule};

use crate::topology::{
    ordered_pair, AtomKey, AtomType, CoordinateSource, ForceField, HydrogenRule,
    ParameterizedAngle, ParameterizedAtom, ParameterizedBond, ParameterizedSystem,
    RebuiltHydrogenSource, ResidueAtom, ResidueTemplate, SelectedMolecule,
};

const HYDROGEN_GB_RADIUS_A: f32 = 1.2;
// Disulfide S-S bonds are about 2.03 A; this cutoff tolerates crystallographic noise.
const DISULFIDE_SG_DISTANCE_CUTOFF_A: f32 = 2.3;
const COMPATIBILITY_ALIAS_REPORT_LIMIT: usize = 6;

#[derive(Debug, Clone)]
pub struct ParameterizeInput {
    pub force_field_path: String,
    pub plugin_dirs: Vec<PathBuf>,
    pub selection: String,
    pub molecules: Vec<SelectedMolecule>,
}

pub fn parameterize(
    force_field: &ForceField,
    source: String,
    selection: String,
    molecules: &[SelectedMolecule],
) -> Result<ParameterizedSystem, String> {
    parameterize_with(force_field, source, selection, molecules, &HashMap::new())
}

/// Like [`parameterize`], but `forced_variants` (keyed by [`residue_key`]) lets
/// the caller override which template a residue resolves to — used by
/// pH-based protonation to select e.g. `HIP`/`HIE`, `ASH`, `LYN`.
pub fn parameterize_with(
    force_field: &ForceField,
    source: String,
    selection: String,
    molecules: &[SelectedMolecule],
    forced_variants: &HashMap<String, String>,
) -> Result<ParameterizedSystem, String> {
    let mut report = String::new();
    let mut atoms = Vec::new();
    let mut bonds = Vec::new();
    let mut key_to_index = HashMap::new();
    let mut residue_atoms: HashMap<String, Vec<usize>> = HashMap::new();
    let mut failures = Vec::new();
    let mut compatibility_aliases = Vec::new();
    let mut restored_hydrogens = 0usize;
    let mut approximated_controls = 0usize;

    for selected in molecules {
        let molecule = &selected.molecule;
        let selected_residues =
            selected_residue_entries(&selected.object, molecule, selected.selection.raw_indices());

        for selected_residue in selected_residues {
            let res_key = &selected_residue.key;
            let resn = forced_variants
                .get(res_key)
                .map(String::as_str)
                .unwrap_or_else(|| selected_residue.resn.as_str());
            let Some((template, residue_alias)) =
                resolve_residue_template(force_field, resn, selected_residue.context)
            else {
                failures.push(format!(
                    "{}: residue {} is not present in force field {}",
                    res_key, resn, force_field.name
                ));
                continue;
            };
            if let Some(alias) = residue_alias {
                compatibility_aliases.push(format!("{res_key}: residue {resn}->{alias}"));
            }

            let residue_match =
                match_residue_atoms(molecule, &selected_residue.atom_indices, template);
            for alias in &residue_match.aliases {
                compatibility_aliases.push(format!("{res_key}: atom {alias}"));
            }
            if !residue_match.extra.is_empty() {
                failures.push(format!(
                    "{res_key}: extra atoms: {}",
                    residue_match.extra.join(", ")
                ));
                continue;
            }

            let mut residue_failures = Vec::new();
            let mut planned_atoms = Vec::new();
            for tmpl_atom in &template.atoms {
                let Some(atom_type) = force_field.atom_type(&tmpl_atom.atom_type).cloned() else {
                    residue_failures.push(format!(
                        "{}:{} uses unknown atom type {}",
                        res_key, tmpl_atom.name, tmpl_atom.atom_type
                    ));
                    continue;
                };
                let coordinate = if let Some(&atom_idx) =
                    residue_match.present.get(tmpl_atom.name.as_str())
                {
                    let atom = molecule.get_atom(AtomIndex(atom_idx as u32)).unwrap();
                    PlannedCoordinate::Loaded {
                        atom_idx,
                        gb_radius_a: gb_radius(atom.effective_vdw()),
                    }
                } else if is_hydrogen_name(&tmpl_atom.name) {
                    match find_hydrogen_rule(force_field, template.name.as_str(), &tmpl_atom.name) {
                        Some(match_rule) => PlannedCoordinate::Rebuilt(match_rule),
                        None => {
                            residue_failures.push(format!(
                                "{res_key}: missing hydrogen {} has no HDB rule",
                                tmpl_atom.name
                            ));
                            continue;
                        }
                    }
                } else {
                    residue_failures.push(format!(
                        "{res_key}: missing non-hydrogen atom {}",
                        tmpl_atom.name
                    ));
                    continue;
                };
                planned_atoms.push(PlannedAtom {
                    template: tmpl_atom,
                    atom_type,
                    coordinate,
                });
            }

            if !residue_failures.is_empty() {
                failures.extend(residue_failures);
                continue;
            }

            let global_start = atoms.len();
            let mut local_indices = HashMap::new();
            for (local_idx, planned) in planned_atoms.iter().enumerate() {
                local_indices.insert(planned.template.name.as_str(), local_idx);
            }

            let mut local_atoms = Vec::with_capacity(planned_atoms.len());
            for (local_idx, planned) in planned_atoms.iter().enumerate() {
                let global_idx = global_start + local_idx;
                let (key, source, synthetic_parent, gb_radius_a) = match &planned.coordinate {
                    PlannedCoordinate::Loaded {
                        atom_idx,
                        gb_radius_a,
                    } => {
                        let key = AtomKey {
                            object: selected.object.clone(),
                            atom_index: *atom_idx,
                        };
                        (
                            key.clone(),
                            CoordinateSource::Loaded(key),
                            None,
                            *gb_radius_a,
                        )
                    }
                    PlannedCoordinate::Rebuilt(_) => (
                        synthetic_atom_key(&selected.object, global_idx),
                        CoordinateSource::RebuiltHydrogen(RebuiltHydrogenSource {
                            controls: Vec::new(),
                            bond_length_a: 1.0,
                            slot: 0,
                            count: 1,
                            method: 0,
                        }),
                        None,
                        gb_radius(HYDROGEN_GB_RADIUS_A),
                    ),
                };
                local_atoms.push(ParameterizedAtom {
                    key,
                    source,
                    synthetic_parent,
                    residue_key: res_key.to_string(),
                    residue_name: resn.to_string(),
                    atom_name: planned.template.name.clone(),
                    atom_type: planned.template.atom_type.clone(),
                    charge: planned.template.charge,
                    mass: planned.atom_type.mass,
                    sigma_a: planned.atom_type.sigma_a,
                    epsilon_kj: planned.atom_type.epsilon_kj,
                    gb_radius_a,
                });
            }

            let mut pending_failures = Vec::new();
            for (local_idx, planned) in planned_atoms.iter().enumerate() {
                let PlannedCoordinate::Rebuilt(match_rule) = &planned.coordinate else {
                    continue;
                };
                let (control_locals, approximated) =
                    resolve_hdb_controls(&match_rule.rule.controls, &local_indices);
                approximated_controls += approximated;
                let Some(&parent_local) = control_locals.first() else {
                    pending_failures.push(format!(
                        "{res_key}: missing hydrogen {} cannot resolve parent control atom from HDB controls {}",
                        planned.template.name,
                        match_rule.rule.controls.join(", ")
                    ));
                    continue;
                };
                let controls = control_locals
                    .into_iter()
                    .map(|idx| global_start + idx)
                    .collect::<Vec<_>>();
                let parent_global = global_start + parent_local;
                let bond_length_a = hydrogen_bond_length_a(
                    force_field,
                    &local_atoms[parent_local],
                    &local_atoms[local_idx],
                );
                local_atoms[local_idx].source =
                    CoordinateSource::RebuiltHydrogen(RebuiltHydrogenSource {
                        controls,
                        bond_length_a,
                        slot: match_rule.slot,
                        count: match_rule.rule.count,
                        method: match_rule.rule.method,
                    });
                local_atoms[local_idx].synthetic_parent = Some(parent_global);
                restored_hydrogens += 1;
            }

            if !pending_failures.is_empty() {
                failures.extend(pending_failures);
                continue;
            }

            for (local_idx, atom) in local_atoms.iter().enumerate() {
                if let CoordinateSource::Loaded(key) = &atom.source {
                    key_to_index.insert(key.clone(), global_start + local_idx);
                }
            }
            atoms.extend(local_atoms);

            let residue_global = residue_atoms.entry(res_key.to_string()).or_default();
            residue_global.extend(global_start..atoms.len());
            for pair in &template.bonds {
                let (Some(&a_local), Some(&b_local)) = (
                    local_indices.get(pair[0].as_str()),
                    local_indices.get(pair[1].as_str()),
                ) else {
                    continue;
                };
                let a = global_start + a_local;
                let b = global_start + b_local;
                let type_a = atoms[a].atom_type.as_str();
                let type_b = atoms[b].atom_type.as_str();
                let bond_type = force_field.bond_type(type_a, type_b);
                bonds.push(ParameterizedBond {
                    a,
                    b,
                    length_nm: bond_type.map(|ty| ty.length_nm),
                    force_kj_mol_nm2: bond_type.map(|ty| ty.force_kj_mol_nm2),
                });
            }
        }
    }

    if !failures.is_empty() {
        report.push_str("Parameterization failed.\n");
        for failure in &failures {
            report.push_str(" - ");
            report.push_str(failure);
            report.push('\n');
        }
        append_compatibility_alias_report(&mut report, &compatibility_aliases);
        return Err(report);
    }

    let exclusions = build_exclusions(&bonds);
    let one_four_pairs = build_one_four_pairs(atoms.len(), &bonds, &exclusions);
    let angles = build_angles(force_field, &atoms, &bonds);

    let total_mass: f64 = atoms.iter().map(|atom| atom.mass).sum();
    let net_charge: f64 = atoms.iter().map(|atom| atom.charge).sum();
    let mut residue_names: Vec<&str> = atoms
        .iter()
        .map(|atom| atom.residue_name.as_str())
        .collect();
    residue_names.sort_unstable();
    residue_names.dedup();
    let first_atom = atoms
        .first()
        .map(|atom| format!("{}:{}", atom.residue_key, atom.atom_name))
        .unwrap_or_else(|| "none".into());
    report.push_str(&format!(
        "Parameterized {} atoms, {} bonds, {} angles from {} residues using {}.\n",
        atoms.len(),
        bonds.len(),
        angles.len(),
        residue_atoms.len(),
        force_field.name
    ));
    report.push_str(&format!(
        "Net charge {:.4} e, total mass {:.3} Da, residue types: {}, first atom: {}.\n",
        net_charge,
        total_mass,
        residue_names.join(", "),
        first_atom
    ));
    if force_field.terminal_patches.n_terminal_files + force_field.terminal_patches.c_terminal_files
        > 0
    {
        report.push_str(&format!(
            "Detected terminal patch databases: N={}, C={}; RTP terminal variants are used when present, .tdb patch chemistry is not applied in v1.\n",
            force_field.terminal_patches.n_terminal_files,
            force_field.terminal_patches.c_terminal_files
        ));
    }
    if restored_hydrogens > 0 {
        report.push_str(&format!(
            "Rebuilt {restored_hydrogens} missing hydrogen atom(s) from GROMACS .hdb rules.\n"
        ));
    }
    if approximated_controls > 0 {
        report.push_str(&format!(
            "Approximated {approximated_controls} external HDB control atom reference(s) with local residue atoms.\n"
        ));
    }
    append_compatibility_alias_report(&mut report, &compatibility_aliases);

    Ok(ParameterizedSystem {
        ff_name: force_field.name.clone(),
        source,
        selection,
        defaults: force_field.defaults,
        atoms,
        bonds,
        angles,
        exclusions,
        one_four_pairs,
        key_to_index,
        report,
    })
}

#[derive(Debug)]
struct ResidueAtomMatch {
    present: HashMap<String, usize>,
    extra: Vec<String>,
    aliases: Vec<String>,
}

#[derive(Debug)]
struct PlannedAtom<'a> {
    template: &'a ResidueAtom,
    atom_type: AtomType,
    coordinate: PlannedCoordinate,
}

#[derive(Debug)]
enum PlannedCoordinate {
    Loaded { atom_idx: usize, gb_radius_a: f64 },
    Rebuilt(HydrogenRuleMatch),
}

#[derive(Debug)]
struct HydrogenRuleMatch {
    rule: HydrogenRule,
    slot: usize,
}

#[derive(Debug, Clone)]
struct SelectedResidue {
    key: String,
    chain: String,
    resn: String,
    resv: i32,
    inscode: char,
    atom_indices: Vec<usize>,
    context: ResidueContext,
}

#[derive(Debug, Clone, Copy, Default)]
struct ResidueContext {
    n_terminal: bool,
    c_terminal: bool,
    disulfide_cys: bool,
}

fn selected_residue_entries(
    object: &str,
    molecule: &ObjectMolecule,
    atom_indices: impl IntoIterator<Item = usize>,
) -> Vec<SelectedResidue> {
    let mut selected_atoms = HashSet::new();
    let mut by_residue: HashMap<String, SelectedResidue> = HashMap::new();
    for atom_idx in atom_indices {
        let Some(atom) = molecule.get_atom(AtomIndex(atom_idx as u32)) else {
            continue;
        };
        selected_atoms.insert(atom_idx);
        let key = residue_key(object, atom);
        by_residue
            .entry(key.clone())
            .or_insert_with(|| SelectedResidue {
                key,
                chain: atom.residue.chain.clone(),
                resn: atom.residue.resn.clone(),
                resv: atom.residue.resv,
                inscode: atom.residue.inscode,
                atom_indices: Vec::new(),
                context: ResidueContext::default(),
            })
            .atom_indices
            .push(atom_idx);
    }

    let disulfide_keys = disulfide_residue_keys(object, molecule, &selected_atoms);
    let mut residues = by_residue.into_values().collect::<Vec<_>>();
    for residue in &mut residues {
        residue.context.disulfide_cys = disulfide_keys.contains(&residue.key);
    }
    mark_chain_terminals(&mut residues);
    residues
}

fn mark_chain_terminals(residues: &mut [SelectedResidue]) {
    let mut by_chain: HashMap<String, Vec<usize>> = HashMap::new();
    for (idx, residue) in residues.iter().enumerate() {
        by_chain.entry(residue.chain.clone()).or_default().push(idx);
    }

    for indices in by_chain.values_mut() {
        indices
            .sort_by(|&left, &right| compare_residue_position(&residues[left], &residues[right]));
        if let Some(&first) = indices.first() {
            residues[first].context.n_terminal = true;
        }
        if let Some(&last) = indices.last() {
            residues[last].context.c_terminal = true;
        }
    }
}

fn compare_residue_position(left: &SelectedResidue, right: &SelectedResidue) -> std::cmp::Ordering {
    left.resv
        .cmp(&right.resv)
        .then_with(|| left.inscode.cmp(&right.inscode))
        .then_with(|| left.key.cmp(&right.key))
}

fn disulfide_residue_keys(
    object: &str,
    molecule: &ObjectMolecule,
    selected_atoms: &HashSet<usize>,
) -> HashSet<String> {
    let mut keys = HashSet::new();
    for bond in molecule.bonds() {
        let left_idx = bond.atom1.as_usize();
        let right_idx = bond.atom2.as_usize();
        if !selected_atoms.contains(&left_idx) || !selected_atoms.contains(&right_idx) {
            continue;
        }
        let (Some(left), Some(right)) =
            (molecule.get_atom(bond.atom1), molecule.get_atom(bond.atom2))
        else {
            continue;
        };
        if is_cys_sg(left) && is_cys_sg(right) {
            mark_disulfide_pair(object, left, right, &mut keys);
        }
    }

    let sg_atoms = selected_atoms
        .iter()
        .filter_map(|&idx| {
            let atom = molecule.get_atom(AtomIndex(idx as u32))?;
            is_cys_sg(atom).then_some((idx, atom))
        })
        .collect::<Vec<_>>();
    for (left_pos, (left_idx, left)) in sg_atoms.iter().enumerate() {
        for (right_idx, right) in sg_atoms.iter().skip(left_pos + 1) {
            if left.residue.key == right.residue.key {
                continue;
            }
            if !sg_atoms_within_disulfide_distance(molecule, *left_idx, *right_idx) {
                continue;
            }
            mark_disulfide_pair(object, left, right, &mut keys);
        }
    }
    keys
}

fn is_cys_sg(atom: &patinae_mol::Atom) -> bool {
    atom.residue.resn == "CYS" && atom.name.as_ref() == "SG"
}

fn sg_atoms_within_disulfide_distance(
    molecule: &ObjectMolecule,
    left_idx: usize,
    right_idx: usize,
) -> bool {
    let (Some(left), Some(right)) = (
        molecule.get_coord(AtomIndex(left_idx as u32), 0),
        molecule.get_coord(AtomIndex(right_idx as u32), 0),
    ) else {
        return false;
    };
    (left - right).magnitude() <= DISULFIDE_SG_DISTANCE_CUTOFF_A
}

fn mark_disulfide_pair(
    object: &str,
    left: &patinae_mol::Atom,
    right: &patinae_mol::Atom,
    keys: &mut HashSet<String>,
) {
    keys.insert(residue_key(object, left));
    keys.insert(residue_key(object, right));
}

fn match_residue_atoms(
    molecule: &ObjectMolecule,
    atom_indices: &[usize],
    template: &ResidueTemplate,
) -> ResidueAtomMatch {
    let expected = template.atom_names();
    let mut present: HashMap<String, usize> = HashMap::new();
    let mut extra = Vec::new();
    let mut aliases = Vec::new();

    for &idx in atom_indices {
        let Some(atom) = molecule.get_atom(AtomIndex(idx as u32)) else {
            continue;
        };
        let name = atom.name.as_ref();
        if expected.contains(name) {
            present.insert(name.to_string(), idx);
        } else if let Some(template_name) = atom_alias_for_template(template.name.as_str(), name) {
            if expected.contains(template_name) {
                present.insert(template_name.to_string(), idx);
                aliases.push(format!("{name}->{template_name}"));
            } else {
                extra.push(name.to_string());
            }
        } else {
            extra.push(name.to_string());
        }
    }

    ResidueAtomMatch {
        present,
        extra,
        aliases,
    }
}

fn resolve_residue_template<'a>(
    force_field: &'a ForceField,
    residue_name: &str,
    context: ResidueContext,
) -> Option<(&'a ResidueTemplate, Option<String>)> {
    residue_template_candidates(residue_name, context)
        .into_iter()
        .find_map(|candidate| {
            force_field.residue(&candidate).map(|template| {
                let alias = (candidate != residue_name).then_some(candidate);
                (template, alias)
            })
        })
}

fn residue_template_candidates(residue_name: &str, context: ResidueContext) -> Vec<String> {
    let mut variants = Vec::new();
    if context.disulfide_cys && residue_name == "CYS" {
        push_unique_candidate(&mut variants, "CYX".to_string());
        push_unique_candidate(&mut variants, "CYS2".to_string());
    }
    push_unique_candidate(&mut variants, residue_name.to_string());
    for alias in residue_alias_candidates(residue_name) {
        push_unique_candidate(&mut variants, (*alias).to_string());
    }

    let mut candidates = Vec::new();
    for variant in variants {
        if context.c_terminal {
            push_unique_candidate(&mut candidates, format!("C{variant}"));
        }
        if context.n_terminal {
            push_unique_candidate(&mut candidates, format!("N{variant}"));
        }
        push_unique_candidate(&mut candidates, variant);
    }
    candidates
}

fn push_unique_candidate(candidates: &mut Vec<String>, candidate: String) {
    if !candidates.iter().any(|existing| existing == &candidate) {
        candidates.push(candidate);
    }
}

fn residue_alias_candidates(residue_name: &str) -> &'static [&'static str] {
    match residue_name {
        "HIS" => &[
            "HIE", "HID", "HIP", "HSE", "HSD", "HSP", "HISE", "HISD", "HISH",
        ],
        "HID" | "HSD" | "HISD" => &["HID", "HSD", "HISD"],
        "HIE" | "HSE" | "HISE" => &["HIE", "HSE", "HISE"],
        "HIP" | "HSP" | "HISH" => &["HIP", "HSP", "HISH"],
        _ => &[],
    }
}

fn atom_alias_for_template(template_residue_name: &str, atom_name: &str) -> Option<&'static str> {
    if template_residue_name.ends_with("ILE") && atom_name == "CD1" {
        return Some("CD");
    }
    None
}

fn append_compatibility_alias_report(report: &mut String, aliases: &[String]) {
    if aliases.is_empty() {
        return;
    }
    report.push_str(&format!(
        "Applied {} PDB/RTP compatibility alias(es).\n",
        aliases.len()
    ));
    for alias in aliases.iter().take(COMPATIBILITY_ALIAS_REPORT_LIMIT) {
        report.push_str(" - ");
        report.push_str(alias);
        report.push('\n');
    }
    if aliases.len() > COMPATIBILITY_ALIAS_REPORT_LIMIT {
        report.push_str(&format!(
            " - ... {} additional alias(es)\n",
            aliases.len() - COMPATIBILITY_ALIAS_REPORT_LIMIT
        ));
    }
}

fn find_hydrogen_rule(
    force_field: &ForceField,
    residue_name: &str,
    atom_name: &str,
) -> Option<HydrogenRuleMatch> {
    force_field
        .hydrogen_rules
        .get(residue_name)?
        .iter()
        .find_map(|rule| {
            rule.expanded_names()
                .into_iter()
                .position(|name| name == atom_name)
                .map(|slot| HydrogenRuleMatch {
                    rule: rule.clone(),
                    slot,
                })
        })
}

fn is_hydrogen_name(name: &str) -> bool {
    name.starts_with('H')
}

fn resolve_hdb_controls(
    controls: &[String],
    local_indices: &HashMap<&str, usize>,
) -> (Vec<usize>, usize) {
    let mut indices = Vec::new();
    let mut approximated = 0usize;
    for control in controls {
        if let Some(&idx) = local_indices.get(control.as_str()) {
            indices.push(idx);
            continue;
        }
        if let Some(stripped) = control
            .strip_prefix('-')
            .or_else(|| control.strip_prefix('+'))
        {
            if let Some(&idx) = local_indices.get(stripped) {
                indices.push(idx);
                approximated += 1;
            }
        }
    }
    (indices, approximated)
}

fn hydrogen_bond_length_a(
    force_field: &ForceField,
    parent: &ParameterizedAtom,
    hydrogen: &ParameterizedAtom,
) -> f32 {
    force_field
        .bond_type(parent.atom_type.as_str(), hydrogen.atom_type.as_str())
        .map(|ty| (ty.length_nm * 10.0) as f32)
        .unwrap_or_else(|| fallback_hydrogen_bond_length_a(parent.atom_name.as_str()))
}

fn fallback_hydrogen_bond_length_a(parent_name: &str) -> f32 {
    match parent_name.chars().next().unwrap_or('C') {
        'O' => 0.96,
        'N' => 1.01,
        'S' => 1.34,
        _ => 1.09,
    }
}

fn synthetic_atom_key(object: &str, global_idx: usize) -> AtomKey {
    AtomKey {
        object: object.to_string(),
        atom_index: usize::MAX - global_idx,
    }
}

fn build_exclusions(bonds: &[ParameterizedBond]) -> HashSet<(usize, usize)> {
    let mut exclusions = HashSet::new();
    let mut adjacency: HashMap<usize, Vec<usize>> = HashMap::new();
    for bond in bonds {
        adjacency.entry(bond.a).or_default().push(bond.b);
        adjacency.entry(bond.b).or_default().push(bond.a);
        exclusions.insert(ordered_pair(bond.a, bond.b));
    }

    for &start in adjacency.keys() {
        let mut queue = VecDeque::from([(start, 0usize)]);
        let mut seen = HashSet::from([start]);
        while let Some((node, depth)) = queue.pop_front() {
            if depth == 2 {
                continue;
            }
            for &next in adjacency.get(&node).into_iter().flatten() {
                if seen.insert(next) {
                    let next_depth = depth + 1;
                    if next_depth <= 2 {
                        exclusions.insert(ordered_pair(start, next));
                    }
                    queue.push_back((next, next_depth));
                }
            }
        }
    }

    exclusions
}

fn build_one_four_pairs(
    atom_count: usize,
    bonds: &[ParameterizedBond],
    exclusions: &HashSet<(usize, usize)>,
) -> HashSet<(usize, usize)> {
    let mut adjacency = vec![Vec::new(); atom_count];
    for bond in bonds {
        adjacency[bond.a].push(bond.b);
        adjacency[bond.b].push(bond.a);
    }

    let mut pairs = HashSet::new();
    for start in 0..atom_count {
        let mut queue = VecDeque::from([(start, 0usize)]);
        let mut seen = HashSet::from([start]);
        while let Some((node, depth)) = queue.pop_front() {
            if depth == 3 {
                pairs.insert(ordered_pair(start, node));
                continue;
            }
            for &next in &adjacency[node] {
                if seen.insert(next) {
                    queue.push_back((next, depth + 1));
                }
            }
        }
    }
    pairs.retain(|pair| pair.0 != pair.1 && !exclusions.contains(pair));
    pairs
}

fn build_angles(
    force_field: &ForceField,
    atoms: &[ParameterizedAtom],
    bonds: &[ParameterizedBond],
) -> Vec<ParameterizedAngle> {
    let mut adjacency: HashMap<usize, Vec<usize>> = HashMap::new();
    for bond in bonds {
        adjacency.entry(bond.a).or_default().push(bond.b);
        adjacency.entry(bond.b).or_default().push(bond.a);
    }

    let mut angles = Vec::new();
    for (&center, neighbors) in &adjacency {
        for i in 0..neighbors.len() {
            for j in (i + 1)..neighbors.len() {
                let a = neighbors[i];
                let c = neighbors[j];
                let angle_type = force_field.angle_type(
                    atoms[a].atom_type.as_str(),
                    atoms[center].atom_type.as_str(),
                    atoms[c].atom_type.as_str(),
                );
                angles.push(ParameterizedAngle {
                    a,
                    b: center,
                    c,
                    angle_deg: angle_type.map(|ty| ty.angle_deg),
                    force_kj_mol_rad2: angle_type.map(|ty| ty.force_kj_mol_rad2),
                });
            }
        }
    }
    angles
}

fn gb_radius(vdw_a: f32) -> f64 {
    let radius = f64::from(vdw_a);
    if radius > 0.0 {
        radius
    } else {
        1.5
    }
}

/// Stable per-residue key: `object/chain/resn/resv+inscode`.
pub fn residue_key(object: &str, atom: &patinae_mol::Atom) -> String {
    format!(
        "{object}/{}/{}/{}{}",
        atom.residue.chain, atom.residue.resn, atom.residue.resv, atom.residue.inscode
    )
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;
    use std::path::PathBuf;
    use std::sync::Arc;

    use lin_alg::f32::Vec3;
    use patinae_mol::{Atom, AtomResidue, CoordSet, Element, ObjectMolecule, ResidueKey};
    use patinae_select::SelectionResult;

    use super::*;
    use crate::topology::{
        AtomType, ForceField, ForceFieldDefaults, HydrogenRule, ResidueAtom, ResidueTemplate,
        TerminalPatchSummary,
    };

    fn tiny_force_field() -> ForceField {
        ForceField {
            name: "tiny.ff".into(),
            defaults: ForceFieldDefaults::default(),
            atom_types: HashMap::from([
                (
                    "C".into(),
                    AtomType {
                        mass: 12.011,
                        sigma_a: 3.4,
                        epsilon_kj: 0.4,
                    },
                ),
                (
                    "H".into(),
                    AtomType {
                        mass: 1.008,
                        sigma_a: 2.5,
                        epsilon_kj: 0.1,
                    },
                ),
            ]),
            bond_types: Vec::new(),
            angle_types: Vec::new(),
            residues: HashMap::from([(
                "LIG".into(),
                ResidueTemplate {
                    name: "LIG".into(),
                    atoms: vec![
                        ResidueAtom {
                            name: "C1".into(),
                            atom_type: "C".into(),
                            charge: 0.1,
                        },
                        ResidueAtom {
                            name: "H1".into(),
                            atom_type: "H".into(),
                            charge: -0.1,
                        },
                    ],
                    bonds: vec![["C1".into(), "H1".into()]],
                    exclusions: Vec::new(),
                    impropers: Vec::new(),
                },
            )]),
            hydrogen_rules: HashMap::from([(
                "LIG".into(),
                vec![HydrogenRule {
                    count: 1,
                    method: 1,
                    name: "H1".into(),
                    controls: vec!["C1".into()],
                }],
            )]),
            terminal_patches: TerminalPatchSummary::default(),
        }
    }

    fn molecule_with_residue(resn: &str, atom_names: &[(&str, Element)]) -> ObjectMolecule {
        let atoms = atom_names
            .iter()
            .enumerate()
            .map(|(idx, (name, element))| (*name, *element, Vec3::new(idx as f32, 0.0, 0.0)))
            .collect::<Vec<_>>();
        molecule_with_residue_atoms("A", resn, 1, &atoms)
    }

    fn molecule_with_residue_atoms(
        chain: &str,
        resn: &str,
        resv: i32,
        atoms: &[(&str, Element, Vec3)],
    ) -> ObjectMolecule {
        let mut mol = ObjectMolecule::new("mol");
        let residue = Arc::new(AtomResidue::new(
            ResidueKey::new(chain, resn, resv, ' '),
            String::new(),
        ));

        for (name, element, _) in atoms {
            let mut atom = Atom::new(*name, *element);
            atom.residue = residue.clone();
            mol.add_atom(atom);
        }
        let coords = atoms.iter().map(|(_, _, coord)| *coord).collect::<Vec<_>>();
        mol.add_coord_set(CoordSet::from_vec3(&coords));
        mol
    }

    fn molecule_with_two_cys_sg(distance_a: f32) -> ObjectMolecule {
        let mut mol = ObjectMolecule::new("mol");
        let first_residue = Arc::new(AtomResidue::new(
            ResidueKey::new("A", "CYS", 1, ' '),
            String::new(),
        ));
        let second_residue = Arc::new(AtomResidue::new(
            ResidueKey::new("A", "CYS", 2, ' '),
            String::new(),
        ));
        let mut first = Atom::new("SG", Element::Sulfur);
        first.residue = first_residue;
        mol.add_atom(first);
        let mut second = Atom::new("SG", Element::Sulfur);
        second.residue = second_residue;
        mol.add_atom(second);
        mol.add_coord_set(CoordSet::from_vec3(&[
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(distance_a, 0.0, 0.0),
        ]));
        mol
    }

    fn tiny_molecule(include_h: bool) -> ObjectMolecule {
        let atom_names = if include_h {
            vec![("C1", Element::Carbon), ("H1", Element::Hydrogen)]
        } else {
            vec![("C1", Element::Carbon)]
        };
        molecule_with_residue("LIG", &atom_names)
    }

    fn test_structure_path(relative: &str) -> Option<PathBuf> {
        let repo_fixture = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../..")
            .join("_tests")
            .join(relative);
        if repo_fixture.exists() {
            return Some(repo_fixture);
        }
        let Some(root) = std::env::var_os("TEST_STRUCTURES_DIR") else {
            eprintln!("skipping fixture test: TEST_STRUCTURES_DIR is not set");
            return None;
        };
        let path = PathBuf::from(root).join(relative);
        if path.exists() {
            Some(path)
        } else {
            eprintln!("skipping fixture test: {} does not exist", path.display());
            None
        }
    }

    #[test]
    fn parameterizes_matching_rtp_residue() {
        let ff = tiny_force_field();
        let mol = tiny_molecule(true);
        let selected = SelectedMolecule {
            object: "mol".into(),
            selection: SelectionResult::all(mol.atom_count()),
            molecule: mol,
        };

        let system = parameterize(&ff, "tiny.ff".into(), "all".into(), &[selected]).unwrap();
        assert_eq!(system.atoms.len(), 2);
        assert_eq!(system.bonds.len(), 1);
        assert!(system.report.contains("Parameterized 2 atoms"));
    }

    #[test]
    fn missing_hydrogen_is_rebuilt_from_hdb() {
        let ff = tiny_force_field();
        let mol = tiny_molecule(false);
        let selected = SelectedMolecule {
            object: "mol".into(),
            selection: SelectionResult::all(mol.atom_count()),
            molecule: mol,
        };

        let system = parameterize(&ff, "tiny.ff".into(), "all".into(), &[selected]).unwrap();
        assert_eq!(system.atoms.len(), 2);
        assert_eq!(system.atoms[1].atom_name, "H1");
        assert_eq!(system.atoms[1].synthetic_parent, Some(0));
        assert!(system.report.contains("Rebuilt 1 missing hydrogen"));
    }

    #[test]
    fn pdb_his_can_match_force_field_histidine_variant() {
        let mut ff = tiny_force_field();
        let mut template = ff.residues.get("LIG").unwrap().clone();
        template.name = "HIE".into();
        ff.residues.insert("HIE".into(), template);
        let mol =
            molecule_with_residue("HIS", &[("C1", Element::Carbon), ("H1", Element::Hydrogen)]);
        let selected = SelectedMolecule {
            object: "mol".into(),
            selection: SelectionResult::all(mol.atom_count()),
            molecule: mol,
        };

        let system = parameterize(&ff, "tiny.ff".into(), "all".into(), &[selected]).unwrap();
        assert_eq!(system.atoms.len(), 2);
        assert_eq!(system.atoms[0].residue_name, "HIS");
        assert!(system.report.contains("residue HIS->HIE"));
    }

    #[test]
    fn pdb_ile_cd1_can_match_gromacs_ile_cd() {
        let mut ff = tiny_force_field();
        ff.residues.insert(
            "ILE".into(),
            ResidueTemplate {
                name: "ILE".into(),
                atoms: vec![ResidueAtom {
                    name: "CD".into(),
                    atom_type: "C".into(),
                    charge: 0.0,
                }],
                bonds: Vec::new(),
                exclusions: Vec::new(),
                impropers: Vec::new(),
            },
        );
        let mol = molecule_with_residue("ILE", &[("CD1", Element::Carbon)]);
        let selected = SelectedMolecule {
            object: "mol".into(),
            selection: SelectionResult::all(mol.atom_count()),
            molecule: mol,
        };

        let system = parameterize(&ff, "tiny.ff".into(), "all".into(), &[selected]).unwrap();
        assert_eq!(system.atoms.len(), 1);
        assert_eq!(system.atoms[0].atom_name, "CD");
        assert!(system.report.contains("atom CD1->CD"));
    }

    #[test]
    fn c_terminal_oxt_selects_prefixed_rtp_template() {
        let mut ff = tiny_force_field();
        ff.residues.insert(
            "CLIG".into(),
            ResidueTemplate {
                name: "CLIG".into(),
                atoms: vec![
                    ResidueAtom {
                        name: "C1".into(),
                        atom_type: "C".into(),
                        charge: 0.0,
                    },
                    ResidueAtom {
                        name: "OXT".into(),
                        atom_type: "C".into(),
                        charge: 0.0,
                    },
                ],
                bonds: vec![["C1".into(), "OXT".into()]],
                exclusions: Vec::new(),
                impropers: Vec::new(),
            },
        );
        let mol = molecule_with_residue_atoms(
            "A",
            "LIG",
            1,
            &[
                ("C1", Element::Carbon, Vec3::new(0.0, 0.0, 0.0)),
                ("OXT", Element::Oxygen, Vec3::new(1.2, 0.0, 0.0)),
            ],
        );
        let selected = SelectedMolecule {
            object: "mol".into(),
            selection: SelectionResult::all(mol.atom_count()),
            molecule: mol,
        };

        let system = parameterize(&ff, "tiny.ff".into(), "all".into(), &[selected]).unwrap();
        assert_eq!(system.atoms.len(), 2);
        assert!(system.report.contains("residue LIG->CLIG"));
    }

    #[test]
    fn disulfide_cys_selects_cyx_rtp_template() {
        let mut ff = tiny_force_field();
        ff.residues.insert(
            "CYS".into(),
            ResidueTemplate {
                name: "CYS".into(),
                atoms: vec![
                    ResidueAtom {
                        name: "SG".into(),
                        atom_type: "C".into(),
                        charge: 0.0,
                    },
                    ResidueAtom {
                        name: "HG".into(),
                        atom_type: "H".into(),
                        charge: 0.0,
                    },
                ],
                bonds: vec![["SG".into(), "HG".into()]],
                exclusions: Vec::new(),
                impropers: Vec::new(),
            },
        );
        ff.residues.insert(
            "CYX".into(),
            ResidueTemplate {
                name: "CYX".into(),
                atoms: vec![ResidueAtom {
                    name: "SG".into(),
                    atom_type: "C".into(),
                    charge: 0.0,
                }],
                bonds: Vec::new(),
                exclusions: Vec::new(),
                impropers: Vec::new(),
            },
        );
        let mol = molecule_with_two_cys_sg(2.03);
        let selected = SelectedMolecule {
            object: "mol".into(),
            selection: SelectionResult::all(mol.atom_count()),
            molecule: mol,
        };

        let system = parameterize(&ff, "tiny.ff".into(), "all".into(), &[selected]).unwrap();
        assert_eq!(system.atoms.len(), 2);
        assert!(system.report.contains("residue CYS->CYX"));
    }

    #[test]
    fn fixture_1fdl_parameterizes_with_bundled_amber() {
        let Some(path) = test_structure_path("1fdl.cif") else {
            return;
        };
        let mol = patinae_io::cif::read_cif(&path).unwrap();
        let ff_path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("forcefields/amber19sb.ff");
        let force_field = crate::gromacs::load_force_field(ff_path).unwrap();
        let input_atom_count = mol.atom_count();
        let selected = SelectedMolecule {
            object: "1fdl".into(),
            selection: SelectionResult::all(mol.atom_count()),
            molecule: mol,
        };

        let system = parameterize(&force_field, "AMBER".into(), "all".into(), &[selected])
            .unwrap_or_else(|err| panic!("{err}"));
        assert!(system.report.contains("Parameterized"));
        assert!(system.atoms.len() > input_atom_count);
    }

    #[test]
    fn missing_heavy_atoms_are_fatal() {
        let mut ff = tiny_force_field();
        ff.residues.get_mut("LIG").unwrap().atoms.push(ResidueAtom {
            name: "O1".into(),
            atom_type: "C".into(),
            charge: 0.0,
        });
        let mol = tiny_molecule(false);
        let selected = SelectedMolecule {
            object: "mol".into(),
            selection: SelectionResult::all(mol.atom_count()),
            molecule: mol,
        };

        let err = parameterize(&ff, "tiny.ff".into(), "all".into(), &[selected]).unwrap_err();
        assert!(err.contains("missing non-hydrogen atom O1"));
    }
}
