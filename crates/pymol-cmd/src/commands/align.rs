//! Structural alignment command
//!
//! Implements `align mobile, target` with three methods:
//! - `kabsch` (default): direct Kabsch fit on corresponding atoms
//! - `sequence`: Needleman-Wunsch sequence alignment → Cα pairs → Kabsch fit
//! - `ce`: Combinatorial Extension structure-based alignment

use lin_alg::f32::{Mat4, Vec3};
use pymol_algos::{
    ce_align, rmsd, CeParams,
    global_align, superpose, AlignedPair, AlignmentScoring, SuperposeParams, SuperposeResult,
    substitution_matrix,
};
use pymol_mol::{residue_to_char, AtomIndex};

use crate::args::ParsedCommand;
use crate::command::{ArgHint, Command, CommandContext, CommandRegistry, ViewerLike};
use crate::commands::selecting::evaluate_selection;
use crate::error::{CmdError, CmdResult};

pub fn register(registry: &mut CommandRegistry) {
    registry.register(AlignCommand);
    registry.register(RmsdCommand);
}

struct AlignCommand;

impl Command for AlignCommand {
    fn name(&self) -> &str {
        "align"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "align" performs structural superposition of one selection onto another.

USAGE

    align mobile, target [, cycles [, cutoff [, method ]]]

ARGUMENTS

    mobile = string: selection for the mobile object (will be moved)

    target = string: selection for the fixed object (stays in place)

    cycles = int: number of outlier rejection cycles {default: 5}

    cutoff = float: outlier rejection cutoff (distance/RMSD ratio) {default: 2.0}

    method = kabsch, sequence, or ce: alignment method {default: kabsch}
        kabsch   — direct fit on matching atoms (selections must have same size)
        sequence — sequence alignment to find Cα correspondences
        ce       — Combinatorial Extension structure-based alignment

    matrix = string: substitution matrix for sequence alignment {default: blosum62}
        Available: blosum62, blosum50, blosum80, pam250, identity

    gap_open = float: gap opening penalty {default: -10.0}

    gap_extend = float: gap extension penalty {default: -1.0}

    win_size = int: CE fragment window size {default: 8}

    gap_max = int: CE maximum gap between aligned fragments {default: 30}

    d0 = float: CE fragment similarity cutoff in Angstroms {default: 3.0}

    d1 = float: CE fragment compatibility cutoff in Angstroms {default: 4.0}

EXAMPLES

    align chain A, chain B
    align mobile, target, cycles=0
    align 1hpx, 1t46, method=sequence
    align 1hpx, 1t46, method=sequence, matrix=blosum50
    align 1hpx, 1t46, method=sequence, gap_open=-12.0, gap_extend=-2.0
    align 1hpx, 1t46, method=ce
    align 1hpx, 1t46, method=ce, win_size=6
"#
    }

    fn execute<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        args: &ParsedCommand,
    ) -> CmdResult {
        let mobile_sel = args
            .get_str(0)
            .ok_or_else(|| CmdError::MissingArgument("mobile selection".into()))?;
        let target_sel = args
            .get_str(1)
            .ok_or_else(|| CmdError::MissingArgument("target selection".into()))?;

        let cycles = args
            .get_int(2)
            .or_else(|| args.get_named_int("cycles"))
            .unwrap_or(5) as u32;

        let cutoff = args
            .get_float(3)
            .or_else(|| args.get_named_float("cutoff"))
            .map(|f| f as f32)
            .unwrap_or(2.0);

        let method = args
            .get_str(4)
            .or_else(|| args.get_named_str("method"))
            .unwrap_or("kabsch");

        let params = SuperposeParams { cycles, cutoff };

        match method {
            "sequence" | "seq" => {
                let matrix_name = args
                    .get_named_str("matrix")
                    .unwrap_or("blosum62");

                let matrix = substitution_matrix::get_matrix(matrix_name)
                    .ok_or_else(|| CmdError::InvalidArgument {
                        name: "matrix".to_string(),
                        reason: format!(
                            "Unknown substitution matrix '{}'. Available: blosum62, blosum50, blosum80, pam250, identity",
                            matrix_name
                        ),
                    })?;

                let scoring = AlignmentScoring {
                    matrix,
                    gap_open: args.get_named_float("gap_open").unwrap_or(-10.0) as f32,
                    gap_extend: args.get_named_float("gap_extend").unwrap_or(-1.0) as f32,
                };

                align_by_sequence(ctx, mobile_sel, target_sel, &params, &scoring)
            }
            "ce" => align_by_ce(ctx, mobile_sel, target_sel, &params, args),
            "kabsch" | _ => align_by_kabsch(ctx, mobile_sel, target_sel, &params),
        }
    }
}

// ============================================================================
// RMSD command (no superposition)
// ============================================================================

struct RmsdCommand;

impl Command for RmsdCommand {
    fn name(&self) -> &str {
        "rmsd"
    }

    fn aliases(&self) -> &[&str] {
        &["rms_cur"]
    }

    fn arg_hints(&self) -> &[ArgHint] {
        &[ArgHint::Selection, ArgHint::Selection]
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "rmsd" computes the root-mean-square deviation between two selections
    without performing any superposition (atoms stay in place).

USAGE

    rmsd sel1, sel2

ARGUMENTS

    sel1 = string: first atom selection

    sel2 = string: second atom selection (must have same number of atoms)

EXAMPLES

    rmsd chain A, chain B
    rmsd mol1 and name CA, mol2 and name CA

SEE ALSO

    align
"#
    }

    fn execute<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        args: &ParsedCommand,
    ) -> CmdResult {
        let sel1 = args
            .get_str(0)
            .ok_or_else(|| CmdError::MissingArgument("first selection".into()))?;
        let sel2 = args
            .get_str(1)
            .ok_or_else(|| CmdError::MissingArgument("second selection".into()))?;

        let results1 = evaluate_selection(ctx.viewer, sel1)?;
        let results2 = evaluate_selection(ctx.viewer, sel2)?;

        let (obj1, indices1) = first_molecule_selection(&results1, sel1)?;
        let (obj2, indices2) = first_molecule_selection(&results2, sel2)?;

        if indices1.len() != indices2.len() {
            return Err(CmdError::Execution(format!(
                "Selections have different atom counts: {} vs {}",
                indices1.len(),
                indices2.len()
            )));
        }

        let n = indices1.len();
        if n == 0 {
            return Err(CmdError::Execution("Selections are empty".into()));
        }

        let coords1 = extract_coords(ctx.viewer, &obj1, &indices1)?;
        let coords2 = extract_coords(ctx.viewer, &obj2, &indices2)?;

        let value = rmsd(&coords1, &coords2);

        ctx.print(&format!(
            " Executive: RMSD = {:8.3}, {} atoms",
            value, n
        ));

        Ok(())
    }
}

/// Direct Kabsch alignment — selections must have equal atom count.
fn align_by_kabsch(
    ctx: &mut CommandContext<'_, '_, dyn ViewerLike + '_>,
    mobile_sel: &str,
    target_sel: &str,
    params: &SuperposeParams,
) -> CmdResult {
    // Evaluate selections
    let mobile_results = evaluate_selection(ctx.viewer, mobile_sel)?;
    let target_results = evaluate_selection(ctx.viewer, target_sel)?;

    // Collect (object_name, atom_indices) for mobile
    let (mobile_obj, mobile_indices) = first_molecule_selection(&mobile_results, mobile_sel)?;
    let (target_obj, target_indices) = first_molecule_selection(&target_results, target_sel)?;

    if mobile_indices.len() != target_indices.len() {
        return Err(CmdError::Execution(format!(
            "Selections have different atom counts: {} vs {} (use method=sequence for unequal sizes)",
            mobile_indices.len(),
            target_indices.len()
        )));
    }

    let n = mobile_indices.len();
    if n < 3 {
        return Err(CmdError::Execution(format!(
            "Need at least 3 atoms for alignment, got {}",
            n
        )));
    }

    // Extract coordinates
    let mobile_coords = extract_coords(ctx.viewer, &mobile_obj, &mobile_indices)?;
    let target_coords = extract_coords(ctx.viewer, &target_obj, &target_indices)?;

    // Build 1:1 pairs
    let pairs: Vec<(usize, usize)> = (0..n).map(|i| (i, i)).collect();

    // Superpose
    let result = superpose(&mobile_coords, &target_coords, &pairs, params)
        .map_err(|e: pymol_algos::AlignError| CmdError::Execution(e.to_string()))?;

    apply_superpose_transform(ctx, &mobile_obj, &result)?;
    print_superpose_result(ctx, &result, params, &[]);

    Ok(())
}

/// Sequence-based alignment — aligns by matching residue sequences, then fits Cα atoms.
fn align_by_sequence(
    ctx: &mut CommandContext<'_, '_, dyn ViewerLike + '_>,
    mobile_sel: &str,
    target_sel: &str,
    params: &SuperposeParams,
    scoring: &AlignmentScoring,
) -> CmdResult {
    // Evaluate selections
    let mobile_results = evaluate_selection(ctx.viewer, mobile_sel)?;
    let target_results = evaluate_selection(ctx.viewer, target_sel)?;

    let (mobile_obj, mobile_indices) = first_molecule_selection(&mobile_results, mobile_sel)?;
    let (target_obj, target_indices) = first_molecule_selection(&target_results, target_sel)?;

    // Extract residue sequences and Cα atom indices from selections
    let (mobile_seq, mobile_ca) =
        extract_residue_sequence(ctx.viewer, &mobile_obj, &mobile_indices)?;
    let (target_seq, target_ca) =
        extract_residue_sequence(ctx.viewer, &target_obj, &target_indices)?;

    if mobile_seq.is_empty() || target_seq.is_empty() {
        return Err(CmdError::Execution(
            "No protein/nucleic acid residues found in one or both selections".into(),
        ));
    }

    // Sequence alignment
    let alignment = global_align(&mobile_seq, &target_seq, scoring);

    // Collect matched Cα atom indices from sequence alignment
    let mut mobile_ca_indices: Vec<AtomIndex> = Vec::new();
    let mut target_ca_indices: Vec<AtomIndex> = Vec::new();

    for pair in &alignment.pairs {
        if let AlignedPair::Match { source, target } = pair {
            if let (Some(&Some(src_ca_idx)), Some(&Some(tgt_ca_idx))) =
                (mobile_ca.get(*source), target_ca.get(*target))
            {
                mobile_ca_indices.push(src_ca_idx);
                target_ca_indices.push(tgt_ca_idx);
            }
        }
    }

    // Extract coordinates in batch
    let mobile_ca_coords = extract_coords(ctx.viewer, &mobile_obj, &mobile_ca_indices)?;
    let target_ca_coords = extract_coords(ctx.viewer, &target_obj, &target_ca_indices)?;
    let ca_pairs: Vec<(usize, usize)> = (0..mobile_ca_coords.len()).map(|i| (i, i)).collect();

    if ca_pairs.len() < 3 {
        return Err(CmdError::Execution(format!(
            "Too few Cα pairs for alignment (need ≥3, got {})",
            ca_pairs.len()
        )));
    }

    // Superpose using Cα pairs
    let result = superpose(&mobile_ca_coords, &target_ca_coords, &ca_pairs, params)
        .map_err(|e: pymol_algos::AlignError| CmdError::Execution(e.to_string()))?;

    apply_superpose_transform(ctx, &mobile_obj, &result)?;
    print_superpose_result(ctx, &result, params, &[
        format!(
            "   Sequence identity: {:.1}% ({} of {} residues)",
            alignment.identity * 100.0,
            alignment.n_matched,
            mobile_seq.len().max(target_seq.len())
        ),
        format!("   Matched Cα pairs:  {}", ca_pairs.len()),
    ]);

    Ok(())
}

/// CE structural alignment — structure-based alignment using Combinatorial Extension.
fn align_by_ce(
    ctx: &mut CommandContext<'_, '_, dyn ViewerLike + '_>,
    mobile_sel: &str,
    target_sel: &str,
    params: &SuperposeParams,
    args: &ParsedCommand,
) -> CmdResult {
    // Evaluate selections
    let mobile_results = evaluate_selection(ctx.viewer, mobile_sel)?;
    let target_results = evaluate_selection(ctx.viewer, target_sel)?;

    let (mobile_obj, mobile_indices) = first_molecule_selection(&mobile_results, mobile_sel)?;
    let (target_obj, target_indices) = first_molecule_selection(&target_results, target_sel)?;

    // Extract residue sequences + Cα atom indices
    let (_, mobile_ca) =
        extract_residue_sequence(ctx.viewer, &mobile_obj, &mobile_indices)?;
    let (_, target_ca) =
        extract_residue_sequence(ctx.viewer, &target_obj, &target_indices)?;

    // Collect Cα atom indices for residues that have Cα atoms
    let mobile_ca_indices: Vec<AtomIndex> = mobile_ca.iter().filter_map(|opt| *opt).collect();
    let target_ca_indices: Vec<AtomIndex> = target_ca.iter().filter_map(|opt| *opt).collect();

    let mobile_ca_coords = extract_coords(ctx.viewer, &mobile_obj, &mobile_ca_indices)?;
    let target_ca_coords = extract_coords(ctx.viewer, &target_obj, &target_ca_indices)?;

    if mobile_ca_coords.is_empty() || target_ca_coords.is_empty() {
        return Err(CmdError::Execution(
            "No Cα atoms found in one or both selections".into(),
        ));
    }

    // Build CE params from named args
    let ce_params = CeParams {
        win_size: args.get_named_int("win_size").unwrap_or(8) as usize,
        gap_max: args.get_named_int("gap_max").unwrap_or(30) as usize,
        d0: args.get_named_float("d0").unwrap_or(3.0) as f32,
        d1: args.get_named_float("d1").unwrap_or(4.0) as f32,
        ..CeParams::default()
    };

    // Run CE alignment
    let ce_result = ce_align(&mobile_ca_coords, &target_ca_coords, &ce_params)
        .map_err(|e: pymol_algos::AlignError| CmdError::Execution(e.to_string()))?;

    if ce_result.pairs.len() < 3 {
        return Err(CmdError::Execution(format!(
            "CE alignment found too few matching residues ({})",
            ce_result.pairs.len()
        )));
    }

    // Extract aligned Cα coordinates for superposition
    let mut aligned_mobile: Vec<[f32; 3]> = Vec::new();
    let mut aligned_target: Vec<[f32; 3]> = Vec::new();
    let mut superpose_pairs: Vec<(usize, usize)> = Vec::new();
    for &(si, ti) in &ce_result.pairs {
        let idx = aligned_mobile.len();
        aligned_mobile.push(mobile_ca_coords[si]);
        aligned_target.push(target_ca_coords[ti]);
        superpose_pairs.push((idx, idx));
    }

    // Superpose with outlier rejection
    let result = superpose(&aligned_mobile, &aligned_target, &superpose_pairs, params)
        .map_err(|e: pymol_algos::AlignError| CmdError::Execution(e.to_string()))?;

    apply_superpose_transform(ctx, &mobile_obj, &result)?;
    print_superpose_result(ctx, &result, params, &[
        format!(
            "   CE alignment:    {} residue pairs (Z-score: {:.1})",
            ce_result.n_aligned, ce_result.z_score
        ),
    ]);

    Ok(())
}

// ============================================================================
// Helper functions
// ============================================================================

/// Get the first molecule object and its selected atom indices from selection results.
fn first_molecule_selection(
    results: &[(String, pymol_select::SelectionResult)],
    sel_name: &str,
) -> CmdResult<(String, Vec<AtomIndex>)> {
    for (obj_name, sel_result) in results {
        let indices: Vec<AtomIndex> = sel_result.indices().collect();
        if !indices.is_empty() {
            return Ok((obj_name.clone(), indices));
        }
    }
    Err(CmdError::Selection(format!(
        "No atoms matching '{}'",
        sel_name
    )))
}

/// Extract coordinates for selected atoms from a molecule object.
fn extract_coords(
    viewer: &dyn ViewerLike,
    obj_name: &str,
    indices: &[AtomIndex],
) -> CmdResult<Vec<[f32; 3]>> {
    let mol_obj = viewer
        .objects()
        .get_molecule(obj_name)
        .ok_or_else(|| CmdError::Execution(format!("Object '{}' not found", obj_name)))?;
    let mol = mol_obj.molecule();
    let cs = mol
        .current_coord_set()
        .ok_or_else(|| CmdError::Execution(format!("No coordinates for '{}'", obj_name)))?;

    let mut coords = Vec::with_capacity(indices.len());
    for &idx in indices {
        let v = cs
            .get_atom_coord(idx)
            .ok_or_else(|| CmdError::Execution(format!("Missing coord for atom {}", idx.0)))?;
        coords.push([v.x, v.y, v.z]);
    }
    Ok(coords)
}

/// Extract residue sequence and Cα atom indices for selected atoms.
///
/// Returns (sequence_chars, ca_indices) where ca_indices[i] is Some(AtomIndex) if
/// residue i has a Cα in the selection, None otherwise.
fn extract_residue_sequence(
    viewer: &dyn ViewerLike,
    obj_name: &str,
    selected_indices: &[AtomIndex],
) -> CmdResult<(Vec<char>, Vec<Option<AtomIndex>>)> {
    let mol_obj = viewer
        .objects()
        .get_molecule(obj_name)
        .ok_or_else(|| CmdError::Execution(format!("Object '{}' not found", obj_name)))?;
    let mol = mol_obj.molecule();

    // Build a set of selected atom indices for fast lookup
    let selected_set: std::collections::HashSet<u32> =
        selected_indices.iter().map(|idx| idx.0).collect();

    let mut sequence = Vec::new();
    let mut ca_indices = Vec::new();

    for residue in mol.residues() {
        if !residue.is_protein() && !residue.is_nucleic() {
            continue;
        }

        // Check if any atom of this residue is in the selection
        let has_selected = residue
            .iter_indexed()
            .any(|(idx, _)| selected_set.contains(&idx.0));

        if !has_selected {
            continue;
        }

        let ch = residue_to_char(residue.resn());
        sequence.push(ch);

        // Find Cα (or C3' for nucleic) in the selection
        let ca = residue.ca().and_then(|(idx, _)| {
            if selected_set.contains(&idx.0) {
                Some(idx)
            } else {
                None
            }
        });
        ca_indices.push(ca);
    }

    Ok((sequence, ca_indices))
}

/// Apply the superposition transform to the mobile object and request a redraw.
fn apply_superpose_transform(
    ctx: &mut CommandContext<'_, '_, dyn ViewerLike + '_>,
    obj_name: &str,
    result: &SuperposeResult,
) -> CmdResult {
    let transform =
        build_transform_mat4(&result.transform.rotation, &result.transform.translation);
    apply_transform_to_object(ctx, obj_name, &transform)?;
    ctx.viewer.request_redraw();
    Ok(())
}

/// Print superposition results. `extra_lines` are method-specific lines
/// printed between the summary and the initial/final RMSD.
fn print_superpose_result(
    ctx: &mut CommandContext<'_, '_, dyn ViewerLike + '_>,
    result: &SuperposeResult,
    params: &SuperposeParams,
    extra_lines: &[String],
) {
    if ctx.quiet {
        return;
    }
    ctx.print(&format!(
        " Executive: RMSD = {:8.3}, {} to {} atoms, {} cycles",
        result.final_rmsd, result.n_aligned, result.n_aligned, result.cycles_performed
    ));
    for line in extra_lines {
        ctx.print(line);
    }
    if params.cycles > 0 {
        ctx.print(&format!(
            "   Initial RMSD:    {:.3}",
            result.initial_rmsd
        ));
        ctx.print(&format!(
            "   Final RMSD:      {:.3} ({} atoms after rejection of {})",
            result.final_rmsd, result.n_aligned, result.n_rejected
        ));
    }
}

/// Build a combined 4×4 transform matrix from rotation + translation.
///
/// The Mat4 from Kabsch has rotation in the 3×3 upper-left.
/// We embed the translation in column 3 (row-major indices 3, 7, 11).
fn build_transform_mat4(rotation: &Mat4, translation: &Vec3) -> Mat4 {
    let r = &rotation.data;
    // Row-major: data[row*4 + col]
    Mat4::new([
        r[0], r[1], r[2], translation.x,
        r[4], r[5], r[6], translation.y,
        r[8], r[9], r[10], translation.z,
        0.0, 0.0, 0.0, 1.0,
    ])
}

/// Apply a rigid-body transform to all atoms in an object (all states).
fn apply_transform_to_object(
    ctx: &mut CommandContext<'_, '_, dyn ViewerLike + '_>,
    obj_name: &str,
    transform: &Mat4,
) -> CmdResult {
    let mol_obj = ctx
        .viewer
        .objects_mut()
        .get_molecule_mut(obj_name)
        .ok_or_else(|| CmdError::Execution(format!("Object '{}' not found", obj_name)))?;
    mol_obj.molecule_mut().transform_all_states(transform);
    Ok(())
}
