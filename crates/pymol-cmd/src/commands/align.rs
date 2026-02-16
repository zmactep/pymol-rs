//! Structural alignment command
//!
//! Implements `align source, target` with two methods:
//! - `kabsch` (default): direct Kabsch fit on corresponding atoms
//! - `sequence`: Needleman-Wunsch sequence alignment → Cα pairs → Kabsch fit

use lin_alg::f32::{Mat4, Vec3};
use pymol_algos::{
    global_align, superpose, AlignedPair, AlignmentScoring, SuperposeParams,
};
use pymol_mol::{residue_to_char, AtomIndex};

use crate::args::ParsedCommand;
use crate::command::{Command, CommandContext, CommandRegistry, ViewerLike};
use crate::commands::selecting::evaluate_selection;
use crate::error::{CmdError, CmdResult};

pub fn register(registry: &mut CommandRegistry) {
    registry.register(AlignCommand);
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

    align source, target [, cycles [, cutoff [, method ]]]

ARGUMENTS

    source = string: selection for the mobile object (will be moved)

    target = string: selection for the fixed object (stays in place)

    cycles = int: number of outlier rejection cycles {default: 5}

    cutoff = float: distance cutoff for outlier rejection in Å {default: 2.0}

    method = kabsch or sequence: alignment method {default: kabsch}
        kabsch   — direct fit on matching atoms (selections must have same size)
        sequence — sequence alignment to find Cα correspondences

EXAMPLES

    align chain A, chain B
    align mobile, target, cycles=0
    align 1hpx, 1t46, method=sequence
"#
    }

    fn execute<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        args: &ParsedCommand,
    ) -> CmdResult {
        let source_sel = args
            .get_str(0)
            .ok_or_else(|| CmdError::MissingArgument("source selection".into()))?;
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
            "sequence" | "seq" => align_by_sequence(ctx, source_sel, target_sel, &params),
            "kabsch" | _ => align_by_kabsch(ctx, source_sel, target_sel, &params),
        }
    }
}

/// Direct Kabsch alignment — selections must have equal atom count.
fn align_by_kabsch(
    ctx: &mut CommandContext<'_, '_, dyn ViewerLike + '_>,
    source_sel: &str,
    target_sel: &str,
    params: &SuperposeParams,
) -> CmdResult {
    // Evaluate selections
    let source_results = evaluate_selection(ctx.viewer, source_sel)?;
    let target_results = evaluate_selection(ctx.viewer, target_sel)?;

    // Collect (object_name, atom_indices) for source
    let (source_obj, source_indices) = first_molecule_selection(&source_results, source_sel)?;
    let (target_obj, target_indices) = first_molecule_selection(&target_results, target_sel)?;

    if source_indices.len() != target_indices.len() {
        return Err(CmdError::Execution(format!(
            "Selections have different atom counts: {} vs {} (use method=sequence for unequal sizes)",
            source_indices.len(),
            target_indices.len()
        )));
    }

    let n = source_indices.len();
    if n < 3 {
        return Err(CmdError::Execution(format!(
            "Need at least 3 atoms for alignment, got {}",
            n
        )));
    }

    // Extract coordinates
    let source_coords = extract_coords(ctx.viewer, &source_obj, &source_indices)?;
    let target_coords = extract_coords(ctx.viewer, &target_obj, &target_indices)?;

    // Build 1:1 pairs
    let pairs: Vec<(usize, usize)> = (0..n).map(|i| (i, i)).collect();

    // Superpose
    let result = superpose(&source_coords, &target_coords, &pairs, params)
        .map_err(|e: pymol_algos::AlignError| CmdError::Execution(e.to_string()))?;

    // Build combined transform matrix (rotation + translation in one Mat4)
    let transform = build_transform_mat4(&result.transform.rotation, &result.transform.translation);

    // Apply to ALL atoms in source object (all states)
    apply_transform_to_object(ctx, &source_obj, &transform)?;

    ctx.viewer.request_redraw();

    // Print results
    if !ctx.quiet {
        ctx.print(&format!(
            " Executive: RMSD = {:8.3}, {} to {} atoms, {} cycles",
            result.final_rmsd,
            result.n_aligned,
            result.n_aligned,
            result.cycles_performed
        ));
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

    Ok(())
}

/// Sequence-based alignment — aligns by matching residue sequences, then fits Cα atoms.
fn align_by_sequence(
    ctx: &mut CommandContext<'_, '_, dyn ViewerLike + '_>,
    source_sel: &str,
    target_sel: &str,
    params: &SuperposeParams,
) -> CmdResult {
    // Evaluate selections
    let source_results = evaluate_selection(ctx.viewer, source_sel)?;
    let target_results = evaluate_selection(ctx.viewer, target_sel)?;

    let (source_obj, source_indices) = first_molecule_selection(&source_results, source_sel)?;
    let (target_obj, target_indices) = first_molecule_selection(&target_results, target_sel)?;

    // Extract residue sequences and Cα atom indices from selections
    let (source_seq, source_ca) =
        extract_residue_sequence(ctx.viewer, &source_obj, &source_indices)?;
    let (target_seq, target_ca) =
        extract_residue_sequence(ctx.viewer, &target_obj, &target_indices)?;

    if source_seq.is_empty() || target_seq.is_empty() {
        return Err(CmdError::Execution(
            "No protein/nucleic acid residues found in one or both selections".into(),
        ));
    }

    // Sequence alignment
    let alignment = global_align(&source_seq, &target_seq, &AlignmentScoring::default());

    // Extract matched Cα pairs
    let mut ca_pairs: Vec<(usize, usize)> = Vec::new();
    let mut source_ca_coords: Vec<[f32; 3]> = Vec::new();
    let mut target_ca_coords: Vec<[f32; 3]> = Vec::new();

    for pair in &alignment.pairs {
        if let AlignedPair::Match { source, target } = pair {
            // Both residues have a position — check if both have Cα
            if let (Some(&src_ca_idx), Some(&tgt_ca_idx)) =
                (source_ca.get(*source), target_ca.get(*target))
            {
                if let (Some(src_ca_idx), Some(tgt_ca_idx)) = (src_ca_idx, tgt_ca_idx) {
                    let src_coord =
                        get_atom_coord(ctx.viewer, &source_obj, src_ca_idx)?;
                    let tgt_coord =
                        get_atom_coord(ctx.viewer, &target_obj, tgt_ca_idx)?;
                    let pair_idx = source_ca_coords.len();
                    source_ca_coords.push(src_coord);
                    target_ca_coords.push(tgt_coord);
                    ca_pairs.push((pair_idx, pair_idx));
                }
            }
        }
    }

    if ca_pairs.len() < 3 {
        return Err(CmdError::Execution(format!(
            "Too few Cα pairs for alignment (need ≥3, got {})",
            ca_pairs.len()
        )));
    }

    // Superpose using Cα pairs
    let result = superpose(&source_ca_coords, &target_ca_coords, &ca_pairs, params)
        .map_err(|e: pymol_algos::AlignError| CmdError::Execution(e.to_string()))?;

    // Build combined transform matrix
    let transform = build_transform_mat4(&result.transform.rotation, &result.transform.translation);

    // Apply to ALL atoms in source object
    apply_transform_to_object(ctx, &source_obj, &transform)?;

    ctx.viewer.request_redraw();

    // Print results
    if !ctx.quiet {
        ctx.print(&format!(
            " Executive: RMSD = {:8.3}, {} to {} atoms, {} cycles",
            result.final_rmsd,
            result.n_aligned,
            result.n_aligned,
            result.cycles_performed
        ));
        ctx.print(&format!(
            "   Sequence identity: {:.1}% ({} of {} residues)",
            alignment.identity * 100.0,
            alignment.n_matched,
            source_seq.len().max(target_seq.len())
        ));
        ctx.print(&format!(
            "   Matched Cα pairs:  {}",
            ca_pairs.len()
        ));
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

/// Get a single atom coordinate.
fn get_atom_coord(
    viewer: &dyn ViewerLike,
    obj_name: &str,
    idx: AtomIndex,
) -> CmdResult<[f32; 3]> {
    let mol_obj = viewer
        .objects()
        .get_molecule(obj_name)
        .ok_or_else(|| CmdError::Execution(format!("Object '{}' not found", obj_name)))?;
    let mol = mol_obj.molecule();
    let cs = mol
        .current_coord_set()
        .ok_or_else(|| CmdError::Execution(format!("No coordinates for '{}'", obj_name)))?;
    let v = cs
        .get_atom_coord(idx)
        .ok_or_else(|| CmdError::Execution(format!("Missing coord for atom {}", idx.0)))?;
    Ok([v.x, v.y, v.z])
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
