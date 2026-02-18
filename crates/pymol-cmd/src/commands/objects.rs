//! Object commands: delete, rename, create, copy, group, ungroup

use pymol_mol::{AtomIndex, ObjectMolecule};
use pymol_scene::{DirtyFlags, MoleculeObject};

use crate::args::ParsedCommand;
use crate::command::{ArgHint, Command, CommandContext, CommandRegistry, ViewerLike};
use crate::commands::selecting::evaluate_selection;
use crate::error::{CmdError, CmdResult};

/// Simple glob matching (prefix* and *suffix patterns)
fn glob_match(pattern: &str, name: &str) -> bool {
    if pattern == "*" || pattern == "all" {
        return true;
    }
    if let Some(prefix) = pattern.strip_suffix('*') {
        return name.starts_with(prefix);
    }
    if let Some(suffix) = pattern.strip_prefix('*') {
        return name.ends_with(suffix);
    }
    pattern == name
}

/// Register object commands
pub fn register(registry: &mut CommandRegistry) {
    registry.register(DeleteCommand);
    registry.register(RenameCommand);
    registry.register(CreateCommand);
    registry.register(CopyCommand);
    registry.register(GroupCommand);
    registry.register(UngroupCommand);
    registry.register(StateCommand);
    registry.register(CountStatesCommand);
    registry.register(SplitStatesCommand);
    registry.register(ExtractCommand);
}

// ============================================================================
// delete command
// ============================================================================

struct DeleteCommand;

impl Command for DeleteCommand {
    fn name(&self) -> &str {
        "delete"
    }

    fn arg_hints(&self) -> &[ArgHint] {
        &[ArgHint::Object]
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "delete" removes objects and/or named selections.

USAGE

    delete name

ARGUMENTS

    name = string: object or selection name (supports glob patterns)

EXAMPLES

    delete protein
    delete obj*
    delete sele
    delete all
"#
    }

    fn execute<'v, 'r>(&self, ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>, args: &ParsedCommand) -> CmdResult {
        let name = args
            .get_str(0)
            .or_else(|| args.get_named_str("name"))
            .ok_or_else(|| CmdError::MissingArgument("name".to_string()))?;

        let mut deleted_objects = 0usize;
        let mut deleted_selections = 0usize;

        // Delete matching objects
        let obj_matches: Vec<String> = ctx
            .viewer
            .objects()
            .matching(name)
            .iter()
            .map(|s| s.to_string())
            .collect();

        for obj_name in &obj_matches {
            ctx.viewer.objects_mut().remove(obj_name);
            deleted_objects += 1;
        }

        // Delete matching named selections
        // "all" deletes all selections too
        if name == "all" {
            let sel_names = ctx.viewer.selections().names();
            deleted_selections = sel_names.len();
            for sel_name in sel_names {
                ctx.viewer.remove_selection(&sel_name);
            }
        } else {
            // Try exact match or glob pattern on selections
            let sel_names = ctx.viewer.selections().names();
            for sel_name in &sel_names {
                if sel_name == name || glob_match(name, sel_name) {
                    ctx.viewer.remove_selection(sel_name);
                    deleted_selections += 1;
                }
            }
        }

        let total = deleted_objects + deleted_selections;
        if total == 0 {
            return Err(CmdError::ObjectNotFound(name.to_string()));
        }

        ctx.viewer.request_redraw();

        if !ctx.quiet {
            match (deleted_objects, deleted_selections) {
                (o, 0) if o == 1 => ctx.print(&format!(" Deleted \"{}\"", name)),
                (o, 0) => ctx.print(&format!(" Deleted {} objects matching \"{}\"", o, name)),
                (0, s) if s == 1 => ctx.print(&format!(" Deleted selection \"{}\"", name)),
                (0, s) => ctx.print(&format!(" Deleted {} selections matching \"{}\"", s, name)),
                (o, s) => ctx.print(&format!(" Deleted {} objects and {} selections matching \"{}\"", o, s, name)),
            }
        }

        Ok(())
    }
}

// ============================================================================
// rename command (also: set_name)
// ============================================================================

struct RenameCommand;

impl Command for RenameCommand {
    fn name(&self) -> &str {
        "set_name"
    }

    fn aliases(&self) -> &[&str] {
        &["rename"]
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "set_name" renames an object or selection.

USAGE

    set_name old_name, new_name

ARGUMENTS

    old_name = string: current object or selection name
    new_name = string: new name

EXAMPLES

    set_name protein, myprotein
    rename ligand, drug
"#
    }

    fn execute<'v, 'r>(&self, ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>, args: &ParsedCommand) -> CmdResult {
        let old_name = args
            .get_str(0)
            .or_else(|| args.get_named_str("old_name"))
            .ok_or_else(|| CmdError::MissingArgument("old_name".to_string()))?;

        let new_name = args
            .get_str(1)
            .or_else(|| args.get_named_str("new_name"))
            .ok_or_else(|| CmdError::MissingArgument("new_name".to_string()))?;

        // Try renaming as an object first, then as a selection
        let renamed = match ctx.viewer.objects_mut().rename(old_name, new_name) {
            Ok(()) => true,
            Err(_) => ctx.viewer.rename_selection(old_name, new_name),
        };

        if !renamed {
            return Err(CmdError::ObjectNotFound(old_name.to_string()));
        }

        if !ctx.quiet {
            ctx.print(&format!(" Renamed \"{}\" to \"{}\"", old_name, new_name));
        }

        Ok(())
    }
}

// ============================================================================
// create command
// ============================================================================

struct CreateCommand;

impl Command for CreateCommand {
    fn name(&self) -> &str {
        "create"
    }

    fn arg_hints(&self) -> &[ArgHint] {
        &[ArgHint::None, ArgHint::Selection]
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "create" creates a new object from a selection.

USAGE

    create name, selection [, source_state [, target_state [, discrete ]]]

ARGUMENTS

    name = string: name for the new object
    selection = string: selection of atoms to include
    source_state = integer: state to copy from (default: 0 = all)
    target_state = integer: state to copy to (default: 0)
    discrete = 0/1: create discrete states (default: 0)

EXAMPLES

    create chainA, chain A
    create ligand, organic
"#
    }

    fn execute<'v, 'r>(&self, ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>, args: &ParsedCommand) -> CmdResult {
        let name = args
            .get_str(0)
            .or_else(|| args.get_named_str("name"))
            .ok_or_else(|| CmdError::MissingArgument("name".to_string()))?;

        let selection = args
            .get_str(1)
            .or_else(|| args.get_named_str("selection"))
            .ok_or_else(|| CmdError::MissingArgument("selection".to_string()))?;

        let source_state = args
            .get_int(2)
            .or_else(|| args.get_named_int("source_state"))
            .unwrap_or(0);

        // Evaluate selection
        let results = evaluate_selection(ctx.viewer, selection)?;

        // Find first object with selected atoms
        let (src_name, sel_result) = results
            .into_iter()
            .find(|(_, sel)| sel.count() > 0)
            .ok_or_else(|| {
                CmdError::Selection(format!("No atoms matching '{}'", selection))
            })?;

        // Convert 1-indexed source_state to 0-indexed Option
        // source_state=0 means all states, N>0 means state N (1-indexed)
        let state_opt = if source_state > 0 {
            Some((source_state - 1) as usize)
        } else {
            None
        };

        // Get source molecule and its representations (read-only borrow)
        let (new_mol, source_reps) = {
            let mol_obj = ctx
                .viewer
                .objects()
                .get_molecule(&src_name)
                .ok_or_else(|| CmdError::ObjectNotFound(src_name.clone()))?;
            let src_mol = mol_obj.molecule();

            (extract_molecule(src_mol, &sel_result, name, state_opt), mol_obj.visible_reps())
        };

        let atom_count = new_mol.atom_count();
        let mut mol_obj = MoleculeObject::with_name(new_mol, name);
        mol_obj.set_visible_reps(source_reps);
        ctx.viewer.objects_mut().add(mol_obj);
        ctx.viewer.request_redraw();

        if !ctx.quiet {
            ctx.print(&format!(
                " Created \"{}\" from \"{}\" with {} atoms",
                name, selection, atom_count
            ));
        }

        Ok(())
    }
}

/// Extract a subset of atoms from a molecule based on a selection result.
///
/// If `source_state` is `None`, all states are copied. If `Some(i)` (0-indexed),
/// only that single state is copied into the new object.
fn extract_molecule(
    src: &ObjectMolecule,
    selection: &pymol_select::SelectionResult,
    name: &str,
    source_state: Option<usize>,
) -> ObjectMolecule {
    use std::collections::HashMap;

    let selected_indices: Vec<AtomIndex> = selection.indices().collect();
    let mut new_mol = ObjectMolecule::with_capacity(name, selected_indices.len(), 0);

    // Build old→new atom index mapping and copy atoms
    let mut index_map: HashMap<u32, AtomIndex> = HashMap::new();
    for &old_idx in &selected_indices {
        let atom = src.get_atom(old_idx).unwrap().clone();
        let new_idx = new_mol.add_atom(atom);
        index_map.insert(old_idx.0, new_idx);
    }

    // Copy bonds where both endpoints are selected
    for bond in src.bonds() {
        if let (Some(&new_a1), Some(&new_a2)) =
            (index_map.get(&bond.atom1.0), index_map.get(&bond.atom2.0))
        {
            let _ = new_mol.add_bond(new_a1, new_a2, bond.order);
        }
    }

    // Determine which states to copy
    let states: Vec<usize> = match source_state {
        Some(s) => vec![s],
        None => (0..src.state_count()).collect(),
    };

    // Copy coordinate sets (extract coords for selected atoms only)
    for state in states {
        if let Some(cs) = src.get_coord_set(state) {
            let mut coords = Vec::with_capacity(selected_indices.len() * 3);
            for &old_idx in &selected_indices {
                if let Some(v) = cs.get_atom_coord(old_idx) {
                    coords.push(v.x);
                    coords.push(v.y);
                    coords.push(v.z);
                } else {
                    coords.extend_from_slice(&[0.0, 0.0, 0.0]);
                }
            }
            let new_cs = pymol_mol::CoordSet::from_coords(coords);
            new_mol.add_coord_set(new_cs);
        }
    }

    // Classify atoms (protein, nucleic, etc.)
    new_mol.classify_atoms();

    new_mol
}

// ============================================================================
// copy command
// ============================================================================

struct CopyCommand;

impl Command for CopyCommand {
    fn name(&self) -> &str {
        "copy"
    }

    fn arg_hints(&self) -> &[ArgHint] {
        &[ArgHint::None, ArgHint::Selection]
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "copy" creates a new object from an object or selection.

USAGE

    copy target, source

ARGUMENTS

    target = string: name for the new object
    source = string: object name or selection expression

EXAMPLES

    copy protein_copy, protein
    copy active_site, byres around 5 ligand
    copy chainA, chain A
"#
    }

    fn execute<'v, 'r>(&self, ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>, args: &ParsedCommand) -> CmdResult {
        let target = args
            .get_str(0)
            .or_else(|| args.get_named_str("target"))
            .ok_or_else(|| CmdError::MissingArgument("target".to_string()))?;

        let source = args
            .get_str(1)
            .or_else(|| args.get_named_str("source"))
            .ok_or_else(|| CmdError::MissingArgument("source".to_string()))?;

        // Try as full object copy first (exact name match, no selection parsing)
        if let Some(mol_obj) = ctx.viewer.objects().get_molecule(source) {
            let cloned_mol = mol_obj.molecule().clone();
            let source_reps = mol_obj.visible_reps();
            let atom_count = cloned_mol.atom_count();

            let mut new_obj = MoleculeObject::with_name(cloned_mol, target);
            new_obj.set_visible_reps(source_reps);
            ctx.viewer.objects_mut().add(new_obj);
            ctx.viewer.request_redraw();

            if !ctx.quiet {
                ctx.print(&format!(
                    " Created \"{}\" as copy of \"{}\" ({} atoms)",
                    target, source, atom_count
                ));
            }
            return Ok(());
        }

        // Otherwise treat as selection expression (like create)
        let results = evaluate_selection(ctx.viewer, source)?;

        let (src_name, sel_result) = results
            .into_iter()
            .find(|(_, sel)| sel.count() > 0)
            .ok_or_else(|| {
                CmdError::Selection(format!("No atoms matching '{}'", source))
            })?;

        let (new_mol, source_reps) = {
            let mol_obj = ctx.viewer.objects().get_molecule(&src_name)
                .ok_or_else(|| CmdError::ObjectNotFound(src_name.clone()))?;
            (extract_molecule(mol_obj.molecule(), &sel_result, target, None), mol_obj.visible_reps())
        };

        let atom_count = new_mol.atom_count();
        let mut new_obj = MoleculeObject::with_name(new_mol, target);
        new_obj.set_visible_reps(source_reps);
        ctx.viewer.objects_mut().add(new_obj);
        ctx.viewer.request_redraw();

        if !ctx.quiet {
            ctx.print(&format!(
                " Created \"{}\" from \"{}\" ({} atoms)",
                target, source, atom_count
            ));
        }

        Ok(())
    }
}

// ============================================================================
// group command
// ============================================================================

struct GroupCommand;

impl Command for GroupCommand {
    fn name(&self) -> &str {
        "group"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "group" creates or modifies a group of objects.

USAGE

    group name [, members [, action ]]

ARGUMENTS

    name = string: group name
    members = string: object names to add (space or comma separated)
    action = add, remove, open, close, toggle, auto, ungroup, empty,
             purge, excise (default: add)

EXAMPLES

    group proteins, obj1 obj2 obj3
    group ligands, lig*, action=add
"#
    }

    fn execute<'v, 'r>(&self, ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>, args: &ParsedCommand) -> CmdResult {
        let name = args
            .get_str(0)
            .or_else(|| args.get_named_str("name"))
            .ok_or_else(|| CmdError::MissingArgument("name".to_string()))?;

        let _members = args
            .get_str(1)
            .or_else(|| args.get_named_str("members"))
            .unwrap_or("");

        let _action = args
            .get_str(2)
            .or_else(|| args.get_named_str("action"))
            .unwrap_or("add");

        // TODO: Implement group functionality
        // This would create/modify GroupObject instances

        if !ctx.quiet {
            ctx.print(&format!(" Group \"{}\" (not yet fully implemented)", name));
        }

        Ok(())
    }
}

// ============================================================================
// ungroup command
// ============================================================================

struct UngroupCommand;

impl Command for UngroupCommand {
    fn name(&self) -> &str {
        "ungroup"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "ungroup" removes objects from a group.

USAGE

    ungroup name

ARGUMENTS

    name = string: group name to ungroup

EXAMPLES

    ungroup proteins
"#
    }

    fn execute<'v, 'r>(&self, ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>, args: &ParsedCommand) -> CmdResult {
        let name = args
            .get_str(0)
            .or_else(|| args.get_named_str("name"))
            .ok_or_else(|| CmdError::MissingArgument("name".to_string()))?;

        // TODO: Implement ungroup functionality

        if !ctx.quiet {
            ctx.print(&format!(" Ungrouped \"{}\" (not yet implemented)", name));
        }

        Ok(())
    }
}

// ============================================================================
// state command
// ============================================================================

struct StateCommand;

impl Command for StateCommand {
    fn name(&self) -> &str {
        "state"
    }

    fn aliases(&self) -> &[&str] {
        &["frame"]
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "state" sets the displayed state for objects.

USAGE

    state N [, object]

ARGUMENTS

    N = integer: state number (1-indexed)
    object = string: object name (default: all objects)

EXAMPLES

    state 1
    state 5, 1nmr
"#
    }

    fn execute<'v, 'r>(&self, ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>, args: &ParsedCommand) -> CmdResult {
        let state_num = args
            .get_int(0)
            .or_else(|| args.get_named_int("N"))
            .ok_or_else(|| CmdError::MissingArgument("N".to_string()))?;

        if state_num < 1 {
            return Err(CmdError::invalid_arg("N", "State number must be >= 1"));
        }

        let state_idx = (state_num - 1) as usize;

        let object = args
            .get_str(1)
            .or_else(|| args.get_named_str("object"));

        if let Some(obj_name) = object {
            // Set state for a specific object
            let mol_obj = ctx.viewer.objects_mut().get_molecule_mut(obj_name)
                .ok_or_else(|| CmdError::ObjectNotFound(obj_name.to_string()))?;
            if !mol_obj.set_display_state(state_idx) {
                return Err(CmdError::execution(format!(
                    "State {} is out of range for \"{}\" ({} states)",
                    state_num, obj_name, mol_obj.molecule().state_count()
                )));
            }
        } else {
            // Set state for all objects
            let names: Vec<String> = ctx.viewer.objects().names().map(|s| s.to_string()).collect();
            for name in &names {
                if let Some(mol_obj) = ctx.viewer.objects_mut().get_molecule_mut(name) {
                    mol_obj.set_display_state(state_idx);
                }
            }
        }

        ctx.viewer.request_redraw();

        if !ctx.quiet {
            if let Some(obj_name) = object {
                ctx.print(&format!(" State {} for \"{}\"", state_num, obj_name));
            } else {
                ctx.print(&format!(" State {}", state_num));
            }
        }

        Ok(())
    }
}

// ============================================================================
// count_states command
// ============================================================================

struct CountStatesCommand;

impl Command for CountStatesCommand {
    fn name(&self) -> &str {
        "count_states"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "count_states" reports the number of states in an object.

USAGE

    count_states [selection]

ARGUMENTS

    selection = string: object or selection (default: all)

EXAMPLES

    count_states
    count_states 1nmr
"#
    }

    fn execute<'v, 'r>(&self, ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>, args: &ParsedCommand) -> CmdResult {
        let selection = args
            .get_str(0)
            .or_else(|| args.get_named_str("selection"));

        let mut max_states = 0usize;

        if let Some(sel) = selection {
            // Try as object name first
            if let Some(mol_obj) = ctx.viewer.objects().get_molecule(sel) {
                let count = mol_obj.molecule().state_count();
                max_states = max_states.max(count);
                if !ctx.quiet {
                    ctx.print(&format!(" \"{}\" has {} states.", sel, count));
                }
            } else {
                // Try as selection — report for all matching objects
                let results = evaluate_selection(ctx.viewer, sel)?;
                for (obj_name, sel_result) in &results {
                    if sel_result.count() > 0 {
                        if let Some(mol_obj) = ctx.viewer.objects().get_molecule(obj_name) {
                            let count = mol_obj.molecule().state_count();
                            max_states = max_states.max(count);
                            if !ctx.quiet {
                                ctx.print(&format!(" \"{}\" has {} states.", obj_name, count));
                            }
                        }
                    }
                }
            }
        } else {
            // Report for all objects
            let names: Vec<String> = ctx.viewer.objects().names().map(|s| s.to_string()).collect();
            for name in &names {
                if let Some(mol_obj) = ctx.viewer.objects().get_molecule(name) {
                    let count = mol_obj.molecule().state_count();
                    max_states = max_states.max(count);
                    if !ctx.quiet {
                        ctx.print(&format!(" \"{}\" has {} states.", name, count));
                    }
                }
            }
        }

        if !ctx.quiet {
            ctx.print(&format!(" cmd.count_states: {} states.", max_states));
        }

        Ok(())
    }
}

// ============================================================================
// split_states command
// ============================================================================

struct SplitStatesCommand;

impl Command for SplitStatesCommand {
    fn name(&self) -> &str {
        "split_states"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "split_states" creates individual objects for each state of a
    multi-state object.

USAGE

    split_states object [, first [, last [, prefix ]]]

ARGUMENTS

    object = string: source object name
    first = integer: first state (1-indexed, default: 1)
    last = integer: last state (1-indexed, default: last)
    prefix = string: name prefix (default: object_)

EXAMPLES

    split_states 1nmr
    split_states 1nmr, first=1, last=5
    split_states 1nmr, prefix=model_
"#
    }

    fn execute<'v, 'r>(&self, ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>, args: &ParsedCommand) -> CmdResult {
        let object = args
            .get_str(0)
            .or_else(|| args.get_named_str("object"))
            .ok_or_else(|| CmdError::MissingArgument("object".to_string()))?;

        let mol_obj = ctx.viewer.objects().get_molecule(object)
            .ok_or_else(|| CmdError::ObjectNotFound(object.to_string()))?;
        let n_states = mol_obj.molecule().state_count();

        if n_states == 0 {
            return Err(CmdError::execution(format!("\"{}\" has no states", object)));
        }

        let first = args
            .get_int(1)
            .or_else(|| args.get_named_int("first"))
            .map(|v| (v.max(1) - 1) as usize)
            .unwrap_or(0);

        let last = args
            .get_int(2)
            .or_else(|| args.get_named_int("last"))
            .map(|v| (v.max(1) as usize).min(n_states))
            .unwrap_or(n_states);

        let prefix = args
            .get_str(3)
            .or_else(|| args.get_named_str("prefix"))
            .map(|s| s.to_string())
            .unwrap_or_else(|| format!("{}_", object));

        // Collect new molecules first (borrow src as read-only)
        let new_mols: Vec<(String, ObjectMolecule)> = {
            let mol_obj = ctx.viewer.objects().get_molecule(object).unwrap();
            let src_mol = mol_obj.molecule();

            // Build an "all atoms" selection
            let all_sel = pymol_select::SelectionResult::all(src_mol.atom_count());

            (first..last)
                .map(|state_idx| {
                    let name = format!("{}{:04}", prefix, state_idx + 1);
                    let mol = extract_molecule(src_mol, &all_sel, &name, Some(state_idx));
                    (name, mol)
                })
                .collect()
        };

        let count = new_mols.len();
        for (name, mol) in new_mols {
            let mol_obj = MoleculeObject::with_name(mol, &name);
            ctx.viewer.objects_mut().add(mol_obj);
        }

        ctx.viewer.request_redraw();

        if !ctx.quiet {
            ctx.print(&format!(
                " Split \"{}\" into {} objects ({}{:04} to {}{:04})",
                object, count, prefix, first + 1, prefix, last
            ));
        }

        Ok(())
    }
}

// ============================================================================
// extract command
// ============================================================================

struct ExtractCommand;

impl Command for ExtractCommand {
    fn name(&self) -> &str {
        "extract"
    }

    fn arg_hints(&self) -> &[ArgHint] {
        &[ArgHint::None, ArgHint::Selection]
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "extract" creates a new object from a selection and removes
    those atoms from the source object.

USAGE

    extract name, selection [, source_state]

ARGUMENTS

    name = string: name for the new object
    selection = string: selection of atoms to extract
    source_state = integer: state to copy from (default: 0 = all)

EXAMPLES

    extract ligand, organic
    extract chainA, chain A
"#
    }

    fn execute<'v, 'r>(&self, ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>, args: &ParsedCommand) -> CmdResult {
        let name = args
            .get_str(0)
            .or_else(|| args.get_named_str("name"))
            .ok_or_else(|| CmdError::MissingArgument("name".to_string()))?;

        let selection = args
            .get_str(1)
            .or_else(|| args.get_named_str("selection"))
            .ok_or_else(|| CmdError::MissingArgument("selection".to_string()))?;

        let source_state = args
            .get_int(2)
            .or_else(|| args.get_named_int("source_state"))
            .unwrap_or(0);

        // Convert 1-indexed source_state to 0-indexed Option
        let state_opt = if source_state > 0 {
            Some((source_state - 1) as usize)
        } else {
            None
        };

        // Evaluate selection
        let results = evaluate_selection(ctx.viewer, selection)?;

        // Find first object with selected atoms
        let (src_name, sel_result) = results
            .into_iter()
            .find(|(_, sel)| sel.count() > 0)
            .ok_or_else(|| {
                CmdError::Selection(format!("No atoms matching '{}'", selection))
            })?;

        // Extract molecule from source (read-only borrow)
        let (new_mol, remove_indices) = {
            let mol_obj = ctx.viewer.objects().get_molecule(&src_name)
                .ok_or_else(|| CmdError::ObjectNotFound(src_name.clone()))?;
            let src_mol = mol_obj.molecule();

            let new_mol = extract_molecule(src_mol, &sel_result, name, state_opt);
            let remove_indices: Vec<AtomIndex> = sel_result.indices().collect();
            (new_mol, remove_indices)
        };

        let atom_count = new_mol.atom_count();

        // Copy visible representations from source object
        let source_reps = ctx.viewer.objects().get_molecule(&src_name)
            .map(|obj| obj.visible_reps())
            .unwrap_or_default();

        // Add the new object with same representations as source
        let mut mol_obj = MoleculeObject::with_name(new_mol, name);
        mol_obj.set_visible_reps(source_reps);
        ctx.viewer.objects_mut().add(mol_obj);

        // Remove atoms from source
        if let Some(src_obj) = ctx.viewer.objects_mut().get_molecule_mut(&src_name) {
            src_obj.molecule_mut().remove_atoms(&remove_indices);
            src_obj.invalidate(DirtyFlags::ALL);
        }

        ctx.viewer.request_redraw();

        if !ctx.quiet {
            ctx.print(&format!(
                " Extracted {} atoms from \"{}\" into \"{}\"",
                atom_count, src_name, name
            ));
        }

        Ok(())
    }
}
