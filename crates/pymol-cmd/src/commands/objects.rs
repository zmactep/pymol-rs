//! Object commands: delete, rename, create, copy, group, ungroup

use pymol_mol::{AtomIndex, ObjectMolecule};
use pymol_scene::MoleculeObject;

use crate::args::ParsedCommand;
use crate::command::{Command, CommandContext, CommandRegistry, ViewerLike};
use crate::commands::selecting::evaluate_selection;
use crate::error::{CmdError, CmdResult};

/// Register object commands
pub fn register(registry: &mut CommandRegistry) {
    registry.register(DeleteCommand);
    registry.register(RenameCommand);
    registry.register(CreateCommand);
    registry.register(CopyCommand);
    registry.register(GroupCommand);
    registry.register(UngroupCommand);
}

// ============================================================================
// delete command
// ============================================================================

struct DeleteCommand;

impl Command for DeleteCommand {
    fn name(&self) -> &str {
        "delete"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "delete" removes objects or selections.

USAGE

    delete name

ARGUMENTS

    name = string: object or selection name pattern

EXAMPLES

    delete protein
    delete obj*
    delete all
"#
    }

    fn execute<'v, 'r>(&self, ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>, args: &ParsedCommand) -> CmdResult {
        let name = args
            .get_str(0)
            .or_else(|| args.get_named_str("name"))
            .ok_or_else(|| CmdError::MissingArgument("name".to_string()))?;

        // Get matching objects
        let matches: Vec<String> = ctx
            .viewer
            .objects()
            .matching(name)
            .iter()
            .map(|s| s.to_string())
            .collect();

        if matches.is_empty() {
            return Err(CmdError::ObjectNotFound(name.to_string()));
        }

        let count = matches.len();
        for obj_name in matches {
            ctx.viewer.objects_mut().remove(&obj_name);
        }

        ctx.viewer.request_redraw();

        if !ctx.quiet {
            if count == 1 {
                ctx.print(&format!(" Deleted \"{}\"", name));
            } else {
                ctx.print(&format!(" Deleted {} objects matching \"{}\"", count, name));
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

        let _source_state = args
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

        // Get source molecule (read-only borrow)
        let new_mol = {
            let mol_obj = ctx
                .viewer
                .objects()
                .get_molecule(&src_name)
                .ok_or_else(|| CmdError::ObjectNotFound(src_name.clone()))?;
            let src_mol = mol_obj.molecule();

            extract_molecule(src_mol, &sel_result, name)
        };

        let atom_count = new_mol.atom_count();
        let mol_obj = MoleculeObject::with_name(new_mol, name);
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
fn extract_molecule(
    src: &ObjectMolecule,
    selection: &pymol_select::SelectionResult,
    name: &str,
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

    // Copy coordinate sets (extract coords for selected atoms only)
    for state in 0..src.state_count() {
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

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "copy" creates a copy of an object.

USAGE

    copy target, source

ARGUMENTS

    target = string: name for the new object
    source = string: object to copy

EXAMPLES

    copy protein_copy, protein
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

        // Clone the source molecule
        let cloned_mol = {
            let mol_obj = ctx
                .viewer
                .objects()
                .get_molecule(source)
                .ok_or_else(|| CmdError::ObjectNotFound(source.to_string()))?;
            mol_obj.molecule().clone()
        };

        let atom_count = cloned_mol.atom_count();
        let mol_obj = MoleculeObject::with_name(cloned_mol, target);
        ctx.viewer.objects_mut().add(mol_obj);
        ctx.viewer.request_redraw();

        if !ctx.quiet {
            ctx.print(&format!(
                " Created \"{}\" as copy of \"{}\" ({} atoms)",
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
