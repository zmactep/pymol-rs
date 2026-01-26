//! Object commands: delete, rename, create, copy, group, ungroup

use crate::args::ParsedCommand;
use crate::command::{Command, CommandContext, CommandRegistry};
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

    fn execute(&self, ctx: &mut CommandContext, args: &ParsedCommand) -> CmdResult {
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

    "set_name" renames an object.

USAGE

    set_name old_name, new_name

ARGUMENTS

    old_name = string: current object name
    new_name = string: new object name

EXAMPLES

    set_name protein, myprotein
    rename ligand, drug
"#
    }

    fn execute(&self, ctx: &mut CommandContext, args: &ParsedCommand) -> CmdResult {
        let old_name = args
            .get_str(0)
            .or_else(|| args.get_named_str("old_name"))
            .ok_or_else(|| CmdError::MissingArgument("old_name".to_string()))?;

        let new_name = args
            .get_str(1)
            .or_else(|| args.get_named_str("new_name"))
            .ok_or_else(|| CmdError::MissingArgument("new_name".to_string()))?;

        ctx.viewer
            .objects_mut()
            .rename(old_name, new_name)
            .map_err(|e| CmdError::Scene(e.to_string()))?;

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

    fn execute(&self, ctx: &mut CommandContext, args: &ParsedCommand) -> CmdResult {
        let name = args
            .get_str(0)
            .or_else(|| args.get_named_str("name"))
            .ok_or_else(|| CmdError::MissingArgument("name".to_string()))?;

        let selection = args
            .get_str(1)
            .or_else(|| args.get_named_str("selection"))
            .ok_or_else(|| CmdError::MissingArgument("selection".to_string()))?;

        // TODO: Implement create with selection support
        // This would:
        // 1. Evaluate the selection
        // 2. Extract atoms matching the selection
        // 3. Create a new ObjectMolecule with those atoms
        // 4. Add to the registry

        if !ctx.quiet {
            ctx.print(&format!(
                " Creating \"{}\" from \"{}\" (not yet implemented)",
                name, selection
            ));
        }

        Ok(())
    }
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

    fn execute(&self, ctx: &mut CommandContext, args: &ParsedCommand) -> CmdResult {
        let target = args
            .get_str(0)
            .or_else(|| args.get_named_str("target"))
            .ok_or_else(|| CmdError::MissingArgument("target".to_string()))?;

        let source = args
            .get_str(1)
            .or_else(|| args.get_named_str("source"))
            .ok_or_else(|| CmdError::MissingArgument("source".to_string()))?;

        // Check that source exists
        if !ctx.viewer.objects().contains(source) {
            return Err(CmdError::ObjectNotFound(source.to_string()));
        }

        // TODO: Implement proper object cloning
        // ObjectMolecule doesn't implement Clone, so we'd need to
        // serialize/deserialize or implement a deep copy method

        if !ctx.quiet {
            ctx.print(&format!(
                " Copy from \"{}\" to \"{}\" (not yet fully implemented)",
                source, target
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

    fn execute(&self, ctx: &mut CommandContext, args: &ParsedCommand) -> CmdResult {
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

    fn execute(&self, ctx: &mut CommandContext, args: &ParsedCommand) -> CmdResult {
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
