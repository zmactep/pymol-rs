//! Control commands: quit, reinitialize, refresh, rebuild

use pymol_scene::DirtyFlags;

use crate::args::ParsedCommand;
use crate::command::{Command, CommandContext, CommandRegistry, ViewerLike};
use crate::error::{CmdError, CmdResult};

/// Register control commands
pub fn register(registry: &mut CommandRegistry) {
    registry.register(QuitCommand);
    registry.register(ReinitializeCommand);
    registry.register(RefreshCommand);
    registry.register(RebuildCommand);
    registry.register(HelpCommand);
}

// ============================================================================
// quit command
// ============================================================================

struct QuitCommand;

impl Command for QuitCommand {
    fn name(&self) -> &str {
        "quit"
    }

    fn aliases(&self) -> &[&str] {
        &["exit"]
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "quit" exits PyMOL.

USAGE

    quit [ code ]

ARGUMENTS

    code = integer: exit code (default: 0)

EXAMPLES

    quit
    quit 1
"#
    }

    fn execute<'v, 'r>(&self, _ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>, args: &ParsedCommand) -> CmdResult {
        let _code = args.get_int(0).or_else(|| args.get_named_int("code")).unwrap_or(0);

        // Signal that we want to quit
        // The actual quit will be handled by the event loop
        Err(CmdError::Aborted)
    }
}

// ============================================================================
// reinitialize command
// ============================================================================

struct ReinitializeCommand;

impl Command for ReinitializeCommand {
    fn name(&self) -> &str {
        "reinitialize"
    }

    fn aliases(&self) -> &[&str] {
        &["reinit"]
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "reinitialize" resets PyMOL to its initial state.

USAGE

    reinitialize [ what [, object ]]

ARGUMENTS

    what = everything, settings, objects, purge_defaults (default: everything)
    object = string: specific object to reinitialize

EXAMPLES

    reinitialize
    reinitialize settings
    reinitialize objects
"#
    }

    fn execute<'v, 'r>(&self, ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>, args: &ParsedCommand) -> CmdResult {
        let what = args
            .get_str(0)
            .or_else(|| args.get_named_str("what"))
            .unwrap_or("everything");

        match what.to_lowercase().as_str() {
            "everything" | "all" => {
                ctx.viewer.objects_mut().clear();
                // TODO: Reset settings to defaults
                ctx.viewer.reset_view();
            }
            "settings" => {
                // TODO: Reset settings to defaults
            }
            "objects" => {
                ctx.viewer.objects_mut().clear();
            }
            _ => {
                return Err(CmdError::invalid_arg(
                    "what",
                    format!("unknown reinitialize target: {}", what),
                ));
            }
        }

        ctx.viewer.request_redraw();

        if !ctx.quiet {
            ctx.print(&format!(" Reinitialized {}", what));
        }

        Ok(())
    }
}

// ============================================================================
// refresh command
// ============================================================================

struct RefreshCommand;

impl Command for RefreshCommand {
    fn name(&self) -> &str {
        "refresh"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "refresh" forces a redraw of the scene.

USAGE

    refresh

EXAMPLES

    refresh
"#
    }

    fn execute<'v, 'r>(&self, ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>, _args: &ParsedCommand) -> CmdResult {
        ctx.viewer.request_redraw();

        if !ctx.quiet {
            ctx.print(" Refreshed");
        }

        Ok(())
    }
}

// ============================================================================
// rebuild command
// ============================================================================

struct RebuildCommand;

impl Command for RebuildCommand {
    fn name(&self) -> &str {
        "rebuild"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "rebuild" forces rebuilding of all representations.

USAGE

    rebuild [ selection ]

ARGUMENTS

    selection = string: selection to rebuild (default: all)

EXAMPLES

    rebuild
    rebuild protein
"#
    }

    fn execute<'v, 'r>(&self, ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>, args: &ParsedCommand) -> CmdResult {
        let selection = args
            .get_str(0)
            .or_else(|| args.get_named_str("selection"))
            .unwrap_or("all");

        // Get objects to rebuild
        let object_names: Vec<String> = if selection == "all" || selection == "*" {
            ctx.viewer.objects().names().map(|s| s.to_string()).collect()
        } else {
            ctx.viewer
                .objects()
                .matching(selection)
                .iter()
                .map(|s| s.to_string())
                .collect()
        };

        // Mark all molecule objects as dirty to force rebuild
        for name in &object_names {
            if let Some(mol) = ctx.viewer.objects_mut().get_molecule_mut(name) {
                mol.invalidate(DirtyFlags::ALL);
            }
        }

        ctx.viewer.request_redraw();

        if !ctx.quiet {
            ctx.print(&format!(" Rebuilt \"{}\"", selection));
        }

        Ok(())
    }
}

// ============================================================================
// help command
// ============================================================================

struct HelpCommand;

impl Command for HelpCommand {
    fn name(&self) -> &str {
        "help"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "help" displays help for commands.

USAGE

    help [ command ]

ARGUMENTS

    command = string: command to get help for (default: list commands)

EXAMPLES

    help
    help load
    help zoom
"#
    }

    fn execute<'v, 'r>(&self, ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>, args: &ParsedCommand) -> CmdResult {
        let command = args.get_str(0).or_else(|| args.get_named_str("command"));

        if let Some(cmd_name) = command {
            // Show help for specific command
            if let Some(registry) = ctx.registry() {
                if let Some(cmd) = registry.get(cmd_name) {
                    ctx.print(cmd.help());
                } else {
                    ctx.print(&format!(" Unknown command: '{}'", cmd_name));
                }
            } else {
                ctx.print(&format!(" Help for '{}' - registry not available", cmd_name));
            }
        } else {
            // List available commands
            ctx.print(" Available commands:");
            ctx.print("   File I/O: load, save, cd, pwd, ls");
            ctx.print("   Viewing:  zoom, center, orient, reset, clip");
            ctx.print("   Display:  show, hide, as, enable, disable, color, bg_color");
            ctx.print("   Objects:  delete, rename, create, copy, group");
            ctx.print("   Settings: set, get, unset");
            ctx.print("   Control:  quit, reinitialize, refresh, rebuild, help");
            ctx.print("");
            ctx.print(" Type 'help <command>' for detailed help");
        }

        Ok(())
    }
}
