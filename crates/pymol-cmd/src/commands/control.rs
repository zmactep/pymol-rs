//! Control commands: quit, reinitialize, refresh, rebuild, run

use pymol_scene::{DirtyFlags, Session};

use crate::args::ParsedCommand;
use crate::command::{ArgHint, Command, CommandContext, CommandRegistry, ViewerLike};
use crate::error::{CmdError, CmdResult};
use crate::executor::CommandExecutor;
use crate::script::ScriptEngine;

use super::io::expand_path;

/// Register control commands
pub fn register(registry: &mut CommandRegistry) {
    registry.register(QuitCommand);
    registry.register(ReinitializeCommand);
    registry.register(RefreshCommand);
    registry.register(RebuildCommand);
    registry.register(HelpCommand);
    registry.register(RunCommand);
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
                ctx.viewer.replace_session(Session::new());
            }
            "settings" => {
                *ctx.viewer.settings_mut() = pymol_settings::Settings::default();
                ctx.viewer.session_mut().apply_default_settings();
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

    fn arg_hints(&self) -> &[ArgHint] {
        &[ArgHint::Command]
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
            let cmd = ctx.registry().and_then(|r| r.get(cmd_name));
            match cmd {
                Some(cmd) => ctx.print(cmd.help()),
                None => ctx.print(&format!(" Unknown command: '{}'", cmd_name)),
            }
        } else {
            // List available commands from registry
            let lines: Vec<String> = ctx.registry().map(|registry| {
                let mut names: Vec<&str> = registry.names().collect();
                names.sort();
                names.chunks(5)
                    .map(|chunk| format!("   {}", chunk.join(", ")))
                    .collect()
            }).unwrap_or_default();
            if !lines.is_empty() {
                ctx.print(" Available commands:");
                for line in &lines {
                    ctx.print(line);
                }
            }
            ctx.print("");
            ctx.print(" Type 'help <command>' for detailed help");
        }

        Ok(())
    }
}

// ============================================================================
// run command
// ============================================================================

struct RunCommand;

impl Command for RunCommand {
    fn name(&self) -> &str {
        "run"
    }

    fn aliases(&self) -> &[&str] {
        &["@"]
    }

    fn arg_hints(&self) -> &[ArgHint] {
        &[ArgHint::Path]
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "run" executes a script file. Supports .pml (PyMOL script) and any
    file type registered by a plugin (e.g., .py via the Python plugin).

USAGE

    run filename

ARGUMENTS

    filename = string: path to script file

EXAMPLES

    run setup.pml
    run ~/scripts/analysis.pml
    run script.py
    @ script.pml
"#
    }

    fn execute<'v, 'r>(&self, ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>, args: &ParsedCommand) -> CmdResult {
        let filename = args
            .get_str(0)
            .or_else(|| args.get_named_str("filename"))
            .ok_or_else(|| CmdError::MissingArgument("filename".to_string()))?;

        let path = expand_path(filename);
        let ext = path
            .extension()
            .and_then(|e| e.to_str())
            .unwrap_or("pml");

        match ext {
            "pml" => {
                let executor = if let Some(registry) = ctx.registry() {
                    CommandExecutor::with_registry(
                        registry.clone(),
                        ctx.script_handlers_map().cloned().unwrap_or_default(),
                        ctx.format_handlers_map().cloned().unwrap_or_default(),
                    )
                } else {
                    CommandExecutor::new()
                };
                let mut engine = ScriptEngine::with_executor(executor);
                engine.run_pml(ctx.viewer, &path)
            }
            other => {
                if let Some(handler) = ctx.script_handler(other) {
                    let path_str = path.to_str().ok_or_else(|| {
                        CmdError::Execution("invalid path encoding".to_string())
                    })?;
                    handler(path_str).map_err(CmdError::Execution)
                } else {
                    Err(CmdError::Execution(format!(
                        "no handler for .{} files. Install a plugin that handles this format.",
                        other
                    )))
                }
            }
        }
    }
}
