//! Settings commands: set, get, unset

use crate::args::ParsedCommand;
use crate::command::{Command, CommandContext, CommandRegistry};
use crate::error::{CmdError, CmdResult};

/// Register settings commands
pub fn register(registry: &mut CommandRegistry) {
    registry.register(SetCommand);
    registry.register(GetCommand);
    registry.register(UnsetCommand);
}

// ============================================================================
// set command
// ============================================================================

struct SetCommand;

impl Command for SetCommand {
    fn name(&self) -> &str {
        "set"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "set" changes a setting value.

USAGE

    set name [, value [, selection [, state ]]]

ARGUMENTS

    name = string: setting name
    value = string: new value (depends on setting type)
    selection = string: apply to specific selection (default: global)
    state = integer: state for state-specific settings (default: 0)

EXAMPLES

    set sphere_scale, 0.5
    set cartoon_color, red, chain A
    set bg_rgb, [1.0, 1.0, 1.0]
"#
    }

    fn execute(&self, ctx: &mut CommandContext, args: &ParsedCommand) -> CmdResult {
        let name = args
            .get_str(0)
            .or_else(|| args.get_named_str("name"))
            .ok_or_else(|| CmdError::MissingArgument("name".to_string()))?;

        let value = args
            .get_str(1)
            .or_else(|| args.get_named_str("value"));

        let _selection = args
            .get_str(2)
            .or_else(|| args.get_named_str("selection"));

        // Try to set the value
        // TODO: Implement setting name to ID lookup
        // For now, we just log the attempt
        if let Some(val) = value {
            // This is a stub - full implementation would look up the setting ID by name
            // and call the appropriate set_* method
            log::debug!("Setting {} = {} (lookup not yet implemented)", name, val);
        }

        ctx.viewer.request_redraw();

        if !ctx.quiet {
            if let Some(val) = value {
                ctx.print(&format!(" Setting {} to {}", name, val));
            } else {
                ctx.print(&format!(" Toggling {}", name));
            }
        }

        Ok(())
    }
}

// ============================================================================
// get command
// ============================================================================

struct GetCommand;

impl Command for GetCommand {
    fn name(&self) -> &str {
        "get"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "get" displays the current value of a setting.

USAGE

    get name [, selection [, state ]]

ARGUMENTS

    name = string: setting name
    selection = string: get from specific selection (default: global)
    state = integer: state for state-specific settings (default: 0)

EXAMPLES

    get sphere_scale
    get bg_rgb
"#
    }

    fn execute(&self, ctx: &mut CommandContext, args: &ParsedCommand) -> CmdResult {
        let name = args
            .get_str(0)
            .or_else(|| args.get_named_str("name"))
            .ok_or_else(|| CmdError::MissingArgument("name".to_string()))?;

        // TODO: Implement setting name to ID lookup
        // For now, we just report that the setting lookup is not implemented
        ctx.print(&format!(" {} = (lookup not yet implemented)", name));

        Ok(())
    }
}

// ============================================================================
// unset command
// ============================================================================

struct UnsetCommand;

impl Command for UnsetCommand {
    fn name(&self) -> &str {
        "unset"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "unset" restores a setting to its default value.

USAGE

    unset name [, selection [, state ]]

ARGUMENTS

    name = string: setting name
    selection = string: apply to specific selection (default: global)
    state = integer: state for state-specific settings (default: 0)

EXAMPLES

    unset sphere_scale
"#
    }

    fn execute(&self, ctx: &mut CommandContext, args: &ParsedCommand) -> CmdResult {
        let name = args
            .get_str(0)
            .or_else(|| args.get_named_str("name"))
            .ok_or_else(|| CmdError::MissingArgument("name".to_string()))?;

        // TODO: Implement unset by restoring default value

        ctx.viewer.request_redraw();

        if !ctx.quiet {
            ctx.print(&format!(" Unset {} (restored to default)", name));
        }

        Ok(())
    }
}
