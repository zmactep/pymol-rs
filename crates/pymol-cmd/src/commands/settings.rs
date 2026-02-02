//! Settings commands: set, get, unset

use crate::args::ParsedCommand;
use crate::command::{Command, CommandContext, CommandRegistry, ViewerLike};
use crate::error::{CmdError, CmdResult};
use pymol_settings::{get_setting, get_setting_id, id as setting_id, SettingType, SettingValue};

/// Register settings commands
pub fn register(registry: &mut CommandRegistry) {
    registry.register(SetCommand);
    registry.register(GetCommand);
    registry.register(UnsetCommand);
}

/// Parse a string value into a SettingValue based on the setting type
fn parse_setting_value(s: &str, setting_type: SettingType) -> CmdResult<SettingValue> {
    match setting_type {
        SettingType::Bool => {
            match s.to_lowercase().as_str() {
                "1" | "on" | "true" | "yes" => Ok(SettingValue::Bool(true)),
                "0" | "off" | "false" | "no" => Ok(SettingValue::Bool(false)),
                _ => Err(CmdError::invalid_arg(
                    "value",
                    format!("Invalid boolean value '{}'. Use 1/0, on/off, true/false, or yes/no", s)
                ))
            }
        }
        SettingType::Int => {
            // Try numeric parse first
            if let Ok(v) = s.parse::<i32>() {
                return Ok(SettingValue::Int(v));
            }
            // Also accept on/off/true/false for integer settings (PyMOL compatibility)
            match s.to_lowercase().as_str() {
                "on" | "true" | "yes" => Ok(SettingValue::Int(1)),
                "off" | "false" | "no" => Ok(SettingValue::Int(0)),
                _ => Err(CmdError::invalid_arg(
                    "value",
                    format!("Invalid integer value '{}'. Expected a number or on/off", s)
                ))
            }
        }
        SettingType::Float => {
            // Try numeric parse first
            if let Ok(v) = s.parse::<f32>() {
                return Ok(SettingValue::Float(v));
            }
            // Also accept on/off for float settings (PyMOL compatibility)
            match s.to_lowercase().as_str() {
                "on" | "true" | "yes" => Ok(SettingValue::Float(1.0)),
                "off" | "false" | "no" => Ok(SettingValue::Float(0.0)),
                _ => Err(CmdError::invalid_arg(
                    "value",
                    format!("Invalid float value '{}'. Expected a number or on/off", s)
                ))
            }
        }
        SettingType::Float3 => {
            // Parse [x, y, z] or x,y,z format
            let cleaned = s.trim_start_matches('[').trim_end_matches(']');
            let parts: Vec<&str> = cleaned.split(',').map(|p| p.trim()).collect();
            if parts.len() != 3 {
                return Err(CmdError::invalid_arg(
                    "value",
                    format!("Invalid float3 value '{}'. Expected [x, y, z] or x,y,z format", s)
                ));
            }
            let x = parts[0].parse::<f32>().map_err(|_| CmdError::invalid_arg("value", "Invalid x component"))?;
            let y = parts[1].parse::<f32>().map_err(|_| CmdError::invalid_arg("value", "Invalid y component"))?;
            let z = parts[2].parse::<f32>().map_err(|_| CmdError::invalid_arg("value", "Invalid z component"))?;
            Ok(SettingValue::Float3([x, y, z]))
        }
        SettingType::Color => {
            // Color can be an integer index or color name (we'll handle names later)
            s.parse::<i32>()
                .map(SettingValue::Color)
                .map_err(|_| CmdError::invalid_arg(
                    "value",
                    format!("Invalid color value '{}'. Expected a color index (integer)", s)
                ))
        }
        SettingType::String => {
            Ok(SettingValue::String(s.to_string()))
        }
        SettingType::Blank => {
            Err(CmdError::invalid_arg("value", "Cannot set a blank/unused setting"))
        }
    }
}

/// Format a SettingValue for display
fn format_setting_value(value: &SettingValue) -> String {
    match value {
        SettingValue::Bool(v) => if *v { "on".to_string() } else { "off".to_string() },
        SettingValue::Int(v) => v.to_string(),
        SettingValue::Float(v) => format!("{:.6}", v),
        SettingValue::Float3(v) => format!("[{:.3}, {:.3}, {:.3}]", v[0], v[1], v[2]),
        SettingValue::Color(v) => v.to_string(),
        SettingValue::String(v) => v.clone(),
    }
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

    fn execute<'v, 'r>(&self, ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>, args: &ParsedCommand) -> CmdResult {
        let name = args
            .get_str(0)
            .or_else(|| args.get_named_str("name"))
            .ok_or_else(|| CmdError::MissingArgument("name".to_string()))?;

        let value_str = args
            .get_arg(1)
            .and_then(|v| v.to_string_repr())
            .or_else(|| args.get_named("value").and_then(|v| v.to_string_repr()));

        let _selection = args
            .get_str(2)
            .or_else(|| args.get_named_str("selection"));

        // Look up the setting by name
        let id = get_setting_id(name)
            .ok_or_else(|| CmdError::invalid_arg("name", format!("Unknown setting: {}", name)))?;

        let setting = get_setting(id)
            .ok_or_else(|| CmdError::invalid_arg("name", "Invalid setting ID"))?;

        // If no value provided, toggle boolean settings or show current value
        let value_str = match value_str {
            Some(v) => v,
            None => {
                // For boolean settings, toggle
                if setting.setting_type == SettingType::Bool {
                    let current = ctx.viewer.settings().get_bool(id);
                    let new_value = !current;
                    ctx.viewer.settings_mut().set_bool(id, new_value)
                        .map_err(|e| CmdError::execution(e.to_string()))?;
                    ctx.viewer.request_redraw();
                    if !ctx.quiet {
                        ctx.print(&format!(" {} = {}", name, if new_value { "on" } else { "off" }));
                    }
                    return Ok(());
                } else {
                    return Err(CmdError::MissingArgument("value".to_string()));
                }
            }
        };

        // Parse the value based on the setting type
        let value = parse_setting_value(&value_str, setting.setting_type)?;

        // Set the value
        ctx.viewer.settings_mut().set(id, value.clone())
            .map_err(|e| CmdError::execution(e.to_string()))?;

        // Special handling for surface_quality - propagate to all molecules
        if id == setting_id::surface_quality {
            if let SettingValue::Int(quality) = &value {
                let names: Vec<_> = ctx.viewer.objects().names().map(|s| s.to_string()).collect();
                for name in names {
                    if let Some(mol) = ctx.viewer.objects_mut().get_molecule_mut(&name) {
                        mol.set_surface_quality(*quality);
                    }
                }
            }
        }

        ctx.viewer.request_redraw();

        if !ctx.quiet {
            ctx.print(&format!(" {} = {}", name, format_setting_value(&value)));
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

    fn execute<'v, 'r>(&self, ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>, args: &ParsedCommand) -> CmdResult {
        let name = args
            .get_str(0)
            .or_else(|| args.get_named_str("name"))
            .ok_or_else(|| CmdError::MissingArgument("name".to_string()))?;

        // Look up the setting by name
        let id = get_setting_id(name)
            .ok_or_else(|| CmdError::invalid_arg("name", format!("Unknown setting: {}", name)))?;

        let setting = get_setting(id)
            .ok_or_else(|| CmdError::invalid_arg("name", "Invalid setting ID"))?;

        // Get the current value
        let value = ctx.viewer.settings().get(id)
            .ok_or_else(|| CmdError::execution("Failed to get setting value"))?;

        // Display the value with type information
        ctx.print(&format!(" {} ({}) = {}", name, setting.setting_type, format_setting_value(&value)));

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
    unset ray_trace_mode
"#
    }

    fn execute<'v, 'r>(&self, ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>, args: &ParsedCommand) -> CmdResult {
        let name = args
            .get_str(0)
            .or_else(|| args.get_named_str("name"))
            .ok_or_else(|| CmdError::MissingArgument("name".to_string()))?;

        // Look up the setting by name
        let id = get_setting_id(name)
            .ok_or_else(|| CmdError::invalid_arg("name", format!("Unknown setting: {}", name)))?;

        let setting = get_setting(id)
            .ok_or_else(|| CmdError::invalid_arg("name", "Invalid setting ID"))?;

        // Reset to default
        ctx.viewer.settings_mut().reset(id)
            .map_err(|e| CmdError::execution(e.to_string()))?;

        ctx.viewer.request_redraw();

        if !ctx.quiet {
            ctx.print(&format!(" {} reset to default: {}", name, format_setting_value(&setting.default)));
        }

        Ok(())
    }
}
