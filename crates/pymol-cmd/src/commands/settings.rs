//! Settings commands: set, get, unset, dss

use crate::args::ParsedCommand;
use crate::command::{Command, CommandContext, CommandRegistry, ViewerLike};
use crate::commands::selecting::evaluate_selection;
use crate::error::{CmdError, CmdResult};
use pymol_mol::dss::{assign_secondary_structure, DssSettings};
use pymol_scene::DirtyFlags;
use pymol_select::AtomIndex;
use pymol_settings::{get_setting, get_setting_id, id as setting_id, SettingType, SettingValue};

/// Register settings commands
pub fn register(registry: &mut CommandRegistry) {
    registry.register(SetCommand);
    registry.register(GetCommand);
    registry.register(UnsetCommand);
    registry.register(DssCommand);
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

        let selection = args
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

        // Check if this is a representation color setting
        let is_rep_color_setting = matches!(id, 
            setting_id::stick_color |
            setting_id::line_color |
            setting_id::cartoon_color |
            setting_id::surface_color |
            setting_id::mesh_color |
            setting_id::sphere_color |
            setting_id::ribbon_color
        );

        // For color settings, try to resolve color names
        let value = if setting.setting_type == SettingType::Color {
            // Try to parse as integer first
            if let Ok(v) = value_str.parse::<i32>() {
                SettingValue::Color(v)
            } else {
                // Check for special color schemes first (these use negative indices)
                // -1: by element (atomic/CPK)
                // -2: by chain
                // -3: by secondary structure
                // -4: by B-factor
                let color_index = match value_str.to_lowercase().as_str() {
                    "atomic" | "cpk" | "element" | "by_element" => -1,
                    "chain" | "by_chain" | "chainbow" => -2,
                    "ss" | "secondary_structure" | "by_ss" | "dssp" => -3,
                    "b" | "b_factor" | "bfactor" | "by_b" => -4,
                    _ => {
                        // Try to resolve as named color
                        ctx.viewer
                            .color_index(&value_str)
                            .map(|idx| idx as i32)
                            .ok_or_else(|| CmdError::invalid_arg("value", format!("Unknown color: {}", value_str)))?
                    }
                };
                SettingValue::Color(color_index)
            }
        } else {
            parse_setting_value(&value_str, setting.setting_type)?
        };

        // Handle representation color settings - default to "all" if no selection provided
        if is_rep_color_setting {
            let selection_str = selection.unwrap_or("all");
            let color_index = match &value {
                SettingValue::Color(idx) => *idx,
                _ => return Err(CmdError::execution("Expected color value")),
            };

            // Evaluate the selection
            let selection_results = evaluate_selection(ctx.viewer, selection_str)?;
            let mut total_affected = 0usize;

            // Apply color to atoms matching the selection
            for (obj_name, selected) in selection_results {
                let count = selected.count();
                if count > 0 {
                    if let Some(mol_obj) = ctx.viewer.objects_mut().get_molecule_mut(&obj_name) {
                        let mol_mut = mol_obj.molecule_mut();
                        for idx in selected.indices() {
                            if let Some(atom) = mol_mut.get_atom_mut(AtomIndex(idx.0)) {
                                match id {
                                    setting_id::stick_color => atom.repr.colors.stick = color_index,
                                    setting_id::line_color => atom.repr.colors.line = color_index,
                                    setting_id::cartoon_color => atom.repr.colors.cartoon = color_index,
                                    setting_id::surface_color => atom.repr.colors.surface = color_index,
                                    setting_id::mesh_color => atom.repr.colors.mesh = color_index,
                                    setting_id::sphere_color => atom.repr.colors.sphere = color_index,
                                    setting_id::ribbon_color => atom.repr.colors.ribbon = color_index,
                                    _ => {}
                                }
                            }
                        }
                        total_affected += count;
                        mol_obj.invalidate(DirtyFlags::COLOR);
                    }
                }
            }

            ctx.viewer.request_redraw();

            if !ctx.quiet {
                ctx.print(&format!(" Set {} = {} for {} atoms", name, value_str, total_affected));
            }

            return Ok(());
        }

        // Handle sphere_scale - can be applied per-atom like in original PyMOL
        if id == setting_id::sphere_scale {
            let scale_value = match &value {
                SettingValue::Float(v) => *v,
                _ => return Err(CmdError::execution("Expected float value for sphere_scale")),
            };

            // If selection is provided, apply per-atom
            if let Some(selection_str) = selection {
                let selection_results = evaluate_selection(ctx.viewer, selection_str)?;
                let mut total_affected = 0usize;

                for (obj_name, selected) in selection_results {
                    let count = selected.count();
                    if count > 0 {
                        if let Some(mol_obj) = ctx.viewer.objects_mut().get_molecule_mut(&obj_name) {
                            let mol_mut = mol_obj.molecule_mut();
                            for idx in selected.indices() {
                                if let Some(atom) = mol_mut.get_atom_mut(AtomIndex(idx.0)) {
                                    atom.repr.sphere_scale = Some(scale_value);
                                }
                            }
                            total_affected += count;
                            mol_obj.invalidate(DirtyFlags::REPS);
                        }
                    }
                }

                ctx.viewer.request_redraw();

                if !ctx.quiet {
                    ctx.print(&format!(" Set sphere_scale = {} for {} atoms", scale_value, total_affected));
                }

                return Ok(());
            }
            // If no selection, fall through to set global value
        }

        // Set the global value
        ctx.viewer.settings_mut().set(id, value.clone())
            .map_err(|e| CmdError::execution(e.to_string()))?;

        // Special handling for surface_quality - propagate to all molecules
        if id == setting_id::surface_quality {
            if let SettingValue::Int(quality) = &value {
                let names: Vec<_> = ctx.viewer.objects().names().map(|s| s.to_string()).collect();
                for obj_name in names {
                    if let Some(mol) = ctx.viewer.objects_mut().get_molecule_mut(&obj_name) {
                        mol.set_surface_quality(*quality);
                    }
                }
            }
        }
        
        // Special handling for transparency - invalidate all molecule representations
        // so surfaces get rebuilt with the new transparency value
        if id == setting_id::transparency {
            let names: Vec<_> = ctx.viewer.objects().names().map(|s| s.to_string()).collect();
            for obj_name in names {
                if let Some(mol) = ctx.viewer.objects_mut().get_molecule_mut(&obj_name) {
                    mol.invalidate_representations();
                }
            }
        }

        // Special handling for sphere_scale - invalidate all molecule representations
        // so spheres get rebuilt with the new scale value
        if id == setting_id::sphere_scale {
            let names: Vec<_> = ctx.viewer.objects().names().map(|s| s.to_string()).collect();
            for obj_name in names {
                if let Some(mol) = ctx.viewer.objects_mut().get_molecule_mut(&obj_name) {
                    mol.invalidate_representations();
                }
            }
        }

        // Special handling for cartoon geometry settings - invalidate all molecule
        // representations so cartoons get rebuilt with the new settings
        if matches!(id, 
            setting_id::cartoon_fancy_helices |
            setting_id::cartoon_fancy_sheets |
            setting_id::cartoon_oval_width |
            setting_id::cartoon_oval_length |
            setting_id::cartoon_rect_width |
            setting_id::cartoon_rect_length |
            setting_id::cartoon_loop_radius |
            setting_id::cartoon_dumbbell_width |
            setting_id::cartoon_dumbbell_length |
            setting_id::cartoon_dumbbell_radius |
            setting_id::cartoon_round_helices |
            setting_id::cartoon_sampling |
            setting_id::cartoon_smooth_loops
        ) {
            let names: Vec<_> = ctx.viewer.objects().names().map(|s| s.to_string()).collect();
            for obj_name in names {
                if let Some(mol) = ctx.viewer.objects_mut().get_molecule_mut(&obj_name) {
                    mol.invalidate_representations();
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

// ============================================================================
// dss command
// ============================================================================

/// Define secondary structure using backbone geometry
struct DssCommand;

impl Command for DssCommand {
    fn name(&self) -> &str {
        "dss"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "dss" calculates secondary structure assignments using backbone
    phi/psi dihedral angles, similar to PyMOL's DSS algorithm.

    This overwrites any existing secondary structure assignments from
    the PDB/CIF file with computed values based on backbone geometry.

USAGE

    dss [selection [, state [, quiet ]]]

ARGUMENTS

    selection = string: atoms to process (default: all)
    state = integer: coordinate state to use (default: 1)
    quiet = 0/1: suppress feedback (default: 0)

NOTES

    The algorithm classifies residues based on phi/psi angles:
    - Alpha helix: phi ~-57°, psi ~-48°
    - Beta strand: phi ~-124°, psi ~124°

    This is automatically applied when loading structures if the
    auto_dss setting is enabled (default: on).

EXAMPLES

    dss
    dss protein
    dss all, 1, quiet=1
"#
    }

    fn execute<'v, 'r>(&self, ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>, args: &ParsedCommand) -> CmdResult {
        let _selection = args
            .get_str(0)
            .or_else(|| args.get_named_str("selection"))
            .unwrap_or("all");

        let state = args
            .get_int(1)
            .or_else(|| args.get_named_int("state"))
            .unwrap_or(1) as usize;
        
        // Convert to 0-indexed state
        let state_idx = if state > 0 { state - 1 } else { 0 };

        let quiet = args
            .get_bool(2)
            .or_else(|| args.get_named_bool("quiet"))
            .unwrap_or(false);

        // Get molecule names first to avoid borrowing conflicts
        let mol_names: Vec<String> = ctx.viewer
            .objects()
            .names()
            .map(|s| s.to_string())
            .collect();

        let settings = DssSettings::default();
        let mut total_updated = 0usize;

        // Apply DSS to all molecule objects
        // TODO: Support selection filtering when selection system is more developed
        for name in mol_names {
            if let Some(mol_obj) = ctx.viewer.objects_mut().get_molecule_mut(&name) {
                let mol = mol_obj.molecule_mut();
                let updated = assign_secondary_structure(mol, state_idx, &settings);
                total_updated += updated;
            }
        }

        ctx.viewer.request_redraw();

        if !quiet {
            ctx.print(&format!(" Assigned secondary structure to {} atoms", total_updated));
        }

        Ok(())
    }
}
