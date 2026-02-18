//! Settings commands: set, get, unset, dss

use crate::args::ParsedCommand;
use crate::command::{Command, CommandContext, CommandRegistry, ViewerLike};
use crate::commands::selecting::evaluate_selection;
use crate::error::{CmdError, CmdResult};
use pymol_mol::dss::{assign_secondary_structure, DssSettings};
use pymol_scene::DirtyFlags;
use pymol_select::AtomIndex;
use pymol_settings::{get_setting, get_setting_id, get_side_effects, id as setting_id, ShadingMode, SideEffectCategory, SettingType, SettingValue};

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

/// Map a rep color setting ID to the corresponding field accessor in AtomReprColors
fn rep_color_field(id: u16) -> Option<fn(&mut pymol_mol::AtomColors) -> &mut i32> {
    match id {
        setting_id::stick_color => Some(|c| &mut c.stick),
        setting_id::line_color => Some(|c| &mut c.line),
        setting_id::cartoon_color => Some(|c| &mut c.cartoon),
        setting_id::surface_color => Some(|c| &mut c.surface),
        setting_id::mesh_color => Some(|c| &mut c.mesh),
        setting_id::sphere_color => Some(|c| &mut c.sphere),
        setting_id::ribbon_color => Some(|c| &mut c.ribbon),
        _ => None,
    }
}

struct SetCommand;

impl SetCommand {
    /// Handle the pseudo-setting "set state, N"
    fn handle_set_state<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        value_str: Option<String>,
    ) -> CmdResult {
        let value_str = value_str.ok_or_else(|| CmdError::MissingArgument("value".to_string()))?;
        let state_num: i64 = value_str.parse()
            .map_err(|_| CmdError::invalid_arg("value", format!("Invalid state number: {}", value_str)))?;
        if state_num < 1 {
            return Err(CmdError::invalid_arg("value", "State number must be >= 1"));
        }
        let state_idx = (state_num - 1) as usize;
        let names: Vec<String> = ctx.viewer.objects().names().map(|s| s.to_string()).collect();
        for obj_name in &names {
            if let Some(mol_obj) = ctx.viewer.objects_mut().get_molecule_mut(obj_name) {
                mol_obj.set_display_state(state_idx);
            }
        }
        ctx.viewer.request_redraw();
        if !ctx.quiet {
            ctx.print(&format!(" state = {}", state_num));
        }
        Ok(())
    }

    /// Parse a value string with per-setting overrides (shading_mode aliases, color names)
    fn parse_value<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        id: u16,
        setting_type: SettingType,
        value_str: &str,
    ) -> CmdResult<SettingValue> {
        // shading_mode: accept string aliases
        if id == setting_id::shading_mode {
            return if let Some(mode) = ShadingMode::from_str_alias(value_str) {
                Ok(SettingValue::Int(i32::from(mode)))
            } else {
                Err(CmdError::invalid_arg(
                    "value",
                    format!("Unknown shading mode '{}'. Use classic or skripkin", value_str),
                ))
            };
        }

        // Color type: resolve names/schemes
        if setting_type == SettingType::Color {
            // Accept [r, g, b] float vector → convert to 0x00RRGGBB integer
            if value_str.starts_with('[') && value_str.ends_with(']') {
                let inner = value_str[1..value_str.len()-1].trim();
                let parts: Vec<&str> = inner.split(',').map(|p| p.trim()).collect();
                if parts.len() == 3 {
                    if let (Ok(r), Ok(g), Ok(b)) = (
                        parts[0].parse::<f32>(),
                        parts[1].parse::<f32>(),
                        parts[2].parse::<f32>(),
                    ) {
                        let ri = (r.clamp(0.0, 1.0) * 255.0) as i32;
                        let gi = (g.clamp(0.0, 1.0) * 255.0) as i32;
                        let bi = (b.clamp(0.0, 1.0) * 255.0) as i32;
                        return Ok(SettingValue::Color((ri << 16) | (gi << 8) | bi));
                    }
                }
            }

            if let Ok(v) = value_str.parse::<i32>() {
                return Ok(SettingValue::Color(v));
            }
            let color_index = match value_str.to_lowercase().as_str() {
                "atomic" | "cpk" | "element" | "by_element" => -1,
                "chain" | "by_chain" | "chainbow" => -2,
                "ss" | "secondary_structure" | "by_ss" | "dssp" => -3,
                "b" | "b_factor" | "bfactor" | "by_b" => -4,
                _ => {
                    ctx.viewer
                        .color_index(value_str)
                        .map(|idx| idx as i32)
                        .ok_or_else(|| CmdError::invalid_arg("value", format!("Unknown color: {}", value_str)))?
                }
            };
            return Ok(SettingValue::Color(color_index));
        }

        // Default type-based parsing
        parse_setting_value(value_str, setting_type)
    }

    /// Apply a rep color per-atom via the dispatch table
    fn apply_per_atom_color<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        selection_str: &str,
        color_index: i32,
        field_accessor: fn(&mut pymol_mol::AtomColors) -> &mut i32,
    ) -> CmdResult<usize> {
        let selection_results = evaluate_selection(ctx.viewer, selection_str)?;
        let mut total_affected = 0usize;

        for (obj_name, selected) in selection_results {
            let count = selected.count();
            if count > 0 {
                if let Some(mol_obj) = ctx.viewer.objects_mut().get_molecule_mut(&obj_name) {
                    let mol_mut = mol_obj.molecule_mut();
                    for idx in selected.indices() {
                        if let Some(atom) = mol_mut.get_atom_mut(AtomIndex(idx.0)) {
                            *field_accessor(&mut atom.repr.colors) = color_index;
                        }
                    }
                    total_affected += count;
                    mol_obj.invalidate(DirtyFlags::COLOR);
                }
            }
        }

        Ok(total_affected)
    }

    /// Apply sphere_scale per-atom
    fn apply_per_atom_sphere_scale<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        selection_str: &str,
        scale_value: f32,
    ) -> CmdResult<usize> {
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

        Ok(total_affected)
    }

    /// Apply side effects based on setting categories from get_side_effects()
    fn apply_side_effects<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        id: u16,
        value: &SettingValue,
    ) {
        let categories = get_side_effects(id);
        for category in categories {
            match category {
                SideEffectCategory::RepresentationRebuild => {
                    // Special case: surface_quality needs propagation
                    if id == setting_id::surface_quality {
                        if let SettingValue::Int(quality) = value {
                            let names: Vec<_> = ctx.viewer.objects().names().map(|s| s.to_string()).collect();
                            for obj_name in names {
                                if let Some(mol) = ctx.viewer.objects_mut().get_molecule_mut(&obj_name) {
                                    mol.set_surface_quality(*quality);
                                }
                            }
                        }
                    }
                    // Invalidate all molecule representations
                    let names: Vec<_> = ctx.viewer.objects().names().map(|s| s.to_string()).collect();
                    for obj_name in names {
                        if let Some(mol) = ctx.viewer.objects_mut().get_molecule_mut(&obj_name) {
                            mol.invalidate_representations();
                        }
                    }
                }
                SideEffectCategory::SceneInvalidate
                | SideEffectCategory::SceneChanged
                | SideEffectCategory::FullRebuild => {
                    // These all need a redraw (handled after this function)
                }
                SideEffectCategory::ViewportUpdate => {
                    // Background color: convert color int to RGB floats and apply
                    if matches!(id, setting_id::bg_rgb | setting_id::bg_rgb_top | setting_id::bg_rgb_bottom) {
                        if let SettingValue::Color(color_int) = value {
                            let r = ((*color_int >> 16) & 0xFF) as f32 / 255.0;
                            let g = ((*color_int >> 8) & 0xFF) as f32 / 255.0;
                            let b = (*color_int & 0xFF) as f32 / 255.0;
                            ctx.viewer.set_background_color(r, g, b);
                        }
                    }
                }
                _ => {
                    // Other categories (ShaderReload, etc.) — not yet wired
                }
            }
        }
    }
}

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
        // === 1. Parse arguments ===
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

        // === 2. Intercept pseudo-settings ===
        if name == "state" {
            return self.handle_set_state(ctx, value_str);
        }

        // === 3. Resolve setting ===
        let id = get_setting_id(name)
            .ok_or_else(|| CmdError::invalid_arg("name", format!("Unknown setting: {}", name)))?;

        let setting = get_setting(id)
            .ok_or_else(|| CmdError::invalid_arg("name", "Invalid setting ID"))?;

        // === 4. Handle missing value (boolean toggle) ===
        let value_str = match value_str {
            Some(v) => v,
            None => {
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

        // === 5. Parse value (with per-setting overrides) ===
        // For Float3 settings: collect 3 comma-separated floats into a single value
        // e.g. "set silhouette_color, 1.0, 1.0, 1.0" → Float3([1.0, 1.0, 1.0])
        // When floats are collected from args 2,3 — selection shifts to arg 4
        // Also collect for Color settings — supports "set bg_rgb, 1.0, 1.0, 1.0"
        let (value_str, selection) = if (setting.setting_type == SettingType::Float3
            || setting.setting_type == SettingType::Color)
            && !value_str.starts_with('[')
            && value_str.parse::<f32>().is_ok()
        {
            if let (Some(y), Some(z)) = (
                args.get_arg(2).and_then(|v| v.to_string_repr()),
                args.get_arg(3).and_then(|v| v.to_string_repr()),
            ) {
                if y.parse::<f32>().is_ok() && z.parse::<f32>().is_ok() {
                    let combined = format!("[{}, {}, {}]", value_str, y, z);
                    let shifted_selection = args.get_str(4).or_else(|| args.get_named_str("selection"));
                    (combined, shifted_selection)
                } else {
                    (value_str, selection)
                }
            } else {
                (value_str, selection)
            }
        } else {
            (value_str, selection)
        };

        let value = self.parse_value(ctx, id, setting.setting_type, &value_str)?;

        // === 6. Apply (per-atom or global) ===
        // Per-atom rep color
        if let Some(field_accessor) = rep_color_field(id) {
            let selection_str = selection.unwrap_or("all");
            let color_index = match &value {
                SettingValue::Color(idx) => *idx,
                _ => return Err(CmdError::execution("Expected color value")),
            };

            let total_affected = self.apply_per_atom_color(ctx, selection_str, color_index, field_accessor)?;
            ctx.viewer.request_redraw();

            if !ctx.quiet {
                ctx.print(&format!(" Set {} = {} for {} atoms", name, value_str, total_affected));
            }
            return Ok(());
        }

        // Per-atom sphere_scale
        if id == setting_id::sphere_scale {
            if let Some(selection_str) = selection {
                let scale_value = match &value {
                    SettingValue::Float(v) => *v,
                    _ => return Err(CmdError::execution("Expected float value for sphere_scale")),
                };

                let total_affected = self.apply_per_atom_sphere_scale(ctx, selection_str, scale_value)?;
                ctx.viewer.request_redraw();

                if !ctx.quiet {
                    ctx.print(&format!(" Set sphere_scale = {} for {} atoms", scale_value, total_affected));
                }
                return Ok(());
            }
            // No selection → fall through to global set
        }

        // Global set
        ctx.viewer.settings_mut().set(id, value.clone())
            .map_err(|e| CmdError::execution(e.to_string()))?;

        // === 7. Side effects (driven by get_side_effects) ===
        self.apply_side_effects(ctx, id, &value);

        ctx.viewer.request_redraw();

        // === 8. Print feedback ===
        if !ctx.quiet {
            if id == setting_id::shading_mode {
                if let SettingValue::Int(v) = &value {
                    let mode = ShadingMode::from(*v);
                    ctx.print(&format!(" {} = {}", name, mode.name()));
                } else {
                    ctx.print(&format!(" {} = {}", name, format_setting_value(&value)));
                }
            } else {
                ctx.print(&format!(" {} = {}", name, format_setting_value(&value)));
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

    fn execute<'v, 'r>(&self, ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>, args: &ParsedCommand) -> CmdResult {
        let name = args
            .get_str(0)
            .or_else(|| args.get_named_str("name"))
            .ok_or_else(|| CmdError::MissingArgument("name".to_string()))?;

        // Intercept "get state" — report displayed state for each object
        if name == "state" {
            let obj_arg = args.get_str(1).or_else(|| args.get_named_str("selection"));
            if let Some(obj_name) = obj_arg {
                if let Some(mol_obj) = ctx.viewer.objects().get_molecule(obj_name) {
                    ctx.print(&format!(" state (int) = {} for \"{}\"", mol_obj.display_state() + 1, obj_name));
                } else {
                    return Err(CmdError::ObjectNotFound(obj_name.to_string()));
                }
            } else {
                let names: Vec<String> = ctx.viewer.objects().names().map(|s| s.to_string()).collect();
                for obj_name in &names {
                    if let Some(mol_obj) = ctx.viewer.objects().get_molecule(obj_name) {
                        ctx.print(&format!(" state (int) = {} for \"{}\"", mol_obj.display_state() + 1, obj_name));
                    }
                }
            }
            return Ok(());
        }

        // Look up the setting by name
        let id = get_setting_id(name)
            .ok_or_else(|| CmdError::invalid_arg("name", format!("Unknown setting: {}", name)))?;

        let setting = get_setting(id)
            .ok_or_else(|| CmdError::invalid_arg("name", "Invalid setting ID"))?;

        // Get the current value
        let value = ctx.viewer.settings().get(id)
            .ok_or_else(|| CmdError::execution("Failed to get setting value"))?;

        // Display the value with type information
        // For shading_mode, show the human-readable name
        if id == setting_id::shading_mode {
            if let SettingValue::Int(v) = &value {
                let mode = ShadingMode::from(*v);
                ctx.print(&format!(" {} = {} ({})", name, v, mode.name()));
            } else {
                ctx.print(&format!(" {} = {}", name, format_setting_value(&value)));
            }
        } else {
            ctx.print(&format!(" {} ({}) = {}", name, setting.setting_type, format_setting_value(&value)));
        }

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
