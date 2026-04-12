//! Settings commands: set, get, unset, dss

use crate::args::ParsedCommand;
use crate::command::{ArgHint, Command, CommandContext, CommandRegistry, DynamicSettingEntry, ViewerLike};
use crate::commands::selecting::evaluate_selection;
use crate::error::{CmdError, CmdResult};
use pymol_scene::{DirtyFlags, Object};
use pymol_select::AtomIndex;
use pymol_settings::{registry, DynamicSettingDescriptor, SideEffectCategory, SettingType, SettingValue};
use pymol_settings::SettingDescriptor;

/// Register settings commands
pub fn register(registry: &mut CommandRegistry) {
    registry.register(SetCommand);
    registry.register(GetCommand);
    registry.register(UnsetCommand);
}

/// Whether a setting stores a color value (packed RGB integer or color index).
///
/// In the typed system, color settings are stored as `i32` (SettingType::Int),
/// but need special parsing (color names, `[r,g,b]` vectors, scheme names).
fn is_color_setting(name: &str) -> bool {
    name.ends_with("_color") || name.starts_with("bg_rgb")
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

/// Format a setting value for display, using hint names when available
fn format_setting_display(desc: &SettingDescriptor, value: &SettingValue) -> String {
    if let Some(name) = desc.hint_name(value) {
        return name.to_string();
    }
    format_setting_value(value)
}

/// Format a dynamic setting value for display, using hint names when available
fn format_dynamic_display(desc: &DynamicSettingDescriptor, value: &SettingValue) -> String {
    if let Some(name) = desc.hint_name(value) {
        return name.to_string();
    }
    format_setting_value(value)
}

/// Apply side effects from a slice of categories.
///
/// Shared between built-in and dynamic settings.
fn apply_side_effects_from_slice<'v, 'r>(
    ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
    side_effects: &[SideEffectCategory],
    setting_name: &str,
    value: &SettingValue,
) {
    for category in side_effects {
        match category {
            SideEffectCategory::RepresentationRebuild => {
                if setting_name == "surface_quality" {
                    if let Some(quality) = value.as_int() {
                        let names: Vec<_> = ctx.viewer.objects().names().map(|s| s.to_string()).collect();
                        for obj_name in names {
                            if let Some(mol) = ctx.viewer.objects_mut().get_molecule_mut(&obj_name) {
                                mol.set_surface_quality(quality);
                            }
                        }
                    }
                }
                let names: Vec<_> = ctx.viewer.objects().names().map(|s| s.to_string()).collect();
                for obj_name in names {
                    if let Some(mol) = ctx.viewer.objects_mut().get_molecule_mut(&obj_name) {
                        mol.invalidate_representations();
                    }
                }
            }
            SideEffectCategory::ColorRebuild => {
                let names: Vec<_> = ctx.viewer.objects().names().map(|s| s.to_string()).collect();
                for obj_name in names {
                    if let Some(mol) = ctx.viewer.objects_mut().get_molecule_mut(&obj_name) {
                        mol.invalidate(DirtyFlags::COLOR);
                    }
                }
            }
            SideEffectCategory::SceneInvalidate
            | SideEffectCategory::SceneChanged
            | SideEffectCategory::FullRebuild => {}
            SideEffectCategory::ViewportUpdate => {
                if setting_name.starts_with("bg_rgb") {
                    if let Some(color_int) = value.as_int() {
                        let [r, g, b] = pymol_color::Color::from_packed_rgb(color_int).to_array();
                        ctx.viewer.set_background_color(r, g, b);
                    }
                }
            }
            _ => {}
        }
    }
}

// ============================================================================
// set command
// ============================================================================

/// Map a rep color setting name to the corresponding mutable field accessor in AtomColors
fn rep_color_field(name: &str) -> Option<fn(&mut pymol_mol::AtomColors) -> &mut i32> {
    match name {
        "stick_color" => Some(|c| &mut c.stick),
        "line_color" => Some(|c| &mut c.line),
        "cartoon_color" => Some(|c| &mut c.cartoon),
        "surface_color" => Some(|c| &mut c.surface),
        "mesh_color" => Some(|c| &mut c.mesh),
        "sphere_color" => Some(|c| &mut c.sphere),
        "ribbon_color" => Some(|c| &mut c.ribbon),
        _ => None,
    }
}

/// Map a rep color setting name to the corresponding read-only field accessor in AtomColors
fn rep_color_field_read(name: &str) -> Option<fn(&pymol_mol::AtomColors) -> i32> {
    match name {
        "stick_color" => Some(|c| c.stick),
        "line_color" => Some(|c| c.line),
        "cartoon_color" => Some(|c| c.cartoon),
        "surface_color" => Some(|c| c.surface),
        "mesh_color" => Some(|c| c.mesh),
        "sphere_color" => Some(|c| c.sphere),
        "ribbon_color" => Some(|c| c.ribbon),
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
        desc: &SettingDescriptor,
        value_str: &str,
    ) -> CmdResult<SettingValue> {
        // Resolve named value hints for settings with named variants
        if desc.has_value_hints() {
            if let Some(val) = desc.resolve_hint(value_str) {
                return Ok(val.clone());
            }
            if value_str.parse::<f64>().is_err() {
                let names: Vec<_> = desc.hint_names().collect();
                return Err(CmdError::invalid_arg(
                    "value",
                    format!("Unknown value '{}' for {}. Use {}", value_str, desc.name, names.join("/")),
                ));
            }
        }

        // Color type: resolve names/schemes
        if is_color_setting(desc.name) {
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
                        let packed = pymol_color::Color::new(
                            r.clamp(0.0, 1.0), g.clamp(0.0, 1.0), b.clamp(0.0, 1.0)
                        ).to_packed_rgb();
                        return Ok(SettingValue::Int(packed));
                    }
                }
            }

            if let Ok(v) = value_str.parse::<i32>() {
                return Ok(SettingValue::Int(v));
            }
            let color_index = if let Some(ci) = pymol_color::ColorIndex::from_scheme_name(value_str) {
                i32::from(ci)
            } else {
                ctx.viewer
                    .color_index(value_str)
                    .map(|idx| idx as i32)
                    .ok_or_else(|| CmdError::invalid_arg("value", format!("Unknown color: {}", value_str)))?
            };
            return Ok(SettingValue::Int(color_index));
        }

        // Default type-based parsing
        parse_setting_value(value_str, desc.setting_type)
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

    /// Execute a `set` for a built-in (static registry) setting.
    fn execute_builtin<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        args: &ParsedCommand,
        desc: &SettingDescriptor,
        name: &str,
        value_str: Option<String>,
        selection: Option<&str>,
    ) -> CmdResult {
        // Handle missing value (boolean toggle)
        let value_str = match value_str {
            Some(v) => v,
            None => {
                if desc.setting_type == SettingType::Bool {
                    let current = (desc.get)(ctx.viewer.settings()).as_bool().unwrap_or(false);
                    let new_value = !current;
                    (desc.set)(ctx.viewer.settings_mut(), SettingValue::Bool(new_value))
                        .map_err(|e| CmdError::execution(e.to_string()))?;
                    self.apply_side_effects(ctx, desc, &SettingValue::Bool(new_value));
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

        // For Float3 settings or color settings: collect 3 comma-separated floats
        let (value_str, selection) = if (desc.setting_type == SettingType::Float3
            || is_color_setting(desc.name))
            && !value_str.starts_with('[')
            && value_str.parse::<f32>().is_ok()
        {
            if let (Some(y), Some(z)) = (
                args.get_arg(2).and_then(|v| v.to_string_repr()),
                args.get_arg(3).and_then(|v| v.to_string_repr()),
            ) {
                if y.parse::<f32>().is_ok() && z.parse::<f32>().is_ok() {
                    let combined = format!("[{}, {}, {}]", value_str, y, z);
                    let shifted_selection = args.str_arg(4, "selection");
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

        let value = self.parse_value(ctx, desc, &value_str)?;

        // Per-atom rep color
        if let Some(field_accessor) = rep_color_field(desc.name) {
            let is_global = selection.is_none();
            let selection_str = selection.unwrap_or("all");
            let color_index = value.as_int()
                .ok_or_else(|| CmdError::execution("Expected color/int value"))?;

            if is_global {
                (desc.set)(ctx.viewer.settings_mut(), value.clone())
                    .map_err(|e| CmdError::execution(e.to_string()))?;
            }

            let total_affected = self.apply_per_atom_color(ctx, selection_str, color_index, field_accessor)?;
            ctx.viewer.request_redraw();

            if !ctx.quiet {
                ctx.print(&format!(" Set {} = {} for {} atoms", name, value_str, total_affected));
            }
            return Ok(());
        }

        // Per-atom sphere_scale
        if desc.name == "sphere_scale" {
            if let Some(selection_str) = selection {
                let scale_value = value.as_float()
                    .ok_or_else(|| CmdError::execution("Expected float value for sphere_scale"))?;

                let total_affected = self.apply_per_atom_sphere_scale(ctx, selection_str, scale_value)?;
                ctx.viewer.request_redraw();

                if !ctx.quiet {
                    ctx.print(&format!(" Set sphere_scale = {} for {} atoms", scale_value, total_affected));
                }
                return Ok(());
            }
        }

        // Per-object override
        if let Some(sel) = selection {
            if let Some(set_override) = desc.set_override {
                if let Some(mol) = ctx.viewer.objects_mut().get_molecule_mut(sel) {
                    set_override(mol.get_or_create_overrides(), value.clone())
                        .map_err(|e| CmdError::execution(e.to_string()))?;
                    mol.invalidate(DirtyFlags::COLOR);
                    ctx.viewer.request_redraw();
                    if !ctx.quiet {
                        let display = format_setting_display(desc, &value);
                        ctx.print(&format!(" {} = {} (object {})", name, display, sel));
                    }
                    return Ok(());
                }
            }
        }

        // Global set
        (desc.set)(ctx.viewer.settings_mut(), value.clone())
            .map_err(|e| CmdError::execution(e.to_string()))?;

        self.apply_side_effects(ctx, desc, &value);
        ctx.viewer.request_redraw();

        if !ctx.quiet {
            let display = format_setting_display(desc, &value);
            ctx.print(&format!(" {} = {}", name, display));
        }

        Ok(())
    }

    /// Execute a `set` for a dynamic (plugin-registered) setting.
    fn execute_dynamic<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        entry: &DynamicSettingEntry,
        name: &str,
        value_str: Option<String>,
        selection: Option<&str>,
    ) -> CmdResult {
        let desc = &entry.descriptor;

        // Handle missing value (boolean toggle)
        let value_str = match value_str {
            Some(v) => v,
            None => {
                if desc.setting_type == SettingType::Bool {
                    let store = entry.store.read().map_err(|e| CmdError::execution(e.to_string()))?;
                    let current = store.get(name)
                        .unwrap_or(&desc.default)
                        .as_bool()
                        .unwrap_or(false);
                    drop(store);
                    let new_value = SettingValue::Bool(!current);
                    entry.store.write().map_err(|e| CmdError::execution(e.to_string()))?
                        .set(name, new_value.clone());
                    apply_side_effects_from_slice(ctx, &desc.side_effects, name, &new_value);
                    ctx.viewer.request_redraw();
                    if !ctx.quiet {
                        ctx.print(&format!(" {} = {}", name, if !current { "on" } else { "off" }));
                    }
                    return Ok(());
                } else {
                    return Err(CmdError::MissingArgument("value".to_string()));
                }
            }
        };

        // Parse value — resolve hints, then type-based parsing
        let value = if desc.has_value_hints() {
            if let Some(val) = desc.resolve_hint(&value_str) {
                val.clone()
            } else if value_str.parse::<f64>().is_err() {
                let names: Vec<_> = desc.hint_names().collect();
                return Err(CmdError::invalid_arg(
                    "value",
                    format!("Unknown value '{}' for {}. Use {}", value_str, name, names.join("/")),
                ));
            } else {
                parse_setting_value(&value_str, desc.setting_type)?
            }
        } else {
            parse_setting_value(&value_str, desc.setting_type)?
        };

        // Validate min/max
        if let (Some(min), Some(max)) = (desc.min, desc.max) {
            if let Some(fval) = value.as_float() {
                if fval < min || fval > max {
                    return Err(CmdError::invalid_arg(
                        "value",
                        format!("Value {} out of range [{}, {}] for {}", fval, min, max, name),
                    ));
                }
            }
        }

        // Per-object override
        if let Some(sel) = selection {
            if desc.object_overridable {
                entry.store.write().map_err(|e| CmdError::execution(e.to_string()))?
                    .set_object(sel, name, value.clone());
                apply_side_effects_from_slice(ctx, &desc.side_effects, name, &value);
                ctx.viewer.request_redraw();
                if !ctx.quiet {
                    let display = format_dynamic_display(desc, &value);
                    ctx.print(&format!(" {} = {} (object {})", name, display, sel));
                }
                return Ok(());
            }
        }

        // Global set
        entry.store.write().map_err(|e| CmdError::execution(e.to_string()))?
            .set(name, value.clone());
        apply_side_effects_from_slice(ctx, &desc.side_effects, name, &value);
        ctx.viewer.request_redraw();

        if !ctx.quiet {
            let display = format_dynamic_display(desc, &value);
            ctx.print(&format!(" {} = {}", name, display));
        }

        Ok(())
    }

    /// Apply side effects based on the descriptor's side_effects list
    fn apply_side_effects<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        desc: &SettingDescriptor,
        value: &SettingValue,
    ) {
        apply_side_effects_from_slice(ctx, desc.side_effects, desc.name, value);
    }
}

impl Command for SetCommand {
    fn name(&self) -> &str {
        "set"
    }

    fn arg_hints(&self) -> &[ArgHint] {
        &[ArgHint::Setting, ArgHint::None, ArgHint::Selection]
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
            .str_arg(0, "name")
            .ok_or_else(|| CmdError::MissingArgument("name".to_string()))?;

        let value_str = args
            .get_arg(1)
            .and_then(|v| v.to_string_repr())
            .or_else(|| args.get_named("value").and_then(|v| v.to_string_repr()));

        let selection = args.str_arg(2, "selection");

        // === 2. Intercept pseudo-settings ===
        if name == "state" {
            return self.handle_set_state(ctx, value_str);
        }

        // === 3. Resolve setting descriptor ===
        if let Some(desc) = registry::lookup_by_name(name) {
            return self.execute_builtin(ctx, args, desc, name, value_str, selection);
        }

        // === 3b. Fall through to dynamic (plugin) settings ===
        if let Some(entry) = ctx.dynamic_setting(name).cloned() {
            return self.execute_dynamic(ctx, &entry, name, value_str, selection);
        }

        Err(CmdError::invalid_arg("name", format!("Unknown setting: {}", name)))
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

    fn arg_hints(&self) -> &[ArgHint] {
        &[ArgHint::Setting, ArgHint::Selection]
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
            .str_arg(0, "name")
            .ok_or_else(|| CmdError::MissingArgument("name".to_string()))?;

        // Intercept "get state" — report displayed state for each object
        if name == "state" {
            let obj_arg = args.str_arg(1, "selection");
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

        // Look up in built-in registry first
        if let Some(desc) = registry::lookup_by_name(name) {
            let selection = args.str_arg(1, "selection");
            if let Some(sel) = selection {
                if let Some(field_reader) = rep_color_field_read(desc.name) {
                    let selection_results = evaluate_selection(ctx.viewer, sel)?;
                    for (obj_name, selected) in &selection_results {
                        if let Some(mol) = ctx.viewer.objects().get_molecule(obj_name) {
                            for idx in selected.indices() {
                                if let Some(atom) = mol.molecule().get_atom(idx) {
                                    let color_val = field_reader(&atom.repr.colors);
                                    let value = SettingValue::Int(color_val);
                                    let display = format_setting_display(desc, &value);
                                    ctx.print(&format!(" {} = {} ({})", name, display, sel));
                                    return Ok(());
                                }
                            }
                        }
                    }
                }

                if let Some(get_override) = desc.get_override {
                    if let Some(mol) = ctx.viewer.objects().get_molecule(sel) {
                        if let Some(overrides) = mol.overrides() {
                            if let Some(value) = get_override(overrides) {
                                let display = format_setting_value(&value);
                                if let Some(hint_name) = desc.hint_name(&value) {
                                    ctx.print(&format!(" {} = {} ({}, object {})", name, display, hint_name, sel));
                                } else {
                                    ctx.print(&format!(" {} ({}) = {} (object {})", name, desc.setting_type, display, sel));
                                }
                                return Ok(());
                            }
                        }
                    }
                }
            }

            let value = (desc.get)(ctx.viewer.settings());
            if let Some(hint_name) = desc.hint_name(&value) {
                ctx.print(&format!(" {} = {} ({})", name, format_setting_value(&value), hint_name));
            } else {
                ctx.print(&format!(" {} ({}) = {}", name, desc.setting_type, format_setting_value(&value)));
            }
            return Ok(());
        }

        // Fall through to dynamic (plugin) settings
        if let Some(entry) = ctx.dynamic_setting(name).cloned() {
            let desc = &entry.descriptor;
            let selection = args.str_arg(1, "selection");
            let store = entry.store.read().map_err(|e| CmdError::execution(e.to_string()))?;

            // Per-object override
            if let Some(sel) = selection {
                if desc.object_overridable {
                    if let Some(value) = store.get_object(sel, name) {
                        let display = format_dynamic_display(desc, value);
                        ctx.print(&format!(" {} ({}) = {} (object {})", name, desc.setting_type, display, sel));
                        return Ok(());
                    }
                }
            }

            // Global value (fall back to default)
            let value = store.get(name).unwrap_or(&desc.default);
            let display = format_dynamic_display(desc, value);
            if let Some(hint_name) = desc.hint_name(value) {
                ctx.print(&format!(" {} = {} ({})", name, format_setting_value(value), hint_name));
            } else {
                ctx.print(&format!(" {} ({}) = {}", name, desc.setting_type, display));
            }
            return Ok(());
        }

        Err(CmdError::invalid_arg("name", format!("Unknown setting: {}", name)))
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

    fn arg_hints(&self) -> &[ArgHint] {
        &[ArgHint::Setting, ArgHint::Selection]
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
            .str_arg(0, "name")
            .ok_or_else(|| CmdError::MissingArgument("name".to_string()))?;

        // Look up in built-in registry first
        if let Some(desc) = registry::lookup_by_name(name) {
            let selection = args.str_arg(1, "selection");

            if let Some(sel) = selection {
                if let Some(field_accessor) = rep_color_field(desc.name) {
                    let selection_results = evaluate_selection(ctx.viewer, sel)?;
                    let mut total = 0usize;
                    for (obj_name, selected) in selection_results {
                        if selected.count() > 0 {
                            if let Some(mol_obj) = ctx.viewer.objects_mut().get_molecule_mut(&obj_name) {
                                for idx in selected.indices() {
                                    if let Some(atom) = mol_obj.molecule_mut().get_atom_mut(AtomIndex(idx.0)) {
                                        *field_accessor(&mut atom.repr.colors) = pymol_mol::COLOR_UNSET;
                                    }
                                }
                                total += selected.count();
                                mol_obj.invalidate(DirtyFlags::COLOR);
                            }
                        }
                    }
                    ctx.viewer.request_redraw();
                    if !ctx.quiet {
                        ctx.print(&format!(" {} reset for {} atoms", name, total));
                    }
                    return Ok(());
                }

                if let Some(unset_override) = desc.unset_override {
                    if let Some(mol) = ctx.viewer.objects_mut().get_molecule_mut(sel) {
                        if let Some(overrides) = mol.overrides_mut() {
                            unset_override(overrides);
                        }
                        mol.invalidate(DirtyFlags::COLOR);
                        ctx.viewer.request_redraw();
                        if !ctx.quiet {
                            ctx.print(&format!(" {} reset to global default (object {})", name, sel));
                        }
                        return Ok(());
                    }
                }
            }

            if let Some(field_accessor) = rep_color_field(desc.name) {
                let names: Vec<_> = ctx.viewer.objects().names().map(|s| s.to_string()).collect();
                for obj_name in names {
                    if let Some(mol_obj) = ctx.viewer.objects_mut().get_molecule_mut(&obj_name) {
                        for atom in mol_obj.molecule_mut().atoms_mut() {
                            *field_accessor(&mut atom.repr.colors) = pymol_mol::COLOR_UNSET;
                        }
                        mol_obj.invalidate(DirtyFlags::COLOR);
                    }
                }
            }

            let default_value = (desc.get)(&pymol_settings::Settings::default());
            (desc.set)(ctx.viewer.settings_mut(), default_value.clone())
                .map_err(|e| CmdError::execution(e.to_string()))?;

            ctx.viewer.request_redraw();

            if !ctx.quiet {
                ctx.print(&format!(" {} reset to default: {}", name, format_setting_value(&default_value)));
            }

            return Ok(());
        }

        // Fall through to dynamic (plugin) settings
        if let Some(entry) = ctx.dynamic_setting(name).cloned() {
            let desc = &entry.descriptor;
            let selection = args.str_arg(1, "selection");

            if let Some(sel) = selection {
                if desc.object_overridable {
                    entry.store.write().map_err(|e| CmdError::execution(e.to_string()))?
                        .remove_object(sel, name);
                    ctx.viewer.request_redraw();
                    if !ctx.quiet {
                        ctx.print(&format!(" {} reset to global default (object {})", name, sel));
                    }
                    return Ok(());
                }
            }

            // Reset global to default
            entry.store.write().map_err(|e| CmdError::execution(e.to_string()))?
                .remove(name);
            ctx.viewer.request_redraw();
            if !ctx.quiet {
                ctx.print(&format!(" {} reset to default: {}", name, format_setting_value(&desc.default)));
            }
            return Ok(());
        }

        Err(CmdError::invalid_arg("name", format!("Unknown setting: {}", name)))
    }
}

