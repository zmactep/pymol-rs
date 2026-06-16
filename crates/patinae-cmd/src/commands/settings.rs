//! Settings commands: set, get, unset, dss

use crate::args::ParsedCommand;
use crate::command::{ArgHint, Command, CommandContext, CommandRegistry, ViewerLike};
use crate::command_help;
use crate::commands::selecting::evaluate_selection;
use crate::error::{CmdError, CmdResult};
use crate::ResolvedSetting;
use patinae_scene::{DirtyFlags, MoleculeObject, Object};
use patinae_select::AtomIndex;
use patinae_settings::{SettingType, SettingValue, SideEffectCategory};

/// Register settings commands
pub fn register(registry: &mut CommandRegistry) {
    registry.register(SetCommand);
    registry.register(GetCommand);
    registry.register(UnsetCommand);
    registry.register(SetupCcdCommand);
}

/// Resolve deprecated setting names to their canonical forms, printing a
/// deprecation warning to `ctx` when an alias is hit. Returns the canonical
/// name; callers should use the returned string for descriptor lookup.
fn resolve_legacy_setting_name<'v, 'r, 'a>(
    ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
    name: &'a str,
) -> std::borrow::Cow<'a, str> {
    match name {
        "transparency" => {
            ctx.print_warning("'transparency' is deprecated, use 'surface_transparency' instead");
            std::borrow::Cow::Borrowed("surface_transparency")
        }
        _ => std::borrow::Cow::Borrowed(name),
    }
}

/// Format a setting source as an unknown-setting error.
fn unknown_setting(name: &str) -> CmdError {
    CmdError::invalid_arg("name", format!("Unknown setting: {}", name))
}

fn format_object_report(
    setting: &ResolvedSetting,
    value: &SettingValue,
    object_name: &str,
) -> String {
    let raw = ResolvedSetting::format_value(value);
    if let Some(hint_name) = setting.hint_name(value) {
        format!(
            "{} = {} ({}, object {})",
            setting.name(),
            raw,
            hint_name,
            object_name
        )
    } else {
        format!(
            "{} ({}) = {} (object {})",
            setting.name(),
            setting.setting_type(),
            raw,
            object_name
        )
    }
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
                        let names: Vec<_> = ctx
                            .viewer
                            .objects()
                            .names()
                            .map(|s| s.to_string())
                            .collect();
                        for obj_name in names {
                            if let Some(mol) = ctx.viewer.objects_mut().get_molecule_mut(&obj_name)
                            {
                                mol.set_surface_quality(quality);
                            }
                        }
                    }
                }
                let names: Vec<_> = ctx
                    .viewer
                    .objects()
                    .names()
                    .map(|s| s.to_string())
                    .collect();
                for obj_name in names {
                    if let Some(mol) = ctx.viewer.objects_mut().get_molecule_mut(&obj_name) {
                        mol.invalidate_representations();
                    }
                }
            }
            SideEffectCategory::ColorRebuild => {
                let names: Vec<_> = ctx
                    .viewer
                    .objects()
                    .names()
                    .map(|s| s.to_string())
                    .collect();
                for obj_name in names {
                    if let Some(mol) = ctx.viewer.objects_mut().get_molecule_mut(&obj_name) {
                        mol.invalidate(DirtyFlags::COLOR);
                    }
                }
            }
            SideEffectCategory::SurfaceTransparency => {
                let names: Vec<_> = ctx
                    .viewer
                    .objects()
                    .names()
                    .map(|s| s.to_string())
                    .collect();
                for obj_name in names {
                    if let Some(mol) = ctx.viewer.objects_mut().get_molecule_mut(&obj_name) {
                        mol.invalidate(DirtyFlags::TRANSPARENCY);
                    }
                }
            }
            SideEffectCategory::SceneInvalidate
            | SideEffectCategory::SceneChanged
            | SideEffectCategory::FullRebuild
            | SideEffectCategory::ViewportUpdate => {}
            _ => {}
        }
    }
}

fn apply_object_side_effects(
    mol: &mut MoleculeObject,
    side_effects: &[SideEffectCategory],
    setting_name: &str,
    value: &SettingValue,
) {
    let mut dirty = DirtyFlags::empty();

    for category in side_effects {
        match category {
            SideEffectCategory::RepresentationRebuild => {
                if setting_name == "surface_quality" {
                    if let Some(quality) = value.as_int() {
                        mol.set_surface_quality(quality);
                    }
                }
                dirty |= DirtyFlags::REPS;
            }
            SideEffectCategory::ColorRebuild => {
                dirty |= DirtyFlags::COLOR;
            }
            SideEffectCategory::SurfaceTransparency => {
                dirty |= DirtyFlags::TRANSPARENCY;
            }
            SideEffectCategory::SceneInvalidate
            | SideEffectCategory::SceneChanged
            | SideEffectCategory::ViewportUpdate => {
                dirty |= DirtyFlags::REPS;
            }
            SideEffectCategory::FullRebuild => {
                dirty |= DirtyFlags::ALL;
            }
            _ => {}
        }
    }

    if dirty.is_empty() {
        dirty = DirtyFlags::REPS;
    }
    mol.invalidate(dirty);
}

// ============================================================================
// set command
// ============================================================================

/// Map a rep color setting name to the corresponding mutable field accessor in AtomColors
fn rep_color_field(name: &str) -> Option<fn(&mut patinae_mol::AtomColors) -> &mut i32> {
    match name {
        "stick_color" => Some(|c| &mut c.stick),
        "line_color" => Some(|c| &mut c.line),
        "cartoon_color" => Some(|c| &mut c.cartoon),
        "surface_color" => Some(|c| &mut c.surface),
        "mesh_color" => Some(|c| &mut c.mesh),
        "sphere_color" => Some(|c| &mut c.sphere),
        "ribbon_color" => Some(|c| &mut c.ribbon),
        "dot_color" => Some(|c| &mut c.dot),
        "ellipsoid_color" => Some(|c| &mut c.ellipsoid),
        _ => None,
    }
}

/// Map a rep color setting name to the corresponding read-only field accessor in AtomColors
fn rep_color_field_read(name: &str) -> Option<fn(&patinae_mol::AtomColors) -> i32> {
    match name {
        "stick_color" => Some(|c| c.stick),
        "line_color" => Some(|c| c.line),
        "cartoon_color" => Some(|c| c.cartoon),
        "surface_color" => Some(|c| c.surface),
        "mesh_color" => Some(|c| c.mesh),
        "sphere_color" => Some(|c| c.sphere),
        "ribbon_color" => Some(|c| c.ribbon),
        "dot_color" => Some(|c| c.dot),
        "ellipsoid_color" => Some(|c| c.ellipsoid),
        _ => None,
    }
}

/// Map a per-atom transparency setting name to the corresponding mutable field accessor
fn per_atom_transparency_field(
    name: &str,
) -> Option<fn(&mut patinae_mol::AtomRepresentation) -> &mut Option<f32>> {
    match name {
        "sphere_transparency" => Some(|r| &mut r.sphere_transparency),
        "stick_transparency" => Some(|r| &mut r.stick_transparency),
        "cartoon_transparency" => Some(|r| &mut r.cartoon_transparency),
        "surface_transparency" => Some(|r| &mut r.surface_transparency),
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
        let value_str = value_str.ok_or_else(|| CmdError::missing_argument("value".to_string()))?;
        let state_num: i64 = value_str.parse().map_err(|_| {
            CmdError::invalid_arg("value", format!("Invalid state number: {}", value_str))
        })?;
        if state_num < 1 {
            return Err(CmdError::invalid_arg("value", "State number must be >= 1"));
        }
        let state_idx = (state_num - 1) as usize;
        let names: Vec<String> = ctx
            .viewer
            .objects()
            .names()
            .map(|s| s.to_string())
            .collect();
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

    /// Parse a value string with per-setting overrides.
    fn parse_value<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        setting: &ResolvedSetting,
        value_str: &str,
    ) -> CmdResult<SettingValue> {
        // Built-in color settings resolve names and color schemes through the viewer.
        if setting.built_in_descriptor().is_some() && setting.setting_type() == SettingType::Color {
            // Accept [r, g, b] float vector and convert it to 0x00RRGGBB.
            if value_str.starts_with('[') && value_str.ends_with(']') {
                let inner = value_str[1..value_str.len() - 1].trim();
                let parts: Vec<&str> = inner.split(',').map(|p| p.trim()).collect();
                if parts.len() == 3 {
                    if let (Ok(r), Ok(g), Ok(b)) = (
                        parts[0].parse::<f32>(),
                        parts[1].parse::<f32>(),
                        parts[2].parse::<f32>(),
                    ) {
                        let packed = patinae_color::Color::new(
                            r.clamp(0.0, 1.0),
                            g.clamp(0.0, 1.0),
                            b.clamp(0.0, 1.0),
                        )
                        .to_packed_rgb();
                        return Ok(SettingValue::Int(packed));
                    }
                }
            }

            if let Ok(v) = value_str.parse::<i32>() {
                return Ok(SettingValue::Int(v));
            }
            let color_index =
                if let Some(ci) = patinae_color::ColorIndex::from_scheme_name(value_str) {
                    i32::from(ci)
                } else if let Some(idx) = ctx.viewer.color_index(value_str) {
                    idx as i32
                } else if let Some(color) = patinae_color::Color::from_hex(value_str) {
                    ctx.viewer.named_palette_mut().set(value_str, color) as i32
                } else {
                    return Err(CmdError::invalid_arg(
                        "value",
                        format!("Unknown color: {}", value_str),
                    ));
                };
            return Ok(SettingValue::Int(color_index));
        }

        setting.parse_value(value_str)
    }

    /// Apply a rep color per-atom via the dispatch table
    fn apply_per_atom_color<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        selection_str: &str,
        color_index: i32,
        field_accessor: fn(&mut patinae_mol::AtomColors) -> &mut i32,
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

    /// Apply per-atom transparency for a selection.
    ///
    /// `dirty` controls how aggressively the affected objects are invalidated:
    /// for surface transparency we use `COLOR`, since `recolor_lut`'s fast path
    /// already respects each atom's `surface_transparency` override. Other reps
    /// (sphere/stick/cartoon) need a full `REPS` rebuild.
    fn apply_per_atom_transparency<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        selection_str: &str,
        transparency: f32,
        field_accessor: fn(&mut patinae_mol::AtomRepresentation) -> &mut Option<f32>,
        dirty: DirtyFlags,
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
                            *field_accessor(&mut atom.repr) = Some(transparency);
                        }
                    }
                    total_affected += count;
                    mol_obj.invalidate(dirty);
                }
            }
        }

        Ok(total_affected)
    }

    /// Execute a `set` for a resolved setting.
    fn execute_setting<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        args: &ParsedCommand,
        setting: ResolvedSetting,
        value_str: Option<String>,
        selection: Option<&str>,
    ) -> CmdResult {
        let name = setting.name().to_string();

        // Handle missing value (boolean toggle)
        let value_str = match value_str {
            Some(v) => v,
            None => {
                if setting.setting_type() == SettingType::Bool {
                    let current = setting
                        .global_value(ctx.viewer)
                        .map_err(CmdError::execution)?
                        .as_bool()
                        .unwrap_or(false);
                    let new_value = !current;
                    let value = SettingValue::Bool(new_value);
                    setting
                        .set_global(ctx.viewer, value.clone())
                        .map_err(CmdError::execution)?;
                    apply_side_effects_from_slice(ctx, setting.side_effects(), &name, &value);
                    ctx.viewer.request_redraw();
                    if !ctx.quiet {
                        ctx.print(&format!(
                            " {} = {}",
                            name,
                            if new_value { "on" } else { "off" }
                        ));
                    }
                    return Ok(());
                } else {
                    return Err(CmdError::missing_argument("value".to_string()));
                }
            }
        };

        // For Float3 settings or color settings: collect 3 comma-separated floats
        let (value_str, selection) = if (setting.setting_type() == SettingType::Float3
            || setting.setting_type() == SettingType::Color)
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

        let value = self.parse_value(ctx, &setting, &value_str)?;
        setting
            .validate_range(&value)
            .map_err(|msg| CmdError::invalid_arg("value", msg))?;

        if let Some(desc) = setting.built_in_descriptor() {
            // Per-atom rep color
            if let Some(field_accessor) = rep_color_field(desc.name) {
                let is_global = selection.is_none();
                let selection_str = selection.unwrap_or("all");
                let color_index = value
                    .as_int()
                    .ok_or_else(|| CmdError::execution("Expected color/int value"))?;

                if is_global {
                    setting
                        .set_global(ctx.viewer, value.clone())
                        .map_err(CmdError::execution)?;
                }

                let total_affected =
                    self.apply_per_atom_color(ctx, selection_str, color_index, field_accessor)?;
                ctx.viewer.request_redraw();

                if !ctx.quiet {
                    ctx.print(&format!(
                        " Set {} = {} for {} atoms",
                        name, value_str, total_affected
                    ));
                }
                return Ok(());
            }

            // Per-atom sphere_scale
            if desc.name == "sphere_scale" {
                if let Some(selection_str) = selection {
                    let scale_value = value.as_float().ok_or_else(|| {
                        CmdError::execution("Expected float value for sphere_scale")
                    })?;

                    let total_affected =
                        self.apply_per_atom_sphere_scale(ctx, selection_str, scale_value)?;
                    ctx.viewer.request_redraw();

                    if !ctx.quiet {
                        ctx.print(&format!(
                            " Set sphere_scale = {} for {} atoms",
                            scale_value, total_affected
                        ));
                    }
                    return Ok(());
                }
            }

            // Per-atom transparency (sphere_transparency, stick_transparency,
            // cartoon_transparency, transparency=surface). Surface uses the
            // recolor_lut fast path (DirtyFlags::COLOR); others need a rebuild.
            if let Some(field_accessor) = per_atom_transparency_field(desc.name) {
                if let Some(selection_str) = selection {
                    let trans_value = value.as_float().ok_or_else(|| {
                        CmdError::execution("Expected float value for transparency")
                    })?;

                    let dirty = if desc.name == "surface_transparency" {
                        DirtyFlags::COLOR
                    } else {
                        DirtyFlags::REPS
                    };
                    let total_affected = self.apply_per_atom_transparency(
                        ctx,
                        selection_str,
                        trans_value,
                        field_accessor,
                        dirty,
                    )?;
                    ctx.viewer.request_redraw();

                    if !ctx.quiet {
                        ctx.print(&format!(
                            " Set {} = {} for {} atoms",
                            name, trans_value, total_affected
                        ));
                    }
                    return Ok(());
                }
            }

            // Built-in object override
            if let Some(sel) = selection {
                if desc.is_object_overridable() {
                    if let Some(mol) = ctx.viewer.objects_mut().get_molecule_mut(sel) {
                        desc.set_override(mol.get_or_create_overrides(), value.clone())
                            .map_err(|e| CmdError::execution(e.to_string()))?;
                        apply_object_side_effects(mol, desc.side_effects, desc.name, &value);
                        ctx.viewer.request_redraw();
                        if !ctx.quiet {
                            let display = setting.format_display(&value);
                            ctx.print(&format!(" {} = {} (object {})", name, display, sel));
                        }
                        return Ok(());
                    }
                }
            }
        }

        // Dynamic object override
        if let Some(sel) = selection {
            if setting.is_object_overridable()
                && setting
                    .set_object_value(sel, value.clone())
                    .map_err(CmdError::execution)?
            {
                apply_side_effects_from_slice(ctx, setting.side_effects(), &name, &value);
                ctx.viewer.request_redraw();
                if !ctx.quiet {
                    let display = setting.format_display(&value);
                    ctx.print(&format!(" {} = {} (object {})", name, display, sel));
                }
                return Ok(());
            }
        }

        // Global set
        setting
            .set_global(ctx.viewer, value.clone())
            .map_err(CmdError::execution)?;
        apply_side_effects_from_slice(ctx, setting.side_effects(), &name, &value);
        ctx.viewer.request_redraw();

        if !ctx.quiet {
            let display = setting.format_display(&value);
            ctx.print(&format!(" {} = {}", name, display));
        }

        Ok(())
    }
}

impl Command for SetCommand {
    fn name(&self) -> &str {
        "set"
    }

    fn arg_hints(&self) -> &[ArgHint] {
        &[ArgHint::Setting, ArgHint::SettingValue, ArgHint::Selection]
    }

    command_help! {
        CMD "set"
        DESCRIPTION [
            "changes a setting value.",
        ]
        USAGE [
            "set name [, value [, selection [, state ]]]",
        ]
        REQUIRED [
            { "name", "string", "setting name" },
        ]
        OPTIONAL [
            { "value", "string", "new value (depends on setting type)", "" },
            { "selection", "string", "apply to specific selection", "global" },
            { "state", "integer", "state for state-specific settings", "0" },
        ]
        EXAMPLES [
            "set sphere_scale, 0.5",
            "set cartoon_color, red, chain A",
            "set cartoon_color, red",
        ]
    }

    fn execute<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        args: &ParsedCommand,
    ) -> CmdResult {
        // === 1. Parse arguments ===
        let raw_name = args
            .str_arg(0, "name")
            .ok_or_else(|| CmdError::missing_argument("name".to_string()))?;
        let resolved = resolve_legacy_setting_name(ctx, raw_name);
        let name = resolved.as_ref();

        let value_str = args
            .get_arg(1)
            .and_then(|v| v.to_string_repr())
            .or_else(|| args.get_named("value").and_then(|v| v.to_string_repr()));

        let selection = args.str_arg(2, "selection");

        // === 2. Intercept pseudo-settings ===
        if name == "state" {
            return self.handle_set_state(ctx, value_str);
        }

        let setting = ctx
            .resolve_setting(name)
            .ok_or_else(|| unknown_setting(name))?;
        self.execute_setting(ctx, args, setting, value_str, selection)
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

    command_help! {
        CMD "get"
        DESCRIPTION [
            "displays the current value of a setting.",
        ]
        USAGE [
            "get name [, selection [, state ]]",
        ]
        REQUIRED [
            { "name", "string", "setting name" },
        ]
        OPTIONAL [
            { "selection", "string", "get from specific selection", "global" },
            { "state", "integer", "state for state-specific settings", "0" },
        ]
        EXAMPLES [
            "get sphere_scale",
            "get antialias",
        ]
    }

    fn execute<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        args: &ParsedCommand,
    ) -> CmdResult {
        let raw_name = args
            .str_arg(0, "name")
            .ok_or_else(|| CmdError::missing_argument("name".to_string()))?;
        let resolved = resolve_legacy_setting_name(ctx, raw_name);
        let name = resolved.as_ref();

        // Intercept "get state" — report displayed state for each object
        if name == "state" {
            let obj_arg = args.str_arg(1, "selection");
            if let Some(obj_name) = obj_arg {
                if let Some(mol_obj) = ctx.viewer.objects().get_molecule(obj_name) {
                    ctx.print(&format!(
                        " state (int) = {} for \"{}\"",
                        mol_obj.display_state() + 1,
                        obj_name
                    ));
                } else {
                    return Err(CmdError::object_not_found(obj_name.to_string()));
                }
            } else {
                let names: Vec<String> = ctx
                    .viewer
                    .objects()
                    .names()
                    .map(|s| s.to_string())
                    .collect();
                for obj_name in &names {
                    if let Some(mol_obj) = ctx.viewer.objects().get_molecule(obj_name) {
                        ctx.print(&format!(
                            " state (int) = {} for \"{}\"",
                            mol_obj.display_state() + 1,
                            obj_name
                        ));
                    }
                }
            }
            return Ok(());
        }

        let setting = ctx
            .resolve_setting(name)
            .ok_or_else(|| unknown_setting(name))?;
        let selection = args.str_arg(1, "selection");

        if let Some(desc) = setting.built_in_descriptor() {
            if let Some(sel) = selection {
                if let Some(field_reader) = rep_color_field_read(desc.name) {
                    let selection_results = evaluate_selection(ctx.viewer, sel)?;
                    for (obj_name, selected) in &selection_results {
                        if let Some(mol) = ctx.viewer.objects().get_molecule(obj_name) {
                            for idx in selected.indices() {
                                if let Some(atom) = mol.molecule().get_atom(idx) {
                                    let color_val = field_reader(&atom.repr.colors);
                                    let value = SettingValue::Int(color_val);
                                    let display = setting.format_display(&value);
                                    ctx.print(&format!(" {} = {} ({})", name, display, sel));
                                    return Ok(());
                                }
                            }
                        }
                    }
                }

                if desc.is_object_overridable() {
                    if let Some(mol) = ctx.viewer.objects().get_molecule(sel) {
                        if let Some(value) = setting.built_in_object_value(mol) {
                            ctx.print(&format!(" {}", format_object_report(&setting, &value, sel)));
                            return Ok(());
                        }
                    }
                }
            }
        }

        if let Some(sel) = selection {
            if setting.is_object_overridable() {
                if let Some(value) = setting.object_value(sel).map_err(CmdError::execution)? {
                    ctx.print(&format!(" {}", format_object_report(&setting, &value, sel)));
                    return Ok(());
                }
            }
        }

        let value = setting
            .global_value(ctx.viewer)
            .map_err(CmdError::execution)?;
        ctx.print(&format!(" {}", setting.format_report(&value)));
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

    fn arg_hints(&self) -> &[ArgHint] {
        &[ArgHint::Setting, ArgHint::Selection]
    }

    command_help! {
        CMD "unset"
        DESCRIPTION [
            "restores a setting to its default value.",
        ]
        USAGE [
            "unset name [, selection [, state ]]",
        ]
        REQUIRED [
            { "name", "string", "setting name" },
        ]
        OPTIONAL [
            { "selection", "string", "apply to specific selection", "global" },
            { "state", "integer", "state for state-specific settings", "0" },
        ]
        EXAMPLES [
            "unset sphere_scale",
            "unset ray_trace_mode",
        ]
    }

    fn execute<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        args: &ParsedCommand,
    ) -> CmdResult {
        let raw_name = args
            .str_arg(0, "name")
            .ok_or_else(|| CmdError::missing_argument("name".to_string()))?;
        let resolved = resolve_legacy_setting_name(ctx, raw_name);
        let name = resolved.as_ref();

        let setting = ctx
            .resolve_setting(name)
            .ok_or_else(|| unknown_setting(name))?;
        let selection = args.str_arg(1, "selection");

        if let Some(desc) = setting.built_in_descriptor() {
            if let Some(sel) = selection {
                if let Some(field_accessor) = rep_color_field(desc.name) {
                    let selection_results = evaluate_selection(ctx.viewer, sel)?;
                    let mut total = 0usize;
                    for (obj_name, selected) in selection_results {
                        if selected.count() > 0 {
                            if let Some(mol_obj) =
                                ctx.viewer.objects_mut().get_molecule_mut(&obj_name)
                            {
                                for idx in selected.indices() {
                                    if let Some(atom) =
                                        mol_obj.molecule_mut().get_atom_mut(AtomIndex(idx.0))
                                    {
                                        *field_accessor(&mut atom.repr.colors) =
                                            patinae_mol::COLOR_UNSET;
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

                // Per-atom transparency unset. Surface uses recolor_lut
                // (COLOR); other reps need full rebuild (REPS).
                if let Some(field_accessor) = per_atom_transparency_field(desc.name) {
                    let unset_dirty = if desc.name == "surface_transparency" {
                        DirtyFlags::COLOR
                    } else {
                        DirtyFlags::REPS
                    };
                    let selection_results = evaluate_selection(ctx.viewer, sel)?;
                    let mut total = 0usize;
                    for (obj_name, selected) in selection_results {
                        if selected.count() > 0 {
                            if let Some(mol_obj) =
                                ctx.viewer.objects_mut().get_molecule_mut(&obj_name)
                            {
                                for idx in selected.indices() {
                                    if let Some(atom) =
                                        mol_obj.molecule_mut().get_atom_mut(AtomIndex(idx.0))
                                    {
                                        *field_accessor(&mut atom.repr) = None;
                                    }
                                }
                                total += selected.count();
                                mol_obj.invalidate(unset_dirty);
                            }
                        }
                    }
                    ctx.viewer.request_redraw();
                    if !ctx.quiet {
                        ctx.print(&format!(" {} reset for {} atoms", name, total));
                    }
                    return Ok(());
                }

                if desc.is_object_overridable() {
                    let inherited_value = desc.get(ctx.viewer.settings());
                    if let Some(mol) = ctx.viewer.objects_mut().get_molecule_mut(sel) {
                        if let Some(overrides) = mol.overrides_mut() {
                            desc.unset_override(overrides);
                        }
                        apply_object_side_effects(
                            mol,
                            desc.side_effects,
                            desc.name,
                            &inherited_value,
                        );
                        ctx.viewer.request_redraw();
                        if !ctx.quiet {
                            ctx.print(&format!(
                                " {} reset to global default (object {})",
                                name, sel
                            ));
                        }
                        return Ok(());
                    }
                }
            }

            if let Some(field_accessor) = rep_color_field(desc.name) {
                let names: Vec<_> = ctx
                    .viewer
                    .objects()
                    .names()
                    .map(|s| s.to_string())
                    .collect();
                for obj_name in names {
                    if let Some(mol_obj) = ctx.viewer.objects_mut().get_molecule_mut(&obj_name) {
                        for atom in mol_obj.molecule_mut().atoms_mut() {
                            *field_accessor(&mut atom.repr.colors) = patinae_mol::COLOR_UNSET;
                        }
                        mol_obj.invalidate(DirtyFlags::COLOR);
                    }
                }
            }
        }

        if let Some(sel) = selection {
            if setting.is_object_overridable()
                && setting
                    .unset_object_value(sel)
                    .map_err(CmdError::execution)?
            {
                ctx.viewer.request_redraw();
                if !ctx.quiet {
                    ctx.print(&format!(
                        " {} reset to global default (object {})",
                        name, sel
                    ));
                }
                return Ok(());
            }
        }

        let default_value = setting
            .unset_global(ctx.viewer)
            .map_err(CmdError::execution)?;
        apply_side_effects_from_slice(ctx, setting.side_effects(), name, &default_value);

        ctx.viewer.request_redraw();

        if !ctx.quiet {
            ctx.print(&format!(
                " {} reset to default: {}",
                name,
                setting.format_display(&default_value)
            ));
        }

        Ok(())
    }
}

// ============================================================================
// setup_ccd command
// ============================================================================

struct SetupCcdCommand;

impl Command for SetupCcdCommand {
    fn name(&self) -> &str {
        "setup_ccd"
    }

    command_help! {
        CMD "setup_ccd"
        DESCRIPTION [
            "loads Chemical Component Dictionary (CCD) bond templates and",
            "re-bonds all objects using template-based connectivity for",
            "known HETATM residues (HEM, ATP, NAD, etc.).",
        ]
        REQUIRED []
        OPTIONAL [
            { "mode", "string", "'scene' to load only templates for residues in the scene", "all" },
        ]
        EXAMPLES [
            "setup_ccd",
            "setup_ccd scene",
        ]
    }

    fn execute<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        args: &ParsedCommand,
    ) -> CmdResult {
        let mode = args.str_arg(0, "mode");

        // Load CCD cache
        let count = if mode == Some("scene") {
            // Collect unique HETATM residue names from all objects
            let mut resnames = std::collections::HashSet::new();
            let names: Vec<String> = ctx
                .viewer
                .objects()
                .names()
                .map(|s| s.to_string())
                .collect();
            for obj_name in &names {
                if let Some(mol_obj) = ctx.viewer.objects().get_molecule(obj_name) {
                    for atom in mol_obj.molecule().atoms_slice() {
                        if atom.state.hetatm {
                            resnames.insert(atom.residue.resn.clone());
                        }
                    }
                }
            }
            patinae_mol::ccd::load_cache_filtered(&resnames).map_err(CmdError::execution)?
        } else {
            patinae_mol::ccd::load_cache().map_err(CmdError::execution)?
        };
        let bond_tolerance = ctx.viewer.settings().behavior.bonding_vdw_cutoff;

        // Re-bond all objects
        let names: Vec<String> = ctx
            .viewer
            .objects()
            .names()
            .map(|s| s.to_string())
            .collect();
        let mut rebonded = 0;
        for obj_name in &names {
            if let Some(mol_obj) = ctx.viewer.objects_mut().get_molecule_mut(obj_name) {
                mol_obj.molecule_mut().rebond(bond_tolerance);
                mol_obj.invalidate(DirtyFlags::ALL);
                rebonded += 1;
            }
        }

        ctx.viewer.request_redraw();
        if !ctx.quiet {
            ctx.print(&format!(
                " CCD: loaded {} templates, re-bonded {} object(s)",
                count, rebonded
            ));
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use patinae_settings::registry;

    #[test]
    fn test_cartoon_smooth_loops_set_on_value_path() {
        let desc = registry::lookup_by_name("cartoon_smooth_loops").unwrap();
        let setting = ResolvedSetting::BuiltIn(desc);
        let value = setting.parse_value("on").unwrap();
        let mut settings = patinae_settings::Settings::default();

        desc.set(&mut settings, value).unwrap();

        assert!(settings.cartoon.smooth_loops);
    }

    #[test]
    fn test_object_side_effects_use_descriptor_dirty_contract() {
        let mol = patinae_mol::MoleculeBuilder::new("obj").build();
        let mut obj = MoleculeObject::with_name(mol, "obj");
        obj.clear_dirty();

        let desc = registry::lookup_by_name("cartoon_smooth_loops").unwrap();
        apply_object_side_effects(
            &mut obj,
            desc.side_effects,
            desc.name,
            &SettingValue::Bool(true),
        );
        assert!(obj.dirty_flags().contains(DirtyFlags::REPS));
        assert!(!obj.dirty_flags().contains(DirtyFlags::COLOR));

        obj.clear_dirty();
        let desc = registry::lookup_by_name("mesh_transparency").unwrap();
        apply_object_side_effects(
            &mut obj,
            desc.side_effects,
            desc.name,
            &SettingValue::Float(0.5),
        );
        assert!(obj.dirty_flags().contains(DirtyFlags::TRANSPARENCY));
    }
}
