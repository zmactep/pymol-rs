//! Display commands: show, hide, enable, disable, color, set_color, bg_color, label

use pymol_mol::{three_to_one, Atom, RepMask};
use pymol_scene::{DirtyFlags, Object};
use pymol_select::AtomIndex;

use crate::args::ParsedCommand;
use crate::command::{ArgHint, Command, CommandContext, CommandRegistry, ViewerLike};
use crate::helpers::{
    for_each_selected_molecule_mut, resolve_object_names, set_enabled_with_group_awareness,
    ResolvedNames,
};
use crate::error::{CmdError, CmdResult};

/// Register display commands
pub fn register(registry: &mut CommandRegistry) {
    registry.register(ShowCommand);
    registry.register(HideCommand);
    registry.register(ShowAsCommand);
    registry.register(EnableCommand);
    registry.register(DisableCommand);
    registry.register(ToggleCommand);
    registry.register(ColorCommand);
    registry.register(SetColorCommand);
    registry.register(BgColorCommand);
    registry.register(LabelCommand);
}

/// Parse a representation name into a RepMask value
fn parse_rep(name: &str) -> Option<RepMask> {
    match name.to_lowercase().as_str() {
        "lines" | "line" => Some(RepMask::LINES),
        "sticks" | "stick" => Some(RepMask::STICKS),
        "spheres" | "sphere" => Some(RepMask::SPHERES),
        "surface" | "surf" => Some(RepMask::SURFACE),
        "mesh" => Some(RepMask::MESH),
        "dots" | "dot" => Some(RepMask::DOTS),
        "cartoon" | "cart" => Some(RepMask::CARTOON),
        "ribbon" | "ribb" => Some(RepMask::RIBBON),
        "labels" | "label" => Some(RepMask::LABELS),
        "nonbonded" | "nb_spheres" => Some(RepMask::NONBONDED),
        "cell" => Some(RepMask::CELL),
        "cgo" => Some(RepMask::CGO),
        "callback" => Some(RepMask::CALLBACK),
        "extent" => Some(RepMask::EXTENT),
        "slice" => Some(RepMask::SLICE),
        "everything" | "all" => Some(RepMask::ALL),
        _ => None,
    }
}

// ============================================================================
// show command
// ============================================================================

struct ShowCommand;

impl Command for ShowCommand {
    fn name(&self) -> &str {
        "show"
    }

    fn arg_hints(&self) -> &[ArgHint] {
        &[ArgHint::Representation, ArgHint::Selection]
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "show" makes representations visible.

USAGE

    show [ representation [, selection ]]

ARGUMENTS

    representation = string: lines, sticks, spheres, surface, mesh, dots,
                            cartoon, ribbon, labels, etc.
    selection = string: atoms to show (default: all)

EXAMPLES

    show
    show cartoon
    show sticks, organic
    show surface, chain A
"#
    }

    fn execute<'v, 'r>(&self, ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>, args: &ParsedCommand) -> CmdResult {
        let rep_name = args.str_arg(0, "representation");
        let selection = args.str_arg_or(1, "selection", "all");

        // If no representation specified, show all
        let rep = if let Some(name) = rep_name {
            parse_rep(name).ok_or_else(|| {
                CmdError::invalid_arg("representation", format!("unknown representation: {}", name))
            })?
        } else {
            RepMask::ALL
        };

        let total_affected = for_each_selected_molecule_mut(
            ctx.viewer, selection, DirtyFlags::REPS,
            |mol_obj, selected| {
                mol_obj.state_mut().visible_reps.set_visible(rep);
                let mol_mut = mol_obj.molecule_mut();
                for idx in selected.indices() {
                    if let Some(atom) = mol_mut.get_atom_mut(AtomIndex(idx.0)) {
                        atom.repr.visible_reps.set_visible(rep);
                    }
                }
            },
        )?;

        ctx.viewer.request_redraw();

        if !ctx.quiet {
            if total_affected == 0 {
                ctx.print_error(&format!(" Show: selection \"{}\" not found", selection));
            } else if let Some(name) = rep_name {
                ctx.print(&format!(" Showing {}", name));
            } else {
                ctx.print(" Showing all representations");
            }
        }

        Ok(())
    }
}

// ============================================================================
// hide command
// ============================================================================

struct HideCommand;

impl Command for HideCommand {
    fn name(&self) -> &str {
        "hide"
    }

    fn arg_hints(&self) -> &[ArgHint] {
        &[ArgHint::Representation, ArgHint::Selection]
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "hide" makes representations invisible.

USAGE

    hide [ representation [, selection ]]

ARGUMENTS

    representation = string: lines, sticks, spheres, surface, etc.
    selection = string: atoms to hide (default: all)

EXAMPLES

    hide
    hide lines
    hide sticks, all
"#
    }

    fn execute<'v, 'r>(&self, ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>, args: &ParsedCommand) -> CmdResult {
        let rep_name = args.str_arg(0, "representation");
        let selection = args.str_arg_or(1, "selection", "all");

        // If no representation specified, hide all
        let rep = if let Some(name) = rep_name {
            parse_rep(name).ok_or_else(|| {
                CmdError::invalid_arg("representation", format!("unknown representation: {}", name))
            })?
        } else {
            RepMask::ALL
        };

        let total_affected = for_each_selected_molecule_mut(
            ctx.viewer, selection, DirtyFlags::REPS,
            |mol_obj, selected| {
                let mol_mut = mol_obj.molecule_mut();
                for idx in selected.indices() {
                    if let Some(atom) = mol_mut.get_atom_mut(AtomIndex(idx.0)) {
                        if rep == RepMask::ALL {
                            atom.repr.visible_reps = RepMask::NONE;
                        } else {
                            atom.repr.visible_reps.set_hidden(rep);
                        }
                    }
                }
            },
        )?;

        ctx.viewer.request_redraw();

        if !ctx.quiet {
            if total_affected == 0 {
                ctx.print_error(&format!(" Hide: selection \"{}\" not found", selection));
            } else if let Some(name) = rep_name {
                ctx.print(&format!(" Hiding {}", name));
            } else {
                ctx.print(" Hiding all representations");
            }
        }

        Ok(())
    }
}

// ============================================================================
// show_as command (alias: as)
// ============================================================================

struct ShowAsCommand;

impl Command for ShowAsCommand {
    fn name(&self) -> &str {
        "show_as"
    }

    fn arg_hints(&self) -> &[ArgHint] {
        &[ArgHint::Representation, ArgHint::Selection]
    }

    fn aliases(&self) -> &[&str] {
        &["as"]
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "show_as" (or "as") hides all representations and shows only the specified one.

USAGE

    as representation [, selection ]

ARGUMENTS

    representation = string: lines, sticks, spheres, surface, cartoon, etc.
    selection = string: atoms to show (default: all)

EXAMPLES

    as cartoon
    as sticks, organic
    as surface, polymer
"#
    }

    fn execute<'v, 'r>(&self, ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>, args: &ParsedCommand) -> CmdResult {
        let rep_name = args.str_arg(0, "representation")
            .ok_or_else(|| CmdError::MissingArgument("representation".to_string()))?;
        let selection = args.str_arg_or(1, "selection", "all");

        let rep = parse_rep(rep_name).ok_or_else(|| {
            CmdError::invalid_arg("representation", format!("unknown representation: {}", rep_name))
        })?;

        let total_affected = for_each_selected_molecule_mut(
            ctx.viewer, selection, DirtyFlags::REPS,
            |mol_obj, selected| {
                mol_obj.state_mut().visible_reps.set_visible(rep);
                let mol_mut = mol_obj.molecule_mut();
                for idx in selected.indices() {
                    if let Some(atom) = mol_mut.get_atom_mut(AtomIndex(idx.0)) {
                        atom.repr.visible_reps = RepMask::NONE;
                        atom.repr.visible_reps.set_visible(rep);
                    }
                }
            },
        )?;

        ctx.viewer.request_redraw();

        if !ctx.quiet {
            if total_affected == 0 {
                ctx.print_error(&format!(" Show as: selection \"{}\" not found", selection));
            } else {
                ctx.print(&format!(" Showing as {}", rep_name));
            }
        }

        Ok(())
    }
}

// ============================================================================
// enable command
// ============================================================================

struct EnableCommand;

impl Command for EnableCommand {
    fn name(&self) -> &str {
        "enable"
    }

    fn arg_hints(&self) -> &[ArgHint] {
        &[ArgHint::Selection]
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "enable" makes objects visible.

USAGE

    enable [ name ]

ARGUMENTS

    name = string: object name pattern (default: all)

EXAMPLES

    enable
    enable protein
    enable obj*
"#
    }

    fn execute<'v, 'r>(&self, ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>, args: &ParsedCommand) -> CmdResult {
        let name = args.str_arg_or(0, "name", "all");
        set_visibility(ctx, name, true)
    }
}

// ============================================================================
// disable command
// ============================================================================

struct DisableCommand;

impl Command for DisableCommand {
    fn name(&self) -> &str {
        "disable"
    }

    fn arg_hints(&self) -> &[ArgHint] {
        &[ArgHint::Selection]
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "disable" makes objects or selections invisible.

USAGE

    disable [ name ]

ARGUMENTS

    name = string: object or selection name pattern (default: all)

EXAMPLES

    disable
    disable protein
    disable sele
"#
    }

    fn execute<'v, 'r>(&self, ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>, args: &ParsedCommand) -> CmdResult {
        let name = args.str_arg_or(0, "name", "all");
        set_visibility(ctx, name, false)
    }
}

/// Shared implementation for enable/disable commands.
fn set_visibility(
    ctx: &mut CommandContext<'_, '_, dyn ViewerLike + '_>,
    name: &str,
    enabled: bool,
) -> CmdResult {
    match resolve_object_names(ctx.viewer.objects(), name) {
        ResolvedNames::All => {
            let all_names: Vec<String> = ctx.viewer.objects().names().map(|s| s.to_string()).collect();
            for obj_name in &all_names {
                set_enabled_with_group_awareness(ctx.viewer.objects_mut(), obj_name, enabled);
            }
        }
        ResolvedNames::Matched(names) => {
            for obj_name in &names {
                set_enabled_with_group_awareness(ctx.viewer.objects_mut(), obj_name, enabled);
            }
        }
        ResolvedNames::Unresolved => {}
    }

    let matching_sels: Vec<String> = ctx.viewer.selections().matching(name)
        .iter().map(|s| s.to_string()).collect();
    for sel_name in &matching_sels {
        ctx.viewer.selections_mut().set_visible(sel_name, enabled);
    }

    ctx.viewer.request_redraw();

    if !ctx.quiet {
        let verb = if enabled { "Enabled" } else { "Disabled" };
        ctx.print(&format!(" {} \"{}\"", verb, name));
    }

    Ok(())
}

// ============================================================================
// toggle command
// ============================================================================

struct ToggleCommand;

impl Command for ToggleCommand {
    fn name(&self) -> &str {
        "toggle"
    }

    fn arg_hints(&self) -> &[ArgHint] {
        &[ArgHint::Selection]
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "toggle" toggles visibility of objects or selections.

USAGE

    toggle name

ARGUMENTS

    name = string: object or selection name

EXAMPLES

    toggle protein
    toggle sele
"#
    }

    fn execute<'v, 'r>(&self, ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>, args: &ParsedCommand) -> CmdResult {
        let name = args.str_arg(0, "name")
            .ok_or_else(|| CmdError::MissingArgument("name".to_string()))?;

        // Try object first
        if let Some(obj) = ctx.viewer.objects().get(name) {
            let currently_enabled = obj.is_enabled();
            set_enabled_with_group_awareness(ctx.viewer.objects_mut(), name, !currently_enabled);

            if !ctx.quiet {
                if currently_enabled {
                    ctx.print(&format!(" Disabled \"{}\"", name));
                } else {
                    ctx.print(&format!(" Enabled \"{}\"", name));
                }
            }
        } else if ctx.viewer.selections().names().contains(&name.to_string()) {
            // Try selection
            let currently_visible = ctx.viewer.selections().is_visible(name);
            ctx.viewer.selections_mut().set_visible(name, !currently_visible);

            if !ctx.quiet {
                if currently_visible {
                    ctx.print(&format!(" Disabled selection \"{}\"", name));
                } else {
                    ctx.print(&format!(" Enabled selection \"{}\"", name));
                }
            }
        } else {
            return Err(CmdError::ObjectNotFound(name.to_string()));
        }

        ctx.viewer.request_redraw();

        Ok(())
    }
}

// ============================================================================
// color command
// ============================================================================

struct ColorCommand;

impl Command for ColorCommand {
    fn name(&self) -> &str {
        "color"
    }

    fn arg_hints(&self) -> &[ArgHint] {
        &[ArgHint::Color, ArgHint::Selection]
    }

    fn aliases(&self) -> &[&str] {
        &["colour"]
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "color" sets the color of atoms or objects.

USAGE

    color color [, selection ]

ARGUMENTS

    color = string: color name or special scheme:
        Named colors: red, green, blue, yellow, cyan, magenta, orange, white, gray, etc.
        Special schemes:
            atomic (cpk, element) - color by element type
            chain (chainbow) - color by chain
            ss (secondary_structure) - color by secondary structure
            b (b_factor, bfactor) - color by B-factor
            residue (residue_type, aa_type) - color by residue type
            index (residue_index, rainbow) - color by residue index
    selection = string: atoms to color (default: all)

EXAMPLES

    color red
    color green, chain A
    color cyan, organic
    color atomic
    color chain
    color ss, polymer
"#
    }

    fn execute<'v, 'r>(&self, ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>, args: &ParsedCommand) -> CmdResult {
        let color_name = args.str_arg(0, "color")
            .ok_or_else(|| CmdError::MissingArgument("color".to_string()))?;
        let selection = args.str_arg_or(1, "selection", "all");

        let color_index: i32 = if let Some(ci) = pymol_color::ColorIndex::from_scheme_name(color_name) {
            i32::from(ci)
        } else {
            ctx.viewer
                .color_index(color_name)
                .map(|idx| idx as i32)
                .ok_or_else(|| CmdError::invalid_arg("color", format!("unknown color: {}", color_name)))?
        };

        let total_colored = for_each_selected_molecule_mut(
            ctx.viewer, selection, DirtyFlags::COLOR,
            |mol_obj, selected| {
                let mol_mut = mol_obj.molecule_mut();
                for idx in selected.indices() {
                    if let Some(atom) = mol_mut.get_atom_mut(AtomIndex(idx.0)) {
                        atom.repr.colors.base = color_index;
                        atom.repr.colors.cartoon = color_index;
                        atom.repr.colors.ribbon = color_index;
                        atom.repr.colors.stick = color_index;
                        atom.repr.colors.line = color_index;
                        atom.repr.colors.sphere = color_index;
                        atom.repr.colors.surface = color_index;
                        atom.repr.colors.mesh = color_index;
                    }
                }
            },
        )?;

        ctx.viewer.request_redraw();

        if !ctx.quiet {
            if total_colored == 0 {
                ctx.print_error(&format!(" Color: 0 atoms colored {} (selection not found)", color_name));
            } else {
                ctx.print(&format!(" Color: {} atoms colored {}", total_colored, color_name));
            }
        }

        Ok(())
    }
}

// ============================================================================
// bg_color command
// ============================================================================

struct BgColorCommand;

impl Command for BgColorCommand {
    fn name(&self) -> &str {
        "bg_color"
    }

    fn arg_hints(&self) -> &[ArgHint] {
        &[ArgHint::Color]
    }

    fn aliases(&self) -> &[&str] {
        &["bg_colour", "background"]
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "bg_color" sets the background color.

USAGE

    bg_color color

ARGUMENTS

    color = string: color name (white, black, gray, etc.)

EXAMPLES

    bg_color white
    bg_color black
    bg_color gray
"#
    }

    fn execute<'v, 'r>(&self, ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>, args: &ParsedCommand) -> CmdResult {
        // Try [r, g, b] vector first (from ArgValue::List)
        if let Some(crate::args::ArgValue::List(items)) = args.get_arg(0) {
            if items.len() == 3 {
                if let (Some(r), Some(g), Some(b)) = (
                    items[0].as_float().map(|v| v as f32),
                    items[1].as_float().map(|v| v as f32),
                    items[2].as_float().map(|v| v as f32),
                ) {
                    ctx.viewer.set_background_color(
                        r.clamp(0.0, 1.0),
                        g.clamp(0.0, 1.0),
                        b.clamp(0.0, 1.0),
                    );
                    if !ctx.quiet {
                        ctx.print(&format!(" Background color set to [{:.2}, {:.2}, {:.2}]", r, g, b));
                    }
                    return Ok(());
                }
            }
        }

        let color_name = args.str_arg(0, "color")
            .ok_or_else(|| CmdError::MissingArgument("color".to_string()))?;

        // Resolve color: named colors registry, then hex
        let color = if let Some((_, color)) = ctx.viewer.named_colors().get_by_name(color_name) {
            color
        } else if let Some(color) = pymol_color::Color::from_hex(color_name) {
            color
        } else {
            return Err(CmdError::invalid_arg(
                "color",
                format!("unknown color: {}", color_name),
            ));
        };

        let [r, g, b] = color.to_array();
        ctx.viewer.set_background_color(r, g, b);

        if !ctx.quiet {
            ctx.print(&format!(" Background color set to {}", color_name));
        }

        Ok(())
    }
}

// ============================================================================
// set_color command
// ============================================================================

struct SetColorCommand;

impl Command for SetColorCommand {
    fn name(&self) -> &str {
        "set_color"
    }

    fn aliases(&self) -> &[&str] {
        &["set_colour"]
    }

    fn arg_hints(&self) -> &[ArgHint] {
        &[ArgHint::Color]
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "set_color" defines a new named color or removes an existing one.

USAGE

    set_color name, [ r, g, b ]
    set_color name

ARGUMENTS

    name = string: the color name to define or remove
    [r, g, b] = list of integers (0-255): the RGB color value

NOTES

    If only a name is provided, the named color is removed.
    If a name and RGB list are provided, the named color is created or updated.

EXAMPLES

    set_color mywhite, [255, 255, 255]
    set_color darkred, [128, 0, 0]
    set_color mywhite
"#
    }

    fn execute<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        args: &ParsedCommand,
    ) -> CmdResult {
        let color_name = args.str_arg(0, "name")
            .ok_or_else(|| CmdError::MissingArgument("name".to_string()))?;

        if let Some(crate::args::ArgValue::List(items)) = args.get_arg(1) {
            if items.len() != 3 {
                return Err(CmdError::invalid_arg(
                    "rgb",
                    format!("expected [r, g, b] (3 values), got {} values", items.len()),
                ));
            }

            let r = items[0]
                .as_int()
                .ok_or_else(|| CmdError::invalid_arg("r", "expected an integer"))? as u8;
            let g = items[1]
                .as_int()
                .ok_or_else(|| CmdError::invalid_arg("g", "expected an integer"))? as u8;
            let b = items[2]
                .as_int()
                .ok_or_else(|| CmdError::invalid_arg("b", "expected an integer"))? as u8;

            let color = pymol_color::Color::from_rgb8(r, g, b);
            let idx = ctx.viewer.named_colors_mut().set(color_name, color);

            ctx.viewer.request_redraw();

            if !ctx.quiet {
                ctx.print(&format!(
                    " Color: \"{}\" defined as [{}, {}, {}] (index {})",
                    color_name, r, g, b, idx
                ));
            }
        } else {
            let removed = ctx.viewer.named_colors_mut().unregister(color_name);

            if !ctx.quiet {
                if removed {
                    ctx.print(&format!(" Color: \"{}\" removed", color_name));
                } else {
                    ctx.print_warning(&format!(
                        " Color: \"{}\" not found (nothing removed)",
                        color_name
                    ));
                }
            }
        }

        Ok(())
    }
}

// ============================================================================
// label command
// ============================================================================

/// Label expression — what property to display as label text
enum LabelExpression {
    Name,
    Resn,
    Resi,
    Chain,
    /// Occupancy (PyMOL convention: q = occupancy)
    Q,
    /// B-factor
    B,
    Segi,
    /// "ATOM" or "HETATM"
    Type,
    FormalCharge,
    PartialCharge,
    /// Element symbol (e.g., "C", "N", "O")
    Elem,
    /// Van der Waals radius
    Vdw,
    /// One-letter amino acid code
    Oneletter,
    StringLiteral(String),
}

/// Parse a label expression string into a LabelExpression
fn parse_label_expr(s: &str) -> Result<LabelExpression, CmdError> {
    let trimmed = s.trim();

    // Check for quoted string literal
    if (trimmed.starts_with('"') && trimmed.ends_with('"'))
        || (trimmed.starts_with('\'') && trimmed.ends_with('\''))
    {
        let inner = &trimmed[1..trimmed.len() - 1];
        return Ok(LabelExpression::StringLiteral(inner.to_string()));
    }

    match trimmed.to_lowercase().as_str() {
        "name" => Ok(LabelExpression::Name),
        "resn" => Ok(LabelExpression::Resn),
        "resi" => Ok(LabelExpression::Resi),
        "chain" => Ok(LabelExpression::Chain),
        "q" => Ok(LabelExpression::Q),
        "b" => Ok(LabelExpression::B),
        "segi" => Ok(LabelExpression::Segi),
        "type" => Ok(LabelExpression::Type),
        "formal_charge" => Ok(LabelExpression::FormalCharge),
        "partial_charge" => Ok(LabelExpression::PartialCharge),
        "elem" | "element" => Ok(LabelExpression::Elem),
        "vdw" => Ok(LabelExpression::Vdw),
        "oneletter" | "one_letter" => Ok(LabelExpression::Oneletter),
        // Anything else is a string literal (quotes are stripped by the command parser)
        _ => Ok(LabelExpression::StringLiteral(trimmed.to_string())),
    }
}

/// Evaluate a label expression for a given atom
fn eval_label_expr(expr: &LabelExpression, atom: &Atom) -> String {
    match expr {
        LabelExpression::Name => atom.name.to_string(),
        LabelExpression::Resn => atom.residue.resn.clone(),
        LabelExpression::Resi => {
            if atom.residue.inscode != ' ' {
                format!("{}{}", atom.residue.resv, atom.residue.inscode)
            } else {
                atom.residue.resv.to_string()
            }
        }
        LabelExpression::Chain => atom.residue.chain.clone(),
        LabelExpression::Q => format!("{:.2}", atom.occupancy),
        LabelExpression::B => format!("{:.2}", atom.b_factor),
        LabelExpression::Segi => atom.residue.segi.clone(),
        LabelExpression::Type => {
            if atom.state.hetatm {
                "HETATM".to_string()
            } else {
                "ATOM".to_string()
            }
        }
        LabelExpression::FormalCharge => atom.formal_charge.to_string(),
        LabelExpression::PartialCharge => format!("{:.4}", atom.partial_charge),
        LabelExpression::Elem => atom.element.symbol().to_string(),
        LabelExpression::Vdw => format!("{:.2}", atom.effective_vdw()),
        LabelExpression::Oneletter => three_to_one(&atom.residue.resn)
            .map(|c| c.to_string())
            .unwrap_or_else(|| atom.residue.resn.clone()),
        LabelExpression::StringLiteral(s) => s.clone(),
    }
}

struct LabelCommand;

impl Command for LabelCommand {
    fn name(&self) -> &str {
        "label"
    }

    fn arg_hints(&self) -> &[ArgHint] {
        &[ArgHint::Selection, ArgHint::LabelProperty]
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "label" sets text labels on atoms based on a property expression.

USAGE

    label [ selection [, expression ]]

ARGUMENTS

    selection = string: atoms to label (default: all)
    expression = string: property to display:
        name           - atom name
        resn           - residue name
        resi           - residue number/identifier
        chain          - chain identifier
        q              - occupancy
        b              - B-factor
        segi           - segment identifier
        type           - ATOM or HETATM
        formal_charge  - formal charge
        partial_charge - partial charge
        "string"       - literal string

    If no expression is given, labels are cleared.

EXAMPLES

    label all, name
    label chain A, resn
    label organic, resi
    label sele, "hello"
    label              # clear all labels
"#
    }

    fn execute<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        args: &ParsedCommand,
    ) -> CmdResult {
        let selection = args.str_arg_or(0, "selection", "all");
        let expr_str = args.str_arg(1, "expression");

        if let Some(expr_str) = expr_str {
            // Set labels
            let expr = parse_label_expr(expr_str)?;

            let total_affected = for_each_selected_molecule_mut(
                ctx.viewer, selection, DirtyFlags::REPS,
                |mol_obj, selected| {
                    mol_obj.state_mut().visible_reps.set_visible(RepMask::LABELS);
                    let mol_mut = mol_obj.molecule_mut();
                    for idx in selected.indices() {
                        if let Some(atom) = mol_mut.get_atom_mut(AtomIndex(idx.0)) {
                            atom.repr.label = eval_label_expr(&expr, atom);
                            atom.repr.visible_reps.set_visible(RepMask::LABELS);
                        }
                    }
                },
            )?;

            if !ctx.quiet {
                ctx.print(&format!(" Label: {} atoms labeled", total_affected));
            }
        } else {
            // Clear labels
            let total_affected = for_each_selected_molecule_mut(
                ctx.viewer, selection, DirtyFlags::REPS,
                |mol_obj, selected| {
                    let mol_mut = mol_obj.molecule_mut();
                    for idx in selected.indices() {
                        if let Some(atom) = mol_mut.get_atom_mut(AtomIndex(idx.0)) {
                            atom.repr.label.clear();
                            atom.repr.visible_reps.set_hidden(RepMask::LABELS);
                        }
                    }
                },
            )?;

            if !ctx.quiet {
                ctx.print(&format!(" Label: {} atoms unlabeled", total_affected));
            }
        }

        ctx.viewer.request_redraw();

        Ok(())
    }
}
