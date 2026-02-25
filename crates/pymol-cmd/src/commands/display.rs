//! Display commands: show, hide, enable, disable, color, bg_color, label

use pymol_mol::{Atom, RepMask};
use pymol_scene::{DirtyFlags, Object};
use pymol_select::AtomIndex;

use crate::args::ParsedCommand;
use crate::command::{ArgHint, Command, CommandContext, CommandRegistry, ViewerLike};
use crate::commands::selecting::evaluate_selection;
use crate::error::{CmdError, CmdResult};

/// Simple glob matching (prefix* and *suffix patterns)
fn glob_match(pattern: &str, name: &str) -> bool {
    if let Some(prefix) = pattern.strip_suffix('*') {
        return name.starts_with(prefix);
    }
    if let Some(suffix) = pattern.strip_prefix('*') {
        return name.ends_with(suffix);
    }
    pattern == name
}

/// Register display commands
pub fn register(registry: &mut CommandRegistry) {
    registry.register(ShowCommand);
    registry.register(HideCommand);
    registry.register(ShowAsCommand);
    registry.register(EnableCommand);
    registry.register(DisableCommand);
    registry.register(ToggleCommand);
    registry.register(ColorCommand);
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
        let rep_name = args
            .get_str(0)
            .or_else(|| args.get_named_str("representation"));

        let selection = args
            .get_str(1)
            .or_else(|| args.get_named_str("selection"))
            .unwrap_or("all");

        // If no representation specified, show all
        let rep = if let Some(name) = rep_name {
            parse_rep(name).ok_or_else(|| {
                CmdError::invalid_arg("representation", format!("unknown representation: {}", name))
            })?
        } else {
            RepMask::ALL
        };

        // Evaluate selection with named selection support
        let selection_results = evaluate_selection(ctx.viewer, selection)?;

        let mut total_affected = 0usize;

        // Apply to atoms matching the selection in each object
        for (obj_name, selected) in selection_results {
            let count = selected.count();
            if count > 0 {
                if let Some(mol_obj) = ctx.viewer.objects_mut().get_molecule_mut(&obj_name) {
                    // Set object-level visibility for this rep
                    mol_obj.state_mut().visible_reps.set_visible(rep);
                    // Get mutable access and show rep on selected atoms
                    let mol_mut = mol_obj.molecule_mut();
                    for idx in selected.indices() {
                        if let Some(atom) = mol_mut.get_atom_mut(AtomIndex(idx.0)) {
                            atom.repr.visible_reps.set_visible(rep);
                        }
                    }
                    mol_obj.invalidate(DirtyFlags::REPS);
                    total_affected += count;
                }
            }
        }

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
        let rep_name = args
            .get_str(0)
            .or_else(|| args.get_named_str("representation"));

        let selection = args
            .get_str(1)
            .or_else(|| args.get_named_str("selection"))
            .unwrap_or("all");

        // If no representation specified, hide all
        let rep = if let Some(name) = rep_name {
            parse_rep(name).ok_or_else(|| {
                CmdError::invalid_arg("representation", format!("unknown representation: {}", name))
            })?
        } else {
            RepMask::ALL
        };

        // Evaluate selection with named selection support
        let selection_results = evaluate_selection(ctx.viewer, selection)?;

        let mut total_affected = 0usize;

        // Apply to atoms matching the selection in each object
        for (obj_name, selected) in selection_results {
            let count = selected.count();
            if count > 0 {
                if let Some(mol_obj) = ctx.viewer.objects_mut().get_molecule_mut(&obj_name) {
                    // Get mutable access and hide rep on selected atoms
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
                    mol_obj.invalidate(DirtyFlags::REPS);
                    total_affected += count;
                }
            }
        }

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
        let rep_name = args
            .get_str(0)
            .or_else(|| args.get_named_str("representation"))
            .ok_or_else(|| CmdError::MissingArgument("representation".to_string()))?;

        let selection = args
            .get_str(1)
            .or_else(|| args.get_named_str("selection"))
            .unwrap_or("all");

        let rep = parse_rep(rep_name).ok_or_else(|| {
            CmdError::invalid_arg("representation", format!("unknown representation: {}", rep_name))
        })?;

        // Evaluate selection with named selection support
        let selection_results = evaluate_selection(ctx.viewer, selection)?;

        let mut total_affected = 0usize;

        // Apply to atoms matching the selection in each object
        for (obj_name, selected) in selection_results {
            let count = selected.count();
            if count > 0 {
                if let Some(mol_obj) = ctx.viewer.objects_mut().get_molecule_mut(&obj_name) {
                    // Set object-level visibility for this rep
                    mol_obj.state_mut().visible_reps.set_visible(rep);
                    // Get mutable access: hide all reps, then show the specified one on selected atoms
                    let mol_mut = mol_obj.molecule_mut();
                    for idx in selected.indices() {
                        if let Some(atom) = mol_mut.get_atom_mut(AtomIndex(idx.0)) {
                            atom.repr.visible_reps = RepMask::NONE;
                            atom.repr.visible_reps.set_visible(rep);
                        }
                    }
                    mol_obj.invalidate(DirtyFlags::REPS);
                    total_affected += count;
                }
            }
        }

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
        let name = args
            .get_str(0)
            .or_else(|| args.get_named_str("name"))
            .unwrap_or("all");

        // Enable matching objects
        let object_names: Vec<String> = ctx
            .viewer
            .objects()
            .matching(name)
            .iter()
            .map(|s| s.to_string())
            .collect();

        for obj_name in &object_names {
            let _ = ctx.viewer.objects_mut().enable(obj_name, true);
        }

        // Enable matching selections
        let sel_names = ctx.viewer.selections().names();
        for sel_name in &sel_names {
            if name == "all" || sel_name == name || glob_match(name, sel_name) {
                ctx.viewer.selections_mut().set_visible(sel_name, true);
            }
        }

        ctx.viewer.request_redraw();

        if !ctx.quiet {
            ctx.print(&format!(" Enabled \"{}\"", name));
        }

        Ok(())
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
        let name = args
            .get_str(0)
            .or_else(|| args.get_named_str("name"))
            .unwrap_or("all");

        // Disable matching objects
        let object_names: Vec<String> = ctx
            .viewer
            .objects()
            .matching(name)
            .iter()
            .map(|s| s.to_string())
            .collect();

        for obj_name in &object_names {
            let _ = ctx.viewer.objects_mut().enable(obj_name, false);
        }

        // Disable matching selections
        let sel_names = ctx.viewer.selections().names();
        for sel_name in &sel_names {
            if name == "all" || sel_name == name || glob_match(name, sel_name) {
                ctx.viewer.selections_mut().set_visible(sel_name, false);
            }
        }

        ctx.viewer.request_redraw();

        if !ctx.quiet {
            ctx.print(&format!(" Disabled \"{}\"", name));
        }

        Ok(())
    }
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
        let name = args
            .get_str(0)
            .or_else(|| args.get_named_str("name"))
            .ok_or_else(|| CmdError::MissingArgument("name".to_string()))?;

        let mut found = false;

        // Try object first
        if let Some(obj) = ctx.viewer.objects().get(name) {
            let currently_enabled = obj.is_enabled();
            let _ = ctx.viewer.objects_mut().enable(name, !currently_enabled);
            found = true;

            if !ctx.quiet {
                if currently_enabled {
                    ctx.print(&format!(" Disabled \"{}\"", name));
                } else {
                    ctx.print(&format!(" Enabled \"{}\"", name));
                }
            }
        }

        // Try selection
        if !found {
            let currently_visible = ctx.viewer.selections().is_visible(name);
            if ctx.viewer.selections().names().contains(&name.to_string()) {
                ctx.viewer.selections_mut().set_visible(name, !currently_visible);
                found = true;

                if !ctx.quiet {
                    if currently_visible {
                        ctx.print(&format!(" Disabled selection \"{}\"", name));
                    } else {
                        ctx.print(&format!(" Enabled selection \"{}\"", name));
                    }
                }
            }
        }

        if !found {
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
        let color_name = args
            .get_str(0)
            .or_else(|| args.get_named_str("color"))
            .ok_or_else(|| CmdError::MissingArgument("color".to_string()))?;

        let selection = args
            .get_str(1)
            .or_else(|| args.get_named_str("selection"))
            .unwrap_or("all");

        let color_index: i32 = if let Some(ci) = pymol_color::ColorIndex::from_scheme_name(color_name) {
            i32::from(ci)
        } else {
            ctx.viewer
                .color_index(color_name)
                .map(|idx| idx as i32)
                .ok_or_else(|| CmdError::invalid_arg("color", format!("unknown color: {}", color_name)))?
        };

        // Evaluate selection with named selection support
        let selection_results = evaluate_selection(ctx.viewer, selection)?;

        let mut total_colored = 0usize;

        // Apply color to atoms matching the selection in each object
        for (obj_name, selected) in selection_results {
            let count = selected.count();
            if count > 0 {
                if let Some(mol_obj) = ctx.viewer.objects_mut().get_molecule_mut(&obj_name) {
                    // Get mutable access and set colors on selected atoms
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
                    total_colored += count;
                    // Mark the molecule as needing color rebuild
                    mol_obj.invalidate(DirtyFlags::COLOR);
                }
            }
        }

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

        let color_name = args
            .get_str(0)
            .or_else(|| args.get_named_str("color"))
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
        let selection = args
            .get_str(0)
            .or_else(|| args.get_named_str("selection"))
            .unwrap_or("all");

        let expr_str = args
            .get_str(1)
            .or_else(|| args.get_named_str("expression"));

        let selection_results = evaluate_selection(ctx.viewer, selection)?;

        let mut total_affected = 0usize;

        if let Some(expr_str) = expr_str {
            // Set labels
            let expr = parse_label_expr(expr_str)?;

            for (obj_name, selected) in selection_results {
                let count = selected.count();
                if count > 0 {
                    if let Some(mol_obj) = ctx.viewer.objects_mut().get_molecule_mut(&obj_name) {
                        mol_obj.state_mut().visible_reps.set_visible(RepMask::LABELS);
                        let mol_mut = mol_obj.molecule_mut();
                        for idx in selected.indices() {
                            if let Some(atom) = mol_mut.get_atom_mut(AtomIndex(idx.0)) {
                                atom.repr.label = eval_label_expr(&expr, atom);
                                atom.repr.visible_reps.set_visible(RepMask::LABELS);
                            }
                        }
                        mol_obj.invalidate(DirtyFlags::REPS);
                        total_affected += count;
                    }
                }
            }

            if !ctx.quiet {
                ctx.print(&format!(" Label: {} atoms labeled", total_affected));
            }
        } else {
            // Clear labels
            for (obj_name, selected) in selection_results {
                let count = selected.count();
                if count > 0 {
                    if let Some(mol_obj) = ctx.viewer.objects_mut().get_molecule_mut(&obj_name) {
                        let mol_mut = mol_obj.molecule_mut();
                        for idx in selected.indices() {
                            if let Some(atom) = mol_mut.get_atom_mut(AtomIndex(idx.0)) {
                                atom.repr.label.clear();
                                atom.repr.visible_reps.set_hidden(RepMask::LABELS);
                            }
                        }
                        mol_obj.invalidate(DirtyFlags::REPS);
                        total_affected += count;
                    }
                }
            }

            if !ctx.quiet {
                ctx.print(&format!(" Label: {} atoms unlabeled", total_affected));
            }
        }

        ctx.viewer.request_redraw();

        Ok(())
    }
}
