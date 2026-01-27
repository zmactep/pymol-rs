//! Display commands: show, hide, enable, disable, color, bg_color

use pymol_mol::RepMask;
use pymol_scene::{DirtyFlags, Object};
use pymol_select::AtomIndex;

use crate::args::ParsedCommand;
use crate::command::{Command, CommandContext, CommandRegistry, ViewerLike};
use crate::commands::selecting::evaluate_selection;
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
    registry.register(BgColorCommand);
}

/// Parse a representation name into a RepMask value
fn parse_rep(name: &str) -> Option<u32> {
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
        "everything" | "all" => Some(RepMask::ALL.0), // Extract inner u32
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

    fn execute<'a>(&self, ctx: &mut CommandContext<'a, dyn ViewerLike + 'a>, args: &ParsedCommand) -> CmdResult {
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
            RepMask::ALL.0
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
                            atom.visible_reps.set_visible(rep);
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

    fn execute<'a>(&self, ctx: &mut CommandContext<'a, dyn ViewerLike + 'a>, args: &ParsedCommand) -> CmdResult {
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
            RepMask::ALL.0
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
                            if rep == RepMask::ALL.0 {
                                atom.visible_reps = RepMask::NONE;
                            } else {
                                atom.visible_reps.set_hidden(rep);
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

    fn execute<'a>(&self, ctx: &mut CommandContext<'a, dyn ViewerLike + 'a>, args: &ParsedCommand) -> CmdResult {
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
                            atom.visible_reps = RepMask::NONE;
                            atom.visible_reps.set_visible(rep);
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

    fn execute<'a>(&self, ctx: &mut CommandContext<'a, dyn ViewerLike + 'a>, args: &ParsedCommand) -> CmdResult {
        let name = args
            .get_str(0)
            .or_else(|| args.get_named_str("name"))
            .unwrap_or("all");

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

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "disable" makes objects invisible.

USAGE

    disable [ name ]

ARGUMENTS

    name = string: object name pattern (default: all)

EXAMPLES

    disable
    disable protein
    disable obj*
"#
    }

    fn execute<'a>(&self, ctx: &mut CommandContext<'a, dyn ViewerLike + 'a>, args: &ParsedCommand) -> CmdResult {
        let name = args
            .get_str(0)
            .or_else(|| args.get_named_str("name"))
            .unwrap_or("all");

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

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "toggle" toggles object visibility.

USAGE

    toggle name

ARGUMENTS

    name = string: object name

EXAMPLES

    toggle protein
"#
    }

    fn execute<'a>(&self, ctx: &mut CommandContext<'a, dyn ViewerLike + 'a>, args: &ParsedCommand) -> CmdResult {
        let name = args
            .get_str(0)
            .or_else(|| args.get_named_str("name"))
            .ok_or_else(|| CmdError::MissingArgument("name".to_string()))?;

        // Get current state
        let currently_enabled = ctx
            .viewer
            .objects()
            .get(name)
            .map(|o| o.is_enabled())
            .unwrap_or(false);

        let _ = ctx.viewer.objects_mut().enable(name, !currently_enabled);

        ctx.viewer.request_redraw();

        if !ctx.quiet {
            if currently_enabled {
                ctx.print(&format!(" Disabled \"{}\"", name));
            } else {
                ctx.print(&format!(" Enabled \"{}\"", name));
            }
        }

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

    fn execute<'a>(&self, ctx: &mut CommandContext<'a, dyn ViewerLike + 'a>, args: &ParsedCommand) -> CmdResult {
        let color_name = args
            .get_str(0)
            .or_else(|| args.get_named_str("color"))
            .ok_or_else(|| CmdError::MissingArgument("color".to_string()))?;

        let selection = args
            .get_str(1)
            .or_else(|| args.get_named_str("selection"))
            .unwrap_or("all");

        // Check for special color schemes first (these use negative indices)
        // -1: by element (atomic/CPK)
        // -2: by chain
        // -3: by secondary structure
        // -4: by B-factor
        let color_index: i32 = match color_name.to_lowercase().as_str() {
            "atomic" | "cpk" | "element" | "by_element" => -1,
            "chain" | "by_chain" | "chainbow" => -2,
            "ss" | "secondary_structure" | "by_ss" | "dssp" => -3,
            "b" | "b_factor" | "bfactor" | "by_b" => -4,
            _ => {
                // Look up as a named color
                ctx.viewer
                    .color_index(color_name)
                    .map(|idx| idx as i32)
                    .ok_or_else(|| CmdError::invalid_arg("color", format!("unknown color: {}", color_name)))?
            }
        };

        // Check if the selection targets a specific non-cartoon/ribbon representation
        // If so, we should preserve cartoon/ribbon colors when changing atom.color
        let selection_lower = selection.to_lowercase();
        let preserve_cartoon_color = selection_lower.contains("rep ")
            && !selection_lower.contains("rep cartoon")
            && !selection_lower.contains("rep cart")
            && !selection_lower.contains("rep ribbon")
            && !selection_lower.contains("rep ribb");

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
                            // If preserving cartoon/ribbon colors and the atom has these reps visible,
                            // save the current color to the rep-specific field before changing atom.color
                            if preserve_cartoon_color {
                                if atom.cartoon_color.is_none() && atom.visible_reps.is_visible(RepMask::CARTOON) {
                                    atom.cartoon_color = Some(atom.color);
                                }
                                if atom.ribbon_color.is_none() && atom.visible_reps.is_visible(RepMask::RIBBON) {
                                    atom.ribbon_color = Some(atom.color);
                                }
                            }
                            atom.color = color_index;
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

    fn execute<'a>(&self, ctx: &mut CommandContext<'a, dyn ViewerLike + 'a>, args: &ParsedCommand) -> CmdResult {
        let color_name = args
            .get_str(0)
            .or_else(|| args.get_named_str("color"))
            .ok_or_else(|| CmdError::MissingArgument("color".to_string()))?;

        // Parse color name to RGB
        let (r, g, b) = match color_name.to_lowercase().as_str() {
            "white" => (1.0, 1.0, 1.0),
            "black" => (0.0, 0.0, 0.0),
            "gray" | "grey" => (0.5, 0.5, 0.5),
            "red" => (1.0, 0.0, 0.0),
            "green" => (0.0, 1.0, 0.0),
            "blue" => (0.0, 0.0, 1.0),
            "yellow" => (1.0, 1.0, 0.0),
            "cyan" => (0.0, 1.0, 1.0),
            "magenta" => (1.0, 0.0, 1.0),
            "orange" => (1.0, 0.5, 0.0),
            _ => {
                // Try to parse as hex color
                if color_name.starts_with("0x") || color_name.starts_with('#') {
                    let hex = color_name.trim_start_matches("0x").trim_start_matches('#');
                    if hex.len() == 6 {
                        let r = u8::from_str_radix(&hex[0..2], 16).unwrap_or(0) as f32 / 255.0;
                        let g = u8::from_str_radix(&hex[2..4], 16).unwrap_or(0) as f32 / 255.0;
                        let b = u8::from_str_radix(&hex[4..6], 16).unwrap_or(0) as f32 / 255.0;
                        (r, g, b)
                    } else {
                        return Err(CmdError::invalid_arg(
                            "color",
                            format!("invalid hex color: {}", color_name),
                        ));
                    }
                } else {
                    return Err(CmdError::invalid_arg(
                        "color",
                        format!("unknown color: {}", color_name),
                    ));
                }
            }
        };

        ctx.viewer.set_background_color(r, g, b);

        if !ctx.quiet {
            ctx.print(&format!(" Background color set to {}", color_name));
        }

        Ok(())
    }
}
