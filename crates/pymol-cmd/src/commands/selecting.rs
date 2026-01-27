//! Selection commands: select, deselect, indicate

use crate::args::ParsedCommand;
use crate::command::{Command, CommandContext, CommandRegistry, ViewerLike};
use crate::error::{CmdError, CmdResult};

use pymol_mol::ObjectMolecule;
use pymol_select::{EvalContext, SelectionResult};

/// Register selection commands
pub fn register(registry: &mut CommandRegistry) {
    registry.register(SelectCommand);
    registry.register(DeselectCommand);
    registry.register(IndicateCommand);
}

// ============================================================================
// Selection evaluation helper
// ============================================================================

/// Evaluate a selection expression with support for named selections from the viewer.
///
/// This function:
/// 1. Gathers all molecules from the viewer's object registry
/// 2. Builds an EvalContext with all molecules
/// 3. Pre-evaluates stored named selections and adds them to the context
/// 4. Evaluates the target selection
/// 5. Returns the total count and a mapping of object name to SelectionResult
///
/// # Arguments
/// * `viewer` - The viewer providing objects and named selections
/// * `selection` - The selection expression to evaluate
///
/// # Returns
/// * Total atom count across all objects
/// * A vector of (object_name, SelectionResult) pairs
pub fn select_with_context(
    viewer: &dyn ViewerLike,
    selection: &str,
) -> CmdResult<(usize, Vec<(String, SelectionResult)>)> {
    // Gather all molecules from the registry
    let object_names: Vec<String> = viewer.objects().names().map(|s| s.to_string()).collect();
    
    // Collect molecule references
    let mut molecules: Vec<(&str, &ObjectMolecule)> = Vec::new();
    for name in &object_names {
        if let Some(mol_obj) = viewer.objects().get_molecule(name) {
            molecules.push((name.as_str(), mol_obj.molecule()));
        }
    }
    
    if molecules.is_empty() {
        // No molecules, return empty result
        return Ok((0, Vec::new()));
    }
    
    // Get the stored named selections
    let selection_names = viewer.selection_names();
    
    let mut total_count = 0;
    let mut results: Vec<(String, SelectionResult)> = Vec::new();
    
    // For each molecule, evaluate the selection in a context that includes named selections
    for (obj_name, mol) in &molecules {
        // Build context for this single molecule
        let mut ctx = EvalContext::single(mol);
        
        // Pre-evaluate and add named selections to the context
        // We evaluate each named selection against this molecule and add it to the context
        for sel_name in &selection_names {
            if let Some(sel_expr) = viewer.get_selection(sel_name) {
                // Parse and evaluate the named selection expression
                if let Ok(expr) = pymol_select::parse(sel_expr) {
                    if let Ok(result) = pymol_select::evaluate(&expr, &ctx) {
                        ctx.add_selection(sel_name.clone(), result);
                    }
                }
            }
        }
        
        // Now evaluate the target selection with named selections available
        match pymol_select::parse(selection) {
            Ok(expr) => {
                match pymol_select::evaluate(&expr, &ctx) {
                    Ok(result) => {
                        let count = result.count();
                        total_count += count;
                        results.push((obj_name.to_string(), result));
                    }
                    Err(e) => {
                        log::debug!("Selection evaluation error for {}: {:?}", obj_name, e);
                        // Continue with other molecules, don't fail completely
                    }
                }
            }
            Err(e) => {
                return Err(CmdError::invalid_arg("selection", format!("parse error: {:?}", e)));
            }
        }
    }
    
    Ok((total_count, results))
}

/// Evaluate a selection expression and return results per object.
///
/// This is a convenience wrapper that evaluates a selection and returns
/// a mapping from object names to their selection results.
///
/// # Arguments
/// * `viewer` - The viewer providing objects and named selections
/// * `selection` - The selection expression to evaluate
///
/// # Returns
/// * A vector of (object_name, SelectionResult) pairs for objects with matches
pub fn evaluate_selection(
    viewer: &dyn ViewerLike,
    selection: &str,
) -> CmdResult<Vec<(String, SelectionResult)>> {
    let (_, results) = select_with_context(viewer, selection)?;
    Ok(results)
}

// ============================================================================
// select command
// ============================================================================

struct SelectCommand;

impl Command for SelectCommand {
    fn name(&self) -> &str {
        "select"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "select" creates a named selection of atoms.

USAGE

    select name, selection [, enable [, quiet [, merge [, state ]]]]

ARGUMENTS

    name = string: selection name
    selection = string: selection expression
    enable = 0/1: enable the selection (default: 1)
    quiet = 0/1: suppress feedback (default: 0)
    merge = 0/1: merge with existing selection (default: 0)
    state = integer: state for selection (default: -1 = all)

EXAMPLES

    select backbone, name C+CA+N+O
    select chain_A, chain A
    select active_site, byres (organic around 4)
"#
    }

    fn execute<'a>(&self, ctx: &mut CommandContext<'a, dyn ViewerLike + 'a>, args: &ParsedCommand) -> CmdResult {
        let name = args
            .get_str(0)
            .or_else(|| args.get_named_str("name"))
            .ok_or_else(|| CmdError::MissingArgument("name".to_string()))?;

        let selection = args
            .get_str(1)
            .or_else(|| args.get_named_str("selection"))
            .ok_or_else(|| CmdError::MissingArgument("selection".to_string()))?;

        // Evaluate the selection to count atoms and validate it
        let (total_count, _results) = select_with_context(ctx.viewer, selection)?;

        // Store the selection expression in the viewer (creates with visible=true by default)
        ctx.viewer.define_selection(name, selection);
        
        // Ensure the selection is visible (show indicators)
        ctx.viewer.set_selection_visible(name, true);

        if !ctx.quiet {
            ctx.print(&format!(" Selector: selection \"{}\" defined with {} atoms.", name, total_count));
        }

        Ok(())
    }
}

// ============================================================================
// deselect command
// ============================================================================

struct DeselectCommand;

impl Command for DeselectCommand {
    fn name(&self) -> &str {
        "deselect"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "deselect" removes a named selection or clears all selections.

USAGE

    deselect [ selection ]

ARGUMENTS

    selection = string: name of selection to remove (default: all selections)

EXAMPLES

    deselect chainA
    deselect
"#
    }

    fn execute<'a>(&self, ctx: &mut CommandContext<'a, dyn ViewerLike + 'a>, args: &ParsedCommand) -> CmdResult {
        // Get optional selection name
        let name = args
            .get_str(0)
            .or_else(|| args.get_named_str("selection"));

        if let Some(name) = name {
            // Remove specific selection
            let existed = ctx.viewer.remove_selection(name);
            if !ctx.quiet {
                if existed {
                    ctx.print(&format!(" Selection \"{}\" deleted.", name));
                } else {
                    ctx.print_error(&format!(" Selection \"{}\" not found.", name));
                }
            }
        } else {
            // Remove all selections
            let names = ctx.viewer.selection_names();
            let count = names.len();
            for name in names {
                ctx.viewer.remove_selection(&name);
            }
            if !ctx.quiet {
                if count > 0 {
                    ctx.print(&format!(" {} selection(s) cleared.", count));
                } else {
                    ctx.print(" No selections to clear.");
                }
            }
        }

        Ok(())
    }
}

// ============================================================================
// indicate command
// ============================================================================

struct IndicateCommand;

impl Command for IndicateCommand {
    fn name(&self) -> &str {
        "indicate"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "indicate" shows a visual indicator (pink/magenta dots) for the specified
    selection in the 3D viewport. This helps visualize which atoms are 
    currently selected.

USAGE

    indicate [ selection ]

ARGUMENTS

    selection = string: selection to indicate (default: sele)
                Use "none" or empty to clear the indication.

EXAMPLES

    indicate chain A
    indicate organic
    indicate sele
    indicate none
"#
    }

    fn execute<'a>(&self, ctx: &mut CommandContext<'a, dyn ViewerLike + 'a>, args: &ParsedCommand) -> CmdResult {
        let selection = args
            .get_str(0)
            .or_else(|| args.get_named_str("selection"));

        match selection {
            Some("none") | Some("") => {
                // Clear indication when "none" or empty is specified
                ctx.viewer.clear_indication();
                if !ctx.quiet {
                    ctx.print(" Indication cleared.");
                }
            }
            None => {
                // No arguments - default to "sele"
                let sel = "sele";
                
                // Check if "sele" is a named selection that exists
                if ctx.viewer.get_selection(sel).is_some() {
                    // Make the existing selection visible
                    ctx.viewer.set_selection_visible(sel, true);
                    let (count, _) = select_with_context(ctx.viewer, sel)?;
                    if !ctx.quiet {
                        ctx.print(&format!(" Indicating \"{}\" ({} atoms)", sel, count));
                    }
                } else {
                    // "sele" doesn't exist - create it with "all" as default
                    let (count, _) = select_with_context(ctx.viewer, "all")?;
                    ctx.viewer.define_selection(sel, "all");
                    ctx.viewer.set_selection_visible(sel, true);
                    if !ctx.quiet {
                        ctx.print(&format!(" Indicating \"{}\" ({} atoms)", sel, count));
                    }
                }
            }
            Some(sel) => {
                // Check if this is an existing named selection
                if ctx.viewer.get_selection(sel).is_some() {
                    // Make the existing selection visible
                    ctx.viewer.set_selection_visible(sel, true);
                    let (count, _) = select_with_context(ctx.viewer, sel)?;
                    if !ctx.quiet {
                        ctx.print(&format!(" Indicating \"{}\" ({} atoms)", sel, count));
                    }
                } else {
                    // It's a selection expression - create/update "indicate" selection
                    let (count, _) = select_with_context(ctx.viewer, sel)?;
                    ctx.viewer.define_selection("indicate", sel);
                    ctx.viewer.set_selection_visible("indicate", true);
                    if !ctx.quiet {
                        ctx.print(&format!(" Indicating \"{}\" ({} atoms)", sel, count));
                    }
                }
            }
        }

        Ok(())
    }
}
