//! Structure transformation commands: translate, rotate, transform_selection
//!
//! These commands modify atomic coordinates (unlike viewing commands which only
//! affect the camera).
//!
//! Full selection expression support:
//! - Property selectors: `name CA`, `resn ALA`, `chain A`, `elem C`
//! - Numeric comparisons: `b > 50`, `resi 1-100`
//! - Logical operators: `name CA and chain A`, `organic or solvent`
//! - Special keywords: `backbone`, `sidechain`, `polymer`, `organic`

use lin_alg::f32::{Mat4, Vec3};
use pymol_mol::{rotation_ttt, ttt_to_mat4, AtomIndex};
use pymol_scene::{DirtyFlags, Object};

use crate::args::ParsedCommand;
use crate::command::{Command, CommandContext, CommandRegistry, ViewerLike};
use crate::helpers::{camera_to_model_vec, state_index_from_user};
use crate::commands::selecting::evaluate_selection;
use crate::error::{CmdError, CmdResult};

/// Register transformation commands
pub fn register(registry: &mut CommandRegistry) {
    registry.register(TranslateCommand);
    registry.register(RotateCommand);
    registry.register(TransformSelectionCommand);
}

// ============================================================================
// translate command
// ============================================================================

struct TranslateCommand;

impl Command for TranslateCommand {
    fn name(&self) -> &str {
        "translate"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "translate" translates the atomic coordinates of atoms in a selection.
    Supports full selection expressions.

USAGE

    translate vector [, selection [, state [, camera ]]]

ARGUMENTS

    vector = float vector: translation vector [x, y, z]

    selection = string: atoms whose coordinates should be modified {default: all}
        Supports full selection expressions:
        - Property selectors: name CA, resn ALA, chain A, elem C
        - Numeric ranges: resi 1-100, b > 50
        - Logical operators: name CA and chain A, organic or solvent
        - Keywords: backbone, sidechain, polymer, organic, solvent

    state > 0: only the indicated state is modified

    state = 0: all states are modified

    state = -1: only the current state is modified {default}

    camera = 0 or 1: is the vector in camera coordinates? {default: 1 (yes)}

EXAMPLES

    translate [1, 0, 0], name CA
    translate [0, 5, 0], chain A
    translate [0, 0, 10], backbone and chain A
    translate [1, 1, 1], resi 50-100
    translate [0, 0, 5], organic
"#
    }

    fn execute<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        args: &ParsedCommand,
    ) -> CmdResult {
        // Parse vector argument (required)
        let vector = parse_vector(args, 0, "vector")?;

        let selection = args.str_arg_or(1, "selection", "all");
        let state = args.int_arg_or(2, "state", -1);
        let camera = args.int_arg_or(3, "camera", 1);

        // Convert vector from camera coordinates if needed
        let shift = if camera != 0 {
            // Transform vector from camera coordinates to model coordinates
            let rotation = &ctx.viewer.camera().current_view().rotation;
            camera_to_model_vec(rotation, vector)
        } else {
            vector
        };

        // When movie is active and selection is a whole object name,
        // modify the object's TTT transform matrix instead of atom coordinates.
        // This allows mview store/interpolation to animate the object.
        let movie_active = !ctx.viewer.movie().is_empty();
        let is_object = ctx.viewer.objects().contains(selection);

        if movie_active && is_object {
            // Build translation matrix (column-major: translation at [12,13,14])
            let translation = Mat4 {
                data: [
                    1.0, 0.0, 0.0, 0.0,
                    0.0, 1.0, 0.0, 0.0,
                    0.0, 0.0, 1.0, 0.0,
                    shift.x, shift.y, shift.z, 1.0,
                ],
            };
            if let Some(mol_obj) = ctx.viewer.objects_mut().get_molecule_mut(selection) {
                let current = mol_obj.state().transform.clone();
                mol_obj.state_mut().set_transform(translation * current);
                mol_obj.invalidate(DirtyFlags::COORDS);
            }
            ctx.viewer.request_redraw();
            if !ctx.quiet {
                ctx.print(&format!(
                    " Translated object \"{}\" by [{:.3}, {:.3}, {:.3}] (object matrix)",
                    selection, shift.x, shift.y, shift.z
                ));
            }
        } else {
            // Apply translation to atom coordinates directly
            let (atom_count, obj_count) =
                apply_selection_transform(ctx, selection, state, |mol, atoms, st| {
                    apply_translation_to_atoms(mol, st, atoms, &shift);
                })?;

            ctx.viewer.request_redraw();

            if !ctx.quiet {
                ctx.print(&format!(
                    " Translated {} atom(s) in {} object(s) by [{:.3}, {:.3}, {:.3}]",
                    atom_count, obj_count, shift.x, shift.y, shift.z
                ));
            }
        }

        Ok(())
    }
}

// ============================================================================
// rotate command
// ============================================================================

struct RotateCommand;

impl Command for RotateCommand {
    fn name(&self) -> &str {
        "rotate"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "rotate" rotates the atomic coordinates of atoms in a selection about
    an axis. Supports full PyMOL selection expressions.

USAGE

    rotate axis, angle [, selection [, state [, camera [, origin ]]]]

ARGUMENTS

    axis = x, y, z, or float vector: axis about which to rotate

    angle = float: degrees of rotation

    selection = string: atoms whose coordinates should be modified {default: all}
        Supports full selection expressions:
        - Property selectors: name CA, resn ALA, chain A, elem C
        - Numeric ranges: resi 1-100, b > 50
        - Logical operators: name CA and chain A
        - Keywords: backbone, sidechain, polymer, organic

    state > 0: only the indicated state is modified

    state = 0: all states are modified

    state = -1: only the current state is modified {default}

    camera = 0 or 1: is the axis in camera coordinates? {default: 1 (yes)}

    origin = float vector: center of rotation {default: view origin}

EXAMPLES

    rotate x, 45, all
    rotate y, 90, chain A
    rotate [1, 1, 1], 10, backbone
    rotate z, 180, resi 100-200, origin=[0, 0, 0]
    rotate x, 45, organic
"#
    }

    fn execute<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        args: &ParsedCommand,
    ) -> CmdResult {
        // Parse axis argument (required)
        let axis = parse_axis(args, 0)?;

        let angle_deg = args.float_arg(1, "angle")
            .ok_or_else(|| CmdError::MissingArgument("angle".to_string()))?;
        let angle = (angle_deg as f32) * std::f32::consts::PI / 180.0;

        let selection = args.str_arg_or(2, "selection", "all");
        let state = args.int_arg_or(3, "state", -1);
        let camera = args.int_arg_or(4, "camera", 1);

        // Parse origin (default: view origin)
        let origin = if let Some(origin_vec) = args
            .get_named("origin")
            .and_then(|v| v.as_vec3())
        {
            origin_vec
        } else {
            // Use view origin as default
            ctx.viewer.camera().current_view().origin
        };

        // Transform axis from camera coordinates if needed
        let rot_axis = if camera != 0 {
            let rotation = &ctx.viewer.camera().current_view().rotation;
            camera_to_model_vec(rotation, axis)
        } else {
            axis
        };

        // Build the TTT matrix for rotation about origin
        let ttt = rotation_ttt(rot_axis, angle, origin);

        // Apply rotation using selection expressions
        let (atom_count, obj_count) =
            apply_selection_transform(ctx, selection, state, |mol, atoms, st| {
                apply_ttt_to_atoms(mol, st, atoms, &ttt);
            })?;

        ctx.viewer.request_redraw();

        if !ctx.quiet {
            ctx.print(&format!(
                " Rotated {} atom(s) in {} object(s) by {:.1} degrees",
                atom_count, obj_count, angle_deg
            ));
        }

        Ok(())
    }
}

// ============================================================================
// transform_selection command
// ============================================================================

struct TransformSelectionCommand;

impl Command for TransformSelectionCommand {
    fn name(&self) -> &str {
        "transform_selection"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "transform_selection" applies a transformation matrix to the atomic
    coordinates of a selection. Supports full PyMOL selection expressions.

USAGE

    transform_selection selection, matrix [, state [, homogenous ]]

ARGUMENTS

    selection = string: atoms to transform {default: all}
        Supports full selection expressions:
        - Property selectors: name CA, resn ALA, chain A
        - Logical operators: name CA and chain A
        - Keywords: backbone, sidechain, polymer, organic

    matrix = list of 16 floats: transformation matrix

    state > 0: only the indicated state is modified

    state = 0: all states are modified

    state = -1: only the current state is modified {default}

    homogenous = 0/1: matrix format {default: 0}
        0 = TTT format (pre-translate, rotate, post-translate)
        1 = Standard homogenous 4x4 matrix

NOTES

    When homogenous=0, the matrix is in TTT format:
    - [0-2, 4-6, 8-10]: 3x3 rotation matrix
    - [3, 7, 11]: post-rotation translation
    - [12, 13, 14]: pre-rotation translation
    - [15]: 1.0

EXAMPLES

    # Translate all atoms by [1, 2, 3]
    transform_selection all, [1,0,0,1, 0,1,0,2, 0,0,1,3, 0,0,0,1], homogenous=1
    
    # Transform only chain A
    transform_selection chain A, [1,0,0,5, 0,1,0,0, 0,0,1,0, 0,0,0,1], homogenous=1
    
    # Transform backbone atoms
    transform_selection backbone, [1,0,0,0, 0,1,0,0, 0,0,1,10, 0,0,0,1], homogenous=1
"#
    }

    fn execute<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        args: &ParsedCommand,
    ) -> CmdResult {
        let selection = args.str_arg_or(0, "selection", "all");
        let matrix = parse_matrix(args, 1)?;
        let state = args.int_arg_or(2, "state", -1);
        let homogenous = args.int_arg_or(3, "homogenous", 0);

        // Apply transformation using selection expressions
        let (atom_count, obj_count) = if homogenous != 0 {
            let mat4 = ttt_to_mat4(&matrix);
            apply_selection_transform(ctx, selection, state, |mol, atoms, st| {
                apply_matrix_to_atoms(mol, st, atoms, &mat4);
            })?
        } else {
            apply_selection_transform(ctx, selection, state, |mol, atoms, st| {
                apply_ttt_to_atoms(mol, st, atoms, &matrix);
            })?
        };

        ctx.viewer.request_redraw();

        if !ctx.quiet {
            ctx.print(&format!(
                " Transformed {} atom(s) in {} object(s)",
                atom_count, obj_count
            ));
        }

        Ok(())
    }
}

// ============================================================================
// Helper functions - Argument parsing
// ============================================================================

/// Parse a vector argument (either [x,y,z] list or x,y,z positional)
fn parse_vector(args: &ParsedCommand, pos: usize, name: &str) -> Result<Vec3, CmdError> {
    args.vec3_arg(pos, name)
        .ok_or_else(|| CmdError::MissingArgument(name.to_string()))
}

/// Parse axis argument (x, y, z, or [ax, ay, az])
fn parse_axis(args: &ParsedCommand, pos: usize) -> Result<Vec3, CmdError> {
    // Check for named axes first
    if let Some(s) = args.get_str(pos) {
        match s.to_lowercase().as_str() {
            "x" => return Ok(Vec3::new(1.0, 0.0, 0.0)),
            "y" => return Ok(Vec3::new(0.0, 1.0, 0.0)),
            "z" => return Ok(Vec3::new(0.0, 0.0, 1.0)),
            _ => {}
        }
    }

    // Try as vector
    parse_vector(args, pos, "axis")
}

/// Parse a 16-element matrix
fn parse_matrix(args: &ParsedCommand, pos: usize) -> Result<[f32; 16], CmdError> {
    let mut matrix = [0.0f32; 16];

    if let Some(arg) = args.get_arg(pos) {
        match arg {
            crate::args::ArgValue::List(items) if items.len() >= 16 => {
                for (i, item) in items.iter().take(16).enumerate() {
                    matrix[i] = item.as_float().unwrap_or(0.0) as f32;
                }
                return Ok(matrix);
            }
            crate::args::ArgValue::String(s) => {
                let s = s.trim().trim_matches(|c| c == '[' || c == ']' || c == '(');
                let parts: Vec<f32> = s
                    .split(',')
                    .filter_map(|p| p.trim().parse::<f32>().ok())
                    .collect();
                if parts.len() >= 16 {
                    for (i, &v) in parts.iter().take(16).enumerate() {
                        matrix[i] = v;
                    }
                    return Ok(matrix);
                }
            }
            _ => {}
        }
    }

    // Try named argument
    if let Some(crate::args::ArgValue::List(items)) = args.get_named("matrix") {
        if items.len() >= 16 {
            for (i, item) in items.iter().take(16).enumerate() {
                matrix[i] = item.as_float().unwrap_or(0.0) as f32;
            }
            return Ok(matrix);
        }
    }

    Err(CmdError::MissingArgument("matrix (16 floats)".to_string()))
}

// ============================================================================
// Helper functions - Selection-based transformations
// ============================================================================

/// Apply a transformation to atoms matching a selection expression.
///
/// Uses `evaluate_selection()` which supports full selection expressions,
/// named selections, object names, and wildcard patterns (e.g. `1fsd_*`).
///
/// Returns (total_atoms_transformed, num_objects_affected)
fn apply_selection_transform<F>(
    ctx: &mut CommandContext<'_, '_, dyn ViewerLike + '_>,
    selection: &str,
    state: i64,
    transform_fn: F,
) -> CmdResult<(usize, usize)>
where
    F: Fn(&mut pymol_mol::ObjectMolecule, &[AtomIndex], i64),
{
    let selection_results = evaluate_selection(ctx.viewer, selection)?;

    let mut total_atoms = 0usize;
    let mut affected_objects = 0usize;

    for (obj_name, selected) in selection_results {
        let atoms: Vec<AtomIndex> = selected.indices().collect();
        if atoms.is_empty() {
            continue;
        }

        if let Some(mol_obj) = ctx.viewer.objects_mut().get_molecule_mut(&obj_name) {
            // Resolve state=-1 (current state) to the object's display state
            let resolved_state = if state == -1 {
                mol_obj.display_state() as i64 + 1 // convert 0-indexed to 1-indexed
            } else {
                state
            };
            transform_fn(mol_obj.molecule_mut(), &atoms, resolved_state);
            total_atoms += atoms.len();
            affected_objects += 1;
        }
    }

    if affected_objects == 0 {
        return Err(CmdError::Selection(format!(
            "No atoms matching '{}'",
            selection
        )));
    }

    Ok((total_atoms, affected_objects))
}

// ============================================================================
// Helper functions - Atom-level transformations
// ============================================================================

/// Collect which 0-indexed states to apply a transformation to.
///
/// - `0` (or negative) → all states
/// - Positive N → just state N-1 (if valid)
fn resolve_state_indices(state: i64, num_states: usize) -> Vec<usize> {
    match state_index_from_user(state) {
        None => (0..num_states).collect(),
        Some(idx) if idx < num_states => vec![idx],
        Some(_) => Vec::new(),
    }
}

/// Apply translation to specific atoms based on state parameter
fn apply_translation_to_atoms(
    mol: &mut pymol_mol::ObjectMolecule,
    state: i64,
    atoms: &[AtomIndex],
    delta: &Vec3,
) {
    for i in resolve_state_indices(state, mol.state_count()) {
        mol.translate_atoms(i, atoms, *delta);
    }
}

/// Apply TTT matrix transform to specific atoms based on state parameter
fn apply_ttt_to_atoms(
    mol: &mut pymol_mol::ObjectMolecule,
    state: i64,
    atoms: &[AtomIndex],
    ttt: &[f32; 16],
) {
    for i in resolve_state_indices(state, mol.state_count()) {
        mol.transform_ttt_atoms(i, atoms, ttt);
    }
}

/// Apply standard 4x4 matrix transform to specific atoms based on state parameter
fn apply_matrix_to_atoms(
    mol: &mut pymol_mol::ObjectMolecule,
    state: i64,
    atoms: &[AtomIndex],
    matrix: &lin_alg::f32::Mat4,
) {
    for i in resolve_state_indices(state, mol.state_count()) {
        mol.transform_atoms(i, atoms, matrix);
    }
}
