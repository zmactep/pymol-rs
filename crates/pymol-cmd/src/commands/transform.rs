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

use lin_alg::f32::Vec3;
use pymol_mol::{rotation_ttt, ttt_to_mat4, AtomIndex};

use crate::args::ParsedCommand;
use crate::command::{Command, CommandContext, CommandRegistry, ViewerLike};
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
    Supports full PyMOL selection expressions.

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

        // Parse selection (default: "all")
        let selection = args
            .get_str(1)
            .or_else(|| args.get_named_str("selection"))
            .unwrap_or("all");

        // Parse state (-1=current, 0=all, >0=specific)
        let state = args
            .get_int(2)
            .or_else(|| args.get_named_int("state"))
            .unwrap_or(-1);

        // Parse camera flag (default: 1 = use camera coordinates)
        let camera = args
            .get_int(3)
            .or_else(|| args.get_named_int("camera"))
            .unwrap_or(1);

        // Convert vector from camera coordinates if needed
        let shift = if camera != 0 {
            // Transform vector from camera coordinates to model coordinates
            let view = ctx.viewer.camera().current_view();
            let r = &view.rotation;
            // Multiply by transpose of rotation matrix (inverse for orthogonal matrix)
            Vec3::new(
                r.data[0] * vector.x + r.data[4] * vector.y + r.data[8] * vector.z,
                r.data[1] * vector.x + r.data[5] * vector.y + r.data[9] * vector.z,
                r.data[2] * vector.x + r.data[6] * vector.y + r.data[10] * vector.z,
            )
        } else {
            vector
        };

        // Apply translation using selection expressions
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

        // Parse angle argument (required)
        let angle_deg = args
            .get_float(1)
            .or_else(|| args.get_named_float("angle"))
            .ok_or_else(|| CmdError::MissingArgument("angle".to_string()))?;
        let angle = (angle_deg as f32) * std::f32::consts::PI / 180.0;

        // Parse selection (default: "all")
        let selection = args
            .get_str(2)
            .or_else(|| args.get_named_str("selection"))
            .unwrap_or("all");

        // Parse state (-1=current, 0=all, >0=specific)
        let state = args
            .get_int(3)
            .or_else(|| args.get_named_int("state"))
            .unwrap_or(-1);

        // Parse camera flag (default: 1 = use camera coordinates)
        let camera = args
            .get_int(4)
            .or_else(|| args.get_named_int("camera"))
            .unwrap_or(1);

        // Parse origin (default: view origin)
        let origin = if let Some(origin_vec) = args
            .get_named("origin")
            .and_then(|v| parse_vector_from_arg(v))
        {
            origin_vec
        } else {
            // Use view origin as default
            ctx.viewer.camera().current_view().origin.clone()
        };

        // Transform axis from camera coordinates if needed
        let rot_axis = if camera != 0 {
            let view = ctx.viewer.camera().current_view();
            let r = &view.rotation;
            Vec3::new(
                r.data[0] * axis.x + r.data[4] * axis.y + r.data[8] * axis.z,
                r.data[1] * axis.x + r.data[5] * axis.y + r.data[9] * axis.z,
                r.data[2] * axis.x + r.data[6] * axis.y + r.data[10] * axis.z,
            )
        } else {
            axis
        };

        // Build the TTT matrix for rotation about origin
        let ttt = rotation_ttt(rot_axis, angle, origin.clone());

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
        0 = PyMOL TTT format (pre-translate, rotate, post-translate)
        1 = Standard homogenous 4x4 matrix

NOTES

    When homogenous=0, the matrix is in PyMOL's TTT format:
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
        // Parse selection
        let selection = args
            .get_str(0)
            .or_else(|| args.get_named_str("selection"))
            .unwrap_or("all");

        // Parse matrix (16 floats)
        let matrix = parse_matrix(args, 1)?;

        // Parse state (-1=current, 0=all, >0=specific)
        let state = args
            .get_int(2)
            .or_else(|| args.get_named_int("state"))
            .unwrap_or(-1);

        // Parse homogenous flag (default: 0 = TTT format)
        let homogenous = args
            .get_int(3)
            .or_else(|| args.get_named_int("homogenous"))
            .unwrap_or(0);

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
    // Try as a list first
    if let Some(arg) = args.get_arg(pos) {
        if let Some(vec) = parse_vector_from_arg(arg) {
            return Ok(vec);
        }
    }

    // Try as named argument
    if let Some(arg) = args.get_named(name) {
        if let Some(vec) = parse_vector_from_arg(arg) {
            return Ok(vec);
        }
    }

    Err(CmdError::MissingArgument(name.to_string()))
}

/// Parse a vector from an argument value
fn parse_vector_from_arg(arg: &crate::args::ArgValue) -> Option<Vec3> {
    use crate::args::ArgValue;

    match arg {
        ArgValue::List(items) if items.len() >= 3 => {
            let x = items.get(0)?.as_float()? as f32;
            let y = items.get(1)?.as_float()? as f32;
            let z = items.get(2)?.as_float()? as f32;
            Some(Vec3::new(x, y, z))
        }
        ArgValue::String(s) => {
            // Try to parse "[x, y, z]" format
            let s = s.trim().trim_matches(|c| c == '[' || c == ']');
            let parts: Vec<&str> = s.split(',').collect();
            if parts.len() >= 3 {
                let x = parts[0].trim().parse::<f32>().ok()?;
                let y = parts[1].trim().parse::<f32>().ok()?;
                let z = parts[2].trim().parse::<f32>().ok()?;
                Some(Vec3::new(x, y, z))
            } else {
                None
            }
        }
        _ => None,
    }
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
    if let Some(arg) = args.get_named("matrix") {
        if let crate::args::ArgValue::List(items) = arg {
            if items.len() >= 16 {
                for (i, item) in items.iter().take(16).enumerate() {
                    matrix[i] = item.as_float().unwrap_or(0.0) as f32;
                }
                return Ok(matrix);
            }
        }
    }

    Err(CmdError::MissingArgument("matrix (16 floats)".to_string()))
}

// ============================================================================
// Helper functions - Selection-based transformations
// ============================================================================

/// Apply a transformation to atoms matching a selection expression
///
/// This function:
/// 1. Iterates over all molecule objects
/// 2. Evaluates the selection expression against each molecule
/// 3. Applies the transform function to matching atoms
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
    let mut total_atoms = 0usize;
    let mut affected_objects = 0usize;

    // Collect molecule names first (to avoid borrow issues)
    let mol_names: Vec<String> = ctx
        .viewer
        .objects()
        .names()
        .filter_map(|name| {
            if ctx.viewer.objects().get_molecule(name).is_some() {
                Some(name.to_string())
            } else {
                None
            }
        })
        .collect();

    // Process each molecule
    for name in mol_names {
        // Get molecule reference for selection evaluation
        let atoms: Vec<AtomIndex> = {
            let mol_obj = match ctx.viewer.objects().get_molecule(&name) {
                Some(m) => m,
                None => continue,
            };
            let mol = mol_obj.molecule();

            // Evaluate selection expression
            match pymol_select::select_atoms(mol, selection) {
                Ok(indices) => indices,
                Err(e) => {
                    // If selection fails, it might be an object name - try pattern matching
                    if is_object_pattern(selection, &name) {
                        // Select all atoms in this object
                        (0..mol.atom_count() as u32).map(AtomIndex).collect()
                    } else {
                        // For complex selections that fail, skip this molecule
                        // (it might match atoms in other molecules)
                        if matches!(e, pymol_select::SelectError::Parse(_)) {
                            return Err(CmdError::Selection(format!(
                                "Invalid selection expression: {}",
                                e
                            )));
                        }
                        continue;
                    }
                }
            }
        };

        if atoms.is_empty() {
            continue;
        }

        // Apply transformation to selected atoms
        if let Some(mol_obj) = ctx.viewer.objects_mut().get_molecule_mut(&name) {
            transform_fn(mol_obj.molecule_mut(), &atoms, state);
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

/// Check if a selection string is an object name pattern matching the given name
fn is_object_pattern(selection: &str, name: &str) -> bool {
    // Exact match
    if selection == name {
        return true;
    }
    // Simple wildcard match
    if selection.contains('*') {
        // Handle simple prefix/suffix wildcards
        if selection == "*" {
            return true;
        }
        if let Some(prefix) = selection.strip_suffix('*') {
            return name.starts_with(prefix);
        }
        if let Some(suffix) = selection.strip_prefix('*') {
            return name.ends_with(suffix);
        }
        // For more complex patterns, check if the non-wildcard parts match
        let parts: Vec<&str> = selection.split('*').collect();
        if parts.len() == 2 {
            return name.starts_with(parts[0]) && name.ends_with(parts[1]);
        }
    }
    false
}

// ============================================================================
// Helper functions - Atom-level transformations
// ============================================================================

/// Apply translation to specific atoms based on state parameter
fn apply_translation_to_atoms(
    mol: &mut pymol_mol::ObjectMolecule,
    state: i64,
    atoms: &[AtomIndex],
    delta: &Vec3,
) {
    let num_states = mol.state_count();
    if num_states == 0 {
        return;
    }

    match state {
        0 => {
            // All states
            for i in 0..num_states {
                mol.translate_atoms(i, atoms, delta.clone());
            }
        }
        -1 => {
            // Current state (use state 0 as current)
            mol.translate_atoms(0, atoms, delta.clone());
        }
        s if s > 0 => {
            // Specific state (1-indexed in PyMOL, 0-indexed internally)
            let state_idx = (s - 1) as usize;
            if state_idx < num_states {
                mol.translate_atoms(state_idx, atoms, delta.clone());
            }
        }
        _ => {}
    }
}

/// Apply TTT matrix transform to specific atoms based on state parameter
fn apply_ttt_to_atoms(
    mol: &mut pymol_mol::ObjectMolecule,
    state: i64,
    atoms: &[AtomIndex],
    ttt: &[f32; 16],
) {
    let num_states = mol.state_count();
    if num_states == 0 {
        return;
    }

    match state {
        0 => {
            // All states
            for i in 0..num_states {
                mol.transform_ttt_atoms(i, atoms, ttt);
            }
        }
        -1 => {
            // Current state
            mol.transform_ttt_atoms(0, atoms, ttt);
        }
        s if s > 0 => {
            let state_idx = (s - 1) as usize;
            if state_idx < num_states {
                mol.transform_ttt_atoms(state_idx, atoms, ttt);
            }
        }
        _ => {}
    }
}

/// Apply standard 4x4 matrix transform to specific atoms based on state parameter
fn apply_matrix_to_atoms(
    mol: &mut pymol_mol::ObjectMolecule,
    state: i64,
    atoms: &[AtomIndex],
    matrix: &lin_alg::f32::Mat4,
) {
    let num_states = mol.state_count();
    if num_states == 0 {
        return;
    }

    match state {
        0 => {
            for i in 0..num_states {
                mol.transform_atoms(i, atoms, matrix);
            }
        }
        -1 => {
            mol.transform_atoms(0, atoms, matrix);
        }
        s if s > 0 => {
            let state_idx = (s - 1) as usize;
            if state_idx < num_states {
                mol.transform_atoms(state_idx, atoms, matrix);
            }
        }
        _ => {}
    }
}
