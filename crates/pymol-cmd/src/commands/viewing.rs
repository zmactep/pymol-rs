//! Viewing commands: zoom, center, orient, reset, clip, move, turn, origin, view

use std::f32::consts::PI;

use lin_alg::f32::{Mat4, Vec3};
use pymol_scene::Object; // For extent() method on MoleculeObject

use crate::args::ParsedCommand;
use crate::command::{Command, CommandContext, CommandRegistry, ViewerLike};
use crate::error::{CmdError, CmdResult};

/// Register viewing commands
pub fn register(registry: &mut CommandRegistry) {
    registry.register(ZoomCommand);
    registry.register(CenterCommand);
    registry.register(OrientCommand);
    registry.register(ResetCommand);
    registry.register(ClipCommand);
    registry.register(GetViewCommand);
    registry.register(SetViewCommand);
    registry.register(MoveCommand);
    registry.register(TurnCommand);
    registry.register(OriginCommand);
    registry.register(ViewCommand);
    registry.register(ViewportCommand);
    registry.register(FullScreenCommand);
}

// ============================================================================
// zoom command
// ============================================================================

struct ZoomCommand;

impl Command for ZoomCommand {
    fn name(&self) -> &str {
        "zoom"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "zoom" scales and translates the camera to cover the indicated selection.

USAGE

    zoom [ selection [, buffer [, state [, complete [, animate ]]]]]

ARGUMENTS

    selection = string: selection to zoom to (default: all)
    buffer = float: extra space around selection in Angstroms (default: 0)
    state = integer: state to zoom to (default: 0 = current)
    complete = 0/1: if 1, always do a complete zoom (default: 0)
    animate = float: animation time in seconds (default: 0)

EXAMPLES

    zoom
    zoom chain A
    zoom organic, buffer=5
"#
    }

    fn execute<'v, 'r>(&self, ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>, args: &ParsedCommand) -> CmdResult {
        let selection = args
            .get_str(0)
            .or_else(|| args.get_named_str("selection"))
            .unwrap_or("all");

        let _buffer = args
            .get_float(1)
            .or_else(|| args.get_named_float("buffer"))
            .unwrap_or(0.0);

        let _state = args
            .get_int(2)
            .or_else(|| args.get_named_int("state"))
            .unwrap_or(0);

        let _animate = args
            .get_float(4)
            .or_else(|| args.get_named_float("animate"))
            .unwrap_or(0.0);

        // For "all", zoom to all objects (preserves rotation)
        if selection == "all" || selection == "*" {
            ctx.viewer.zoom_all();
            if !ctx.quiet {
                ctx.print(" Zoomed to all objects");
            }
            return Ok(());
        }

        // Try to zoom to a specific object by name
        if ctx.viewer.objects().contains(selection) {
            ctx.viewer.zoom_on(selection);
            if !ctx.quiet {
                ctx.print(&format!(" Zoomed to \"{}\"", selection));
            }
            return Ok(());
        }

        // Try pattern matching
        let matches: Vec<String> = ctx.viewer.objects().matching(selection)
            .iter().map(|s| s.to_string()).collect();
        if !matches.is_empty() {
            // For now, zoom on first match
            if let Some(name) = matches.first() {
                ctx.viewer.zoom_on(name);
                if !ctx.quiet {
                    ctx.print(&format!(" Zoomed to \"{}\"", name));
                }
            }
            return Ok(());
        }

        // TODO: Use pymol-select for atom selections

        Err(CmdError::Selection(format!(
            "No objects matching '{}'",
            selection
        )))
    }
}

// ============================================================================
// center command
// ============================================================================

struct CenterCommand;

impl Command for CenterCommand {
    fn name(&self) -> &str {
        "center"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "center" translates the camera origin to the center of the selection.

USAGE

    center [ selection [, state [, origin [, animate ]]]]

ARGUMENTS

    selection = string: selection to center on (default: all)
    state = integer: state to center on (default: 0)
    origin = 0/1: also set origin (default: 1)
    animate = float: animation time in seconds (default: 0)

EXAMPLES

    center
    center chain A
"#
    }

    fn execute<'v, 'r>(&self, ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>, args: &ParsedCommand) -> CmdResult {
        let selection = args
            .get_str(0)
            .or_else(|| args.get_named_str("selection"))
            .unwrap_or("all");

        // Delegate to zoom for now
        if selection == "all" || selection == "*" {
            ctx.viewer.center_all();
        } else if ctx.viewer.objects().contains(selection) {
            ctx.viewer.center_on(selection);
        } else {
            // Try pattern matching
            let matches: Vec<String> = ctx.viewer.objects().matching(selection)
                .iter().map(|s| s.to_string()).collect();
            if let Some(name) = matches.first() {
                ctx.viewer.center_on(name);
            }
        }

        if !ctx.quiet {
            ctx.print(&format!(" Centered on \"{}\"", selection));
        }

        Ok(())
    }
}

// ============================================================================
// orient command
// ============================================================================

struct OrientCommand;

impl Command for OrientCommand {
    fn name(&self) -> &str {
        "orient"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "orient" aligns the principal axes of the selection with the camera.

USAGE

    orient [ selection [, state [, animate ]]]

ARGUMENTS

    selection = string: selection to orient (default: all)
    state = integer: state to orient (default: 0)
    animate = float: animation time in seconds (default: 0)

EXAMPLES

    orient
    orient polymer
"#
    }

    fn execute<'v, 'r>(&self, ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>, args: &ParsedCommand) -> CmdResult {
        let selection = args
            .get_str(0)
            .or_else(|| args.get_named_str("selection"))
            .unwrap_or("all");

        // Reset view which includes resetting rotation
        // TODO: implement proper principal axis orientation
        ctx.viewer.reset_view();

        if !ctx.quiet {
            ctx.print(&format!(" Oriented to \"{}\"", selection));
        }

        Ok(())
    }
}

// ============================================================================
// reset command
// ============================================================================

struct ResetCommand;

impl Command for ResetCommand {
    fn name(&self) -> &str {
        "reset"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "reset" resets the camera to the default view.

USAGE

    reset [ object ]

ARGUMENTS

    object = string: object to reset (default: all)

EXAMPLES

    reset
"#
    }

    fn execute<'v, 'r>(&self, ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>, _args: &ParsedCommand) -> CmdResult {
        ctx.viewer.reset_view();

        if !ctx.quiet {
            ctx.print(" Reset view");
        }

        Ok(())
    }
}

// ============================================================================
// clip command
// ============================================================================

struct ClipCommand;

impl Command for ClipCommand {
    fn name(&self) -> &str {
        "clip"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "clip" alters the near and far clipping planes.

USAGE

    clip mode, distance [, selection [, state ]]

ARGUMENTS

    mode = near, far, move, slab, or atoms
    distance = float: clipping distance adjustment
    selection = string: selection for reference (default: all)

MODES

    near: adjust near clipping plane
    far: adjust far clipping plane
    move: move both planes
    slab: set slab thickness
    atoms: clip to atoms

EXAMPLES

    clip near, -5
    clip far, 10
    clip slab, 20
"#
    }

    fn execute<'v, 'r>(&self, ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>, args: &ParsedCommand) -> CmdResult {
        let mode = args
            .get_str(0)
            .or_else(|| args.get_named_str("mode"))
            .ok_or_else(|| CmdError::MissingArgument("mode".to_string()))?;

        let distance = args
            .get_float(1)
            .or_else(|| args.get_named_float("distance"))
            .ok_or_else(|| CmdError::MissingArgument("distance".to_string()))?
            as f32;

        let (clip_front, clip_back) = {
            let view = ctx.viewer.camera_mut().view_mut();

            match mode.to_lowercase().as_str() {
                "near" => {
                    view.clip_front = (view.clip_front + distance).max(0.01);
                }
                "far" => {
                    view.clip_back = (view.clip_back + distance).max(view.clip_front + 0.01);
                }
                "move" => {
                    view.clip_front = (view.clip_front + distance).max(0.01);
                    view.clip_back = (view.clip_back + distance).max(view.clip_front + 0.01);
                }
                "slab" => {
                    // Set slab thickness centered on current position
                    let center = (view.clip_front + view.clip_back) / 2.0;
                    view.clip_front = (center - distance / 2.0).max(0.01);
                    view.clip_back = center + distance / 2.0;
                }
                _ => {
                    return Err(CmdError::invalid_arg(
                        "mode",
                        format!("unknown clip mode: {}", mode),
                    ));
                }
            }
            (view.clip_front, view.clip_back)
        };

        ctx.viewer.request_redraw();

        if !ctx.quiet {
            ctx.print(&format!(
                " Clip: near={:.1}, far={:.1}",
                clip_front, clip_back
            ));
        }

        Ok(())
    }
}

// ============================================================================
// get_view command
// ============================================================================

struct GetViewCommand;

impl Command for GetViewCommand {
    fn name(&self) -> &str {
        "get_view"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "get_view" returns the current view matrix as a tuple of 18 values.

USAGE

    get_view [ output ]

ARGUMENTS

    output = 0/1/2: output mode (default: 0)

EXAMPLES

    get_view
"#
    }

    fn execute<'v, 'r>(&self, ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>, _args: &ParsedCommand) -> CmdResult {
        let view = ctx.viewer.camera().current_view();

        // Format view as PyMOL-style output
        let mut output = String::from("### cut below here and paste into PyMOL ###\n");
        output.push_str("set_view (\\\n");

        // Rotation matrix (4x4, but we output 3x3 part for PyMOL compatibility)
        let r = &view.rotation;
        output.push_str(&format!(
            "  {:12.6}, {:12.6}, {:12.6},\\\n",
            r.data[0], r.data[1], r.data[2]
        ));
        output.push_str(&format!(
            "  {:12.6}, {:12.6}, {:12.6},\\\n",
            r.data[4], r.data[5], r.data[6]
        ));
        output.push_str(&format!(
            "  {:12.6}, {:12.6}, {:12.6},\\\n",
            r.data[8], r.data[9], r.data[10]
        ));

        // Camera position
        output.push_str(&format!(
            "  {:12.6}, {:12.6}, {:12.6},\\\n",
            view.position.x, view.position.y, view.position.z
        ));

        // Origin
        output.push_str(&format!(
            "  {:12.6}, {:12.6}, {:12.6},\\\n",
            view.origin.x, view.origin.y, view.origin.z
        ));

        // Clip planes and FOV
        output.push_str(&format!(
            "  {:12.6}, {:12.6}, {:12.6} )\n",
            view.clip_front, view.clip_back, view.fov
        ));

        output.push_str("### cut above here and paste into PyMOL ###");

        ctx.print(&output);

        Ok(())
    }
}

// ============================================================================
// set_view command
// ============================================================================

struct SetViewCommand;

impl Command for SetViewCommand {
    fn name(&self) -> &str {
        "set_view"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "set_view" sets the camera view matrix from a tuple of 18 values.

USAGE

    set_view (v1, v2, v3, ...)

ARGUMENTS

    view = tuple of 18 floats defining the view matrix

EXAMPLES

    set_view (\
      1.0, 0.0, 0.0,\
      0.0, 1.0, 0.0,\
      0.0, 0.0, 1.0,\
      0.0, 0.0, -50.0,\
      0.0, 0.0, 0.0,\
      10.0, 100.0, 0.0)
"#
    }

    fn execute<'v, 'r>(&self, ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>, args: &ParsedCommand) -> CmdResult {
        // Parse the view values
        // This is complex because PyMOL accepts a tuple as the argument
        let arg0 = args.get_arg(0).ok_or_else(|| CmdError::MissingArgument("view".to_string()))?;

        let values: Vec<f32> = match arg0 {
            crate::args::ArgValue::List(items) => {
                items.iter().filter_map(|v| v.as_float().map(|f| f as f32)).collect()
            }
            crate::args::ArgValue::String(s) => {
                // Try to parse as comma-separated values
                let s = s.trim_matches(|c| c == '(' || c == ')');
                s.split(',')
                    .filter_map(|v| v.trim().parse::<f32>().ok())
                    .collect()
            }
            _ => return Err(CmdError::invalid_arg("view", "expected tuple of 18 values")),
        };

        if values.len() < 18 {
            return Err(CmdError::invalid_arg(
                "view",
                format!("expected 18 values, got {}", values.len()),
            ));
        }

        let view = ctx.viewer.camera_mut().view_mut();

        // Set rotation matrix (3x3 -> 4x4)
        // PyMOL uses row-major 3x3, we use column-major 4x4
        let mut rotation_data = [[0.0f32; 4]; 4];
        rotation_data[0] = [values[0], values[1], values[2], 0.0];
        rotation_data[1] = [values[3], values[4], values[5], 0.0];
        rotation_data[2] = [values[6], values[7], values[8], 0.0];
        rotation_data[3] = [0.0, 0.0, 0.0, 1.0];
        view.rotation = Mat4::from(rotation_data);

        // Set camera position
        view.position = Vec3::new(values[9], values[10], values[11]);

        // Set origin
        view.origin = Vec3::new(values[12], values[13], values[14]);

        // Set clip planes and FOV
        view.clip_front = values[15];
        view.clip_back = values[16];
        view.fov = values[17];

        ctx.viewer.request_redraw();

        if !ctx.quiet {
            ctx.print(" View set");
        }

        Ok(())
    }
}

// ============================================================================
// move command
// ============================================================================

struct MoveCommand;

impl Command for MoveCommand {
    fn name(&self) -> &str {
        "move"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "move" translates the camera about one of the three primary axes.

USAGE

    move axis, distance

ARGUMENTS

    axis = x, y, or z: axis along which to translate

    distance = float: distance to move in Angstroms

EXAMPLES

    move x, 3
    move y, -1
    move z, 10

NOTES

    Positive x moves right, positive y moves up, positive z moves
    toward the viewer.

SEE ALSO

    turn, rotate, translate, zoom, center, clip
"#
    }

    fn execute<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        args: &ParsedCommand,
    ) -> CmdResult {
        // Parse axis argument (required)
        let axis = args
            .get_str(0)
            .or_else(|| args.get_named_str("axis"))
            .ok_or_else(|| CmdError::MissingArgument("axis".to_string()))?;

        // Parse distance argument (required)
        let distance = args
            .get_float(1)
            .or_else(|| args.get_named_float("distance"))
            .ok_or_else(|| CmdError::MissingArgument("distance".to_string()))?
            as f32;

        // Create translation vector based on axis
        let delta = match axis.to_lowercase().as_str() {
            "x" => Vec3::new(distance, 0.0, 0.0),
            "y" => Vec3::new(0.0, distance, 0.0),
            "z" => Vec3::new(0.0, 0.0, distance),
            _ => {
                return Err(CmdError::invalid_arg(
                    "axis",
                    format!("unknown axis '{}'. Use x, y, or z.", axis),
                ));
            }
        };

        // Apply translation to camera
        ctx.viewer.camera_mut().translate(delta);
        ctx.viewer.request_redraw();

        if !ctx.quiet {
            ctx.print(&format!(" Moved {} by {:.3}", axis, distance));
        }

        Ok(())
    }
}

// ============================================================================
// turn command
// ============================================================================

struct TurnCommand;

impl Command for TurnCommand {
    fn name(&self) -> &str {
        "turn"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "turn" rotates the camera about one of the three primary axes,
    centered at the origin.

USAGE

    turn axis, angle

ARGUMENTS

    axis = x, y, or z: axis about which to rotate

    angle = float: degrees of rotation

EXAMPLES

    turn x, 90
    turn y, 45
    turn z, -30

NOTES

    Rotations follow the right-hand rule. For example, a positive
    rotation about the y-axis will rotate the view to the right.

SEE ALSO

    move, rotate, translate, zoom, center, clip
"#
    }

    fn execute<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        args: &ParsedCommand,
    ) -> CmdResult {
        // Parse axis argument (required)
        let axis = args
            .get_str(0)
            .or_else(|| args.get_named_str("axis"))
            .ok_or_else(|| CmdError::MissingArgument("axis".to_string()))?;

        // Parse angle argument (required)
        let angle_deg = args
            .get_float(1)
            .or_else(|| args.get_named_float("angle"))
            .ok_or_else(|| CmdError::MissingArgument("angle".to_string()))?
            as f32;

        // Convert to radians
        let angle_rad = angle_deg * PI / 180.0;

        // Apply rotation based on axis
        match axis.to_lowercase().as_str() {
            "x" => ctx.viewer.camera_mut().rotate_x(angle_rad),
            "y" => ctx.viewer.camera_mut().rotate_y(angle_rad),
            "z" => ctx.viewer.camera_mut().rotate_z(angle_rad),
            _ => {
                return Err(CmdError::invalid_arg(
                    "axis",
                    format!("unknown axis '{}'. Use x, y, or z.", axis),
                ));
            }
        };

        ctx.viewer.request_redraw();

        if !ctx.quiet {
            ctx.print(&format!(" Turned {} by {:.1} degrees", axis, angle_deg));
        }

        Ok(())
    }
}

// ============================================================================
// origin command
// ============================================================================

struct OriginCommand;

impl Command for OriginCommand {
    fn name(&self) -> &str {
        "origin"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "origin" sets the center of rotation about a selection or position.
    If an object name is specified, it can be used to set the center of
    rotation for the object (for use in animation and editing).

USAGE

    origin [ selection [, object [, position [, state ]]]]

ARGUMENTS

    selection = string: selection-expression or name-list {default: (all)}

    object = string: object name (optional)

    position = [x, y, z]: explicit position coordinates {default: None}

    state = integer: state to use for coordinates {default: 0}

EXAMPLES

    origin chain A
    origin position=[1.0, 2.0, 3.0]
    origin resi 100

SEE ALSO

    zoom, orient, center, reset
"#
    }

    fn execute<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        args: &ParsedCommand,
    ) -> CmdResult {
        // Check for explicit position first
        if let Some(position_arg) = args.get_named("position") {
            let position = parse_position_from_arg(position_arg)?;
            ctx.viewer.camera_mut().view_mut().origin = position.clone();
            ctx.viewer.request_redraw();

            if !ctx.quiet {
                ctx.print(&format!(
                    " Origin set to [{:.3}, {:.3}, {:.3}]",
                    position.x, position.y, position.z
                ));
            }
            return Ok(());
        }

        // Get selection (default: "all")
        let selection = args
            .get_str(0)
            .or_else(|| args.get_named_str("selection"))
            .unwrap_or("all");

        // For "all", set origin to center of all objects
        if selection == "all" || selection == "*" {
            // Get bounding box of all objects
            if let Some((min, max)) = ctx.viewer.objects().extent() {
                let center = Vec3::new(
                    (min.x + max.x) * 0.5,
                    (min.y + max.y) * 0.5,
                    (min.z + max.z) * 0.5,
                );
                ctx.viewer.camera_mut().view_mut().origin = center.clone();
                ctx.viewer.request_redraw();

                if !ctx.quiet {
                    ctx.print(&format!(
                        " Origin set to [{:.3}, {:.3}, {:.3}]",
                        center.x, center.y, center.z
                    ));
                }
            } else if !ctx.quiet {
                ctx.print(" No objects to set origin from");
            }
            return Ok(());
        }

        // Try to set origin to a specific object's center
        if ctx.viewer.objects().contains(selection) {
            if let Some(mol) = ctx.viewer.objects().get_molecule(selection) {
                if let Some((min, max)) = mol.extent() {
                    let center = Vec3::new(
                        (min.x + max.x) * 0.5,
                        (min.y + max.y) * 0.5,
                        (min.z + max.z) * 0.5,
                    );
                    ctx.viewer.camera_mut().view_mut().origin = center.clone();
                    ctx.viewer.request_redraw();

                    if !ctx.quiet {
                        ctx.print(&format!(
                            " Origin set to \"{}\" at [{:.3}, {:.3}, {:.3}]",
                            selection, center.x, center.y, center.z
                        ));
                    }
                    return Ok(());
                }
            }
        }

        // Try pattern matching
        let matches: Vec<String> = ctx
            .viewer
            .objects()
            .matching(selection)
            .iter()
            .map(|s| s.to_string())
            .collect();
        if !matches.is_empty() {
            // Use first match
            if let Some(name) = matches.first() {
                if let Some(mol) = ctx.viewer.objects().get_molecule(name) {
                    if let Some((min, max)) = mol.extent() {
                        let center = Vec3::new(
                            (min.x + max.x) * 0.5,
                            (min.y + max.y) * 0.5,
                            (min.z + max.z) * 0.5,
                        );
                        ctx.viewer.camera_mut().view_mut().origin = center.clone();
                        ctx.viewer.request_redraw();

                        if !ctx.quiet {
                            ctx.print(&format!(
                                " Origin set to \"{}\" at [{:.3}, {:.3}, {:.3}]",
                                name, center.x, center.y, center.z
                            ));
                        }
                        return Ok(());
                    }
                }
            }
        }

        // TODO: Use pymol-select for atom selections

        Err(CmdError::Selection(format!(
            "No objects matching '{}'",
            selection
        )))
    }
}

/// Parse a position vector from an argument value
fn parse_position_from_arg(arg: &crate::args::ArgValue) -> Result<Vec3, CmdError> {
    use crate::args::ArgValue;

    match arg {
        ArgValue::List(items) if items.len() >= 3 => {
            let x = items
                .get(0)
                .and_then(|v| v.as_float())
                .ok_or_else(|| CmdError::invalid_arg("position", "invalid x coordinate"))?
                as f32;
            let y = items
                .get(1)
                .and_then(|v| v.as_float())
                .ok_or_else(|| CmdError::invalid_arg("position", "invalid y coordinate"))?
                as f32;
            let z = items
                .get(2)
                .and_then(|v| v.as_float())
                .ok_or_else(|| CmdError::invalid_arg("position", "invalid z coordinate"))?
                as f32;
            Ok(Vec3::new(x, y, z))
        }
        ArgValue::String(s) => {
            // Try to parse "[x, y, z]" format
            let s = s.trim().trim_matches(|c| c == '[' || c == ']');
            let parts: Vec<&str> = s.split(',').collect();
            if parts.len() >= 3 {
                let x = parts[0]
                    .trim()
                    .parse::<f32>()
                    .map_err(|_| CmdError::invalid_arg("position", "invalid x coordinate"))?;
                let y = parts[1]
                    .trim()
                    .parse::<f32>()
                    .map_err(|_| CmdError::invalid_arg("position", "invalid y coordinate"))?;
                let z = parts[2]
                    .trim()
                    .parse::<f32>()
                    .map_err(|_| CmdError::invalid_arg("position", "invalid z coordinate"))?;
                Ok(Vec3::new(x, y, z))
            } else {
                Err(CmdError::invalid_arg(
                    "position",
                    "expected [x, y, z] format",
                ))
            }
        }
        _ => Err(CmdError::invalid_arg(
            "position",
            "expected [x, y, z] format",
        )),
    }
}

// ============================================================================
// view command
// ============================================================================

struct ViewCommand;

impl Command for ViewCommand {
    fn name(&self) -> &str {
        "view"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "view" saves and restores named camera views.

    Unlike "scene", "view" only stores camera state (rotation, position,
    origin, clipping planes, FOV) without colors, representations, or
    frame state.

USAGE

    view key [, action [, animate ]]

ARGUMENTS

    key = string or *: view name, or * for all views

    action = store, recall, clear: {default: recall}

    animate = float: animation duration in seconds {default: 0}

NOTES

    Views F1 through F12 are automatically bound to function keys
    provided that "set_key" has not been used to redefine the
    behavior of the respective key, and that a "scene" has not been
    defined for that key.

EXAMPLES

    view F1, store          # Store current view as "F1"
    view F1                 # Recall view "F1"
    view F1, recall, 1.0    # Recall with 1 second animation
    view F1, clear          # Delete view "F1"
    view *, clear           # Delete all views

SEE ALSO

    scene, get_view, set_view
"#
    }

    fn execute<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        args: &ParsedCommand,
    ) -> CmdResult {
        let key = args
            .get_str(0)
            .or_else(|| args.get_named_str("key"))
            .unwrap_or("");

        let action = args
            .get_str(1)
            .or_else(|| args.get_named_str("action"))
            .unwrap_or("recall");

        let animate = args
            .get_float(2)
            .or_else(|| args.get_named_float("animate"))
            .unwrap_or(0.0) as f32;

        // Handle wildcard key for listing/clearing
        if key == "*" {
            match action.to_lowercase().as_str() {
                "clear" => {
                    ctx.viewer.view_clear();
                    if !ctx.quiet {
                        ctx.print(" All views cleared.");
                    }
                }
                _ => {
                    // List all views
                    let views = ctx.viewer.view_list();
                    if views.is_empty() {
                        ctx.print(" No views stored.");
                    } else {
                        ctx.print(&format!(" Stored views: {}", views.join(", ")));
                    }
                }
            }
            return Ok(());
        }

        match action.to_lowercase().as_str() {
            "store" => {
                if key.is_empty() {
                    return Err(CmdError::MissingArgument("key".to_string()));
                }
                ctx.viewer.view_store(key);
                if !ctx.quiet {
                    ctx.print(&format!(" View \"{}\" stored.", key));
                }
            }

            "recall" | "" => {
                if key.is_empty() {
                    return Err(CmdError::MissingArgument("key".to_string()));
                }
                ctx.viewer
                    .view_recall(key, animate)
                    .map_err(|e| CmdError::InvalidArgument {
                        name: "key".to_string(),
                        reason: e,
                    })?;
                if !ctx.quiet {
                    ctx.print(&format!(" View \"{}\" recalled.", key));
                }
            }

            "clear" | "delete" => {
                if key.is_empty() {
                    ctx.viewer.view_clear();
                    if !ctx.quiet {
                        ctx.print(" All views cleared.");
                    }
                } else if ctx.viewer.view_delete(key) {
                    if !ctx.quiet {
                        ctx.print(&format!(" View \"{}\" deleted.", key));
                    }
                } else {
                    return Err(CmdError::InvalidArgument {
                        name: "key".to_string(),
                        reason: format!("View '{}' not found.", key),
                    });
                }
            }

            _ => {
                return Err(CmdError::invalid_arg(
                    "action",
                    format!("unknown action '{}'. Use store, recall, or clear.", action),
                ));
            }
        }

        Ok(())
    }
}

// ============================================================================
// viewport command
// ============================================================================

struct ViewportCommand;

impl Command for ViewportCommand {
    fn name(&self) -> &str {
        "viewport"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "viewport" changes the size of the graphics display area.

USAGE

    viewport [ width [, height ]]

ARGUMENTS

    width = integer: width in pixels {default: current}

    height = integer: height in pixels {default: current}

NOTES

    If only width is specified, height will be scaled to maintain
    the current aspect ratio.

    If no arguments are given, the current viewport size is displayed.

EXAMPLES

    viewport                # Show current size
    viewport 800, 600       # Set to 800x600
    viewport 1920, 1080     # Set to 1920x1080

SEE ALSO

    full_screen, get_view, set_view
"#
    }

    fn execute<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        args: &ParsedCommand,
    ) -> CmdResult {
        let width = args
            .get_int(0)
            .or_else(|| args.get_named_int("width"));

        let height = args
            .get_int(1)
            .or_else(|| args.get_named_int("height"));

        match (width, height) {
            (None, None) => {
                // Display current viewport size
                let (w, h) = ctx.viewer.viewport_size();
                if w > 0 && h > 0 {
                    ctx.print(&format!(" Viewport: {} x {}", w, h));
                } else {
                    ctx.print(" Viewport size not available.");
                }
            }
            (Some(w), None) => {
                // Width only - maintain aspect ratio
                let (current_w, current_h) = ctx.viewer.viewport_size();
                let aspect = if current_w > 0 {
                    current_h as f64 / current_w as f64
                } else {
                    0.75 // Default 4:3 aspect
                };
                let h = (w as f64 * aspect).round() as u32;
                ctx.viewer.viewport_set_size(w as u32, h);
                if !ctx.quiet {
                    ctx.print(&format!(" Viewport set to {} x {}", w, h));
                }
            }
            (Some(w), Some(h)) => {
                // Both specified
                ctx.viewer.viewport_set_size(w as u32, h as u32);
                if !ctx.quiet {
                    ctx.print(&format!(" Viewport set to {} x {}", w, h));
                }
            }
            (None, Some(h)) => {
                // Height only - maintain aspect ratio
                let (current_w, current_h) = ctx.viewer.viewport_size();
                let aspect = if current_h > 0 {
                    current_w as f64 / current_h as f64
                } else {
                    4.0 / 3.0 // Default 4:3 aspect
                };
                let w = (h as f64 * aspect).round() as u32;
                ctx.viewer.viewport_set_size(w, h as u32);
                if !ctx.quiet {
                    ctx.print(&format!(" Viewport set to {} x {}", w, h));
                }
            }
        }

        Ok(())
    }
}

// ============================================================================
// full_screen command
// ============================================================================

struct FullScreenCommand;

impl Command for FullScreenCommand {
    fn name(&self) -> &str {
        "full_screen"
    }

    fn aliases(&self) -> &[&str] {
        &["fullscreen"]
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "full_screen" enables or disables full screen mode.

USAGE

    full_screen [ toggle ]

ARGUMENTS

    toggle = on, off, or toggle {default: toggle}
        "on" or "1" enables fullscreen
        "off" or "0" disables fullscreen
        "toggle" or "-1" toggles the current state (default)

EXAMPLES

    full_screen             # Toggle fullscreen
    full_screen on          # Enable fullscreen
    full_screen off         # Disable fullscreen

NOTES

    This does not work correctly on all platforms. If you encounter
    trouble, try using the maximize button on the viewer window
    instead.

SEE ALSO

    viewport
"#
    }

    fn execute<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        args: &ParsedCommand,
    ) -> CmdResult {
        let toggle = args
            .get_str(0)
            .or_else(|| args.get_named_str("toggle"));

        match toggle {
            None | Some("toggle") | Some("-1") => {
                // Toggle
                ctx.viewer.toggle_fullscreen();
                let state = if ctx.viewer.is_fullscreen() { "on" } else { "off" };
                if !ctx.quiet {
                    ctx.print(&format!(" Fullscreen {}.", state));
                }
            }
            Some("on") | Some("1") | Some("true") | Some("yes") => {
                ctx.viewer.set_fullscreen(true);
                if !ctx.quiet {
                    ctx.print(" Fullscreen on.");
                }
            }
            Some("off") | Some("0") | Some("false") | Some("no") => {
                ctx.viewer.set_fullscreen(false);
                if !ctx.quiet {
                    ctx.print(" Fullscreen off.");
                }
            }
            Some(other) => {
                return Err(CmdError::invalid_arg(
                    "toggle",
                    format!("unknown value '{}'. Use on, off, or toggle.", other),
                ));
            }
        }

        Ok(())
    }
}
