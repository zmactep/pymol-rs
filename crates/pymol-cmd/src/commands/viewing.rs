//! Viewing commands: zoom, center, orient, reset, clip, view

use lin_alg::f32::{Mat4, Vec3};

use crate::args::ParsedCommand;
use crate::command::{Command, CommandContext, CommandRegistry};
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

    fn execute(&self, ctx: &mut CommandContext, args: &ParsedCommand) -> CmdResult {
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

    fn execute(&self, ctx: &mut CommandContext, args: &ParsedCommand) -> CmdResult {
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

    fn execute(&self, ctx: &mut CommandContext, args: &ParsedCommand) -> CmdResult {
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

    fn execute(&self, ctx: &mut CommandContext, _args: &ParsedCommand) -> CmdResult {
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

    fn execute(&self, ctx: &mut CommandContext, args: &ParsedCommand) -> CmdResult {
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

    fn execute(&self, ctx: &mut CommandContext, _args: &ParsedCommand) -> CmdResult {
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

    fn execute(&self, ctx: &mut CommandContext, args: &ParsedCommand) -> CmdResult {
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
