//! Movie commands: mplay, mstop, mpause, mset, mview, mpng, madd, frame, forward, backward, rewind, ending, rock

use crate::args::{ArgValue, ParsedCommand};
use crate::command::{Command, CommandContext, CommandRegistry, ViewerLike};
use crate::error::{CmdError, CmdResult};

/// Register movie commands
pub fn register(registry: &mut CommandRegistry) {
    registry.register(MplayCommand);
    registry.register(MstopCommand);
    registry.register(MpauseCommand);
    registry.register(MtoggleCommand);
    registry.register(MsetCommand);
    registry.register(MviewCommand);
    registry.register(MpngCommand);
    registry.register(MaddCommand);
    registry.register(FrameCommand);
    registry.register(ForwardCommand);
    registry.register(BackwardCommand);
    registry.register(RewindCommand);
    registry.register(MiddleCommand);
    registry.register(EndingCommand);
    registry.register(RockCommand);
}

// ============================================================================
// mplay command
// ============================================================================

struct MplayCommand;

impl Command for MplayCommand {
    fn name(&self) -> &str {
        "mplay"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "mplay" starts movie playback.

USAGE

    mplay

EXAMPLES

    mplay

SEE ALSO

    mstop, mpause, mset, frame
"#
    }

    fn execute<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        _args: &ParsedCommand,
    ) -> CmdResult {
        ctx.viewer.movie_play();
        if !ctx.quiet {
            ctx.print(" Movie playing.");
        }
        Ok(())
    }
}

// ============================================================================
// mstop command
// ============================================================================

struct MstopCommand;

impl Command for MstopCommand {
    fn name(&self) -> &str {
        "mstop"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "mstop" stops movie playback and resets to frame 1.

USAGE

    mstop

EXAMPLES

    mstop

SEE ALSO

    mplay, mpause, frame
"#
    }

    fn execute<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        _args: &ParsedCommand,
    ) -> CmdResult {
        ctx.viewer.movie_stop();
        if !ctx.quiet {
            ctx.print(" Movie stopped.");
        }
        Ok(())
    }
}

// ============================================================================
// mpause command
// ============================================================================

struct MpauseCommand;

impl Command for MpauseCommand {
    fn name(&self) -> &str {
        "mpause"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "mpause" pauses movie playback at the current frame.

USAGE

    mpause

EXAMPLES

    mpause

SEE ALSO

    mplay, mstop, frame
"#
    }

    fn execute<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        _args: &ParsedCommand,
    ) -> CmdResult {
        ctx.viewer.movie_pause();
        if !ctx.quiet {
            ctx.print(" Movie paused.");
        }
        Ok(())
    }
}

// ============================================================================
// mtoggle command
// ============================================================================

struct MtoggleCommand;

impl Command for MtoggleCommand {
    fn name(&self) -> &str {
        "mtoggle"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "mtoggle" toggles movie playback between play and pause.

USAGE

    mtoggle

EXAMPLES

    mtoggle

SEE ALSO

    mplay, mpause, mstop
"#
    }

    fn execute<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        _args: &ParsedCommand,
    ) -> CmdResult {
        ctx.viewer.movie_toggle();
        let state = if ctx.viewer.movie_is_playing() {
            "playing"
        } else {
            "paused"
        };
        if !ctx.quiet {
            ctx.print(&format!(" Movie {}.", state));
        }
        Ok(())
    }
}

// ============================================================================
// mset command
// ============================================================================

struct MsetCommand;

/// Parse an mset specification string into a Vec of state indices.
///
/// Syntax:
///   `1 x60`         → 60 frames, all state 1
///   `1 -30`         → 30 frames for states 1..30
///   `1 x30 1 -60`   → 30 frames of state 1, then 30 frames ramping 1..60
///   `1 x10 2 x10`   → 10 frames of state 1, then 10 of state 2
///
/// A simple integer N is shorthand for `1 xN`.
fn parse_mset_spec(args: &ParsedCommand) -> Result<Vec<usize>, CmdError> {
    // Reconstruct raw spec string from parsed args
    let mut spec = String::new();
    for (_, val) in &args.args {
        if !spec.is_empty() {
            spec.push(' ');
        }
        match val {
            ArgValue::Int(i) => spec.push_str(&i.to_string()),
            ArgValue::Float(f) => spec.push_str(&f.to_string()),
            ArgValue::String(s) => spec.push_str(s),
            _ => {}
        }
    }
    let spec = spec.trim();

    if spec.is_empty() {
        return Err(CmdError::MissingArgument("specification".to_string()));
    }

    // Tokenize: split by whitespace
    let tokens: Vec<&str> = spec.split_whitespace().collect();

    // Simple case: single integer with no operators
    if tokens.len() == 1 {
        if let Ok(n) = tokens[0].parse::<usize>() {
            // Just a number: create n frames of state 1
            return Ok(vec![1; n]);
        }
    }

    let mut states: Vec<usize> = Vec::new();
    let mut i = 0;

    while i < tokens.len() {
        let token = tokens[i];

        // Parse a state number
        let state: usize = token.parse().map_err(|_| CmdError::InvalidArgument {
            name: "specification".to_string(),
            reason: format!("expected state number, got '{}'", token),
        })?;

        i += 1;

        if i >= tokens.len() {
            // Trailing state number with no operator: single frame
            states.push(state);
            break;
        }

        let next = tokens[i];

        if let Some(count_str) = next.strip_prefix('x').or_else(|| next.strip_prefix('X')) {
            // Repeat: state xN
            let count: usize = count_str.parse().map_err(|_| CmdError::InvalidArgument {
                name: "specification".to_string(),
                reason: format!("expected count after 'x', got '{}'", next),
            })?;
            for _ in 0..count {
                states.push(state);
            }
            i += 1;
        } else if next.starts_with('-') || next.starts_with("−") {
            // Range: state_start -state_end (or negative number)
            let end_str = next.trim_start_matches('-').trim_start_matches('−');
            let end_state: usize = end_str.parse().map_err(|_| CmdError::InvalidArgument {
                name: "specification".to_string(),
                reason: format!("expected end state after '-', got '{}'", next),
            })?;
            if end_state >= state {
                for s in state..=end_state {
                    states.push(s);
                }
            } else {
                for s in (end_state..=state).rev() {
                    states.push(s);
                }
            }
            i += 1;
        } else {
            // No operator following this state: single frame of this state,
            // and the next token starts a new clause
            states.push(state);
        }
    }

    if states.is_empty() {
        return Err(CmdError::InvalidArgument {
            name: "specification".to_string(),
            reason: "empty specification".to_string(),
        });
    }

    Ok(states)
}

impl Command for MsetCommand {
    fn name(&self) -> &str {
        "mset"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "mset" sets the number of frames in the movie.

USAGE

    mset specification

ARGUMENTS

    specification = frame specification string
        "1 x30" creates 30 copies of state 1
        "1 -30" creates frames for states 1 through 30
        "1 x30 1 -60" 30 frames of state 1, then ramp 1..60
        Simple number N is shorthand for "1 xN"

EXAMPLES

    mset 100            # 100 frames of state 1
    mset 1 x60          # 60 frames of state 1
    mset 1 -60          # Frames for states 1-60
    mset 1 x30 1 -60    # 30 of state 1, then ramp states 1-60

SEE ALSO

    mplay, madd, frame, mview
"#
    }

    fn execute<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        args: &ParsedCommand,
    ) -> CmdResult {
        let states = parse_mset_spec(args)?;
        let count = states.len();
        ctx.viewer.movie_set_from_spec(states);

        if !ctx.quiet {
            ctx.print(&format!(" Movie: {} frames.", count));
        }
        Ok(())
    }
}

// ============================================================================
// frame command
// ============================================================================

struct FrameCommand;

impl Command for FrameCommand {
    fn name(&self) -> &str {
        "frame"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "frame" goes to a specific frame in the movie.

USAGE

    frame frame_number

ARGUMENTS

    frame_number = int: frame number (1-indexed)

EXAMPLES

    frame 1             # Go to first frame
    frame 50            # Go to frame 50

SEE ALSO

    mplay, mset, forward, backward
"#
    }

    fn execute<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        args: &ParsedCommand,
    ) -> CmdResult {
        let frame = args
            .get_int(0)
            .or_else(|| args.get_named_int("frame"))
            .unwrap_or(1) as usize;

        // PyMOL uses 1-based frame numbers
        ctx.viewer.movie_goto(frame.saturating_sub(1));

        if !ctx.quiet {
            ctx.print(&format!(" Frame {}.", frame));
        }
        Ok(())
    }
}

// ============================================================================
// forward command
// ============================================================================

struct ForwardCommand;

impl Command for ForwardCommand {
    fn name(&self) -> &str {
        "forward"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "forward" advances to the next frame.

USAGE

    forward

SEE ALSO

    backward, frame, mplay
"#
    }

    fn execute<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        _args: &ParsedCommand,
    ) -> CmdResult {
        ctx.viewer.movie_next();
        let frame = ctx.viewer.movie_current_frame() + 1; // 1-indexed for display
        if !ctx.quiet {
            ctx.print(&format!(" Frame {}.", frame));
        }
        Ok(())
    }
}

// ============================================================================
// backward command
// ============================================================================

struct BackwardCommand;

impl Command for BackwardCommand {
    fn name(&self) -> &str {
        "backward"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "backward" goes back to the previous frame.

USAGE

    backward

SEE ALSO

    forward, frame, mplay
"#
    }

    fn execute<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        _args: &ParsedCommand,
    ) -> CmdResult {
        ctx.viewer.movie_prev();
        let frame = ctx.viewer.movie_current_frame() + 1; // 1-indexed for display
        if !ctx.quiet {
            ctx.print(&format!(" Frame {}.", frame));
        }
        Ok(())
    }
}

// ============================================================================
// rewind command
// ============================================================================

struct RewindCommand;

impl Command for RewindCommand {
    fn name(&self) -> &str {
        "rewind"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "rewind" goes to the first frame of the movie.

USAGE

    rewind

SEE ALSO

    ending, middle, frame
"#
    }

    fn execute<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        _args: &ParsedCommand,
    ) -> CmdResult {
        ctx.viewer.movie_goto(0);
        if !ctx.quiet {
            ctx.print(" Frame 1.");
        }
        Ok(())
    }
}

// ============================================================================
// middle command
// ============================================================================

struct MiddleCommand;

impl Command for MiddleCommand {
    fn name(&self) -> &str {
        "middle"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "middle" goes to the middle frame of the movie.

USAGE

    middle

SEE ALSO

    rewind, ending, frame
"#
    }

    fn execute<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        _args: &ParsedCommand,
    ) -> CmdResult {
        let count = ctx.viewer.movie_frame_count();
        if count > 0 {
            let middle = count / 2;
            ctx.viewer.movie_goto(middle);
            if !ctx.quiet {
                ctx.print(&format!(" Frame {}.", middle + 1));
            }
        } else if !ctx.quiet {
            ctx.print(" No movie frames.");
        }
        Ok(())
    }
}

// ============================================================================
// ending command
// ============================================================================

struct EndingCommand;

impl Command for EndingCommand {
    fn name(&self) -> &str {
        "ending"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "ending" goes to the last frame of the movie.

USAGE

    ending

SEE ALSO

    rewind, middle, frame
"#
    }

    fn execute<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        _args: &ParsedCommand,
    ) -> CmdResult {
        let count = ctx.viewer.movie_frame_count();
        if count > 0 {
            ctx.viewer.movie_goto(count - 1);
            if !ctx.quiet {
                ctx.print(&format!(" Frame {}.", count));
            }
        } else if !ctx.quiet {
            ctx.print(" No movie frames.");
        }
        Ok(())
    }
}

// ============================================================================
// rock command
// ============================================================================

struct RockCommand;

impl Command for RockCommand {
    fn name(&self) -> &str {
        "rock"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "rock" toggles Y-axis rocking (oscillation) animation.

    When enabled, the view continuously rocks back and forth around
    the Y-axis. This is useful for presentations and visual inspection.

USAGE

    rock [mode]

ARGUMENTS

    mode = on/off/1/0 (optional)
        If omitted, toggles the current state.
        "on" or "1" enables rock mode.
        "off" or "0" disables rock mode.

EXAMPLES

    rock            # Toggle rock mode
    rock on         # Enable rock mode
    rock off        # Disable rock mode
    rock 1          # Enable rock mode
    rock 0          # Disable rock mode

SEE ALSO

    mplay, turn, rotate
"#
    }

    fn execute<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        args: &ParsedCommand,
    ) -> CmdResult {
        // Check for optional mode argument
        if let Some(mode_str) = args.get_str(0) {
            let enabled = match mode_str.to_lowercase().as_str() {
                "on" | "1" | "true" | "yes" => true,
                "off" | "0" | "false" | "no" => false,
                _ => {
                    return Err(CmdError::InvalidArgument {
                        name: "mode".to_string(),
                        reason: format!(
                            "Invalid mode '{}'. Use on/off, 1/0, true/false, or yes/no.",
                            mode_str
                        ),
                    });
                }
            };
            ctx.viewer.rock_set(enabled);
            if !ctx.quiet {
                let state = if enabled { "on" } else { "off" };
                ctx.print(&format!(" Rock {}.", state));
            }
        } else {
            // Toggle mode
            ctx.viewer.rock_toggle();
            if !ctx.quiet {
                let state = if ctx.viewer.rock_is_enabled() {
                    "on"
                } else {
                    "off"
                };
                ctx.print(&format!(" Rock {}.", state));
            }
        }
        Ok(())
    }
}

// ============================================================================
// mview command
// ============================================================================

struct MviewCommand;

impl Command for MviewCommand {
    fn name(&self) -> &str {
        "mview"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "mview" stores, recalls, or clears camera keyframes in the movie.

USAGE

    mview action [, frame [, scene=name [, object=name [, state=N]]]]

ARGUMENTS

    action = store, recall, clear, interpolate, reinterpolate
    frame = int: frame number (1-indexed, default: current frame)
    scene = string: scene name for scene keyframes
    object = string: object name for object keyframes
    state = int: coordinate set state

EXAMPLES

    mview store              # Store keyframe at current frame
    mview store, 30          # Store keyframe at frame 30
    mview store, 1, scene=001   # Store scene keyframe at frame 1
    mview interpolate        # Recompute interpolation
    mview clear              # Clear all keyframes

SEE ALSO

    mset, mplay, mpng, scene
"#
    }

    fn execute<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        args: &ParsedCommand,
    ) -> CmdResult {
        let action = args
            .get_str(0)
            .ok_or_else(|| CmdError::MissingArgument("action".to_string()))?;

        match action.to_lowercase().as_str() {
            "store" => {
                // Get frame (1-indexed), default to current frame
                let frame = args
                    .get_int(1)
                    .or_else(|| args.get_named_int("frame"))
                    .map(|f| (f as usize).saturating_sub(1))
                    .unwrap_or_else(|| ctx.viewer.movie_current_frame());

                let scene = args
                    .get_named("scene")
                    .and_then(|v| v.to_string_repr());
                let object = args
                    .get_named("object")
                    .and_then(|v| v.to_string_repr());
                let state = args.get_named_int("state").map(|s| s as usize);

                if let Some(ref scene_name) = scene {
                    ctx.viewer.movie_store_scene(frame, scene_name);
                    if !ctx.quiet {
                        ctx.print(&format!(
                            " Scene \"{}\" keyframe stored at frame {}.",
                            scene_name,
                            frame + 1
                        ));
                    }
                } else if let Some(ref obj_name) = object {
                    ctx.viewer.movie_store_object(frame, obj_name, state);
                    if !ctx.quiet {
                        ctx.print(&format!(
                            " Object \"{}\" keyframe stored at frame {}.",
                            obj_name,
                            frame + 1
                        ));
                    }
                } else {
                    ctx.viewer.movie_store_view(frame);
                    if !ctx.quiet {
                        ctx.print(&format!(
                            " Camera keyframe stored at frame {}.",
                            frame + 1
                        ));
                    }
                }
            }

            "recall" => {
                let frame = args
                    .get_int(1)
                    .or_else(|| args.get_named_int("frame"))
                    .map(|f| (f as usize).saturating_sub(1))
                    .unwrap_or_else(|| ctx.viewer.movie_current_frame());

                ctx.viewer.movie_goto(frame);
                if !ctx.quiet {
                    ctx.print(&format!(" Recalled keyframe at frame {}.", frame + 1));
                }
            }

            "clear" => {
                let frame = args
                    .get_int(1)
                    .or_else(|| args.get_named_int("frame"))
                    .map(|f| Some((f as usize).saturating_sub(1)))
                    .unwrap_or(None);

                ctx.viewer.movie_clear_view(frame);
                if !ctx.quiet {
                    match frame {
                        Some(f) => ctx.print(&format!(" Cleared keyframe at frame {}.", f + 1)),
                        None => ctx.print(" All keyframes cleared."),
                    }
                }
            }

            "interpolate" | "reinterpolate" => {
                ctx.viewer.movie_interpolate();
                if !ctx.quiet {
                    ctx.print(" Keyframes interpolated.");
                }
            }

            _ => {
                return Err(CmdError::InvalidArgument {
                    name: "action".to_string(),
                    reason: format!(
                        "Unknown action '{}'. Use store, recall, clear, or interpolate.",
                        action
                    ),
                });
            }
        }

        Ok(())
    }
}

// ============================================================================
// mpng command
// ============================================================================

struct MpngCommand;

impl Command for MpngCommand {
    fn name(&self) -> &str {
        "mpng"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "mpng" renders movie frames to a PNG image sequence.

USAGE

    mpng prefix [, first [, last [, preserve [, width [, height]]]]]

ARGUMENTS

    prefix = string: output path prefix (e.g., /tmp/movie/frame)
    first = int: first frame to render (1-indexed, default: 1)
    last = int: last frame to render (1-indexed, default: last frame)
    preserve = int: unused (PyMOL compat)
    width = int: image width in pixels (default: viewport width)
    height = int: image height in pixels (default: viewport height)

EXAMPLES

    mpng /tmp/movie/frame           # Render all frames
    mpng /tmp/movie/frame, 1, 30    # Render frames 1-30
    mpng /tmp/movie/frame, width=1920, height=1080

SEE ALSO

    mset, mview, mplay
"#
    }

    fn execute<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        args: &ParsedCommand,
    ) -> CmdResult {
        let prefix = args
            .get_str(0)
            .ok_or_else(|| CmdError::MissingArgument("prefix".to_string()))?;

        let total = ctx.viewer.movie_frame_count();
        if total == 0 {
            return Err(CmdError::InvalidArgument {
                name: "movie".to_string(),
                reason: "No movie frames. Use mset first.".to_string(),
            });
        }

        let first = args
            .get_int(1)
            .or_else(|| args.get_named_int("first"))
            .map(|f| (f as usize).saturating_sub(1))
            .unwrap_or(0);

        let last = args
            .get_int(2)
            .or_else(|| args.get_named_int("last"))
            .map(|f| (f as usize).saturating_sub(1))
            .unwrap_or(total - 1);

        // skip arg index 3 (preserve)
        let width = args
            .get_int(4)
            .or_else(|| args.get_named_int("width"))
            .map(|w| w as u32);

        let height = args
            .get_int(5)
            .or_else(|| args.get_named_int("height"))
            .map(|h| h as u32);

        // Create output directory if needed
        let parent = std::path::Path::new(prefix).parent();
        if let Some(dir) = parent {
            if !dir.as_os_str().is_empty() {
                std::fs::create_dir_all(dir).map_err(|e| CmdError::InvalidArgument {
                    name: "prefix".to_string(),
                    reason: format!("Cannot create directory: {}", e),
                })?;
            }
        }

        let frame_count = last - first + 1;
        if !ctx.quiet {
            ctx.print(&format!(
                " Rendering {} frames to {}NNNN.png...",
                frame_count, prefix
            ));
        }

        for frame in first..=last {
            let path = format!("{}{:04}.png", prefix, frame + 1);
            ctx.viewer
                .capture_frame_png(frame, std::path::Path::new(&path), width, height)
                .map_err(|e| CmdError::InvalidArgument {
                    name: "render".to_string(),
                    reason: format!("Frame {}: {}", frame + 1, e),
                })?;

            if !ctx.quiet {
                ctx.print(&format!(
                    " Ray: frame {} of {} ({}).",
                    frame - first + 1,
                    frame_count,
                    path
                ));
            }
        }

        if !ctx.quiet {
            ctx.print(&format!(
                " {} frames rendered to {}NNNN.png.",
                frame_count, prefix
            ));
        }
        Ok(())
    }
}

// ============================================================================
// madd command
// ============================================================================

struct MaddCommand;

impl Command for MaddCommand {
    fn name(&self) -> &str {
        "madd"
    }

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "madd" appends frames to the movie using the same specification
    syntax as mset.

USAGE

    madd specification

ARGUMENTS

    specification = frame specification (same syntax as mset)

EXAMPLES

    madd 1 x30          # Append 30 frames of state 1
    madd 1 -60          # Append frames for states 1-60

SEE ALSO

    mset, mview, mplay
"#
    }

    fn execute<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        args: &ParsedCommand,
    ) -> CmdResult {
        let states = parse_mset_spec(args)?;
        let count = states.len();
        ctx.viewer.movie_append_from_spec(states);

        if !ctx.quiet {
            ctx.print(&format!(
                " Movie: appended {} frames ({} total).",
                count,
                ctx.viewer.movie_frame_count()
            ));
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::args::ParsedCommand;
    use crate::parser::parse_command;

    #[test]
    fn test_mset_simple_count() {
        let cmd = ParsedCommand::new("mset").with_arg(ArgValue::Int(60));
        let states = parse_mset_spec(&cmd).unwrap();
        assert_eq!(states.len(), 60);
        assert!(states.iter().all(|&s| s == 1));
    }

    #[test]
    fn test_mset_1_x60() {
        // "mset 1 x60" → parser produces Int(1), String("x60")
        let cmd = ParsedCommand::new("mset")
            .with_arg(ArgValue::Int(1))
            .with_arg(ArgValue::String("x60".to_string()));
        let states = parse_mset_spec(&cmd).unwrap();
        assert_eq!(states.len(), 60);
        assert!(states.iter().all(|&s| s == 1));
    }

    #[test]
    fn test_mset_range() {
        // "mset 1 -30" → parser may produce Int(1), Int(-30) or Int(1), String("-30")
        let cmd = ParsedCommand::new("mset")
            .with_arg(ArgValue::Int(1))
            .with_arg(ArgValue::Int(-30));
        let states = parse_mset_spec(&cmd).unwrap();
        assert_eq!(states.len(), 30);
        assert_eq!(states[0], 1);
        assert_eq!(states[29], 30);
    }

    #[test]
    fn test_mset_combined() {
        // "mset 1 x10 2 x10" → 10 of state 1, then 10 of state 2
        let cmd = ParsedCommand::new("mset")
            .with_arg(ArgValue::Int(1))
            .with_arg(ArgValue::String("x10".to_string()))
            .with_arg(ArgValue::Int(2))
            .with_arg(ArgValue::String("x10".to_string()));
        let states = parse_mset_spec(&cmd).unwrap();
        assert_eq!(states.len(), 20);
        assert!(states[..10].iter().all(|&s| s == 1));
        assert!(states[10..].iter().all(|&s| s == 2));
    }

    #[test]
    fn test_mset_mixed() {
        // "mset 1 x30 1 -60" → 30 of state 1, then states 1..60
        let cmd = ParsedCommand::new("mset")
            .with_arg(ArgValue::Int(1))
            .with_arg(ArgValue::String("x30".to_string()))
            .with_arg(ArgValue::Int(1))
            .with_arg(ArgValue::Int(-60));
        let states = parse_mset_spec(&cmd).unwrap();
        assert_eq!(states.len(), 90); // 30 + 60
        assert!(states[..30].iter().all(|&s| s == 1));
        assert_eq!(states[30], 1);
        assert_eq!(states[89], 60);
    }

    #[test]
    fn test_parse_mset_1_x60_from_parser() {
        let cmd = parse_command("mset 1 x60").unwrap();
        assert_eq!(cmd.name, "mset");
        let states = parse_mset_spec(&cmd).unwrap();
        assert_eq!(states.len(), 60);
        assert!(states.iter().all(|&s| s == 1));
    }

    #[test]
    fn test_parse_mview_store() {
        let cmd = parse_command("mview store").unwrap();
        assert_eq!(cmd.name, "mview");
        assert_eq!(cmd.get_str(0), Some("store"));
    }

    #[test]
    fn test_parse_mview_store_frame() {
        let cmd = parse_command("mview store, 30").unwrap();
        assert_eq!(cmd.name, "mview");
        assert_eq!(cmd.get_str(0), Some("store"));
        assert_eq!(cmd.get_int(1), Some(30));
    }

    #[test]
    fn test_parse_mview_store_scene() {
        let cmd = parse_command("mview store, 1, scene=001").unwrap();
        assert_eq!(cmd.name, "mview");
        assert_eq!(cmd.get_str(0), Some("store"));
        assert_eq!(cmd.get_int(1), Some(1));
        // scene=001 may parse as Int(1) or String("001") depending on parser
        let scene = cmd.get_named("scene").and_then(|v| v.to_string_repr());
        assert_eq!(scene.as_deref(), Some("1"));
    }

    #[test]
    fn test_parse_mpng() {
        let cmd = parse_command("mpng /tmp/movie_test/frame").unwrap();
        assert_eq!(cmd.name, "mpng");
        assert_eq!(cmd.get_str(0), Some("/tmp/movie_test/frame"));
    }

    #[test]
    fn test_parse_mview_interpolate() {
        let cmd = parse_command("mview interpolate").unwrap();
        assert_eq!(cmd.name, "mview");
        assert_eq!(cmd.get_str(0), Some("interpolate"));
    }

    #[test]
    fn test_parse_mview_store_scene_002() {
        let cmd = parse_command("mview store, 30, scene=002").unwrap();
        assert_eq!(cmd.name, "mview");
        assert_eq!(cmd.get_str(0), Some("store"));
        assert_eq!(cmd.get_int(1), Some(30));
        let scene = cmd.get_named("scene").and_then(|v| v.to_string_repr());
        assert_eq!(scene.as_deref(), Some("2"));
    }
}
