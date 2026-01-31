//! Movie commands: mplay, mstop, mpause, mset, frame, forward, backward, rewind, ending, rock

use crate::args::ParsedCommand;
use crate::command::{Command, CommandContext, CommandRegistry, ViewerLike};
use crate::error::{CmdError, CmdResult};

/// Register movie commands
pub fn register(registry: &mut CommandRegistry) {
    registry.register(MplayCommand);
    registry.register(MstopCommand);
    registry.register(MpauseCommand);
    registry.register(MtoggleCommand);
    registry.register(MsetCommand);
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

    specification = string/int: number of frames or frame specification
        Simple number creates that many frames
        "1 x30" creates 30 copies of state 1
        "1 -30" creates frames for states 1 through 30

EXAMPLES

    mset 100            # Create 100 frames
    mset 1 x30          # 30 frames of state 1
    mset 1 -60          # Frames for states 1-60

SEE ALSO

    mplay, frame, mdo
"#
    }

    fn execute<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        args: &ParsedCommand,
    ) -> CmdResult {
        // Simple case: just a number of frames
        let count = args
            .get_int(0)
            .ok_or_else(|| CmdError::MissingArgument("specification".to_string()))?
            as usize;

        ctx.viewer.movie_set_frame_count(count);

        if !ctx.quiet {
            ctx.print(&format!(" Movie set to {} frames.", count));
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
