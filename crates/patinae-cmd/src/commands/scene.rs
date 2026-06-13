//! Scene commands: scene (store/recall/delete/list/rename/clear)

use crate::args::ParsedCommand;
use crate::command::{ArgHint, Command, CommandContext, CommandRegistry, ViewerLike};
use crate::command_help;
use crate::error::{CmdError, CmdResult};

/// Register scene commands
pub fn register(registry: &mut CommandRegistry) {
    registry.register(SceneCommand);
}

// ============================================================================
// scene command
// ============================================================================

struct SceneCommand;

impl Command for SceneCommand {
    fn name(&self) -> &str {
        "scene"
    }

    command_help! {
        CMD "scene"
        DESCRIPTION [
            "saves and recalls named scenes.",
            "",
            "A scene stores the camera view, object states, colors, representations,",
            "and current frame state.",
        ]
        USAGE [
            "scene key [, action [, message [, view [, color [, rep [, frame [, animate ]]]]]]]",
        ]
        REQUIRED [
            { "key", "string", "scene name" },
        ]
        OPTIONAL [
            { "action", "string", "store, recall, delete, clear, list, rename", "recall" },
            { "message", "string", "message to display when recalled", "\"\"" },
            { "view", "0/1", "store/recall view", "1" },
            { "color", "0/1", "store/recall colors", "1" },
            { "rep", "0/1", "store/recall representations", "1" },
            { "frame", "0/1", "store/recall frame state", "1" },
            { "animate", "float", "animation duration in seconds", "0.5" },
        ]
        NOTES("ACTIONS") [
            "store   - save current state as named scene",
            "recall  - recall a named scene (default)",
            "delete  - delete a named scene (use key=\"*\" for all)",
            "clear   - delete all scenes",
            "list    - list all stored scenes",
            "rename  - rename a scene (new_name as second positional arg)",
        ]
        EXAMPLES [
            "scene F1, store         # Store current view as \"F1\"",
            "scene F1                # Recall scene \"F1\"",
            "scene F1, recall        # Same as above",
            "scene F1, delete        # Delete scene \"F1\"",
            "scene *, delete         # Delete all scenes",
            "scene , clear           # Clear all scenes",
            "scene , list            # List all scenes",
            "scene F1, rename, F2    # Rename \"F1\" to \"F2\"",
        ]
        SEE ALSO [
            "get_view, set_view, rock, mplay",
        ]
    }

    fn arg_hints(&self) -> &[ArgHint] {
        &[
            ArgHint::None,
            ArgHint::Keywords(&["store", "recall", "delete", "clear", "list", "rename"]),
        ]
    }

    fn execute<'v, 'r>(
        &self,
        ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>,
        args: &ParsedCommand,
    ) -> CmdResult {
        let key = args.str_arg_or(0, "key", "");

        let action = args.str_arg_or(1, "action", "recall");

        // Get optional parameters
        let _view = args.int_arg_or(3, "view", 1);
        let _color = args.int_arg_or(4, "color", 1);
        let _rep = args.int_arg_or(5, "rep", 1);
        let _frame = args.int_arg_or(6, "frame", 1);
        let animate = args.float_arg_or(7, "animate", 0.5) as f32;

        // Build storemask from flags (0x3F = ALL)
        // VIEW=0x01, ACTIVE=0x02, COLOR=0x04, REP=0x08, FRAME=0x10
        let storemask = 0x3F_u32; // For now, always store everything

        match action.to_lowercase().as_str() {
            "store" => {
                if key.is_empty() {
                    return Err(CmdError::missing_argument("key".to_string()));
                }
                ctx.viewer.scene_store(key, storemask);
                if !ctx.quiet {
                    ctx.print(&format!(" Scene \"{}\" stored.", key));
                }
            }

            "recall" | "" => {
                if key.is_empty() {
                    return Err(CmdError::missing_argument("key".to_string()));
                }
                ctx.viewer
                    .scene_recall(key, true, animate)
                    .map_err(|e| CmdError::invalid_arg("key", e))?;
                if !ctx.quiet {
                    ctx.print(&format!(" Scene \"{}\" recalled.", key));
                }
            }

            "delete" => {
                if key == "*" || key.is_empty() {
                    ctx.viewer.scene_clear();
                    if !ctx.quiet {
                        ctx.print(" All scenes deleted.");
                    }
                } else if ctx.viewer.scene_delete(key) {
                    if !ctx.quiet {
                        ctx.print(&format!(" Scene \"{}\" deleted.", key));
                    }
                } else {
                    return Err(CmdError::invalid_arg(
                        "key",
                        format!("Scene '{}' not found.", key),
                    ));
                }
            }

            "clear" => {
                ctx.viewer.scene_clear();
                if !ctx.quiet {
                    ctx.print(" All scenes cleared.");
                }
            }

            "list" => {
                let scenes = ctx.viewer.scene_list();
                if scenes.is_empty() {
                    ctx.print(" No scenes stored.");
                } else {
                    ctx.print(&format!(" Scenes: {}", scenes.join(", ")));
                }
            }

            "rename" => {
                if key.is_empty() {
                    return Err(CmdError::missing_argument("key".to_string()));
                }
                let new_key = args
                    .str_arg(2, "new_key")
                    .ok_or_else(|| CmdError::missing_argument("new_key".to_string()))?;

                ctx.viewer
                    .scene_rename(key, new_key)
                    .map_err(|e| CmdError::invalid_arg("key", e))?;

                if !ctx.quiet {
                    ctx.print(&format!(" Scene \"{}\" renamed to \"{}\".", key, new_key));
                }
            }

            _ => {
                return Err(CmdError::invalid_arg(
                    "action",
                    format!("Unknown action: '{}'", action),
                ));
            }
        }

        Ok(())
    }
}
