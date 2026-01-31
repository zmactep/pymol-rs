//! Scene commands: scene (store/recall/delete/list/rename/clear)

use crate::args::ParsedCommand;
use crate::command::{Command, CommandContext, CommandRegistry, ViewerLike};
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

    fn help(&self) -> &str {
        r#"
DESCRIPTION

    "scene" saves and recalls named scenes.

    A scene stores the camera view, object states, colors, representations,
    and current frame state.

USAGE

    scene key [, action [, message [, view [, color [, rep [, frame [, animate ]]]]]]]

ARGUMENTS

    key = string: scene name (required)
    action = store, recall, delete, clear, list, rename (default: recall)
    message = string: message to display when recalled (default: "")
    view = 1/0: store/recall view (default: 1)
    color = 1/0: store/recall colors (default: 1)
    rep = 1/0: store/recall representations (default: 1)
    frame = 1/0: store/recall frame state (default: 1)
    animate = float: animation duration in seconds (default: 0.5)

ACTIONS

    store   - save current state as named scene
    recall  - recall a named scene (default)
    delete  - delete a named scene (use key="*" for all)
    clear   - delete all scenes
    list    - list all stored scenes
    rename  - rename a scene (new_name as second positional arg)

EXAMPLES

    scene F1, store         # Store current view as "F1"
    scene F1                # Recall scene "F1"
    scene F1, recall        # Same as above
    scene F1, delete        # Delete scene "F1"
    scene *, delete         # Delete all scenes
    scene , clear           # Clear all scenes
    scene , list            # List all scenes
    scene F1, rename, F2    # Rename "F1" to "F2"

SEE ALSO

    get_view, set_view, rock, mplay
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

        // Get optional parameters
        let _view = args
            .get_int(3)
            .or_else(|| args.get_named_int("view"))
            .unwrap_or(1);
        let _color = args
            .get_int(4)
            .or_else(|| args.get_named_int("color"))
            .unwrap_or(1);
        let _rep = args
            .get_int(5)
            .or_else(|| args.get_named_int("rep"))
            .unwrap_or(1);
        let _frame = args
            .get_int(6)
            .or_else(|| args.get_named_int("frame"))
            .unwrap_or(1);
        let animate = args
            .get_float(7)
            .or_else(|| args.get_named_float("animate"))
            .unwrap_or(0.5) as f32;

        // Build storemask from flags (0x3F = ALL)
        // VIEW=0x01, ACTIVE=0x02, COLOR=0x04, REP=0x08, FRAME=0x10
        let storemask = 0x3F_u32; // For now, always store everything

        match action.to_lowercase().as_str() {
            "store" => {
                if key.is_empty() {
                    return Err(CmdError::MissingArgument("key".to_string()));
                }
                ctx.viewer.scene_store(key, storemask);
                if !ctx.quiet {
                    ctx.print(&format!(" Scene \"{}\" stored.", key));
                }
            }

            "recall" | "" => {
                if key.is_empty() {
                    return Err(CmdError::MissingArgument("key".to_string()));
                }
                ctx.viewer
                    .scene_recall(key, true, animate)
                    .map_err(|e| CmdError::InvalidArgument {
                        name: "key".to_string(),
                        reason: e,
                    })?;
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
                    return Err(CmdError::InvalidArgument {
                        name: "key".to_string(),
                        reason: format!("Scene '{}' not found.", key),
                    });
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
                    return Err(CmdError::MissingArgument("key".to_string()));
                }
                let new_key = args
                    .get_str(2)
                    .or_else(|| args.get_named_str("new_key"))
                    .ok_or_else(|| CmdError::MissingArgument("new_key".to_string()))?;

                ctx.viewer
                    .scene_rename(key, new_key)
                    .map_err(|e| CmdError::InvalidArgument {
                        name: "key".to_string(),
                        reason: e,
                    })?;

                if !ctx.quiet {
                    ctx.print(&format!(" Scene \"{}\" renamed to \"{}\".", key, new_key));
                }
            }

            _ => {
                return Err(CmdError::InvalidArgument {
                    name: "action".to_string(),
                    reason: format!("Unknown action: '{}'", action),
                });
            }
        }

        Ok(())
    }
}
