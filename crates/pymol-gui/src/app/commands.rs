//! Command Execution
//!
//! Parsing, dispatching, and executing PyMOL commands.
//! Also handles async task result processing.

use pymol_cmd::{CommandAction, MessageKind};
use pymol_scene::SessionAdapter;

use pymol_framework::message::AppMessage;

use super::App;

impl App {
    /// Execute a command using the CommandExecutor with SessionAdapter
    ///
    /// Returns Ok(()) on success, Err(message) on failure.
    /// When `silent` is true, no output is written to the GUI output buffer.
    pub(crate) fn execute_command(&mut self, cmd: &str, silent: bool) -> Result<(), String> {
        let cmd = cmd.trim();
        if cmd.is_empty() || cmd.starts_with('#') {
            return Ok(());
        }

        // Echo the command to output
        if !silent {
            self.bus.send(AppMessage::PrintCommand(format!("PyMOL> {}", cmd)));
        }

        self.execute_builtin_command_internal(cmd, silent)
    }

    /// Execute a built-in command (internal helper)
    ///
    /// Returns Ok(()) on success, Err(message) on failure.
    /// When `silent` is true, no output is written to the GUI output buffer.
    pub(crate) fn execute_builtin_command_internal(&mut self, cmd: &str, silent: bool) -> Result<(), String> {
        // Calculate default size for PNG capture from viewport or window
        let default_size = self.view.viewport_rect
            .map(|r| (r.width().max(1.0) as u32, r.height().max(1.0) as u32))
            .unwrap_or((1024, 768));

        // Create a SessionAdapter that wraps our state.
        // Use a proxy bool because SessionAdapter expects `&mut bool` for needs_redraw,
        // but we store scene_dirty on App (which is already borrowed via &mut self.state).
        let mut dirty_proxy = self.scene_dirty;
        let task_runner = &self.task_runner;
        let mut adapter = SessionAdapter {
            session: &mut self.state,
            render_context: self.view.gpu.render_context.as_ref(),
            default_size,
            needs_redraw: &mut dirty_proxy,
            async_fetch_fn: Some(Box::new(|code: &str, name: &str, format: u8| {
                let fmt = match format {
                    1 => pymol_io::FetchFormat::Pdb,
                    0 => pymol_io::FetchFormat::Cif,
                    _ => pymol_io::FetchFormat::Bcif,
                };
                task_runner.spawn(crate::fetch::FetchTask::new(
                    code.to_string(),
                    name.to_string(),
                    fmt,
                ));
                true
            })),
        };

        // Temporarily take the executor out to avoid a borrow conflict
        // (adapter already borrows self.state mutably).
        let mut executor = std::mem::replace(&mut self.executor, pymol_cmd::CommandExecutor::new());
        let result = executor.do_with_options(&mut adapter, cmd, false);
        self.executor = executor;

        // Drop the adapter so we can read dirty_proxy
        drop(adapter);

        // Write back the proxy — if the command dirtied the scene, propagate it
        if dirty_proxy {
            self.scene_dirty = true;
        }

        match result {
            Ok(output) => {
                // Display any output messages from the command with appropriate styling
                if !silent {
                    for msg in output.messages {
                        match msg.kind {
                            MessageKind::Info => self.bus.print_info(msg.text),
                            MessageKind::Warning => self.bus.print_warning(msg.text),
                            MessageKind::Error => self.bus.print_error(msg.text),
                        }
                    }
                }
                // Dispatch side-effect actions requested by the command
                for action in output.actions {
                    match action {
                        CommandAction::ShowPanel(name) => self.bus.send(AppMessage::ShowPanel(name)),
                        CommandAction::HidePanel(name) => self.bus.send(AppMessage::HidePanel(name)),
                        CommandAction::Quit => self.frame.quit_requested = true,
                    }
                }
                Ok(())
            }
            Err(e) => {
                // Print the error to the GUI output
                let msg = format!("{}", e);
                if !silent {
                    self.bus.print_error(msg.clone());
                }
                Err(msg)
            }
        }
    }

    /// Process completed async tasks
    ///
    /// This should be called each frame to handle results from background operations
    /// like fetch, save, etc.
    pub(crate) fn process_async_tasks(&mut self) {
        let mut had_results = false;
        while let Some(result) = self.task_runner.poll() {
            // Each result knows how to apply itself via the TaskResult trait
            result.apply(self);
            had_results = true;
        }
        // Request window redraw if any tasks completed
        if had_results {
            self.mark_dirty();
        }
    }
}
