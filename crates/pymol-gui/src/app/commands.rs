//! Command Execution
//!
//! Parsing, dispatching, and executing PyMOL commands (built-in and external).
//! Also handles async task result processing.

use pymol_cmd::{CmdError, MessageKind, parse_command};
use pymol_scene::SessionAdapter;

use crate::ipc::{IpcCallbackTask, IpcResponse};

use super::App;

impl App {
    /// Execute a command using the CommandExecutor with SessionAdapter
    ///
    /// This method is IPC-aware - if the command is an external command
    /// registered via IPC, it will send a callback request to the client.
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
            self.output.print_command(format!("PyMOL> {}", cmd));
        }

        // Parse using pymol-cmd parser for proper handling of quotes, commas, etc.
        let parsed = match parse_command(cmd) {
            Ok(p) => p,
            Err(e) => {
                let msg = format!("Parse error: {}", e);
                if !silent {
                    self.output.print_error(msg.clone());
                }
                return Err(msg);
            }
        };

        // Check if it's an external command (registered via IPC)
        if self.external_commands.contains(&parsed.name) {
            self.execute_external_command(&parsed);
            // External commands are async, errors are handled separately via callback
            return Ok(());
        }

        // Execute as built-in command using CommandExecutor
        self.execute_builtin_command_internal(cmd, silent)
    }

    /// Execute an external command via IPC callback
    ///
    /// This spawns an async task that waits for the callback response.
    fn execute_external_command(&mut self, parsed: &pymol_cmd::ParsedCommand) {
        let server = match self.ipc_server.as_mut() {
            Some(s) => s,
            None => {
                self.output.print_error("No IPC server available");
                return;
            }
        };

        let id = server.next_callback_id();
        // Convert parsed arguments to strings using ArgValue's Display impl
        let args: Vec<String> = parsed.args
            .iter()
            .map(|(_, v)| v.to_string())
            .collect();

        // Create oneshot channel for the callback response
        let (tx, rx) = tokio::sync::oneshot::channel();

        // Register the callback with the server
        server.register_callback(id, tx);

        // Send the CallbackRequest to the client
        let response = IpcResponse::CallbackRequest {
            id,
            name: parsed.name.clone(),
            args,
        };

        if let Err(e) = server.send(response) {
            self.output.print_error(format!("IPC error: {}", e));
            return;
        }

        // Spawn the async task that waits for the response
        let task = IpcCallbackTask::new(id, parsed.name.clone(), rx);
        self.task_runner.spawn(task);
    }

    /// Execute a built-in command (internal helper)
    ///
    /// Returns Ok(()) on success, Err(message) on failure.
    /// When `silent` is true, no output is written to the GUI output buffer.
    fn execute_builtin_command_internal(&mut self, cmd: &str, silent: bool) -> Result<(), String> {
        // Calculate default size for PNG capture from viewport or window
        let default_size = self.view.viewport_rect
            .map(|r| (r.width().max(1.0) as u32, r.height().max(1.0) as u32))
            .unwrap_or((1024, 768));

        // Create a SessionAdapter that wraps our state
        let task_runner = &self.task_runner;
        let mut adapter = SessionAdapter {
            session: &mut self.state,
            render_context: self.view.render_context.as_ref(),
            default_size,
            needs_redraw: &mut self.needs_redraw,
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
        // (adapter already borrows self.state and self.needs_redraw mutably).
        let mut executor = std::mem::replace(&mut self.executor, pymol_cmd::CommandExecutor::new());
        let result = executor.do_with_options(&mut adapter, cmd, true, false);
        self.executor = executor;

        match result {
            Ok(output) => {
                // Display any output messages from the command with appropriate styling
                if !silent {
                    for msg in output.messages {
                        match msg.kind {
                            MessageKind::Info => self.output.print_info(msg.text),
                            MessageKind::Warning => self.output.print_warning(msg.text),
                            MessageKind::Error => self.output.print_error(msg.text),
                        }
                    }
                }
                Ok(())
            }
            Err(CmdError::Aborted) => {
                // Quit/exit command was issued - signal application to close
                self.quit_requested = true;
                Ok(())  // Not an error to the caller
            }
            Err(e) => {
                // Print the error to the GUI output
                let msg = format!("{}", e);
                if !silent {
                    self.output.print_error(msg.clone());
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
            self.request_redraw();
        }
    }
}
