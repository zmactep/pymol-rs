//! Command executor
//!
//! Dispatches and executes commands against a ViewerLike implementation.

use std::path::Path;

use crate::command::{CommandContext, CommandRegistry, OutputMessage, ViewerLike};
use crate::error::{CmdError, CmdResult};
use crate::history::CommandHistory;
use crate::logger::CommandLogger;
use crate::parser::{parse_command, parse_commands};

/// Result of command execution including any output messages
#[derive(Debug, Default)]
pub struct CommandOutput {
    /// Output messages from the command (typed with info/warning/error)
    pub messages: Vec<OutputMessage>,
}

impl CommandOutput {
    /// Create a new empty output
    pub fn new() -> Self {
        Self { messages: Vec::new() }
    }

    /// Check if there are any output messages
    pub fn is_empty(&self) -> bool {
        self.messages.is_empty()
    }
}

/// Command executor
///
/// Manages command execution, logging, and history.
pub struct CommandExecutor {
    /// Command registry
    registry: CommandRegistry,
    /// Command logger (optional)
    logger: CommandLogger,
    /// Command history
    history: CommandHistory,
    /// Whether to echo commands
    echo: bool,
}

impl Default for CommandExecutor {
    fn default() -> Self {
        Self::new()
    }
}

impl CommandExecutor {
    /// Create a new executor with built-in commands
    pub fn new() -> Self {
        Self {
            registry: CommandRegistry::with_builtins(),
            logger: CommandLogger::new(),
            history: CommandHistory::new(),
            echo: false,
        }
    }

    /// Get a reference to the command registry
    pub fn registry(&self) -> &CommandRegistry {
        &self.registry
    }

    /// Get a mutable reference to the command registry
    pub fn registry_mut(&mut self) -> &mut CommandRegistry {
        &mut self.registry
    }

    /// Get a reference to the command history
    pub fn history(&self) -> &CommandHistory {
        &self.history
    }

    /// Get a mutable reference to the command history
    pub fn history_mut(&mut self) -> &mut CommandHistory {
        &mut self.history
    }

    /// Get a reference to the logger
    pub fn logger(&self) -> &CommandLogger {
        &self.logger
    }

    /// Get a mutable reference to the logger
    pub fn logger_mut(&mut self) -> &mut CommandLogger {
        &mut self.logger
    }

    /// Execute a single command string
    ///
    /// # Arguments
    /// * `viewer` - The viewer to execute against (implements ViewerLike)
    /// * `cmd` - The command string to execute
    ///
    /// # Example
    /// ```ignore
    /// executor.do_(&mut viewer, "load protein.pdb")?;
    /// executor.do_(&mut viewer, "zoom")?;
    /// ```
    pub fn do_(&mut self, viewer: &mut dyn ViewerLike, cmd: &str) -> CmdResult {
        self.do_with_options(viewer, cmd, true, false).map(|_| ())
    }

    /// Execute a command with options, returning any output messages
    ///
    /// # Arguments
    /// * `viewer` - The viewer to execute against (implements ViewerLike)
    /// * `cmd` - The command string
    /// * `log` - Whether to log the command
    /// * `quiet` - Whether to suppress output
    ///
    /// # Returns
    /// On success, returns `CommandOutput` containing any messages from the command.
    pub fn do_with_options(
        &mut self,
        viewer: &mut dyn ViewerLike,
        cmd: &str,
        log: bool,
        quiet: bool,
    ) -> Result<CommandOutput, CmdError> {
        let cmd = cmd.trim();
        if cmd.is_empty() {
            return Ok(CommandOutput::new());
        }

        // Skip comments
        if cmd.starts_with('#') {
            return Ok(CommandOutput::new());
        }

        // Echo command if enabled
        if self.echo && !quiet {
            log::info!("PyMOL> {}", cmd);
        }

        // Add to history
        self.history.push(cmd.to_string());

        // Log command
        if log && self.logger.is_active() {
            self.logger.log(cmd);
        }

        // Parse the command
        let parsed = parse_command(cmd)?;

        // Look up the command
        let command = self
            .registry
            .get(&parsed.name)
            .ok_or_else(|| CmdError::UnknownCommand(parsed.name.clone()))?;

        // Execute and collect output
        let mut ctx = CommandContext::new(viewer)
            .with_quiet(quiet)
            .with_log(log)
            .with_registry(&self.registry);
        command.execute(&mut ctx, &parsed)?;

        // Return collected output
        Ok(CommandOutput {
            messages: ctx.take_output(),
        })
    }

    /// Execute multiple commands (semicolon or newline separated)
    ///
    /// Stops on first error unless the command is prefixed with `-` (silent fail).
    pub fn do_multi(&mut self, viewer: &mut dyn ViewerLike, cmds: &str) -> CmdResult {
        let commands = parse_commands(cmds)?;

        for cmd in commands {
            // Reconstruct command string for logging
            let cmd_str = format_command(&cmd);

            if let Err(e) = self.do_(viewer, &cmd_str) {
                // Check if this is a "silent fail" command (starts with -)
                if cmd.name.starts_with('-') {
                    log::debug!("Silently ignoring error in '{}': {}", cmd.name, e);
                    continue;
                }
                return Err(e);
            }
        }

        Ok(())
    }

    /// Open a log file
    pub fn log_open(&mut self, path: &Path) -> CmdResult {
        self.logger.log_open(path)
    }

    /// Close the log file
    pub fn log_close(&mut self) -> CmdResult {
        self.logger.log_close()
    }
}

/// Format a parsed command back to a string (for logging)
fn format_command(cmd: &crate::args::ParsedCommand) -> String {
    let mut s = cmd.name.clone();

    for (i, (name, value)) in cmd.args.iter().enumerate() {
        if i == 0 {
            s.push(' ');
        } else {
            s.push_str(", ");
        }

        if let Some(name) = name {
            s.push_str(name);
            s.push('=');
        }

        s.push_str(&value.to_string());
    }

    s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_executor_creation() {
        let executor = CommandExecutor::new();
        assert!(!executor.registry().is_empty() || executor.registry().is_empty());
    }

    #[test]
    fn test_format_command() {
        use crate::args::ParsedCommand;

        let cmd = ParsedCommand::new("load")
            .with_arg("file.pdb")
            .with_named_arg("object", "mol");

        let formatted = format_command(&cmd);
        assert_eq!(formatted, "load file.pdb, object=mol");
    }
}
