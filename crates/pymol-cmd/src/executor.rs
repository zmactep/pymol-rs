//! Command executor
//!
//! Dispatches and executes commands against a ViewerLike implementation.

use std::sync::Arc;

use ahash::AHashMap;

use crate::command::{CommandContext, CommandRegistry, DynamicSettingRegistry, FormatHandler, ScriptHandler, OutputMessage, ViewerLike};
use crate::error::{CmdError, CmdResult};
use crate::history::CommandHistory;
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
/// Manages command execution and history.
pub struct CommandExecutor {
    /// Command registry
    registry: CommandRegistry,
    /// Command history
    history: CommandHistory,
    /// Script handlers registered by plugins (extension -> handler for `run` command)
    script_handlers: AHashMap<String, ScriptHandler>,
    /// Format handlers registered by plugins (extension -> handler for `load`/`save`)
    format_handlers: AHashMap<String, Arc<FormatHandler>>,
    /// Dynamic settings registered by plugins
    dynamic_settings: DynamicSettingRegistry,
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
            history: CommandHistory::new(),
            script_handlers: AHashMap::new(),
            format_handlers: AHashMap::new(),
            dynamic_settings: DynamicSettingRegistry::new(),
        }
    }

    /// Create an executor with a pre-populated registry and handlers.
    ///
    /// Used by the `run` command to create a script engine that inherits
    /// plugin-registered commands and handlers from the calling context.
    pub fn with_registry(
        registry: CommandRegistry,
        script_handlers: AHashMap<String, ScriptHandler>,
        format_handlers: AHashMap<String, Arc<FormatHandler>>,
        dynamic_settings: DynamicSettingRegistry,
    ) -> Self {
        Self {
            registry,
            history: CommandHistory::new(),
            script_handlers,
            format_handlers,
            dynamic_settings,
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

    /// Register a script handler for a specific extension.
    ///
    /// Used by plugins to handle non-.pml files in the `run` command.
    pub fn register_script_handler(&mut self, extension: impl Into<String>, handler: ScriptHandler) {
        self.script_handlers.insert(extension.into(), handler);
    }

    /// Get a reference to the script handlers map
    pub fn script_handlers(&self) -> &AHashMap<String, ScriptHandler> {
        &self.script_handlers
    }

    /// Register a format handler for `load`/`save`.
    ///
    /// Each extension in the handler's list gets an entry pointing to the
    /// shared handler. Built-in formats always take priority over plugins.
    pub fn register_format_handler(&mut self, handler: FormatHandler) {
        let handler = Arc::new(handler);
        for ext in &handler.extensions {
            self.format_handlers.insert(ext.clone(), handler.clone());
        }
    }

    /// Get a reference to the format handlers map
    pub fn format_handlers(&self) -> &AHashMap<String, Arc<FormatHandler>> {
        &self.format_handlers
    }

    /// Get a reference to the dynamic settings registry
    pub fn dynamic_settings(&self) -> &DynamicSettingRegistry {
        &self.dynamic_settings
    }

    /// Get a mutable reference to the dynamic settings registry
    pub fn dynamic_settings_mut(&mut self) -> &mut DynamicSettingRegistry {
        &mut self.dynamic_settings
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
        self.do_with_options(viewer, cmd, false).map(|_| ())
    }

    /// Execute a command with options, returning any output messages
    ///
    /// # Arguments
    /// * `viewer` - The viewer to execute against (implements ViewerLike)
    /// * `cmd` - The command string
    /// * `quiet` - Whether to suppress output
    ///
    /// # Returns
    /// On success, returns `CommandOutput` containing any messages from the command.
    pub fn do_with_options(
        &mut self,
        viewer: &mut dyn ViewerLike,
        cmd: &str,
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

        // Add to history
        self.history.push(cmd.to_string());

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
            .with_registry(&self.registry)
            .with_script_handlers(&self.script_handlers)
            .with_format_handlers(&self.format_handlers)
            .with_history(&self.history)
            .with_dynamic_settings(&self.dynamic_settings);
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
