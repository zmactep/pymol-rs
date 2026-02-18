//! Command trait and registry
//!
//! Defines the interface for commands and the registry that maps names to implementations.

use std::sync::Arc;

use ahash::AHashMap;

// Re-export ViewerLike from pymol-scene
pub use pymol_scene::ViewerLike;

use crate::args::ParsedCommand;
use crate::error::CmdResult;

// ============================================================================
// Argument hints for completion
// ============================================================================

/// Hint about what type of argument a command expects at a given position.
/// Used by the completion system to provide context-aware suggestions.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ArgHint {
    /// No specific hint (generic argument)
    #[default]
    None,
    /// File or directory path
    Path,
    /// Selection expression
    Selection,
    /// Object name
    Object,
    /// Representation type (cartoon, sticks, spheres, etc.)
    Representation,
    /// Color name or value
    Color,
    /// Setting name
    Setting,
    /// Named selection only (no objects)
    NamedSelection,
}

// ============================================================================
// Output message types
// ============================================================================

/// Kind of command output message
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum MessageKind {
    /// Informational message (default)
    #[default]
    Info,
    /// Warning message
    Warning,
    /// Error message
    Error,
}

/// A typed output message from command execution
#[derive(Debug, Clone)]
pub struct OutputMessage {
    /// The message text
    pub text: String,
    /// The message kind (info, warning, error)
    pub kind: MessageKind,
}

impl OutputMessage {
    /// Create an info message
    pub fn info(text: impl Into<String>) -> Self {
        Self {
            text: text.into(),
            kind: MessageKind::Info,
        }
    }

    /// Create a warning message
    pub fn warning(text: impl Into<String>) -> Self {
        Self {
            text: text.into(),
            kind: MessageKind::Warning,
        }
    }

    /// Create an error message
    pub fn error(text: impl Into<String>) -> Self {
        Self {
            text: text.into(),
            kind: MessageKind::Error,
        }
    }
}

/// Command execution context
///
/// Provides access to the viewer state and execution options.
/// Generic over `V: ViewerLike` to support different viewer implementations.
///
/// The two lifetime parameters are:
/// - `'v` - lifetime of the viewer reference
/// - `'r` - lifetime of the registry reference (may differ from viewer)
pub struct CommandContext<'v, 'r, V: ViewerLike + ?Sized> {
    /// Reference to the viewer (scene, objects, camera, settings)
    pub viewer: &'v mut V,
    /// Whether to suppress output messages
    pub quiet: bool,
    /// Whether to log this command
    pub log: bool,
    /// Collected output messages (for GUI display)
    output_buffer: Vec<OutputMessage>,
    /// Optional reference to command registry (for help lookups)
    registry: Option<&'r CommandRegistry>,
}

impl<'v, 'r, V: ViewerLike + ?Sized> CommandContext<'v, 'r, V> {
    /// Create a new command context
    pub fn new(viewer: &'v mut V) -> Self {
        Self {
            viewer,
            quiet: false,
            log: true,
            output_buffer: Vec::new(),
            registry: None,
        }
    }

    /// Set the quiet flag
    pub fn with_quiet(mut self, quiet: bool) -> Self {
        self.quiet = quiet;
        self
    }

    /// Set the log flag
    pub fn with_log(mut self, log: bool) -> Self {
        self.log = log;
        self
    }

    /// Set the command registry reference
    pub fn with_registry(mut self, registry: &'r CommandRegistry) -> Self {
        self.registry = Some(registry);
        self
    }

    /// Get a reference to the command registry (if available)
    pub fn registry(&self) -> Option<&CommandRegistry> {
        self.registry
    }

    /// Print an info message (unless quiet mode is enabled)
    /// 
    /// The message is both logged and collected in the output buffer
    /// for retrieval by the GUI.
    pub fn print(&mut self, msg: &str) {
        if !self.quiet {
            log::info!("{}", msg);
            self.output_buffer.push(OutputMessage::info(msg));
        }
    }

    /// Print a warning message (unless quiet mode is enabled)
    pub fn print_warning(&mut self, msg: &str) {
        if !self.quiet {
            log::warn!("{}", msg);
            self.output_buffer.push(OutputMessage::warning(msg));
        }
    }

    /// Print an error message (even in quiet mode)
    /// 
    /// Error messages are always shown regardless of quiet mode.
    pub fn print_error(&mut self, msg: &str) {
        log::error!("{}", msg);
        self.output_buffer.push(OutputMessage::error(msg));
    }

    /// Take the collected output messages, clearing the buffer
    pub fn take_output(&mut self) -> Vec<OutputMessage> {
        std::mem::take(&mut self.output_buffer)
    }

    /// Get a reference to the collected output messages
    pub fn output(&self) -> &[OutputMessage] {
        &self.output_buffer
    }
}

/// Trait for command implementations
///
/// Commands receive a context with access to the viewer state and
/// parsed arguments, and return a result indicating success or failure.
pub trait Command: Send + Sync {
    /// Get the command name
    fn name(&self) -> &str;

    /// Execute the command
    /// 
    /// Uses lifetime parameters to allow the context to borrow from non-'static viewers.
    /// - `'v` - lifetime of the viewer reference
    /// - `'r` - lifetime of the registry reference
    fn execute<'v, 'r>(&self, ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>, args: &ParsedCommand) -> CmdResult;

    /// Get help text for this command
    fn help(&self) -> &str {
        "No help available."
    }

    /// Get list of command aliases
    fn aliases(&self) -> &[&str] {
        &[]
    }

    /// Get argument hints for completion
    /// 
    /// Returns a slice of hints indicating what type of argument is expected
    /// at each position. Used by the completion system to provide context-aware
    /// suggestions (e.g., file paths for "load", selections for "select").
    fn arg_hints(&self) -> &[ArgHint] {
        &[]
    }
}

/// Registry mapping command names to implementations
pub struct CommandRegistry {
    /// Commands indexed by name
    commands: AHashMap<String, Arc<dyn Command>>,
    /// Aliases mapping alias -> command name
    aliases: AHashMap<String, String>,
}

impl Default for CommandRegistry {
    fn default() -> Self {
        Self::new()
    }
}

impl CommandRegistry {
    /// Create a new empty registry
    pub fn new() -> Self {
        Self {
            commands: AHashMap::new(),
            aliases: AHashMap::new(),
        }
    }

    /// Create a registry with all built-in commands registered
    pub fn with_builtins() -> Self {
        let mut registry = Self::new();
        crate::commands::register_all(&mut registry);
        registry
    }

    /// Register a command
    ///
    /// Also registers any aliases defined by the command.
    pub fn register<C: Command + 'static>(&mut self, cmd: C) {
        let name = cmd.name().to_string();
        let aliases: Vec<String> = cmd.aliases().iter().map(|s| s.to_string()).collect();
        let cmd = Arc::new(cmd);

        // Register aliases
        for alias in aliases {
            self.aliases.insert(alias, name.clone());
        }

        self.commands.insert(name, cmd);
    }

    /// Add an alias for an existing command
    pub fn add_alias(&mut self, alias: impl Into<String>, command: impl Into<String>) {
        self.aliases.insert(alias.into(), command.into());
    }

    /// Look up a command by name or alias
    pub fn get(&self, name: &str) -> Option<Arc<dyn Command>> {
        // Try direct lookup first
        if let Some(cmd) = self.commands.get(name) {
            return Some(cmd.clone());
        }

        // Try alias lookup
        if let Some(real_name) = self.aliases.get(name) {
            return self.commands.get(real_name).cloned();
        }

        None
    }

    /// Check if a command exists
    pub fn contains(&self, name: &str) -> bool {
        self.commands.contains_key(name) || self.aliases.contains_key(name)
    }

    /// Get all command names (not including aliases)
    pub fn names(&self) -> impl Iterator<Item = &str> {
        self.commands.keys().map(|s| s.as_str())
    }

    /// Get all command names and aliases
    pub fn all_names(&self) -> impl Iterator<Item = &str> {
        self.commands
            .keys()
            .map(|s| s.as_str())
            .chain(self.aliases.keys().map(|s| s.as_str()))
    }

    /// Get the number of registered commands
    pub fn len(&self) -> usize {
        self.commands.len()
    }

    /// Check if the registry is empty
    pub fn is_empty(&self) -> bool {
        self.commands.is_empty()
    }

    /// Remove a command
    pub fn remove(&mut self, name: &str) -> Option<Arc<dyn Command>> {
        // Remove aliases pointing to this command
        self.aliases.retain(|_, v| v != name);
        self.commands.remove(name)
    }

    /// Clear all commands
    pub fn clear(&mut self) {
        self.commands.clear();
        self.aliases.clear();
    }

    /// Get all command names (including aliases) that expect a specific argument hint
    /// at the given position (0-indexed).
    /// 
    /// This is used by the completion system to determine which commands should
    /// trigger specific completion types (e.g., file path completion).
    pub fn commands_with_hint(&self, hint: ArgHint, position: usize) -> Vec<&str> {
        let mut result = Vec::new();
        
        for (name, cmd) in &self.commands {
            let hints = cmd.arg_hints();
            if hints.get(position) == Some(&hint) {
                result.push(name.as_str());
                // Also include aliases for this command
                for (alias, target) in &self.aliases {
                    if target == name {
                        result.push(alias.as_str());
                    }
                }
            }
        }
        
        result
    }

    /// Check if a command expects a specific argument hint at the given position
    pub fn command_has_hint(&self, name: &str, hint: ArgHint, position: usize) -> bool {
        self.get(name)
            .map(|cmd| cmd.arg_hints().get(position) == Some(&hint))
            .unwrap_or(false)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    struct TestCommand {
        name: String,
    }

    impl Command for TestCommand {
        fn name(&self) -> &str {
            &self.name
        }

        fn execute<'v, 'r>(&self, _ctx: &mut CommandContext<'v, 'r, dyn ViewerLike + 'v>, _args: &ParsedCommand) -> CmdResult {
            Ok(())
        }

        fn help(&self) -> &str {
            "Test command"
        }

        fn aliases(&self) -> &[&str] {
            &["test_alias"]
        }
    }

    #[test]
    fn test_registry() {
        let mut registry = CommandRegistry::new();

        registry.register(TestCommand {
            name: "test".to_string(),
        });

        assert!(registry.contains("test"));
        assert!(registry.contains("test_alias"));
        assert!(!registry.contains("unknown"));

        let cmd = registry.get("test").unwrap();
        assert_eq!(cmd.name(), "test");

        let cmd = registry.get("test_alias").unwrap();
        assert_eq!(cmd.name(), "test");
    }
}
