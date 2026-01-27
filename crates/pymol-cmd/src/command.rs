//! Command trait and registry
//!
//! Defines the interface for commands and the registry that maps names to implementations.

use std::sync::Arc;

use ahash::AHashMap;

// Re-export ViewerLike from pymol-scene
pub use pymol_scene::ViewerLike;

use crate::args::{ArgDef, ParsedCommand};
use crate::error::CmdResult;

/// Command execution context
///
/// Provides access to the viewer state and execution options.
/// Generic over `V: ViewerLike` to support different viewer implementations.
pub struct CommandContext<'a, V: ViewerLike + ?Sized> {
    /// Reference to the viewer (scene, objects, camera, settings)
    pub viewer: &'a mut V,
    /// Whether to suppress output messages
    pub quiet: bool,
    /// Whether to log this command
    pub log: bool,
}

impl<'a, V: ViewerLike + ?Sized> CommandContext<'a, V> {
    /// Create a new command context
    pub fn new(viewer: &'a mut V) -> Self {
        Self {
            viewer,
            quiet: false,
            log: true,
        }
    }

    /// Create a quiet context (no output)
    pub fn quiet(viewer: &'a mut V) -> Self {
        Self {
            viewer,
            quiet: true,
            log: true,
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

    /// Print a message (unless quiet mode is enabled)
    pub fn print(&self, msg: &str) {
        if !self.quiet {
            log::info!("{}", msg);
        }
    }

    /// Print an error message (even in quiet mode)
    pub fn error(&self, msg: &str) {
        log::error!("{}", msg);
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
    /// Uses a lifetime parameter to allow the context to borrow from non-'static viewers.
    fn execute<'a>(&self, ctx: &mut CommandContext<'a, dyn ViewerLike + 'a>, args: &ParsedCommand) -> CmdResult;

    /// Get help text for this command
    fn help(&self) -> &str {
        "No help available."
    }

    /// Get argument definitions for validation and help
    fn args(&self) -> &[ArgDef] {
        &[]
    }

    /// Get list of command aliases
    fn aliases(&self) -> &[&str] {
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

    /// Register a command with an Arc
    pub fn register_arc(&mut self, cmd: Arc<dyn Command>) {
        let name = cmd.name().to_string();

        // Register aliases
        for alias in cmd.aliases() {
            self.aliases.insert(alias.to_string(), name.clone());
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
}

/// Helper macro to create simple commands from closures
#[macro_export]
macro_rules! simple_command {
    ($name:expr, $help:expr, |$ctx:ident, $args:ident| $body:expr) => {{
        struct SimpleCommand;
        impl $crate::Command for SimpleCommand {
            fn name(&self) -> &str {
                $name
            }
            fn help(&self) -> &str {
                $help
            }
            fn execute(
                &self,
                $ctx: &mut $crate::CommandContext,
                $args: &$crate::ParsedCommand,
            ) -> $crate::CmdResult {
                $body
            }
        }
        SimpleCommand
    }};
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

        fn execute<'a>(&self, _ctx: &mut CommandContext<'a, dyn ViewerLike + 'a>, _args: &ParsedCommand) -> CmdResult {
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
