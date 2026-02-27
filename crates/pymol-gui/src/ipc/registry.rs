//! External Command Registry
//!
//! Tracks commands registered via IPC that should be forwarded back to the client
//! for execution (callbacks). Help text is stored separately in the CommandRegistry.

use std::collections::HashSet;

/// Registry for external commands registered via IPC
///
/// When a command in this registry is invoked from the GUI command line,
/// the GUI sends a CallbackRequest to the IPC client instead of executing
/// it internally.
#[derive(Debug, Default)]
pub struct ExternalCommandRegistry {
    /// Registered command names
    commands: HashSet<String>,
}

impl ExternalCommandRegistry {
    /// Create a new empty registry
    pub fn new() -> Self {
        Self {
            commands: HashSet::new(),
        }
    }

    /// Register an external command
    ///
    /// The command will appear in autocomplete and when invoked,
    /// will trigger a callback to the IPC client.
    pub fn register(&mut self, name: String) {
        self.commands.insert(name);
    }

    /// Unregister an external command
    ///
    /// Returns true if the command was registered, false otherwise.
    pub fn unregister(&mut self, name: &str) -> bool {
        self.commands.remove(name)
    }

    /// Check if a command is registered
    pub fn contains(&self, name: &str) -> bool {
        self.commands.contains(name)
    }

    /// Get all registered command names
    pub fn names(&self) -> impl Iterator<Item = &str> {
        self.commands.iter().map(|s| s.as_str())
    }

    /// Get the number of registered commands
    pub fn len(&self) -> usize {
        self.commands.len()
    }

    /// Check if the registry is empty
    pub fn is_empty(&self) -> bool {
        self.commands.is_empty()
    }

    /// Clear all registered commands
    pub fn clear(&mut self) {
        self.commands.clear();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_register_unregister() {
        let mut registry = ExternalCommandRegistry::new();

        registry.register("highlight".to_string());
        assert!(registry.contains("highlight"));

        assert!(registry.unregister("highlight"));
        assert!(!registry.contains("highlight"));
        assert!(!registry.unregister("highlight")); // Already removed
    }

    #[test]
    fn test_names_iterator() {
        let mut registry = ExternalCommandRegistry::new();
        registry.register("cmd1".to_string());
        registry.register("cmd2".to_string());
        registry.register("cmd3".to_string());

        let mut names: Vec<_> = registry.names().collect();
        names.sort();
        assert_eq!(names, vec!["cmd1", "cmd2", "cmd3"]);
    }
}
