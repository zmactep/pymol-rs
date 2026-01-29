//! External Command Registry
//!
//! Tracks commands registered via IPC that should be forwarded back to the client
//! for execution (callbacks).

use std::collections::HashMap;

/// Registry for external commands registered via IPC
///
/// When a command in this registry is invoked from the GUI command line,
/// the GUI sends a CallbackRequest to the IPC client instead of executing
/// it internally.
#[derive(Debug, Default)]
pub struct ExternalCommandRegistry {
    /// command name -> help text
    commands: HashMap<String, Option<String>>,
}

impl ExternalCommandRegistry {
    /// Create a new empty registry
    pub fn new() -> Self {
        Self {
            commands: HashMap::new(),
        }
    }

    /// Register an external command
    ///
    /// The command will appear in autocomplete and when invoked,
    /// will trigger a callback to the IPC client.
    pub fn register(&mut self, name: String, help: Option<String>) {
        self.commands.insert(name, help);
    }

    /// Unregister an external command
    ///
    /// Returns true if the command was registered, false otherwise.
    pub fn unregister(&mut self, name: &str) -> bool {
        self.commands.remove(name).is_some()
    }

    /// Check if a command is registered
    pub fn contains(&self, name: &str) -> bool {
        self.commands.contains_key(name)
    }

    /// Get all registered command names
    pub fn names(&self) -> impl Iterator<Item = &str> {
        self.commands.keys().map(|s| s.as_str())
    }

    /// Get help text for a command
    pub fn help(&self, name: &str) -> Option<&str> {
        self.commands.get(name).and_then(|h| h.as_deref())
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
        
        registry.register("highlight".to_string(), Some("Highlight atoms".to_string()));
        assert!(registry.contains("highlight"));
        assert_eq!(registry.help("highlight"), Some("Highlight atoms"));
        
        assert!(registry.unregister("highlight"));
        assert!(!registry.contains("highlight"));
        assert!(!registry.unregister("highlight")); // Already removed
    }

    #[test]
    fn test_names_iterator() {
        let mut registry = ExternalCommandRegistry::new();
        registry.register("cmd1".to_string(), None);
        registry.register("cmd2".to_string(), None);
        registry.register("cmd3".to_string(), None);
        
        let mut names: Vec<_> = registry.names().collect();
        names.sort();
        assert_eq!(names, vec!["cmd1", "cmd2", "cmd3"]);
    }
}
