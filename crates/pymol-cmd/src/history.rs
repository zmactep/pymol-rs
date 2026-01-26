//! Command history tracking
//!
//! Stores executed commands for recall and navigation.

use std::collections::VecDeque;

/// Maximum number of commands to store in history
const DEFAULT_MAX_HISTORY: usize = 1000;

/// Command history for recall and navigation
#[derive(Debug)]
pub struct CommandHistory {
    /// History of executed commands (most recent at back)
    commands: VecDeque<String>,
    /// Maximum number of commands to store
    max_size: usize,
    /// Current position for navigation (None = at end)
    position: Option<usize>,
}

impl Default for CommandHistory {
    fn default() -> Self {
        Self::new()
    }
}

impl CommandHistory {
    /// Create a new empty history
    pub fn new() -> Self {
        Self::with_capacity(DEFAULT_MAX_HISTORY)
    }

    /// Create a new history with specified capacity
    pub fn with_capacity(max_size: usize) -> Self {
        Self {
            commands: VecDeque::with_capacity(max_size.min(1024)),
            max_size,
            position: None,
        }
    }

    /// Add a command to history
    ///
    /// Resets navigation position to the end.
    pub fn push(&mut self, command: String) {
        // Don't add empty commands or duplicates of the last command
        if command.is_empty() {
            return;
        }
        if self.commands.back().map(|s| s.as_str()) == Some(&command) {
            return;
        }

        // Remove oldest if at capacity
        while self.commands.len() >= self.max_size {
            self.commands.pop_front();
        }

        self.commands.push_back(command);
        self.position = None;
    }

    /// Get the previous command in history (going back)
    pub fn previous(&mut self) -> Option<&str> {
        if self.commands.is_empty() {
            return None;
        }

        let new_pos = match self.position {
            None => self.commands.len().saturating_sub(1),
            Some(0) => 0,
            Some(p) => p - 1,
        };

        self.position = Some(new_pos);
        self.commands.get(new_pos).map(|s| s.as_str())
    }

    /// Get the next command in history (going forward)
    pub fn next(&mut self) -> Option<&str> {
        match self.position {
            None => None,
            Some(p) if p + 1 >= self.commands.len() => {
                self.position = None;
                None
            }
            Some(p) => {
                self.position = Some(p + 1);
                self.commands.get(p + 1).map(|s| s.as_str())
            }
        }
    }

    /// Reset navigation position to end
    pub fn reset_position(&mut self) {
        self.position = None;
    }

    /// Get the current navigation position
    pub fn position(&self) -> Option<usize> {
        self.position
    }

    /// Get the number of commands in history
    pub fn len(&self) -> usize {
        self.commands.len()
    }

    /// Check if history is empty
    pub fn is_empty(&self) -> bool {
        self.commands.is_empty()
    }

    /// Clear all history
    pub fn clear(&mut self) {
        self.commands.clear();
        self.position = None;
    }

    /// Get command at index
    pub fn get(&self, index: usize) -> Option<&str> {
        self.commands.get(index).map(|s| s.as_str())
    }

    /// Iterate over all commands (oldest first)
    pub fn iter(&self) -> impl Iterator<Item = &str> {
        self.commands.iter().map(|s| s.as_str())
    }

    /// Get the last N commands
    pub fn last_n(&self, n: usize) -> impl Iterator<Item = &str> {
        let start = self.commands.len().saturating_sub(n);
        self.commands.iter().skip(start).map(|s| s.as_str())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_push_and_navigate() {
        let mut history = CommandHistory::new();

        history.push("cmd1".to_string());
        history.push("cmd2".to_string());
        history.push("cmd3".to_string());

        assert_eq!(history.len(), 3);

        // Navigate back
        assert_eq!(history.previous(), Some("cmd3"));
        assert_eq!(history.previous(), Some("cmd2"));
        assert_eq!(history.previous(), Some("cmd1"));
        assert_eq!(history.previous(), Some("cmd1")); // Stay at beginning

        // Navigate forward
        assert_eq!(history.next(), Some("cmd2"));
        assert_eq!(history.next(), Some("cmd3"));
        assert_eq!(history.next(), None); // At end
    }

    #[test]
    fn test_no_duplicates() {
        let mut history = CommandHistory::new();

        history.push("cmd1".to_string());
        history.push("cmd1".to_string()); // Duplicate, should not be added

        assert_eq!(history.len(), 1);
    }

    #[test]
    fn test_max_capacity() {
        let mut history = CommandHistory::with_capacity(3);

        history.push("cmd1".to_string());
        history.push("cmd2".to_string());
        history.push("cmd3".to_string());
        history.push("cmd4".to_string()); // Should remove cmd1

        assert_eq!(history.len(), 3);
        assert_eq!(history.get(0), Some("cmd2"));
    }
}
