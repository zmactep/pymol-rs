//! Command Line Model
//!
//! Pure domain model for command line input and history navigation.
//! No egui dependency — testable, serializable, headless-compatible.

/// Command line input model (history, current input, browsing state)
#[derive(Debug)]
pub struct CommandLineModel {
    /// Current command input text
    pub input: String,
    /// Command history
    pub history: Vec<String>,
    /// Current position in command history (None = not browsing)
    pub history_index: Option<usize>,
    /// Saved input when browsing history
    pub saved_input: String,
}

impl Default for CommandLineModel {
    fn default() -> Self {
        Self::new()
    }
}

impl CommandLineModel {
    /// Create a new command line model with default values
    pub fn new() -> Self {
        Self {
            input: String::new(),
            history: Vec::new(),
            history_index: None,
            saved_input: String::new(),
        }
    }

    /// Add command to history
    pub fn add_to_history(&mut self, command: String) {
        // Don't add empty commands or duplicates of last command
        if command.is_empty() {
            return;
        }
        if self.history.last() == Some(&command) {
            return;
        }
        self.history.push(command);
    }

    /// Navigate to previous command in history
    pub fn history_previous(&mut self) {
        if self.history.is_empty() {
            return;
        }

        match self.history_index {
            None => {
                // Start browsing - save current input
                self.saved_input = self.input.clone();
                self.history_index = Some(self.history.len() - 1);
                self.input = self.history.last().unwrap().clone();
            }
            Some(0) => {
                // Already at oldest - do nothing
            }
            Some(idx) => {
                self.history_index = Some(idx - 1);
                self.input = self.history[idx - 1].clone();
            }
        }
    }

    /// Navigate to next command in history
    pub fn history_next(&mut self) {
        match self.history_index {
            None => {
                // Not browsing - do nothing
            }
            Some(idx) => {
                if idx + 1 >= self.history.len() {
                    // At end - restore saved input
                    self.history_index = None;
                    self.input = std::mem::take(&mut self.saved_input);
                } else {
                    self.history_index = Some(idx + 1);
                    self.input = self.history[idx + 1].clone();
                }
            }
        }
    }

    /// Reset history browsing state
    pub fn reset_history_browse(&mut self) {
        self.history_index = None;
        self.saved_input.clear();
    }

    /// Take the current command input and clear it
    pub fn take_command(&mut self) -> String {
        self.reset_history_browse();
        std::mem::take(&mut self.input)
    }
}
