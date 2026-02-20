//! Autocomplete/Completion State
//!
//! Manages the state for command line autocompletion including suggestions,
//! selection tracking, and popup visibility.

/// Autocomplete/completion state for command line
#[derive(Debug, Clone, Default)]
pub struct CompletionState {
    /// List of current completion suggestions
    pub suggestions: Vec<String>,
    /// Currently selected suggestion index
    pub selected: usize,
    /// Whether the completion popup is visible
    pub visible: bool,
    /// Position in input where completion starts (byte offset)
    pub start_pos: usize,
    /// Whether to scroll to the selected item (set on arrow key navigation, cleared after scroll)
    pub scroll_to_selected: bool,
}

impl CompletionState {
    /// Create a new empty completion state
    pub fn new() -> Self {
        Self::default()
    }

    /// Reset the completion state (hide popup and clear suggestions)
    pub fn reset(&mut self) {
        self.suggestions.clear();
        self.selected = 0;
        self.visible = false;
        self.start_pos = 0;
    }

    /// Set new suggestions and show the popup
    pub fn show(&mut self, start_pos: usize, suggestions: Vec<String>) {
        if suggestions.is_empty() {
            self.reset();
        } else {
            self.start_pos = start_pos;
            self.suggestions = suggestions;
            self.selected = 0;
            self.visible = true;
        }
    }

    /// Move selection to next item (wrapping)
    pub fn select_next(&mut self) {
        if !self.suggestions.is_empty() {
            self.selected = (self.selected + 1) % self.suggestions.len();
            self.scroll_to_selected = true;
        }
    }

    /// Move selection to previous item (wrapping)
    pub fn select_previous(&mut self) {
        if !self.suggestions.is_empty() {
            self.selected = self.selected
                .checked_sub(1)
                .unwrap_or(self.suggestions.len() - 1);
            self.scroll_to_selected = true;
        }
    }

    /// Get the currently selected suggestion
    pub fn selected_suggestion(&self) -> Option<&str> {
        self.suggestions.get(self.selected).map(|s| s.as_str())
    }

    /// Apply the selected suggestion to the input string
    /// Returns true if a completion was applied
    pub fn apply_to_input(&mut self, input: &mut String) -> bool {
        if let Some(suggestion) = self.suggestions.get(self.selected).cloned() {
            input.truncate(self.start_pos);
            input.push_str(&suggestion);
            // Add space after command name completion only (not for arguments or paths)
            if self.start_pos == 0 && !suggestion.ends_with('/') && !suggestion.contains('/') {
                input.push(' ');
            }
            self.reset();
            true
        } else {
            false
        }
    }
}
