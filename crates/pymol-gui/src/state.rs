//! GUI State Management
//!
//! Maintains UI state separate from the viewer state, including output buffer,
//! command history, and UI visibility settings.

use std::collections::VecDeque;

/// Maximum number of lines in the output buffer
const MAX_OUTPUT_LINES: usize = 1000;

/// Mouse interaction mode
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum MouseMode {
    /// 3-button viewing mode (default)
    #[default]
    ThreeButtonViewing,
    /// 3-button editing mode
    ThreeButtonEditing,
    /// 2-button viewing mode
    TwoButtonViewing,
    /// 1-button viewing mode
    OneButtonViewing,
}

impl MouseMode {
    /// Get display name for the mouse mode
    pub fn display_name(&self) -> &'static str {
        match self {
            MouseMode::ThreeButtonViewing => "3-Button Viewing",
            MouseMode::ThreeButtonEditing => "3-Button Editing",
            MouseMode::TwoButtonViewing => "2-Button Viewing",
            MouseMode::OneButtonViewing => "1-Button Viewing",
        }
    }
}

/// Selection granularity mode
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SelectingMode {
    /// Select atoms
    Atoms,
    /// Select residues (default)
    #[default]
    Residues,
    /// Select chains
    Chains,
    /// Select segments
    Segments,
    /// Select objects
    Objects,
    /// Select molecules
    Molecules,
}

impl SelectingMode {
    /// Get display name for the selecting mode
    pub fn display_name(&self) -> &'static str {
        match self {
            SelectingMode::Atoms => "Atoms",
            SelectingMode::Residues => "Residues",
            SelectingMode::Chains => "Chains",
            SelectingMode::Segments => "Segments",
            SelectingMode::Objects => "Objects",
            SelectingMode::Molecules => "Molecules",
        }
    }

    /// Cycle to the next selecting mode
    pub fn next(&self) -> Self {
        match self {
            SelectingMode::Atoms => SelectingMode::Residues,
            SelectingMode::Residues => SelectingMode::Chains,
            SelectingMode::Chains => SelectingMode::Segments,
            SelectingMode::Segments => SelectingMode::Objects,
            SelectingMode::Objects => SelectingMode::Molecules,
            SelectingMode::Molecules => SelectingMode::Atoms,
        }
    }
}

// ============================================================================
// Autocomplete State
// ============================================================================

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
        }
    }

    /// Move selection to previous item (wrapping)
    pub fn select_previous(&mut self) {
        if !self.suggestions.is_empty() {
            self.selected = self.selected
                .checked_sub(1)
                .unwrap_or(self.suggestions.len() - 1);
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
            // Add space after command completion (not for paths ending in /)
            if !suggestion.ends_with('/') && !suggestion.contains('/') {
                input.push(' ');
            }
            self.reset();
            true
        } else {
            false
        }
    }
}

// ============================================================================
// Output Messages
// ============================================================================

/// Output message type for colored display
#[derive(Debug, Clone)]
pub struct OutputMessage {
    /// The message text
    pub text: String,
    /// Message type for coloring
    pub kind: OutputKind,
}

/// Kind of output message
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputKind {
    /// Normal output
    Normal,
    /// Info/status message
    Info,
    /// Warning message
    Warning,
    /// Error message
    Error,
    /// Command echo
    Command,
}

impl OutputMessage {
    /// Create a normal message
    pub fn normal(text: impl Into<String>) -> Self {
        Self {
            text: text.into(),
            kind: OutputKind::Normal,
        }
    }

    /// Create an info message
    pub fn info(text: impl Into<String>) -> Self {
        Self {
            text: text.into(),
            kind: OutputKind::Info,
        }
    }

    /// Create a warning message
    pub fn warning(text: impl Into<String>) -> Self {
        Self {
            text: text.into(),
            kind: OutputKind::Warning,
        }
    }

    /// Create an error message
    pub fn error(text: impl Into<String>) -> Self {
        Self {
            text: text.into(),
            kind: OutputKind::Error,
        }
    }

    /// Create a command echo message
    pub fn command(text: impl Into<String>) -> Self {
        Self {
            text: text.into(),
            kind: OutputKind::Command,
        }
    }
}

/// GUI state maintained separately from viewer state
#[derive(Debug)]
pub struct GuiState {
    // =========================================================================
    // Output/Log Buffer
    // =========================================================================
    /// Output messages buffer (circular)
    pub output_buffer: VecDeque<OutputMessage>,
    /// Whether to auto-scroll output to bottom
    pub output_auto_scroll: bool,

    // =========================================================================
    // Command Line State
    // =========================================================================
    /// Current command input text
    pub command_input: String,
    /// Command history
    pub command_history: Vec<String>,
    /// Current position in command history (None = not browsing)
    pub history_index: Option<usize>,
    /// Saved input when browsing history
    pub saved_input: String,

    // =========================================================================
    // UI Visibility
    // =========================================================================
    /// Show the output panel at top
    pub show_output_panel: bool,
    /// Show the right-side control panel
    pub show_control_panel: bool,
    /// Output panel height in pixels
    pub output_panel_height: f32,
    /// Right panel width in pixels
    pub right_panel_width: f32,

    // =========================================================================
    // Mouse/Selection Mode
    // =========================================================================
    /// Current mouse interaction mode
    pub mouse_mode: MouseMode,
    /// Current selection granularity
    pub selecting_mode: SelectingMode,

    // =========================================================================
    // Playback State
    // =========================================================================
    /// Whether animation is playing
    pub is_playing: bool,
    /// Current state index (1-based for display)
    pub current_state: usize,
    /// Total number of states
    pub total_states: usize,
    /// Whether rocking animation is enabled
    pub is_rocking: bool,

    // =========================================================================
    // Pending Actions
    // =========================================================================
    /// File path to load on next frame
    pub pending_load_file: Option<String>,

    // =========================================================================
    // Input State Tracking
    // =========================================================================
    /// Whether the command input should request focus (set after command execution)
    pub command_wants_focus: bool,
    /// Whether the command input currently has focus
    pub command_has_focus: bool,
    /// Whether the mouse is currently over the 3D viewport
    pub viewport_hovered: bool,

    // =========================================================================
    // Autocomplete State
    // =========================================================================
    /// Current autocomplete/completion state
    pub completion: CompletionState,
    /// Cached list of command names for autocomplete (populated on startup)
    pub command_names: Vec<String>,

    // =========================================================================
    // Application Control
    // =========================================================================
    /// Whether the quit command was issued
    pub quit_requested: bool,
}

impl Default for GuiState {
    fn default() -> Self {
        Self::new()
    }
}

impl GuiState {
    /// Create a new GUI state with default values
    pub fn new() -> Self {
        let mut state = Self {
            output_buffer: VecDeque::with_capacity(MAX_OUTPUT_LINES),
            output_auto_scroll: true,
            command_input: String::new(),
            command_history: Vec::new(),
            history_index: None,
            saved_input: String::new(),
            show_output_panel: true,
            show_control_panel: true,
            output_panel_height: 150.0,
            right_panel_width: 200.0,
            mouse_mode: MouseMode::default(),
            selecting_mode: SelectingMode::default(),
            is_playing: false,
            current_state: 1,
            total_states: 1,
            is_rocking: false,
            pending_load_file: None,
            command_wants_focus: true, // Focus on startup
            command_has_focus: false,
            viewport_hovered: false,
            completion: CompletionState::new(),
            command_names: Vec::new(), // Populated by App on startup
            quit_requested: false,
        };

        // Add initial welcome messages
        state.add_output(OutputMessage::info("PyMOL-RS - Molecular Visualization"));
        state.add_output(OutputMessage::info("Type commands at the prompt below."));

        state
    }

    /// Add an output message to the buffer
    pub fn add_output(&mut self, message: OutputMessage) {
        self.output_buffer.push_back(message);
        // Trim if over capacity
        while self.output_buffer.len() > MAX_OUTPUT_LINES {
            self.output_buffer.pop_front();
        }
    }

    /// Add a normal output line
    pub fn print(&mut self, text: impl Into<String>) {
        self.add_output(OutputMessage::normal(text));
    }

    /// Add an info output line
    pub fn print_info(&mut self, text: impl Into<String>) {
        self.add_output(OutputMessage::info(text));
    }

    /// Add a warning output line
    pub fn print_warning(&mut self, text: impl Into<String>) {
        self.add_output(OutputMessage::warning(text));
    }

    /// Add an error output line
    pub fn print_error(&mut self, text: impl Into<String>) {
        self.add_output(OutputMessage::error(text));
    }

    /// Add a command echo line
    pub fn print_command(&mut self, text: impl Into<String>) {
        self.add_output(OutputMessage::command(text));
    }

    /// Clear the output buffer
    pub fn clear_output(&mut self) {
        self.output_buffer.clear();
    }

    /// Add command to history
    pub fn add_to_history(&mut self, command: String) {
        // Don't add empty commands or duplicates of last command
        if command.is_empty() {
            return;
        }
        if self.command_history.last() == Some(&command) {
            return;
        }
        self.command_history.push(command);
    }

    /// Navigate to previous command in history
    pub fn history_previous(&mut self) {
        if self.command_history.is_empty() {
            return;
        }

        match self.history_index {
            None => {
                // Start browsing - save current input
                self.saved_input = self.command_input.clone();
                self.history_index = Some(self.command_history.len() - 1);
                self.command_input = self.command_history.last().unwrap().clone();
            }
            Some(0) => {
                // Already at oldest - do nothing
            }
            Some(idx) => {
                self.history_index = Some(idx - 1);
                self.command_input = self.command_history[idx - 1].clone();
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
                if idx + 1 >= self.command_history.len() {
                    // At end - restore saved input
                    self.history_index = None;
                    self.command_input = std::mem::take(&mut self.saved_input);
                } else {
                    self.history_index = Some(idx + 1);
                    self.command_input = self.command_history[idx + 1].clone();
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
        std::mem::take(&mut self.command_input)
    }

    /// Request focus on the command input for the next frame
    pub fn request_command_focus(&mut self) {
        self.command_wants_focus = true;
    }
}
