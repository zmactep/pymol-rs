//! Output Buffer State
//!
//! Manages the output/log buffer for displaying messages to the user.

use std::collections::VecDeque;

/// Maximum number of lines in the output buffer
const MAX_OUTPUT_LINES: usize = 1000;

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

/// Output message type for colored display
#[derive(Debug, Clone)]
pub struct OutputMessage {
    /// The message text
    pub text: String,
    /// Message type for coloring
    pub kind: OutputKind,
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

/// Output buffer state for log/message display
#[derive(Debug)]
pub struct OutputBufferState {
    /// Output messages buffer (circular)
    pub buffer: VecDeque<OutputMessage>,
    /// Whether to auto-scroll output to bottom
    pub auto_scroll: bool,
}

impl Default for OutputBufferState {
    fn default() -> Self {
        Self::new()
    }
}

impl OutputBufferState {
    /// Create a new output buffer state with default values
    pub fn new() -> Self {
        let mut state = Self {
            buffer: VecDeque::with_capacity(MAX_OUTPUT_LINES),
            auto_scroll: true,
        };

        // Add initial welcome messages
        state.add(OutputMessage::info("PyMOL-RS - Molecular Visualization"));
        state.add(OutputMessage::info("Type commands at the prompt below."));

        state
    }

    /// Add an output message to the buffer
    pub fn add(&mut self, message: OutputMessage) {
        self.buffer.push_back(message);
        // Trim if over capacity
        while self.buffer.len() > MAX_OUTPUT_LINES {
            self.buffer.pop_front();
        }
    }

    /// Add a normal output line
    pub fn print(&mut self, text: impl Into<String>) {
        self.add(OutputMessage::normal(text));
    }

    /// Add an info output line
    pub fn print_info(&mut self, text: impl Into<String>) {
        self.add(OutputMessage::info(text));
    }

    /// Add a warning output line
    pub fn print_warning(&mut self, text: impl Into<String>) {
        self.add(OutputMessage::warning(text));
    }

    /// Add an error output line
    pub fn print_error(&mut self, text: impl Into<String>) {
        self.add(OutputMessage::error(text));
    }

    /// Add a command echo line
    pub fn print_command(&mut self, text: impl Into<String>) {
        self.add(OutputMessage::command(text));
    }

    /// Clear the output buffer
    pub fn clear(&mut self) {
        self.buffer.clear();
    }
}
