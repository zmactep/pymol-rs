//! IPC Protocol definitions
//!
//! These definitions mirror the protocol in pymol-gui/src/ipc/protocol.rs
//! to ensure compatibility between Python and GUI.

use serde::{Deserialize, Serialize};

/// Output message kind for categorizing output
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum OutputKind {
    /// Informational message (default)
    Info,
    /// Warning message
    Warning,
    /// Error message
    Error,
}

impl Default for OutputKind {
    fn default() -> Self {
        Self::Info
    }
}

/// Output message for GUI output view
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OutputMessage {
    /// The message text
    pub text: String,
    /// The message kind
    #[serde(default)]
    pub kind: OutputKind,
}

impl OutputMessage {
    /// Create a new info message
    pub fn info(text: impl Into<String>) -> Self {
        Self {
            text: text.into(),
            kind: OutputKind::Info,
        }
    }

    /// Create a new warning message
    pub fn warning(text: impl Into<String>) -> Self {
        Self {
            text: text.into(),
            kind: OutputKind::Warning,
        }
    }

    /// Create a new error message
    pub fn error(text: impl Into<String>) -> Self {
        Self {
            text: text.into(),
            kind: OutputKind::Error,
        }
    }
}

/// Message FROM client (Python) TO GUI
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "type")]
pub enum IpcRequest {
    /// Execute a command string
    Execute {
        /// Request ID for matching responses
        id: u64,
        /// Command string to execute
        command: String,
    },

    /// Register an external command (appears in GUI autocomplete)
    RegisterCommand {
        /// Command name
        name: String,
        /// Optional help text
        help: Option<String>,
    },

    /// Unregister an external command
    UnregisterCommand {
        /// Command name to unregister
        name: String,
    },

    /// Response to a CallbackRequest from GUI
    CallbackResponse {
        /// Request ID this is responding to
        id: u64,
        /// Whether execution succeeded
        success: bool,
        /// Error message if failed
        error: Option<String>,
        /// Captured output (print statements, warnings, errors)
        output: Vec<OutputMessage>,
    },

    /// Get current state
    GetState {
        /// Request ID
        id: u64,
    },

    /// Get object names
    GetNames {
        /// Request ID
        id: u64,
    },

    /// Count atoms matching a selection
    CountAtoms {
        /// Request ID
        id: u64,
        /// Selection expression
        selection: String,
    },

    /// Close the GUI application
    Quit,

    /// Health check / keepalive
    Ping {
        /// Request ID
        id: u64,
    },

    /// Show the application window (make it visible)
    ShowWindow {
        /// Request ID
        id: u64,
    },

    /// Hide the application window (make it invisible)
    HideWindow {
        /// Request ID
        id: u64,
    },
}

/// Message FROM GUI TO client (Python)
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "type")]
pub enum IpcResponse {
    /// Command executed successfully
    Ok {
        /// Request ID this is responding to
        id: u64,
    },

    /// Command failed
    Error {
        /// Request ID this is responding to
        id: u64,
        /// Error message
        message: String,
    },

    /// Return value (for queries)
    Value {
        /// Request ID this is responding to
        id: u64,
        /// The value as JSON
        value: serde_json::Value,
    },

    /// GUI requests client to execute an external command or script
    CallbackRequest {
        /// Request ID for matching responses
        id: u64,
        /// Command name
        name: String,
        /// Command arguments
        args: Vec<String>,
    },

    /// Pong response to Ping
    Pong {
        /// Request ID this is responding to
        id: u64,
    },

    /// GUI is closing
    Closing,
}

impl IpcRequest {
    /// Get the request ID if this request has one
    pub fn id(&self) -> Option<u64> {
        match self {
            Self::Execute { id, .. } => Some(*id),
            Self::RegisterCommand { .. } => None,
            Self::UnregisterCommand { .. } => None,
            Self::CallbackResponse { id, .. } => Some(*id),
            Self::GetState { id } => Some(*id),
            Self::GetNames { id } => Some(*id),
            Self::CountAtoms { id, .. } => Some(*id),
            Self::Quit => None,
            Self::Ping { id } => Some(*id),
            Self::ShowWindow { id } => Some(*id),
            Self::HideWindow { id } => Some(*id),
        }
    }
}

impl IpcResponse {
    /// Get the request ID if this response has one
    pub fn id(&self) -> Option<u64> {
        match self {
            Self::Ok { id } => Some(*id),
            Self::Error { id, .. } => Some(*id),
            Self::Value { id, .. } => Some(*id),
            Self::CallbackRequest { id, .. } => Some(*id),
            Self::Pong { id } => Some(*id),
            Self::Closing => None,
        }
    }
}
