//! IPC Protocol definitions
//!
//! Defines the messages exchanged between the GUI and external clients (like pymol-python).
//! The protocol is bidirectional:
//! - Client → GUI: Send commands, register external commands, respond to callbacks
//! - GUI → Client: Send responses, request callback execution for external commands

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

/// Output message from client (for GUI output view)
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

/// Message FROM client TO GUI
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "type")]
pub enum IpcRequest {
    /// Execute a command string (parsed by CommandExecutor)
    Execute {
        /// Request ID for matching responses
        id: u64,
        /// Command string to execute
        command: String,
    },

    /// Register an external command (appears in autocomplete)
    /// When invoked from GUI command line, sends CallbackRequest to client
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
    /// Includes captured output from execution
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

    /// Get current state (objects, settings, etc.)
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

/// Message FROM GUI TO client
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
    /// Sent when user invokes a registered command from GUI command line
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_request_serialization() {
        let req = IpcRequest::Execute {
            id: 1,
            command: "load protein.pdb".to_string(),
        };
        let json = serde_json::to_string(&req).unwrap();
        assert!(json.contains("Execute"));
        assert!(json.contains("load protein.pdb"));

        let parsed: IpcRequest = serde_json::from_str(&json).unwrap();
        if let IpcRequest::Execute { id, command } = parsed {
            assert_eq!(id, 1);
            assert_eq!(command, "load protein.pdb");
        } else {
            panic!("Wrong variant");
        }
    }

    #[test]
    fn test_response_serialization() {
        let resp = IpcResponse::CallbackRequest {
            id: 42,
            name: "highlight".to_string(),
            args: vec!["chain".to_string(), "A".to_string()],
        };
        let json = serde_json::to_string(&resp).unwrap();
        assert!(json.contains("CallbackRequest"));

        let parsed: IpcResponse = serde_json::from_str(&json).unwrap();
        if let IpcResponse::CallbackRequest { id, name, args } = parsed {
            assert_eq!(id, 42);
            assert_eq!(name, "highlight");
            assert_eq!(args, vec!["chain", "A"]);
        } else {
            panic!("Wrong variant");
        }
    }

    #[test]
    fn test_output_message() {
        let msg = OutputMessage::error("Something went wrong");
        let json = serde_json::to_string(&msg).unwrap();
        let parsed: OutputMessage = serde_json::from_str(&json).unwrap();
        assert_eq!(parsed.text, "Something went wrong");
        assert_eq!(parsed.kind, OutputKind::Error);
    }
}
