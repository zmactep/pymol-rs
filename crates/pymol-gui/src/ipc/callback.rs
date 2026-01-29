//! IPC callback task and result types
//!
//! This module provides async task support for IPC callback operations,
//! allowing the GUI to wait asynchronously for responses from external
//! clients (like pymol-python).

use std::future::Future;
use std::pin::Pin;
use std::time::Duration;

use tokio::sync::oneshot;

use super::protocol::{IpcRequest, OutputKind, OutputMessage};
use crate::async_tasks::{AsyncTask, TaskContext, TaskResult};

/// Timeout for waiting on IPC callback response (30 seconds)
const CALLBACK_TIMEOUT_SECS: u64 = 30;

// ============================================================================
// IpcCallbackResult
// ============================================================================

/// Result of an IPC callback operation
pub struct IpcCallbackResult {
    /// Callback ID for matching
    pub id: u64,
    /// Command that was executed
    pub command: String,
    /// Whether execution succeeded
    pub success: bool,
    /// Error message if failed
    pub error: Option<String>,
    /// Output lines from client (print statements, etc.)
    pub output: Vec<OutputMessage>,
}

impl TaskResult for IpcCallbackResult {
    fn apply(self: Box<Self>, ctx: &mut dyn TaskContext) {
        // Display captured output
        for msg in self.output {
            match msg.kind {
                OutputKind::Info => ctx.print_info(msg.text),
                OutputKind::Warning => ctx.print_warning(msg.text),
                OutputKind::Error => ctx.print_error(msg.text),
            }
        }

        // Display error if failed
        if !self.success {
            if let Some(err) = self.error {
                ctx.print_error(format!("{}: {}", self.command, err));
            }
        }

        ctx.request_redraw();
    }
}

// ============================================================================
// IpcCallbackTask
// ============================================================================

/// Task that waits for an IPC callback response
///
/// This task is created after:
/// 1. Registering the callback with IpcServer
/// 2. Sending the CallbackRequest to the client
///
/// The task then waits on the oneshot receiver for the response.
pub struct IpcCallbackTask {
    /// Callback ID for matching
    pub id: u64,
    /// Command name being executed
    pub command: String,
    /// Receiver for the callback response
    pub receiver: oneshot::Receiver<IpcRequest>,
}

impl IpcCallbackTask {
    /// Create a new IPC callback task
    pub fn new(id: u64, command: String, receiver: oneshot::Receiver<IpcRequest>) -> Self {
        Self {
            id,
            command,
            receiver,
        }
    }
}

impl AsyncTask for IpcCallbackTask {
    fn notification_message(&self) -> String {
        format!("Running {}...", self.command)
    }

    fn execute(self) -> Pin<Box<dyn Future<Output = Box<dyn TaskResult>> + Send>> {
        Box::pin(async move {
            let timeout = Duration::from_secs(CALLBACK_TIMEOUT_SECS);

            let result = match tokio::time::timeout(timeout, self.receiver).await {
                Ok(Ok(IpcRequest::CallbackResponse {
                    id: _,
                    success,
                    error,
                    output,
                })) => IpcCallbackResult {
                    id: self.id,
                    command: self.command,
                    success,
                    error,
                    output,
                },
                Ok(Ok(_)) => {
                    // Unexpected message type
                    IpcCallbackResult {
                        id: self.id,
                        command: self.command,
                        success: false,
                        error: Some("Unexpected IPC response type".to_string()),
                        output: vec![],
                    }
                }
                Ok(Err(_)) => {
                    // Channel closed (client disconnected)
                    IpcCallbackResult {
                        id: self.id,
                        command: self.command,
                        success: false,
                        error: Some("IPC client disconnected".to_string()),
                        output: vec![],
                    }
                }
                Err(_) => {
                    // Timeout
                    IpcCallbackResult {
                        id: self.id,
                        command: self.command,
                        success: false,
                        error: Some(format!(
                            "Callback timed out after {} seconds",
                            CALLBACK_TIMEOUT_SECS
                        )),
                        output: vec![],
                    }
                }
            };

            Box::new(result) as Box<dyn TaskResult>
        })
    }
}
