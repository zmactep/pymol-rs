//! IPC (Inter-Process Communication) module
//!
//! Provides an optional IPC server for external control of the GUI application.
//! This allows external processes (like pymol-python) to send commands to the GUI.
//!
//! The IPC server is completely Python-agnostic - it just processes JSON messages
//! according to the defined protocol.

mod callback;
mod protocol;
mod registry;
mod server;

pub use callback::{IpcCallbackResult, IpcCallbackTask};
pub use protocol::{IpcRequest, IpcResponse, OutputKind, OutputMessage};
pub use registry::ExternalCommandRegistry;
pub use server::IpcServer;
