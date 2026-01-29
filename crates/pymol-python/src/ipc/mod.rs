//! IPC (Inter-Process Communication) module for pymol-python
//!
//! Provides an IPC client that communicates with the pymol-rs GUI application
//! via Unix domain sockets.

mod client;
mod listener;
mod protocol;

pub use client::IpcClient;
pub use listener::{CallbackListener, ExtendedCommands};
pub use protocol::{IpcRequest, IpcResponse, OutputKind, OutputMessage};
