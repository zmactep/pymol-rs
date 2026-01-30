//! PyMOL-RS GUI Application
//!
//! This crate provides a PyMOL-like graphical user interface for molecular visualization.
//! It combines the rendering capabilities of `pymol-scene` with an egui-based UI overlay.
//!
//! ## Features
//!
//! - 3D molecular viewer with wgpu rendering
//! - Command line interface for PyMOL commands
//! - Output/log panel showing command results
//! - Object list panel with visibility controls
//! - Control buttons (Reset, Zoom, Orient, etc.)
//! - Mouse mode information panel
//! - State/playback controls for multi-state molecules

pub mod app;
pub mod async_tasks;
pub mod fetch;
pub mod ipc;
pub mod state;
pub mod ui;
pub mod view;
pub mod viewer_adapter;

pub use app::App;
pub use viewer_adapter::ViewerAdapter;
pub use async_tasks::{AsyncTask, TaskContext, TaskResult, TaskRunner};
pub use fetch::{FetchResult, FetchTask};
pub use ipc::{
    ExternalCommandRegistry, IpcCallbackResult, IpcCallbackTask, IpcRequest, IpcResponse,
    IpcServer, OutputKind as IpcOutputKind, OutputMessage as IpcOutputMessage,
};

// Re-export state types
pub use state::{
    AppState, CommandLineState, CompletionState, OutputBufferState, OutputKind, OutputMessage,
};

// Re-export view types
pub use view::{AppView, UiConfig};
