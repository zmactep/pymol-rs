//! PyMOL-RS GUI Application
//!
//! This crate provides a PyMOL-like graphical user interface for molecular visualization.
//! It combines the rendering capabilities of `pymol-scene` with an egui-based UI overlay.
//!
//! ## Architecture
//!
//! The GUI is built around a component system:
//! - [`Component`] — trait for self-contained UI panels (from `pymol-framework`)
//! - [`SharedContext`] — read-only application state shared with components
//! - [`ComponentStore`] — holds all components with typed + dynamic access
//! - [`layout::Layout`] — configurable panel placement around the 3D viewport
//! - [`MessageBus`](pymol_framework::message::MessageBus) — per-frame message queue for inter-component communication

pub mod app;
pub mod async_tasks;
pub mod components;
pub mod fetch;
pub mod ipc;
pub mod layout;
pub mod model;
pub mod plugin_manager;
pub mod ui;
pub mod view;

pub use app::App;
pub use async_tasks::{AsyncTask, TaskContext, TaskResult, TaskRunner};
pub use pymol_framework::component::{Component, SharedContext};
pub use pymol_framework::component_store::ComponentStore;
pub use components::{
    MovieComponent, ObjectListComponent, ReplComponent, SequenceComponent,
};
pub use fetch::{FetchResult, FetchTask};
pub use ipc::{
    ExternalCommandRegistry, IpcCallbackResult, IpcCallbackTask, IpcRequest, IpcResponse,
    IpcServer, OutputKind as IpcOutputKind, OutputMessage as IpcOutputMessage,
};
pub use layout::Layout;
pub use plugin_manager::PluginManager;

// Re-export model types
pub use model::{CommandLineModel, OutputModel, OutputKind, OutputMessage};

// Re-export UI state types
pub use ui::{CommandLineUiState, ObjectListUiState, SequenceUiState};

// Re-export view types
pub use view::AppView;
