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
pub mod state;
pub mod ui;

pub use app::App;
pub use async_tasks::{AsyncTask, AsyncTaskResult, AsyncTaskRunner, FetchTask};
pub use state::GuiState;
