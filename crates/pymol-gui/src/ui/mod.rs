//! UI Components Module
//!
//! This module contains all the egui-based UI panels and widgets for the PyMOL GUI.
//! Each panel is a view function that takes model + UI state + message bus.

pub mod command;
pub mod completion;
pub mod drag_drop_overlay;
pub mod movie;
pub mod notification;
pub mod objects;
pub mod output;
pub mod sequence;

pub use command::{CommandLinePanel, CommandLineUiState};
pub use completion::{generate_completions, CompletionResult};
pub use drag_drop_overlay::DragDropOverlay;
pub use movie::MoviePanel;
pub use notification::NotificationOverlay;
pub use objects::{ObjectListPanel, ObjectListUiState};
pub use output::OutputPanel;
pub use sequence::{SequencePanel, SequenceUiState};
