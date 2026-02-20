//! UI Components Module
//!
//! This module contains all the egui-based UI panels and widgets for the PyMOL GUI.

pub mod command;
pub mod completion;
pub mod drag_drop_overlay;
pub mod movie;
pub mod notification;
pub mod objects;
pub mod output;
pub mod sequence;

pub use command::CommandLinePanel;
pub use completion::{generate_completions, CompletionResult};
pub use drag_drop_overlay::DragDropOverlay;
pub use movie::{MovieAction, MoviePanel};
pub use notification::NotificationOverlay;
pub use objects::ObjectListPanel;
pub use output::OutputPanel;
pub use sequence::SequencePanel;
