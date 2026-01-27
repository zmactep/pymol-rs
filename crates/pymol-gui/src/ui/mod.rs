//! UI Components Module
//!
//! This module contains all the egui-based UI panels and widgets for the PyMOL GUI.

pub mod command;
pub mod controls;
pub mod mouse_info;
pub mod objects;
pub mod output;
pub mod state_bar;

pub use command::CommandLinePanel;
pub use controls::ControlsPanel;
pub use mouse_info::MouseInfoPanel;
pub use objects::ObjectListPanel;
pub use output::OutputPanel;
pub use state_bar::StateBar;
