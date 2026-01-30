//! State Management
//!
//! This module contains all application state types, organized into focused components:
//! - `AppState` - Scene/domain data (registry, camera, selections, settings, colors, executor)
//! - `OutputBufferState` - Log/output message buffer
//! - `CommandLineState` - Command line input, history, and autocomplete
//! - `CompletionState` - Autocomplete popup state

mod app_state;
mod command_line;
mod completion;
mod output;

pub use app_state::AppState;
pub use command_line::CommandLineState;
pub use completion::CompletionState;
pub use output::{OutputBufferState, OutputKind, OutputMessage};
