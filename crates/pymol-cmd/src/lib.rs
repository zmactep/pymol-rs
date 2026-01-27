//! PyMOL Command System
//!
//! This crate provides the command parsing, registration, and execution layer for PyMOL-RS.
//!
//! # Overview
//!
//! The command system allows executing PyMOL commands from text strings, supporting:
//! - Positional and named arguments
//! - PyMOL selection syntax
//! - Script (.pml) file execution
//! - Command logging and history
//!
//! # Example
//!
//! ```rust,ignore
//! use pymol_cmd::CommandExecutor;
//! use pymol_scene::Viewer;
//!
//! let mut viewer = Viewer::new();
//! let mut executor = CommandExecutor::new();
//!
//! // Execute commands
//! executor.do_(&mut viewer, "load protein.pdb")?;
//! executor.do_(&mut viewer, "show cartoon")?;
//! executor.do_(&mut viewer, "color green, chain A")?;
//! executor.do_(&mut viewer, "zoom")?;
//! ```
//!
//! # Architecture
//!
//! The command system consists of several components:
//!
//! - **Parser**: Parses command strings into structured `ParsedCommand` objects
//! - **Command trait**: Interface for implementing commands
//! - **CommandRegistry**: Maps command names to implementations
//! - **CommandExecutor**: Dispatches and executes commands
//! - **ScriptEngine**: Executes .pml script files
//! - **CommandLogger**: Records commands for replay

mod args;
mod command;
pub mod commands;
mod error;
mod executor;
mod history;
mod logger;
mod parser;
mod script;

// Re-export main types
pub use args::{ArgDef, ArgType, ArgValue, ParsedCommand};
pub use command::{Command, CommandContext, CommandRegistry, ViewerLike};
pub use error::{CmdError, CmdResult, ParseError};
pub use executor::CommandExecutor;
pub use history::CommandHistory;
pub use logger::{CommandLogger, LogFormat};
pub use parser::{parse_command, parse_commands};
pub use script::ScriptEngine;

/// Prelude for convenient imports
pub mod prelude {
    pub use crate::args::{ArgDef, ArgValue, ParsedCommand};
    pub use crate::command::{Command, CommandContext, CommandRegistry, ViewerLike};
    pub use crate::error::{CmdError, CmdResult};
    pub use crate::executor::CommandExecutor;
    pub use crate::parser::{parse_command, parse_commands};
}
