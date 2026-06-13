//! Command system.
//!
//! This crate provides the command parsing, registration, and execution layer.
//!
//! # Overview
//!
//! The command system allows executing text commands, supporting:
//! - Positional and named arguments
//! - Selection expressions
//! - Script (.pml) file execution
//! - Command history
//!
//! # Example
//!
//! ```rust,ignore
//! use patinae_cmd::CommandExecutor;
//! use patinae_scene::Viewer;
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

mod args;
mod command;
pub mod commands;
mod dynamic;
mod error;
mod executor;
pub mod helpers;
mod history;
mod parser;
mod script;

// Re-export main types
pub use args::{ArgValue, ParsedCommand};
pub use command::{
    ArgHint, AsyncCommandRequest, AsyncCommandSink, Command, CommandAction, CommandContext,
    CommandRegistry, CommandRuntimeRequirements, CommandSource, DynamicSettingEntry,
    DynamicSettingRegistry, FetchFormatCode, FetchRequest, FormatHandler, MessageKind,
    OutputMessage, PluginReaderFn, PluginWriterFn, ScriptHandler, ViewerLike,
};
pub use dynamic::{DynamicCommand, DynamicCommandInvocation};
pub use error::{CmdError, CmdResult, ParseError};
pub use executor::{CommandExecutor, CommandOutput};
pub use history::CommandHistory;
pub use parser::{join_continued_lines, parse_command, parse_commands};
pub use script::ScriptEngine;

/// Prelude for convenient imports
pub mod prelude {
    pub use crate::args::{ArgValue, ParsedCommand};
    pub use crate::command::{
        ArgHint, AsyncCommandRequest, AsyncCommandSink, Command, CommandAction, CommandContext,
        CommandRegistry, CommandRuntimeRequirements, CommandSource, FetchFormatCode, FetchRequest,
        FormatHandler, MessageKind, OutputMessage, PluginReaderFn, PluginWriterFn, ScriptHandler,
        ViewerLike,
    };
    pub use crate::error::{CmdError, CmdResult};
    pub use crate::executor::CommandExecutor;
    pub use crate::parser::{parse_command, parse_commands};
}
