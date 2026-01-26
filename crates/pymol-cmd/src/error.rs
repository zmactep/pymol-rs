//! Error types for the command system
//!
//! This module defines error types for command parsing and execution.

use thiserror::Error;

/// Result type for command operations
pub type CmdResult<T = ()> = Result<T, CmdError>;

/// Errors that can occur during command execution
#[derive(Debug, Error)]
pub enum CmdError {
    /// Command parsing failed
    #[error("parse error: {0}")]
    Parse(#[from] ParseError),

    /// Command not found in registry
    #[error("unknown command: {0}")]
    UnknownCommand(String),

    /// Invalid argument provided
    #[error("invalid argument '{name}': {reason}")]
    InvalidArgument { name: String, reason: String },

    /// Missing required argument
    #[error("missing required argument: {0}")]
    MissingArgument(String),

    /// Too many arguments provided
    #[error("too many arguments: expected at most {expected}, got {got}")]
    TooManyArguments { expected: usize, got: usize },

    /// Object not found
    #[error("object not found: {0}")]
    ObjectNotFound(String),

    /// Selection error
    #[error("selection error: {0}")]
    Selection(String),

    /// File I/O error
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    /// File format error
    #[error("file format error: {0}")]
    FileFormat(String),

    /// Setting error
    #[error("setting error: {0}")]
    Setting(String),

    /// Scene error (from pymol-scene)
    #[error("scene error: {0}")]
    Scene(String),

    /// Script execution error
    #[error("script error at line {line}: {message}")]
    Script { line: usize, message: String },

    /// Command execution was aborted
    #[error("command aborted")]
    Aborted,

    /// Generic execution error
    #[error("{0}")]
    Execution(String),
}

/// Errors that can occur during command parsing
#[derive(Debug, Error, Clone, PartialEq)]
pub enum ParseError {
    /// Unexpected end of input
    #[error("unexpected end of input")]
    UnexpectedEof,

    /// Unexpected character
    #[error("unexpected character '{0}' at position {1}")]
    UnexpectedChar(char, usize),

    /// Unterminated string
    #[error("unterminated string starting at position {0}")]
    UnterminatedString(usize),

    /// Unbalanced parentheses
    #[error("unbalanced parentheses")]
    UnbalancedParens,

    /// Invalid number format
    #[error("invalid number: {0}")]
    InvalidNumber(String),

    /// Empty command
    #[error("empty command")]
    EmptyCommand,

    /// Invalid argument name
    #[error("invalid argument name: {0}")]
    InvalidArgName(String),

    /// Duplicate argument name
    #[error("duplicate argument: {0}")]
    DuplicateArgument(String),

    /// Generic parse error with message
    #[error("parse error: {0}")]
    Generic(String),
}

impl From<nom::Err<nom::error::Error<&str>>> for ParseError {
    fn from(err: nom::Err<nom::error::Error<&str>>) -> Self {
        match err {
            nom::Err::Incomplete(_) => ParseError::UnexpectedEof,
            nom::Err::Error(e) | nom::Err::Failure(e) => {
                ParseError::Generic(format!("at '{}...'", &e.input[..e.input.len().min(20)]))
            }
        }
    }
}

impl CmdError {
    /// Create an invalid argument error
    pub fn invalid_arg(name: impl Into<String>, reason: impl Into<String>) -> Self {
        CmdError::InvalidArgument {
            name: name.into(),
            reason: reason.into(),
        }
    }

    /// Create a selection error
    pub fn selection(msg: impl Into<String>) -> Self {
        CmdError::Selection(msg.into())
    }

    /// Create an execution error
    pub fn execution(msg: impl Into<String>) -> Self {
        CmdError::Execution(msg.into())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_error_display() {
        let err = CmdError::UnknownCommand("foo".to_string());
        assert_eq!(format!("{}", err), "unknown command: foo");

        let err = CmdError::InvalidArgument {
            name: "count".to_string(),
            reason: "must be positive".to_string(),
        };
        assert_eq!(
            format!("{}", err),
            "invalid argument 'count': must be positive"
        );
    }

    #[test]
    fn test_parse_error() {
        let err = ParseError::UnterminatedString(5);
        assert_eq!(
            format!("{}", err),
            "unterminated string starting at position 5"
        );
    }
}
