//! Error types for selection parsing and evaluation
//!
//! Provides error types for all failure modes in the selection system.

use thiserror::Error;

/// Errors that can occur during selection parsing
#[derive(Debug, Error, Clone, PartialEq)]
pub enum ParseError {
    /// Unexpected end of input
    #[error("unexpected end of input")]
    UnexpectedEof,

    /// Unexpected token encountered
    #[error("unexpected token: {0}")]
    UnexpectedToken(String),

    /// Expected a specific token
    #[error("expected {expected}, found {found}")]
    Expected { expected: String, found: String },

    /// Invalid number format
    #[error("invalid number: {0}")]
    InvalidNumber(String),

    /// Invalid residue specification
    #[error("invalid residue specification: {0}")]
    InvalidResidue(String),

    /// Unknown keyword
    #[error("unknown keyword: {0}")]
    UnknownKeyword(String),

    /// Missing required argument
    #[error("missing argument for {0}")]
    MissingArgument(String),

    /// Invalid pattern
    #[error("invalid pattern: {0}")]
    InvalidPattern(String),

    /// Unmatched parenthesis
    #[error("unmatched parenthesis")]
    UnmatchedParen,

    /// Invalid macro syntax
    #[error("invalid macro syntax: {0}")]
    InvalidMacro(String),

    /// General parse error with context
    #[error("parse error at position {position}: {message}")]
    General { position: usize, message: String },
}

/// Errors that can occur during selection evaluation
#[derive(Debug, Error, Clone, PartialEq)]
pub enum EvalError {
    /// Named selection not found
    #[error("selection not found: {0}")]
    SelectionNotFound(String),

    /// Model/object not found
    #[error("model not found: {0}")]
    ModelNotFound(String),

    /// No molecules in context
    #[error("no molecules in context")]
    NoMolecules,

    /// State index out of bounds
    #[error("state {0} out of bounds")]
    StateOutOfBounds(usize),

    /// Invalid property name
    #[error("invalid property: {0}")]
    InvalidProperty(String),

    /// Type mismatch in comparison
    #[error("type mismatch: expected {expected}, got {got}")]
    TypeMismatch { expected: String, got: String },

    /// No coordinates available for distance calculation
    #[error("no coordinates available for state {0}")]
    NoCoordinates(usize),

    /// Internal error
    #[error("internal error: {0}")]
    Internal(String),
}

/// Combined error type for selection operations
#[derive(Debug, Error, Clone, PartialEq)]
pub enum SelectError {
    /// Parse error
    #[error("parse error: {0}")]
    Parse(#[from] ParseError),

    /// Evaluation error
    #[error("evaluation error: {0}")]
    Eval(#[from] EvalError),
}

/// Result type for parsing operations
#[allow(dead_code)]
pub type ParseResult<T> = Result<T, ParseError>;

/// Result type for evaluation operations
pub type EvalResult<T> = Result<T, EvalError>;

/// Result type for general selection operations
pub type SelectResult<T> = Result<T, SelectError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_error_display() {
        let err = ParseError::UnknownKeyword("foobar".to_string());
        assert_eq!(format!("{}", err), "unknown keyword: foobar");
    }

    #[test]
    fn test_eval_error_display() {
        let err = EvalError::SelectionNotFound("sele1".to_string());
        assert_eq!(format!("{}", err), "selection not found: sele1");
    }

    #[test]
    fn test_select_error_from() {
        let parse_err = ParseError::UnexpectedEof;
        let select_err: SelectError = parse_err.into();
        assert!(matches!(select_err, SelectError::Parse(_)));

        let eval_err = EvalError::NoMolecules;
        let select_err: SelectError = eval_err.into();
        assert!(matches!(select_err, SelectError::Eval(_)));
    }
}
