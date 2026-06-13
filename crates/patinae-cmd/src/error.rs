//! Error types for the command system.
//!
//! This module defines stable inspection APIs for parsing and execution errors.

use std::backtrace::Backtrace;
use std::error::Error as StdError;
use std::fmt;

/// Result type for command operations.
pub type CmdResult<T = ()> = Result<T, CmdError>;

/// Errors that can occur during command execution.
#[derive(Debug)]
pub struct CmdError {
    kind: CmdErrorKind,
    backtrace: Backtrace,
}

#[derive(Debug)]
enum CmdErrorKind {
    Parse(ParseError),
    UnknownCommand(String),
    InvalidArgument { name: String, reason: String },
    MissingArgument(String),
    TooManyArguments { expected: usize, got: usize },
    ObjectNotFound(String),
    Selection(String),
    Io(std::io::Error),
    FileFormat(String),
    Setting(String),
    Scene(String),
    Script { line: usize, message: String },
    Execution(String),
}

/// Errors that can occur during command parsing.
#[derive(Debug, Clone, PartialEq)]
pub struct ParseError {
    kind: ParseErrorKind,
}

#[derive(Debug, Clone, PartialEq)]
enum ParseErrorKind {
    UnexpectedEof,
    UnexpectedChar(char, usize),
    UnterminatedString(usize),
    UnbalancedParens,
    InvalidNumber(String),
    EmptyCommand,
    InvalidArgName(String),
    DuplicateArgument(String),
    Generic(String),
}

impl CmdError {
    fn new(kind: CmdErrorKind) -> Self {
        Self {
            kind,
            backtrace: Backtrace::capture(),
        }
    }

    /// Create a command parse error.
    pub fn parse(error: ParseError) -> Self {
        Self::new(CmdErrorKind::Parse(error))
    }

    /// Create an unknown command error.
    pub fn unknown_command(name: impl Into<String>) -> Self {
        Self::new(CmdErrorKind::UnknownCommand(name.into()))
    }

    /// Create an invalid argument error.
    pub fn invalid_arg(name: impl Into<String>, reason: impl Into<String>) -> Self {
        Self::new(CmdErrorKind::InvalidArgument {
            name: name.into(),
            reason: reason.into(),
        })
    }

    /// Create an invalid argument error.
    pub fn invalid_argument(name: impl Into<String>, reason: impl Into<String>) -> Self {
        Self::invalid_arg(name, reason)
    }

    /// Create a missing argument error.
    pub fn missing_argument(name: impl Into<String>) -> Self {
        Self::new(CmdErrorKind::MissingArgument(name.into()))
    }

    /// Create a too-many-arguments error.
    pub fn too_many_arguments(expected: usize, got: usize) -> Self {
        Self::new(CmdErrorKind::TooManyArguments { expected, got })
    }

    /// Create an object-not-found error.
    pub fn object_not_found(name: impl Into<String>) -> Self {
        Self::new(CmdErrorKind::ObjectNotFound(name.into()))
    }

    /// Create a selection error.
    pub fn selection(message: impl Into<String>) -> Self {
        Self::new(CmdErrorKind::Selection(message.into()))
    }

    /// Create a standard I/O error.
    pub fn io(error: std::io::Error) -> Self {
        Self::new(CmdErrorKind::Io(error))
    }

    /// Create a file format error.
    pub fn file_format(message: impl Into<String>) -> Self {
        Self::new(CmdErrorKind::FileFormat(message.into()))
    }

    /// Create a setting error.
    pub fn setting(message: impl Into<String>) -> Self {
        Self::new(CmdErrorKind::Setting(message.into()))
    }

    /// Create a scene error.
    pub fn scene(message: impl Into<String>) -> Self {
        Self::new(CmdErrorKind::Scene(message.into()))
    }

    /// Create a script execution error.
    pub fn script(line: usize, message: impl Into<String>) -> Self {
        Self::new(CmdErrorKind::Script {
            line,
            message: message.into(),
        })
    }

    /// Create a generic execution error.
    pub fn execution(message: impl Into<String>) -> Self {
        Self::new(CmdErrorKind::Execution(message.into()))
    }

    /// Return the captured backtrace.
    pub fn backtrace(&self) -> &Backtrace {
        &self.backtrace
    }

    /// Return the upstream source error, if available.
    pub fn source(&self) -> Option<&(dyn StdError + 'static)> {
        match &self.kind {
            CmdErrorKind::Parse(error) => Some(error),
            CmdErrorKind::Io(error) => Some(error),
            _ => None,
        }
    }

    /// Return `true` if command parsing failed.
    pub fn is_parse(&self) -> bool {
        matches!(self.kind, CmdErrorKind::Parse(_))
    }

    /// Return `true` if the command name was unknown.
    pub fn is_unknown_command(&self) -> bool {
        matches!(self.kind, CmdErrorKind::UnknownCommand(_))
    }

    /// Return `true` if an argument was invalid.
    pub fn is_invalid_argument(&self) -> bool {
        matches!(self.kind, CmdErrorKind::InvalidArgument { .. })
    }

    /// Return `true` if a required argument was missing.
    pub fn is_missing_argument(&self) -> bool {
        matches!(self.kind, CmdErrorKind::MissingArgument(_))
    }

    /// Return `true` if too many arguments were provided.
    pub fn is_too_many_arguments(&self) -> bool {
        matches!(self.kind, CmdErrorKind::TooManyArguments { .. })
    }

    /// Return `true` if an object lookup failed.
    pub fn is_object_not_found(&self) -> bool {
        matches!(self.kind, CmdErrorKind::ObjectNotFound(_))
    }

    /// Return `true` if selection parsing or evaluation failed.
    pub fn is_selection(&self) -> bool {
        matches!(self.kind, CmdErrorKind::Selection(_))
    }

    /// Return `true` if this error wraps standard I/O.
    pub fn is_io(&self) -> bool {
        matches!(self.kind, CmdErrorKind::Io(_))
    }

    /// Return `true` if file format handling failed.
    pub fn is_file_format(&self) -> bool {
        matches!(self.kind, CmdErrorKind::FileFormat(_))
    }

    /// Return `true` if setting handling failed.
    pub fn is_setting(&self) -> bool {
        matches!(self.kind, CmdErrorKind::Setting(_))
    }

    /// Return `true` if scene handling failed.
    pub fn is_scene(&self) -> bool {
        matches!(self.kind, CmdErrorKind::Scene(_))
    }

    /// Return `true` if script execution failed.
    pub fn is_script(&self) -> bool {
        matches!(self.kind, CmdErrorKind::Script { .. })
    }

    /// Return `true` if generic command execution failed.
    pub fn is_execution(&self) -> bool {
        matches!(self.kind, CmdErrorKind::Execution(_))
    }

    /// Return the argument name for argument errors.
    pub fn argument_name(&self) -> Option<&str> {
        match &self.kind {
            CmdErrorKind::InvalidArgument { name, .. } | CmdErrorKind::MissingArgument(name) => {
                Some(name)
            }
            _ => None,
        }
    }

    /// Return the textual payload for string-backed errors.
    pub fn message(&self) -> Option<&str> {
        match &self.kind {
            CmdErrorKind::UnknownCommand(message)
            | CmdErrorKind::MissingArgument(message)
            | CmdErrorKind::ObjectNotFound(message)
            | CmdErrorKind::Selection(message)
            | CmdErrorKind::FileFormat(message)
            | CmdErrorKind::Setting(message)
            | CmdErrorKind::Scene(message)
            | CmdErrorKind::Execution(message) => Some(message),
            CmdErrorKind::InvalidArgument { reason, .. }
            | CmdErrorKind::Script {
                message: reason, ..
            } => Some(reason),
            CmdErrorKind::Parse(_)
            | CmdErrorKind::Io(_)
            | CmdErrorKind::TooManyArguments { .. } => None,
        }
    }
}

impl fmt::Display for CmdError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self.kind {
            CmdErrorKind::Parse(error) => write!(f, "parse error: {error}"),
            CmdErrorKind::UnknownCommand(command) => write!(f, "unknown command: {command}"),
            CmdErrorKind::InvalidArgument { name, reason } => {
                write!(f, "invalid argument '{name}': {reason}")
            }
            CmdErrorKind::MissingArgument(argument) => {
                write!(f, "missing required argument: {argument}")
            }
            CmdErrorKind::TooManyArguments { expected, got } => {
                write!(
                    f,
                    "too many arguments: expected at most {expected}, got {got}"
                )
            }
            CmdErrorKind::ObjectNotFound(object) => write!(f, "object not found: {object}"),
            CmdErrorKind::Selection(message) => write!(f, "selection error: {message}"),
            CmdErrorKind::Io(error) => write!(f, "I/O error: {error}"),
            CmdErrorKind::FileFormat(message) => write!(f, "file format error: {message}"),
            CmdErrorKind::Setting(message) => write!(f, "setting error: {message}"),
            CmdErrorKind::Scene(message) => write!(f, "scene error: {message}"),
            CmdErrorKind::Script { line, message } => {
                write!(f, "script error at line {line}: {message}")
            }
            CmdErrorKind::Execution(message) => write!(f, "{message}"),
        }
    }
}

impl StdError for CmdError {
    fn source(&self) -> Option<&(dyn StdError + 'static)> {
        CmdError::source(self)
    }
}

impl From<ParseError> for CmdError {
    fn from(error: ParseError) -> Self {
        Self::parse(error)
    }
}

impl From<std::io::Error> for CmdError {
    fn from(error: std::io::Error) -> Self {
        Self::io(error)
    }
}

impl ParseError {
    fn new(kind: ParseErrorKind) -> Self {
        Self { kind }
    }

    /// Create an unexpected end-of-input error.
    pub fn unexpected_eof() -> Self {
        Self::new(ParseErrorKind::UnexpectedEof)
    }

    /// Create an unexpected character error.
    pub fn unexpected_char(ch: char, position: usize) -> Self {
        Self::new(ParseErrorKind::UnexpectedChar(ch, position))
    }

    /// Create an unterminated string error.
    pub fn unterminated_string(position: usize) -> Self {
        Self::new(ParseErrorKind::UnterminatedString(position))
    }

    /// Create an unbalanced parentheses error.
    pub fn unbalanced_parens() -> Self {
        Self::new(ParseErrorKind::UnbalancedParens)
    }

    /// Create an invalid number error.
    pub fn invalid_number(message: impl Into<String>) -> Self {
        Self::new(ParseErrorKind::InvalidNumber(message.into()))
    }

    /// Create an empty command error.
    pub fn empty_command() -> Self {
        Self::new(ParseErrorKind::EmptyCommand)
    }

    /// Create an invalid argument name error.
    pub fn invalid_arg_name(name: impl Into<String>) -> Self {
        Self::new(ParseErrorKind::InvalidArgName(name.into()))
    }

    /// Create a duplicate argument error.
    pub fn duplicate_argument(name: impl Into<String>) -> Self {
        Self::new(ParseErrorKind::DuplicateArgument(name.into()))
    }

    /// Create a generic parse error.
    pub fn generic(message: impl Into<String>) -> Self {
        Self::new(ParseErrorKind::Generic(message.into()))
    }

    /// Return `true` if input ended unexpectedly.
    pub fn is_unexpected_eof(&self) -> bool {
        matches!(self.kind, ParseErrorKind::UnexpectedEof)
    }

    /// Return `true` if an unexpected character was found.
    pub fn is_unexpected_char(&self) -> bool {
        matches!(self.kind, ParseErrorKind::UnexpectedChar(_, _))
    }

    /// Return `true` if a string was unterminated.
    pub fn is_unterminated_string(&self) -> bool {
        matches!(self.kind, ParseErrorKind::UnterminatedString(_))
    }

    /// Return `true` if parentheses were unbalanced.
    pub fn is_unbalanced_parens(&self) -> bool {
        matches!(self.kind, ParseErrorKind::UnbalancedParens)
    }

    /// Return `true` if a number could not be parsed.
    pub fn is_invalid_number(&self) -> bool {
        matches!(self.kind, ParseErrorKind::InvalidNumber(_))
    }

    /// Return `true` if the command was empty.
    pub fn is_empty_command(&self) -> bool {
        matches!(self.kind, ParseErrorKind::EmptyCommand)
    }

    /// Return `true` if an argument name was invalid.
    pub fn is_invalid_arg_name(&self) -> bool {
        matches!(self.kind, ParseErrorKind::InvalidArgName(_))
    }

    /// Return `true` if an argument name was duplicated.
    pub fn is_duplicate_argument(&self) -> bool {
        matches!(self.kind, ParseErrorKind::DuplicateArgument(_))
    }

    /// Return `true` if this is a generic parse error.
    pub fn is_generic(&self) -> bool {
        matches!(self.kind, ParseErrorKind::Generic(_))
    }

    /// Return the position associated with this parse error.
    pub fn position(&self) -> Option<usize> {
        match self.kind {
            ParseErrorKind::UnexpectedChar(_, position)
            | ParseErrorKind::UnterminatedString(position) => Some(position),
            _ => None,
        }
    }

    /// Return the unexpected character if available.
    pub fn character(&self) -> Option<char> {
        match self.kind {
            ParseErrorKind::UnexpectedChar(ch, _) => Some(ch),
            _ => None,
        }
    }

    /// Return the textual payload for string-backed parse errors.
    pub fn message(&self) -> Option<&str> {
        match &self.kind {
            ParseErrorKind::InvalidNumber(message)
            | ParseErrorKind::InvalidArgName(message)
            | ParseErrorKind::DuplicateArgument(message)
            | ParseErrorKind::Generic(message) => Some(message),
            ParseErrorKind::UnexpectedEof
            | ParseErrorKind::UnexpectedChar(_, _)
            | ParseErrorKind::UnterminatedString(_)
            | ParseErrorKind::UnbalancedParens
            | ParseErrorKind::EmptyCommand => None,
        }
    }

    /// Return the upstream source error, if available.
    pub fn source(&self) -> Option<&(dyn StdError + 'static)> {
        None
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self.kind {
            ParseErrorKind::UnexpectedEof => write!(f, "unexpected end of input"),
            ParseErrorKind::UnexpectedChar(ch, position) => {
                write!(f, "unexpected character '{ch}' at position {position}")
            }
            ParseErrorKind::UnterminatedString(position) => {
                write!(f, "unterminated string starting at position {position}")
            }
            ParseErrorKind::UnbalancedParens => write!(f, "unbalanced parentheses"),
            ParseErrorKind::InvalidNumber(message) => write!(f, "invalid number: {message}"),
            ParseErrorKind::EmptyCommand => write!(f, "empty command"),
            ParseErrorKind::InvalidArgName(name) => write!(f, "invalid argument name: {name}"),
            ParseErrorKind::DuplicateArgument(name) => write!(f, "duplicate argument: {name}"),
            ParseErrorKind::Generic(message) => write!(f, "parse error: {message}"),
        }
    }
}

impl StdError for ParseError {
    fn source(&self) -> Option<&(dyn StdError + 'static)> {
        ParseError::source(self)
    }
}

impl From<nom::Err<nom::error::Error<&str>>> for ParseError {
    fn from(err: nom::Err<nom::error::Error<&str>>) -> Self {
        match err {
            nom::Err::Incomplete(_) => ParseError::unexpected_eof(),
            nom::Err::Error(e) | nom::Err::Failure(e) => {
                ParseError::generic(format!("at '{}...'", &e.input[..e.input.len().min(20)]))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_error_display() {
        let err = CmdError::unknown_command("foo");
        assert_eq!(format!("{}", err), "unknown command: foo");

        let err = CmdError::invalid_arg("count", "must be positive");
        assert_eq!(
            format!("{}", err),
            "invalid argument 'count': must be positive"
        );
    }

    #[test]
    fn test_parse_error() {
        let err = ParseError::unterminated_string(5);
        assert_eq!(
            format!("{}", err),
            "unterminated string starting at position 5"
        );
        assert!(err.is_unterminated_string());
        assert_eq!(err.position(), Some(5));
    }

    #[test]
    fn test_cmd_error_preserves_parse_source() {
        let err = CmdError::from(ParseError::empty_command());

        assert_eq!(err.to_string(), "parse error: empty command");
        assert!(err.is_parse());
        assert!(StdError::source(&err).is_some());
    }

    #[test]
    fn test_cmd_error_preserves_io_source() {
        let err = CmdError::from(std::io::Error::new(
            std::io::ErrorKind::PermissionDenied,
            "script denied",
        ));

        assert_eq!(err.to_string(), "I/O error: script denied");
        assert!(err.is_io());
        assert!(StdError::source(&err).is_some());
    }
}
