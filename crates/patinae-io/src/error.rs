//! Error types for molecular file I/O.
//!
//! Provides stable inspection APIs for parsing and writing molecular formats.

use std::backtrace::Backtrace;
use std::error::Error as StdError;
use std::fmt;

/// Errors that can occur during molecular file I/O.
#[derive(Debug)]
pub struct IoError {
    kind: IoErrorKind,
    backtrace: Backtrace,
}

#[derive(Debug)]
enum IoErrorKind {
    Io(std::io::Error),
    Parse {
        line: usize,
        message: String,
    },
    UnknownFormat(String),
    InvalidRecord(String),
    MissingField(String),
    Unsupported(String),
    InvalidElement(String),
    InvalidBond(String),
    InvalidCoordinate(String),
    EmptyFile,
    Decompression(String),
    #[cfg(any(feature = "fetch", feature = "fetch-async"))]
    Fetch(String),
    #[cfg(any(feature = "fetch", feature = "fetch-async"))]
    InvalidPdbId(String),
}

impl IoError {
    fn new(kind: IoErrorKind) -> Self {
        Self {
            kind,
            backtrace: Backtrace::capture(),
        }
    }

    /// Create a standard I/O error.
    pub fn io(error: std::io::Error) -> Self {
        Self::new(IoErrorKind::Io(error))
    }

    /// Create a parse error at a specific line.
    pub fn parse(line: usize, message: impl Into<String>) -> Self {
        Self::new(IoErrorKind::Parse {
            line,
            message: message.into(),
        })
    }

    /// Create a parse error without line information.
    pub fn parse_msg(message: impl Into<String>) -> Self {
        Self::parse(0, message)
    }

    /// Create an unknown format error.
    pub fn unknown_format(format: impl Into<String>) -> Self {
        Self::new(IoErrorKind::UnknownFormat(format.into()))
    }

    /// Create an invalid record error.
    pub fn invalid_record(record: impl Into<String>) -> Self {
        Self::new(IoErrorKind::InvalidRecord(record.into()))
    }

    /// Create a missing field error.
    pub fn missing_field(field: impl Into<String>) -> Self {
        Self::new(IoErrorKind::MissingField(field.into()))
    }

    /// Create an unsupported feature error.
    pub fn unsupported(feature: impl Into<String>) -> Self {
        Self::new(IoErrorKind::Unsupported(feature.into()))
    }

    /// Create an invalid element symbol error.
    pub fn invalid_element(element: impl Into<String>) -> Self {
        Self::new(IoErrorKind::InvalidElement(element.into()))
    }

    /// Create an invalid bond specification error.
    pub fn invalid_bond(bond: impl Into<String>) -> Self {
        Self::new(IoErrorKind::InvalidBond(bond.into()))
    }

    /// Create an invalid coordinate error.
    pub fn invalid_coordinate(coordinate: impl Into<String>) -> Self {
        Self::new(IoErrorKind::InvalidCoordinate(coordinate.into()))
    }

    /// Create an empty file error.
    pub fn empty_file() -> Self {
        Self::new(IoErrorKind::EmptyFile)
    }

    /// Create a decompression error.
    pub fn decompression(message: impl Into<String>) -> Self {
        Self::new(IoErrorKind::Decompression(message.into()))
    }

    /// Create a fetch error.
    #[cfg(any(feature = "fetch", feature = "fetch-async"))]
    pub fn fetch(message: impl Into<String>) -> Self {
        Self::new(IoErrorKind::Fetch(message.into()))
    }

    /// Create an invalid PDB ID error.
    #[cfg(any(feature = "fetch", feature = "fetch-async"))]
    pub fn invalid_pdb_id(pdb_id: impl Into<String>) -> Self {
        Self::new(IoErrorKind::InvalidPdbId(pdb_id.into()))
    }

    /// Return the captured backtrace.
    pub fn backtrace(&self) -> &Backtrace {
        &self.backtrace
    }

    /// Return the upstream source error, if available.
    pub fn source(&self) -> Option<&(dyn StdError + 'static)> {
        match &self.kind {
            IoErrorKind::Io(error) => Some(error),
            _ => None,
        }
    }

    /// Return `true` if this error wraps a standard I/O error.
    pub fn is_io(&self) -> bool {
        matches!(self.kind, IoErrorKind::Io(_))
    }

    /// Return `true` if this error describes a parse failure.
    pub fn is_parse(&self) -> bool {
        matches!(self.kind, IoErrorKind::Parse { .. })
    }

    /// Return `true` if this error describes an unknown format.
    pub fn is_unknown_format(&self) -> bool {
        matches!(self.kind, IoErrorKind::UnknownFormat(_))
    }

    /// Return `true` if this error describes an invalid record.
    pub fn is_invalid_record(&self) -> bool {
        matches!(self.kind, IoErrorKind::InvalidRecord(_))
    }

    /// Return `true` if this error describes a missing field.
    pub fn is_missing_field(&self) -> bool {
        matches!(self.kind, IoErrorKind::MissingField(_))
    }

    /// Return `true` if this error describes unsupported functionality.
    pub fn is_unsupported(&self) -> bool {
        matches!(self.kind, IoErrorKind::Unsupported(_))
    }

    /// Return `true` if this error describes an invalid element symbol.
    pub fn is_invalid_element(&self) -> bool {
        matches!(self.kind, IoErrorKind::InvalidElement(_))
    }

    /// Return `true` if this error describes an invalid bond.
    pub fn is_invalid_bond(&self) -> bool {
        matches!(self.kind, IoErrorKind::InvalidBond(_))
    }

    /// Return `true` if this error describes an invalid coordinate.
    pub fn is_invalid_coordinate(&self) -> bool {
        matches!(self.kind, IoErrorKind::InvalidCoordinate(_))
    }

    /// Return `true` if this error describes an empty input.
    pub fn is_empty_file(&self) -> bool {
        matches!(self.kind, IoErrorKind::EmptyFile)
    }

    /// Return `true` if this error describes decompression failure.
    pub fn is_decompression(&self) -> bool {
        matches!(self.kind, IoErrorKind::Decompression(_))
    }

    /// Return `true` if this error describes fetch failure.
    #[cfg(any(feature = "fetch", feature = "fetch-async"))]
    pub fn is_fetch(&self) -> bool {
        matches!(self.kind, IoErrorKind::Fetch(_))
    }

    /// Return `true` if this error describes an invalid PDB ID.
    #[cfg(any(feature = "fetch", feature = "fetch-async"))]
    pub fn is_invalid_pdb_id(&self) -> bool {
        matches!(self.kind, IoErrorKind::InvalidPdbId(_))
    }

    /// Return the parse line for parse errors.
    pub fn parse_line(&self) -> Option<usize> {
        match self.kind {
            IoErrorKind::Parse { line, .. } => Some(line),
            _ => None,
        }
    }

    /// Return the textual payload for string-backed errors.
    pub fn message(&self) -> Option<&str> {
        match &self.kind {
            IoErrorKind::Parse { message, .. }
            | IoErrorKind::UnknownFormat(message)
            | IoErrorKind::InvalidRecord(message)
            | IoErrorKind::MissingField(message)
            | IoErrorKind::Unsupported(message)
            | IoErrorKind::InvalidElement(message)
            | IoErrorKind::InvalidBond(message)
            | IoErrorKind::InvalidCoordinate(message)
            | IoErrorKind::Decompression(message) => Some(message),
            #[cfg(any(feature = "fetch", feature = "fetch-async"))]
            IoErrorKind::Fetch(message) | IoErrorKind::InvalidPdbId(message) => Some(message),
            IoErrorKind::Io(_) | IoErrorKind::EmptyFile => None,
        }
    }
}

impl fmt::Display for IoError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self.kind {
            IoErrorKind::Io(error) => write!(f, "I/O error: {error}"),
            IoErrorKind::Parse { line, message } => {
                write!(f, "Parse error at line {line}: {message}")
            }
            IoErrorKind::UnknownFormat(format) => write!(f, "Unknown format: {format}"),
            IoErrorKind::InvalidRecord(record) => write!(f, "Invalid record: {record}"),
            IoErrorKind::MissingField(field) => write!(f, "Missing required field: {field}"),
            IoErrorKind::Unsupported(feature) => write!(f, "Unsupported feature: {feature}"),
            IoErrorKind::InvalidElement(element) => {
                write!(f, "Invalid element symbol: {element}")
            }
            IoErrorKind::InvalidBond(bond) => write!(f, "Invalid bond: {bond}"),
            IoErrorKind::InvalidCoordinate(coordinate) => {
                write!(f, "Invalid coordinate: {coordinate}")
            }
            IoErrorKind::EmptyFile => write!(f, "Empty file or no molecules found"),
            IoErrorKind::Decompression(message) => write!(f, "Decompression error: {message}"),
            #[cfg(any(feature = "fetch", feature = "fetch-async"))]
            IoErrorKind::Fetch(message) => write!(f, "Fetch error: {message}"),
            #[cfg(any(feature = "fetch", feature = "fetch-async"))]
            IoErrorKind::InvalidPdbId(pdb_id) => write!(f, "Invalid PDB ID: {pdb_id}"),
        }
    }
}

impl StdError for IoError {
    fn source(&self) -> Option<&(dyn StdError + 'static)> {
        IoError::source(self)
    }
}

impl From<std::io::Error> for IoError {
    fn from(error: std::io::Error) -> Self {
        Self::io(error)
    }
}

/// Result type for molecular file I/O operations.
pub type IoResult<T> = Result<T, IoError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn display_and_helpers_match_legacy_messages() {
        let err = IoError::parse(42, "bad atom serial");
        assert_eq!(err.to_string(), "Parse error at line 42: bad atom serial");
        assert!(err.is_parse());
        assert_eq!(err.parse_line(), Some(42));
        assert_eq!(err.message(), Some("bad atom serial"));

        let err = IoError::empty_file();
        assert_eq!(err.to_string(), "Empty file or no molecules found");
        assert!(err.is_empty_file());
    }

    #[test]
    fn io_error_preserves_source_chain() {
        let err = IoError::from(std::io::Error::new(
            std::io::ErrorKind::PermissionDenied,
            "cannot read molecule",
        ));

        assert_eq!(err.to_string(), "I/O error: cannot read molecule");
        assert!(err.is_io());
        assert!(StdError::source(&err).is_some());
    }
}
