//! Error types for molecular file I/O
//!
//! Provides error types for parsing and writing molecular file formats.

use thiserror::Error;

/// Errors that can occur during molecular file I/O
#[derive(Error, Debug)]
pub enum IoError {
    /// Standard I/O error
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    /// Parse error with location information
    #[error("Parse error at line {line}: {message}")]
    Parse {
        /// Line number where the error occurred (1-based)
        line: usize,
        /// Error message
        message: String,
    },

    /// Unknown or unsupported file format
    #[error("Unknown format: {0}")]
    UnknownFormat(String),

    /// Invalid record in the file
    #[error("Invalid record: {0}")]
    InvalidRecord(String),

    /// Missing required field
    #[error("Missing required field: {0}")]
    MissingField(String),

    /// Unsupported feature in the file format
    #[error("Unsupported feature: {0}")]
    Unsupported(String),

    /// Invalid element symbol
    #[error("Invalid element symbol: {0}")]
    InvalidElement(String),

    /// Invalid bond specification
    #[error("Invalid bond: {0}")]
    InvalidBond(String),

    /// Coordinate parsing error
    #[error("Invalid coordinate: {0}")]
    InvalidCoordinate(String),

    /// File is empty or contains no molecules
    #[error("Empty file or no molecules found")]
    EmptyFile,

    /// Decompression error
    #[error("Decompression error: {0}")]
    Decompression(String),

    /// Network/HTTP error during fetch
    #[cfg(any(feature = "fetch", feature = "fetch-async"))]
    #[error("Fetch error: {0}")]
    Fetch(String),

    /// Invalid PDB ID format
    #[cfg(any(feature = "fetch", feature = "fetch-async"))]
    #[error("Invalid PDB ID: {0}")]
    InvalidPdbId(String),
}

impl IoError {
    /// Create a parse error at a specific line
    pub fn parse(line: usize, message: impl Into<String>) -> Self {
        IoError::Parse {
            line,
            message: message.into(),
        }
    }

    /// Create a parse error without line information (line 0)
    pub fn parse_msg(message: impl Into<String>) -> Self {
        IoError::Parse {
            line: 0,
            message: message.into(),
        }
    }

    /// Create an invalid record error
    pub fn invalid_record(record: impl Into<String>) -> Self {
        IoError::InvalidRecord(record.into())
    }

    /// Create a missing field error
    pub fn missing_field(field: impl Into<String>) -> Self {
        IoError::MissingField(field.into())
    }

    /// Create an unsupported feature error
    pub fn unsupported(feature: impl Into<String>) -> Self {
        IoError::Unsupported(feature.into())
    }

    /// Create a fetch error
    #[cfg(any(feature = "fetch", feature = "fetch-async"))]
    pub fn fetch(message: impl Into<String>) -> Self {
        IoError::Fetch(message.into())
    }

    /// Create an invalid PDB ID error
    #[cfg(any(feature = "fetch", feature = "fetch-async"))]
    pub fn invalid_pdb_id(pdb_id: impl Into<String>) -> Self {
        IoError::InvalidPdbId(pdb_id.into())
    }
}

/// Result type for molecular file I/O operations
pub type IoResult<T> = Result<T, IoError>;
