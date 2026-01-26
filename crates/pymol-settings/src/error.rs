//! Error types for the settings system

use thiserror::Error;

/// Errors that can occur when working with settings
#[derive(Error, Debug, Clone)]
pub enum SettingError {
    /// Setting not found
    #[error("Setting not found: {0}")]
    NotFound(String),

    /// Type mismatch when getting or setting a value
    #[error("Type mismatch: expected {expected}, got {actual}")]
    TypeMismatch {
        expected: &'static str,
        actual: &'static str,
    },

    /// Invalid value for the setting
    #[error("Invalid value for setting '{name}': {reason}")]
    InvalidValue { name: String, reason: String },

    /// Setting is read-only
    #[error("Setting '{0}' is read-only")]
    ReadOnly(String),

    /// Setting level mismatch (trying to apply setting at wrong level)
    #[error("Setting '{name}' cannot be applied at level '{level}'")]
    LevelMismatch { name: String, level: String },

    /// Serialization error
    #[error("Serialization error: {0}")]
    Serialization(String),

    /// Deserialization error
    #[error("Deserialization error: {0}")]
    Deserialization(String),
}

impl SettingError {
    /// Create a not found error from a setting ID
    pub fn not_found_id(id: u16) -> Self {
        SettingError::NotFound(format!("id:{}", id))
    }

    /// Create a not found error from a setting name
    pub fn not_found_name(name: impl Into<String>) -> Self {
        SettingError::NotFound(name.into())
    }
}
