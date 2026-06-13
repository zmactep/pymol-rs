//! Error types for the color system

use thiserror::Error;

/// Errors that can occur when working with colors
#[derive(Error, Debug)]
pub enum ColorError {
    /// Color not found by name or index
    #[error("Color not found: {0}")]
    NotFound(String),

    /// Invalid color value
    #[error("Invalid color value: {0}")]
    InvalidValue(String),

    /// Invalid ramp definition
    #[error("Invalid color ramp: {0}")]
    InvalidRamp(String),
}
