//! Error types for the scene crate

use thiserror::Error;

/// Scene-related errors
#[derive(Debug, Error)]
pub enum SceneError {
    /// Object not found in registry
    #[error("Object not found: {0}")]
    ObjectNotFound(String),

    /// Object already exists with this name
    #[error("Object already exists: {0}")]
    ObjectExists(String),

    /// Invalid object type for operation
    #[error("Invalid object type: expected {expected}, got {actual}")]
    InvalidObjectType { expected: String, actual: String },

    /// Scene not found
    #[error("Scene not found: {0}")]
    SceneNotFound(String),

    /// Invalid state index
    #[error("Invalid state index: {index} (object has {count} states)")]
    InvalidState { index: usize, count: usize },

    /// Camera animation error
    #[error("Camera animation error: {0}")]
    AnimationError(String),

    /// Window error
    #[error("Window error: {0}")]
    WindowError(String),

    /// Render error
    #[error("Render error: {0}")]
    RenderError(String),

    /// wgpu surface error
    #[error("Surface error: {0}")]
    SurfaceError(#[from] wgpu::SurfaceError),

    /// wgpu request device error
    #[error("GPU device error: {0}")]
    DeviceError(#[from] wgpu::RequestDeviceError),

    /// IO error
    #[error("IO error: {0}")]
    IoError(#[from] std::io::Error),
}

/// Result type for scene operations
pub type SceneResult<T> = Result<T, SceneError>;

/// Window-specific errors
#[derive(Debug, Error)]
pub enum WindowError {
    /// Failed to create window
    #[error("Failed to create window: {0}")]
    CreationFailed(String),

    /// Failed to create surface
    #[error("Failed to create surface: {0}")]
    SurfaceFailed(String),

    /// Surface lost
    #[error("Surface lost")]
    SurfaceLost,

    /// Window not available
    #[error("Window not available")]
    NotAvailable,

    /// winit OS error
    #[error("OS error: {0}")]
    OsError(#[from] winit::error::OsError),
}

/// Viewer-level errors
#[derive(Debug, Error)]
pub enum ViewerError {
    /// Scene error
    #[error(transparent)]
    Scene(#[from] SceneError),

    /// Window error
    #[error(transparent)]
    Window(#[from] WindowError),

    /// Failed to initialize GPU
    #[error("Failed to initialize GPU: {0}")]
    GpuInitFailed(String),

    /// File loading error
    #[error("Failed to load file: {0}")]
    LoadError(String),

    /// Render error
    #[error("Render error: {0}")]
    RenderError(#[from] pymol_render::RenderError),

    /// IO error (for pymol-io integration)
    #[error("IO error: {0}")]
    IoError(String),
}
