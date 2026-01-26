//! Error types for the render crate

use thiserror::Error;

/// Errors that can occur during rendering operations
#[derive(Error, Debug)]
pub enum RenderError {
    /// Failed to create wgpu adapter
    #[error("Failed to request wgpu adapter")]
    AdapterRequestFailed,

    /// Failed to create wgpu device
    #[error("Failed to request wgpu device: {0}")]
    DeviceRequestFailed(#[from] wgpu::RequestDeviceError),

    /// Shader compilation failed
    #[error("Shader compilation failed: {0}")]
    ShaderCompilationFailed(String),

    /// Pipeline creation failed
    #[error("Pipeline creation failed: {0}")]
    PipelineCreationFailed(String),

    /// Buffer creation failed
    #[error("Buffer creation failed: {0}")]
    BufferCreationFailed(String),

    /// Invalid vertex data
    #[error("Invalid vertex data: {0}")]
    InvalidVertexData(String),

    /// Surface error
    #[error("Surface error: {0}")]
    SurfaceError(#[from] wgpu::SurfaceError),
}
