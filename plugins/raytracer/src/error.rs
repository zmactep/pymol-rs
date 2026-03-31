//! Error types for the raytracer

use thiserror::Error;

/// Raytracing errors
#[derive(Error, Debug)]
pub enum RaytraceError {
    /// GPU initialization failed
    #[error("GPU initialization failed: {0}")]
    GpuInit(String),

    /// Shader compilation failed
    #[error("Shader compilation failed: {0}")]
    ShaderCompilation(String),

    /// Buffer creation failed
    #[error("Buffer creation failed: {0}")]
    BufferCreation(String),

    /// Pipeline creation failed
    #[error("Pipeline creation failed: {0}")]
    PipelineCreation(String),

    /// Invalid parameters
    #[error("Invalid parameters: {0}")]
    InvalidParams(String),

    /// BVH construction failed
    #[error("BVH construction failed: {0}")]
    BvhConstruction(String),

    /// Render failed
    #[error("Render failed: {0}")]
    RenderFailed(String),

    /// GPU timeout
    #[error("GPU operation timed out")]
    Timeout,

    /// No primitives to render
    #[error("No primitives to render")]
    NoPrimitives,
}

/// Result type for raytracing operations
pub type RaytraceResult<T> = Result<T, RaytraceError>;
