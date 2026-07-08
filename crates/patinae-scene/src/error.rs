//! Error types for the scene crate.
//!
//! This module defines stable inspection APIs for scene and viewer errors.

use std::backtrace::Backtrace;
use std::error::Error as StdError;
use std::fmt;

/// Scene-related errors.
#[derive(Debug)]
pub struct SceneError {
    kind: SceneErrorKind,
    backtrace: Backtrace,
}

#[derive(Debug)]
enum SceneErrorKind {
    ObjectNotFound(String),
    ObjectExists(String),
    InvalidObjectType { expected: String, actual: String },
    SceneNotFound(String),
    ViewNotFound(String),
    ViewExists(String),
    InvalidState { index: usize, count: usize },
    AnimationError(String),
    WindowError(String),
    RenderError(String),
    SurfaceError(String),
    DeviceError(wgpu::RequestDeviceError),
    IoError(std::io::Error),
}

/// Result type for scene operations.
pub type SceneResult<T> = Result<T, SceneError>;

/// Viewer-level errors.
#[derive(Debug)]
pub struct ViewerError {
    kind: ViewerErrorKind,
    backtrace: Backtrace,
}

#[derive(Debug)]
enum ViewerErrorKind {
    Scene(Box<SceneError>),
    GpuInitFailed(String),
    LoadError(String),
    IoError(String),
    CaptureError(String),
}

impl SceneError {
    fn new(kind: SceneErrorKind) -> Self {
        Self {
            kind,
            backtrace: Backtrace::capture(),
        }
    }

    /// Create an object-not-found error.
    pub fn object_not_found(name: impl Into<String>) -> Self {
        Self::new(SceneErrorKind::ObjectNotFound(name.into()))
    }

    /// Create an object-already-exists error.
    pub fn object_exists(name: impl Into<String>) -> Self {
        Self::new(SceneErrorKind::ObjectExists(name.into()))
    }

    /// Create an invalid-object-type error.
    pub fn invalid_object_type(expected: impl Into<String>, actual: impl Into<String>) -> Self {
        Self::new(SceneErrorKind::InvalidObjectType {
            expected: expected.into(),
            actual: actual.into(),
        })
    }

    /// Create a scene-not-found error.
    pub fn scene_not_found(name: impl Into<String>) -> Self {
        Self::new(SceneErrorKind::SceneNotFound(name.into()))
    }

    /// Create a view-not-found error.
    pub fn view_not_found(name: impl Into<String>) -> Self {
        Self::new(SceneErrorKind::ViewNotFound(name.into()))
    }

    /// Create a view-already-exists error.
    pub fn view_exists(name: impl Into<String>) -> Self {
        Self::new(SceneErrorKind::ViewExists(name.into()))
    }

    /// Create an invalid-state-index error.
    pub fn invalid_state(index: usize, count: usize) -> Self {
        Self::new(SceneErrorKind::InvalidState { index, count })
    }

    /// Create a camera-animation error.
    pub fn animation_error(message: impl Into<String>) -> Self {
        Self::new(SceneErrorKind::AnimationError(message.into()))
    }

    /// Create a window error.
    pub fn window_error(message: impl Into<String>) -> Self {
        Self::new(SceneErrorKind::WindowError(message.into()))
    }

    /// Create a render error.
    pub fn render_error(message: impl Into<String>) -> Self {
        Self::new(SceneErrorKind::RenderError(message.into()))
    }

    /// Create a surface error.
    pub fn surface_error(error: impl fmt::Display) -> Self {
        Self::new(SceneErrorKind::SurfaceError(error.to_string()))
    }

    /// Create a `wgpu` device-request error.
    pub fn device_error(error: wgpu::RequestDeviceError) -> Self {
        Self::new(SceneErrorKind::DeviceError(error))
    }

    /// Create a standard I/O error.
    pub fn io_error(error: std::io::Error) -> Self {
        Self::new(SceneErrorKind::IoError(error))
    }

    /// Return the captured backtrace.
    pub fn backtrace(&self) -> &Backtrace {
        &self.backtrace
    }

    /// Return the upstream source error, if available.
    pub fn source(&self) -> Option<&(dyn StdError + 'static)> {
        match &self.kind {
            SceneErrorKind::DeviceError(error) => Some(error),
            SceneErrorKind::IoError(error) => Some(error),
            _ => None,
        }
    }

    /// Return `true` if an object lookup failed.
    pub fn is_object_not_found(&self) -> bool {
        matches!(self.kind, SceneErrorKind::ObjectNotFound(_))
    }

    /// Return `true` if an object name already exists.
    pub fn is_object_exists(&self) -> bool {
        matches!(self.kind, SceneErrorKind::ObjectExists(_))
    }

    /// Return `true` if an object lookup or creation failed.
    pub fn is_object_error(&self) -> bool {
        matches!(
            self.kind,
            SceneErrorKind::ObjectNotFound(_)
                | SceneErrorKind::ObjectExists(_)
                | SceneErrorKind::InvalidObjectType { .. }
        )
    }

    /// Return `true` if an object had an unexpected type.
    pub fn is_invalid_object_type(&self) -> bool {
        matches!(self.kind, SceneErrorKind::InvalidObjectType { .. })
    }

    /// Return `true` if a scene lookup failed.
    pub fn is_scene_not_found(&self) -> bool {
        matches!(self.kind, SceneErrorKind::SceneNotFound(_))
    }

    /// Return `true` if a view lookup failed.
    pub fn is_view_not_found(&self) -> bool {
        matches!(self.kind, SceneErrorKind::ViewNotFound(_))
    }

    /// Return `true` if a view name already exists.
    pub fn is_view_exists(&self) -> bool {
        matches!(self.kind, SceneErrorKind::ViewExists(_))
    }

    /// Return `true` if a view lookup or creation failed.
    pub fn is_view_error(&self) -> bool {
        matches!(
            self.kind,
            SceneErrorKind::ViewNotFound(_) | SceneErrorKind::ViewExists(_)
        )
    }

    /// Return `true` if a state index was invalid.
    pub fn is_invalid_state(&self) -> bool {
        matches!(self.kind, SceneErrorKind::InvalidState { .. })
    }

    /// Return `true` if camera animation failed.
    pub fn is_animation_error(&self) -> bool {
        matches!(self.kind, SceneErrorKind::AnimationError(_))
    }

    /// Return `true` if window handling failed.
    pub fn is_window_error(&self) -> bool {
        matches!(self.kind, SceneErrorKind::WindowError(_))
    }

    /// Return `true` if rendering failed.
    pub fn is_render_error(&self) -> bool {
        matches!(self.kind, SceneErrorKind::RenderError(_))
    }

    /// Return `true` if GPU surface handling failed.
    pub fn is_surface_error(&self) -> bool {
        matches!(self.kind, SceneErrorKind::SurfaceError(_))
    }

    /// Return `true` if GPU device creation failed.
    pub fn is_device_error(&self) -> bool {
        matches!(self.kind, SceneErrorKind::DeviceError(_))
    }

    /// Return `true` if any GPU or render operation failed.
    pub fn is_render_or_gpu_error(&self) -> bool {
        matches!(
            self.kind,
            SceneErrorKind::RenderError(_)
                | SceneErrorKind::SurfaceError(_)
                | SceneErrorKind::DeviceError(_)
        )
    }

    /// Return `true` if this error wraps standard I/O.
    pub fn is_io_error(&self) -> bool {
        matches!(self.kind, SceneErrorKind::IoError(_))
    }

    /// Return the object name for object lookup errors.
    pub fn object_name(&self) -> Option<&str> {
        match &self.kind {
            SceneErrorKind::ObjectNotFound(name) | SceneErrorKind::ObjectExists(name) => Some(name),
            _ => None,
        }
    }

    /// Return the scene name for scene lookup errors.
    pub fn scene_name(&self) -> Option<&str> {
        match &self.kind {
            SceneErrorKind::SceneNotFound(name) => Some(name),
            _ => None,
        }
    }

    /// Return the view name for view lookup errors.
    pub fn view_name(&self) -> Option<&str> {
        match &self.kind {
            SceneErrorKind::ViewNotFound(name) | SceneErrorKind::ViewExists(name) => Some(name),
            _ => None,
        }
    }

    /// Return the expected and actual object types for type errors.
    pub fn object_types(&self) -> Option<(&str, &str)> {
        match &self.kind {
            SceneErrorKind::InvalidObjectType { expected, actual } => {
                Some((expected.as_str(), actual.as_str()))
            }
            _ => None,
        }
    }

    /// Return the invalid state index and available state count.
    pub fn state_index_and_count(&self) -> Option<(usize, usize)> {
        match self.kind {
            SceneErrorKind::InvalidState { index, count } => Some((index, count)),
            _ => None,
        }
    }

    /// Return the textual payload for string-backed errors.
    pub fn message(&self) -> Option<&str> {
        match &self.kind {
            SceneErrorKind::ObjectNotFound(message)
            | SceneErrorKind::ObjectExists(message)
            | SceneErrorKind::SceneNotFound(message)
            | SceneErrorKind::ViewNotFound(message)
            | SceneErrorKind::ViewExists(message)
            | SceneErrorKind::AnimationError(message)
            | SceneErrorKind::WindowError(message)
            | SceneErrorKind::RenderError(message) => Some(message),
            SceneErrorKind::InvalidObjectType { .. }
            | SceneErrorKind::InvalidState { .. }
            | SceneErrorKind::SurfaceError(_)
            | SceneErrorKind::DeviceError(_)
            | SceneErrorKind::IoError(_) => None,
        }
    }
}

impl fmt::Display for SceneError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self.kind {
            SceneErrorKind::ObjectNotFound(name) => write!(f, "Object not found: {name}"),
            SceneErrorKind::ObjectExists(name) => write!(f, "Object already exists: {name}"),
            SceneErrorKind::InvalidObjectType { expected, actual } => {
                write!(f, "Invalid object type: expected {expected}, got {actual}")
            }
            SceneErrorKind::SceneNotFound(name) => write!(f, "Scene not found: {name}"),
            SceneErrorKind::ViewNotFound(name) => write!(f, "View not found: {name}"),
            SceneErrorKind::ViewExists(name) => write!(f, "View already exists: {name}"),
            SceneErrorKind::InvalidState { index, count } => {
                write!(
                    f,
                    "Invalid state index: {index} (object has {count} states)"
                )
            }
            SceneErrorKind::AnimationError(message) => {
                write!(f, "Camera animation error: {message}")
            }
            SceneErrorKind::WindowError(message) => write!(f, "Window error: {message}"),
            SceneErrorKind::RenderError(message) => write!(f, "Render error: {message}"),
            SceneErrorKind::SurfaceError(error) => write!(f, "Surface error: {error}"),
            SceneErrorKind::DeviceError(error) => write!(f, "GPU device error: {error}"),
            SceneErrorKind::IoError(error) => write!(f, "IO error: {error}"),
        }
    }
}

impl StdError for SceneError {
    fn source(&self) -> Option<&(dyn StdError + 'static)> {
        SceneError::source(self)
    }
}

impl From<wgpu::RequestDeviceError> for SceneError {
    fn from(error: wgpu::RequestDeviceError) -> Self {
        Self::device_error(error)
    }
}

impl From<std::io::Error> for SceneError {
    fn from(error: std::io::Error) -> Self {
        Self::io_error(error)
    }
}

impl ViewerError {
    fn new(kind: ViewerErrorKind) -> Self {
        Self {
            kind,
            backtrace: Backtrace::capture(),
        }
    }

    /// Create a viewer error from a scene error.
    pub fn scene(error: SceneError) -> Self {
        Self::new(ViewerErrorKind::Scene(Box::new(error)))
    }

    /// Create a GPU initialization error.
    pub fn gpu_init_failed(message: impl Into<String>) -> Self {
        Self::new(ViewerErrorKind::GpuInitFailed(message.into()))
    }

    /// Create a file load error.
    pub fn load_error(message: impl Into<String>) -> Self {
        Self::new(ViewerErrorKind::LoadError(message.into()))
    }

    /// Create a viewer I/O integration error.
    pub fn io_error(message: impl Into<String>) -> Self {
        Self::new(ViewerErrorKind::IoError(message.into()))
    }

    /// Create a capture error.
    pub fn capture_error(message: impl Into<String>) -> Self {
        Self::new(ViewerErrorKind::CaptureError(message.into()))
    }

    /// Return the captured backtrace.
    pub fn backtrace(&self) -> &Backtrace {
        &self.backtrace
    }

    /// Return the upstream source error, if available.
    pub fn source(&self) -> Option<&(dyn StdError + 'static)> {
        match &self.kind {
            ViewerErrorKind::Scene(error) => Some(error.as_ref()),
            _ => None,
        }
    }

    /// Return `true` if this wraps a scene error.
    pub fn is_scene(&self) -> bool {
        matches!(self.kind, ViewerErrorKind::Scene(_))
    }

    /// Return `true` if GPU initialization failed.
    pub fn is_gpu_init_failed(&self) -> bool {
        matches!(self.kind, ViewerErrorKind::GpuInitFailed(_))
    }

    /// Return `true` if file loading failed.
    pub fn is_load_error(&self) -> bool {
        matches!(self.kind, ViewerErrorKind::LoadError(_))
    }

    /// Return `true` if viewer I/O integration failed.
    pub fn is_io_error(&self) -> bool {
        matches!(self.kind, ViewerErrorKind::IoError(_))
    }

    /// Return `true` if PNG or movie capture failed.
    pub fn is_capture_error(&self) -> bool {
        matches!(self.kind, ViewerErrorKind::CaptureError(_))
    }

    /// Return the wrapped scene error, if available.
    pub fn scene_error(&self) -> Option<&SceneError> {
        match &self.kind {
            ViewerErrorKind::Scene(error) => Some(error.as_ref()),
            _ => None,
        }
    }

    /// Return the textual payload for string-backed errors.
    pub fn message(&self) -> Option<&str> {
        match &self.kind {
            ViewerErrorKind::GpuInitFailed(message)
            | ViewerErrorKind::LoadError(message)
            | ViewerErrorKind::IoError(message)
            | ViewerErrorKind::CaptureError(message) => Some(message),
            ViewerErrorKind::Scene(_) => None,
        }
    }
}

impl fmt::Display for ViewerError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self.kind {
            ViewerErrorKind::Scene(error) => write!(f, "{error}"),
            ViewerErrorKind::GpuInitFailed(message) => {
                write!(f, "Failed to initialize GPU: {message}")
            }
            ViewerErrorKind::LoadError(message) => write!(f, "Failed to load file: {message}"),
            ViewerErrorKind::IoError(message) => write!(f, "IO error: {message}"),
            ViewerErrorKind::CaptureError(message) => write!(f, "Capture error: {message}"),
        }
    }
}

impl StdError for ViewerError {
    fn source(&self) -> Option<&(dyn StdError + 'static)> {
        ViewerError::source(self)
    }
}

impl From<SceneError> for ViewerError {
    fn from(error: SceneError) -> Self {
        Self::scene(error)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scene_error_display_matches_legacy_messages() {
        let err = SceneError::object_not_found("mol1");
        assert_eq!(err.to_string(), "Object not found: mol1");
        assert!(err.is_object_not_found());
        assert_eq!(err.object_name(), Some("mol1"));

        let err = SceneError::invalid_state(2, 1);
        assert_eq!(
            err.to_string(),
            "Invalid state index: 2 (object has 1 states)"
        );
        assert!(err.is_invalid_state());
        assert_eq!(err.state_index_and_count(), Some((2, 1)));
    }

    #[test]
    fn scene_error_preserves_io_source() {
        let err = SceneError::io_error(std::io::Error::new(
            std::io::ErrorKind::NotFound,
            "missing scene file",
        ));

        assert_eq!(err.to_string(), "IO error: missing scene file");
        assert!(err.is_io_error());
        assert!(StdError::source(&err).is_some());
    }

    #[test]
    fn viewer_error_display_matches_legacy_messages() {
        let err = ViewerError::capture_error("frame failed");
        assert_eq!(err.to_string(), "Capture error: frame failed");
        assert!(err.is_capture_error());
        assert_eq!(err.message(), Some("frame failed"));

        let err = ViewerError::gpu_init_failed("adapter unavailable");
        assert_eq!(
            err.to_string(),
            "Failed to initialize GPU: adapter unavailable"
        );
        assert!(err.is_gpu_init_failed());
    }

    #[test]
    fn viewer_error_preserves_scene_source() {
        let err = ViewerError::from(SceneError::scene_not_found("scene1"));

        assert_eq!(err.to_string(), "Scene not found: scene1");
        assert!(err.is_scene());
        assert!(err.scene_error().is_some());
        assert!(StdError::source(&err).is_some());
    }
}
