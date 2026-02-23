//! PyMOL-RS Scene Management
//!
//! This crate provides the scene management layer for PyMOL-RS, including:
//!
//! - [`Camera`] - Interactive camera with rotation, translation, zoom
//! - [`SceneView`] - PyMOL-compatible 25-value view state
//! - [`Object`] trait and [`ObjectRegistry`] - Named object management
//! - [`MoleculeObject`] - Molecular object with cached representations
//! - [`Scene`] and [`SceneManager`] - Named scene snapshots
//! - [`ViewManager`] - Named view storage (camera state only)
//! - [`Session`] - Pure scene state (objects, camera, settings, colors)
//! - [`SessionAdapter`] - Bridges `Session` to `ViewerLike` for command execution
//! - [`ViewerLike`] - Trait for command dispatch
//!
//! # Architecture
//!
//! The scene crate serves as the central coordination layer between:
//! - `pymol-render`: GPU rendering (wgpu pipelines, representations)
//! - `pymol-mol`: Molecular data structures
//! - `pymol-color`: Color system
//! - `pymol-settings`: Configuration
//! - `winit`: Window and event handling

mod camera;
mod capture;
pub mod serde_helpers;
mod error;
mod input;
mod keybindings;
mod movie;
mod object;
pub mod quat;
mod pick;
mod raytrace;
mod scene;
mod selection;
mod session;
mod session_adapter;
mod uniform;
mod view;
mod viewer_trait;
mod window;

// Re-export main types
pub use camera::{Camera, CameraAnimation, Projection, SceneView, normalize_matrix};
pub use capture::capture_png_to_file;
pub use error::{SceneError, SceneResult, ViewerError, WindowError};
pub use input::{CameraDelta, InputState};
pub use keybindings::{KeyBinding, KeyBindings, KeyCode};
pub use movie::{LoopMode, Movie, MovieFrame, ObjectKeyframe, PlaybackState, PlayDirection};
pub use object::{DirtyFlags, Label, LabelAnchor, LabelObject, MoleculeObject, MoleculeObjectSnapshot, Object, ObjectRegistry, ObjectRegistrySnapshot, ObjectState, ObjectType};
pub use pick::{expand_pick_to_selection, pick_expression_for_hit, PickHit, Picker};
pub use scene::{Scene, SceneAtomData, SceneManager, SceneObjectData, ScenePerAtomData, SceneStoreMask};
pub use selection::{SelectionEntry, SelectionManager};
pub use session::Session;
pub use session_adapter::SessionAdapter;
pub use raytrace::{raytrace_scene, raytrace_to_file, RaytraceError, RaytraceInput};
pub use uniform::setup_uniforms;
pub use view::ViewManager;
pub use viewer_trait::{RaytracedImage, ViewerLike};
pub use window::Window;

// Re-export types from dependencies that are part of the public API
pub use pymol_color::NamedColors;

/// Re-export commonly used types for convenience
pub mod prelude {
    pub use crate::camera::{Camera, Projection, SceneView};
    pub use crate::error::{SceneError, SceneResult};
    pub use crate::object::{MoleculeObject, Object, ObjectRegistry, ObjectState, ObjectType};
    pub use crate::session::Session;
    pub use crate::session_adapter::SessionAdapter;
    pub use crate::viewer_trait::{RaytracedImage, ViewerLike};
    pub use crate::scene::{Scene, SceneManager};
    pub use crate::view::ViewManager;
}
