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
//! - [`Window`] - winit window with wgpu surface
//! - [`Viewer`] - Main application coordinating rendering
//!
//! # Architecture
//!
//! The scene crate serves as the central coordination layer between:
//! - `pymol-render`: GPU rendering (wgpu pipelines, representations)
//! - `pymol-mol`: Molecular data structures
//! - `pymol-color`: Color system
//! - `pymol-settings`: Configuration
//! - `winit`: Window and event handling
//!
//! # Example
//!
//! ```ignore
//! use pymol_scene::{Viewer, MoleculeObject};
//! use pymol_mol::ObjectMolecule;
//!
//! // Create viewer
//! let mut viewer = pollster::block_on(Viewer::new())?;
//!
//! // Load a molecule
//! viewer.load("protein.pdb", None)?;
//!
//! // Center camera on molecule
//! viewer.center_all();
//!
//! // Run the event loop
//! viewer.run()?;
//! ```

mod camera;
mod capture;
mod error;
mod input;
mod keybindings;
mod movie;
mod object;
mod pick;
mod scene;
mod selection;
mod uniform;
mod view;
mod viewer;
mod viewer_trait;
mod window;

// Re-export main types
pub use camera::{Camera, CameraAnimation, Projection, SceneView};
pub use capture::capture_png_to_file;
pub use error::{SceneError, SceneResult, ViewerError, WindowError};
pub use input::{CameraDelta, InputState};
pub use keybindings::{KeyBinding, KeyBindings, KeyCode};
pub use movie::{LoopMode, Movie, MovieFrame, PlaybackState, PlayDirection};
pub use object::{DirtyFlags, Label, LabelAnchor, LabelObject, MoleculeObject, Object, ObjectRegistry, ObjectState, ObjectType};
pub use pick::{PickHit, Picker};
pub use scene::{Scene, SceneAtomData, SceneManager, SceneObjectData, ScenePerAtomData, SceneStoreMask};
pub use selection::{SelectionEntry, SelectionManager};
pub use uniform::setup_uniforms;
pub use view::ViewManager;
pub use viewer::{run, Viewer};
pub use viewer_trait::ViewerLike;
pub use window::Window;

// Re-export types from dependencies that are part of the public API
pub use pymol_color::NamedColors;

/// Re-export commonly used types for convenience
pub mod prelude {
    pub use crate::camera::{Camera, Projection, SceneView};
    pub use crate::error::{SceneError, SceneResult, ViewerError};
    pub use crate::object::{MoleculeObject, Object, ObjectRegistry, ObjectState, ObjectType};
    pub use crate::viewer::Viewer;
    pub use crate::viewer_trait::ViewerLike;
    pub use crate::scene::{Scene, SceneManager};
    pub use crate::view::ViewManager;
}
