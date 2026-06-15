//! Scene management.
//!
//! This crate provides the scene management layer, including:
//!
//! - [`Camera`] - Interactive camera with rotation, translation, zoom
//! - [`SceneView`] - 25-value view state
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
//! The scene crate coordinates molecular data, color resolution, settings,
//! camera state, object registries, selections, and command-facing adapters.

#[cfg(feature = "render-bridge")]
pub mod bridge;
mod camera;
mod error;
mod gpu_runtime;
pub mod highlight_state;
mod input;
mod keybindings;
mod movie;
mod object;
mod pick;
pub mod quat;
mod render_target;
mod scene;
mod selection;
pub mod serde_helpers;
mod session;
mod session_adapter;
mod view;
mod viewer_trait;
// Re-export main types
pub use camera::{normalize_matrix, Camera, CameraAnimation, Projection, SceneView};
pub use error::{SceneError, SceneResult, ViewerError};
pub use gpu_runtime::{
    GpuAddressMode, GpuBatchCommand, GpuBatchResult, GpuBindGroupDescriptor, GpuBindGroupEntry,
    GpuBindGroupLayoutDescriptor, GpuBindGroupLayoutEntry, GpuBindingResource, GpuBindingType,
    GpuBufferBinding, GpuBufferBindingType, GpuBufferDescriptor, GpuBufferUsage, GpuCacheStats,
    GpuCacheStatus, GpuCachedHandle, GpuCompareFunction, GpuComputePipelineDescriptor,
    GpuDeviceLimits, GpuExtent3d, GpuFilterMode, GpuHandle, GpuHandleKind, GpuOrigin3d,
    GpuPipelineLayoutDescriptor, GpuSamplerBindingType, GpuSamplerDescriptor,
    GpuShaderModuleDescriptor, GpuShaderStages, GpuStorageTextureAccess, GpuSubmitBatch,
    GpuTexelCopyBufferLayout, GpuTexelCopyTextureInfo, GpuTextureAspect, GpuTextureDescriptor,
    GpuTextureDimension, GpuTextureFormat, GpuTextureSampleType, GpuTextureUsage,
    GpuTextureViewDescriptor, GpuTextureViewDimension, RenderArtifactBufferDescriptor,
    RenderArtifactBufferRole, RenderArtifactPrimitiveTopology, RenderArtifactRepDescriptor,
    RenderArtifactRepKind, RenderArtifactSnapshotDescriptor,
};
pub use highlight_state::HighlightState;
pub use input::{ButtonState, CameraDelta, InputState, Modifiers, MouseButton, ScrollDelta};
pub use keybindings::{parse_key_string, KeyBinding, KeyBindings, KeyCode};
pub use movie::{LoopMode, Movie, MovieFrame, ObjectKeyframe, PlayDirection, PlaybackState};
pub use object::{
    DirtyFlags, GroupObject, Label, LabelAnchor, LabelObject, MapData, MapDisplayMode, MapObject,
    Measurement, MeasurementObject, MeasurementType, MoleculeObject, MoleculeObjectSnapshot,
    Object, ObjectRegistry, ObjectRegistrySnapshot, ObjectState, ObjectType, RenderObjectId,
};
pub use pick::{expand_pick_to_selection, pick_expression_for_hit, PickHit};
pub use render_target::CaptureRenderer;
pub use scene::{
    Scene, SceneAtomData, SceneManager, SceneObjectData, ScenePerAtomData, SceneStoreMask,
};
pub use selection::{SelectionEntry, SelectionManager};
pub use session::{AnimationUpdate, HoverTarget, MovieStateSnapshot, Session};
pub use session_adapter::SessionAdapter;
pub use view::ViewManager;
pub use viewer_trait::{ViewerLike, ViewportImage};
// Re-export types from dependencies that are part of the public API
pub use patinae_color::NamedPalette;

/// Re-export commonly used types for convenience
pub mod prelude {
    pub use crate::camera::{Camera, Projection, SceneView};
    pub use crate::error::{SceneError, SceneResult};
    pub use crate::object::{
        MoleculeObject, Object, ObjectRegistry, ObjectState, ObjectType, RenderObjectId,
    };
    pub use crate::scene::{Scene, SceneManager};
    pub use crate::session::Session;
    pub use crate::session_adapter::SessionAdapter;
    pub use crate::view::ViewManager;
    pub use crate::viewer_trait::{ViewerLike, ViewportImage};
}
