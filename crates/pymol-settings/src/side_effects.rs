//! Side effect categories for setting changes.
//!
//! When settings change, the system reacts based on which categories
//! are associated with the setting (via `SettingDescriptor::side_effects`).

/// Categories of side effects for batch processing
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SideEffectCategory {
    /// Scene needs to be invalidated/redrawn
    SceneInvalidate,
    /// Scene has changed (more significant than invalidate)
    SceneChanged,
    /// Shader variables need to be reloaded
    ShaderReload,
    /// Shader computation (lighting) needs to be recalculated
    ShaderComputeLighting,
    /// Ortho/UI needs to be marked dirty
    OrthoDirty,
    /// Sequence view needs update
    SeqChanged,
    /// Stereo mode update
    StereoUpdate,
    /// Representations need to be rebuilt
    RepresentationRebuild,
    /// Representation colors need updating (no geometry rebuild)
    ColorRebuild,
    /// Full rebuild of all objects
    FullRebuild,
    /// Viewport update
    ViewportUpdate,
}
