//! Shading pipeline abstraction.
//!
//! Each shading mode (Classic, Skripkin) is encapsulated as a struct that
//! implements [`ShadingPipeline`]. [`ShadingManager`] owns all mode structs and
//! delegates to whichever is active, handling mode transitions automatically.

pub mod classic;
pub mod skripkin;

pub use classic::ClassicPipeline;
pub use skripkin::SkripkinPipeline;

use crate::multishadow::{MultishadowAtlas, ShadowPipelines};
use crate::RenderContext;
use pymol_settings::GlobalSettings;
use pymol_settings::ShadingMode;

/// Resources required to drive shadow depth render passes.
/// Returned by `ShadingPipeline::shadow_pass_state()` when passes are needed.
pub struct ShadowPassState<'a> {
    pub atlas: &'a MultishadowAtlas,
    pub pipelines: &'a ShadowPipelines,
    pub matrices: &'a [[[f32; 4]; 4]],
}

/// Per-mode rendering pipeline.
///
/// `prepare()` performs all non-rendering setup for the frame (atlas allocation,
/// matrix computation, params upload). It does NOT create render passes — the
/// caller is responsible for driving shadow passes via the pipeline's public
/// state (see [`SkripkinPipeline`]). Returns `true` if shadow render passes
/// need to be executed this frame.
pub trait ShadingPipeline {
    /// Perform per-frame setup. Returns `true` if shadow passes must be rendered.
    fn prepare(
        &mut self,
        context: &mut RenderContext,
        settings: &GlobalSettings,
    ) -> bool;

    /// Clean up when switching away from this mode.
    fn deactivate(&mut self, context: &mut RenderContext);

    /// Whether shadow maps need rebuilding.
    fn needs_shadow_update(&self) -> bool;

    /// Mark shadow maps as stale (call when geometry changes).
    fn invalidate_shadows(&mut self);

    /// Update the scene bounding sphere used for pre-pass computations
    /// (e.g. shadow matrix coverage). No-op for modes that don't need it.
    fn set_scene_bounds(&mut self, _center: [f32; 3], _radius: f32) {}

    /// Called after the caller has rendered shadow passes to mark them complete.
    /// No-op for modes that don't produce shadow passes.
    fn mark_shadow_done(&mut self) {}

    /// Return the resources needed to drive shadow render passes, or `None` if
    /// this mode does not require them. Only valid after a `prepare()` call that
    /// returned `true`.
    fn shadow_pass_state(&self) -> Option<ShadowPassState<'_>> { None }
}

/// Owns all shading mode structs and routes calls to the active one.
pub struct ShadingManager {
    pub classic: ClassicPipeline,
    pub skripkin: SkripkinPipeline,
    pub active_mode: ShadingMode,
}

impl ShadingManager {
    pub fn new() -> Self {
        Self {
            classic: ClassicPipeline,
            skripkin: SkripkinPipeline::new(),
            active_mode: ShadingMode::Classic,
        }
    }

    /// Switch to a new mode, calling `deactivate()` on the old one if needed.
    pub fn set_mode(&mut self, mode: ShadingMode, context: &mut RenderContext) {
        if mode != self.active_mode {
            match self.active_mode {
                ShadingMode::Classic => self.classic.deactivate(context),
                ShadingMode::Skripkin => self.skripkin.deactivate(context),
            }
            self.active_mode = mode;
            context.set_active_shading_mode(mode);
        }
    }

    /// Run the active pipeline's prepare step. Returns `true` if shadow passes
    /// need to be rendered this frame (only relevant for Skripkin mode).
    pub fn prepare(
        &mut self,
        context: &mut RenderContext,
        settings: &GlobalSettings,
    ) -> bool {
        match self.active_mode {
            ShadingMode::Classic => self.classic.prepare(context, settings),
            ShadingMode::Skripkin => self.skripkin.prepare(context, settings),
        }
    }

    /// Invalidate shadows in all modes (only relevant for Skripkin).
    pub fn invalidate_shadows(&mut self) {
        self.skripkin.invalidate_shadows();
    }

    /// Forward scene bounds to the active pipeline.
    pub fn set_scene_bounds(&mut self, center: [f32; 3], radius: f32) {
        match self.active_mode {
            ShadingMode::Classic => self.classic.set_scene_bounds(center, radius),
            ShadingMode::Skripkin => self.skripkin.set_scene_bounds(center, radius),
        }
    }

    /// Notify the active pipeline that shadow passes have been rendered.
    pub fn finish_shadow_passes(&mut self) {
        match self.active_mode {
            ShadingMode::Classic => self.classic.mark_shadow_done(),
            ShadingMode::Skripkin => self.skripkin.mark_shadow_done(),
        }
    }

    /// Return shadow pass resources for the active pipeline, if any.
    pub fn shadow_pass_state(&self) -> Option<ShadowPassState<'_>> {
        match self.active_mode {
            ShadingMode::Classic => self.classic.shadow_pass_state(),
            ShadingMode::Skripkin => self.skripkin.shadow_pass_state(),
        }
    }
}
