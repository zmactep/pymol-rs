//! Shading pipeline abstraction.
//!
//! Each shading mode (Classic, Skripkin) is encapsulated as a struct that
//! implements [`ShadingPipeline`]. [`ShadingManager`] owns all mode structs and
//! delegates to whichever is active, handling mode transitions automatically.

pub mod classic;
pub mod full;
pub mod skripkin;

pub use classic::ClassicPipeline;
pub use full::FullPipeline;
pub use skripkin::SkripkinPipeline;

use crate::multishadow::{
    MultishadowAtlas, ShadowParams, ShadowPipelines, create_disabled_shadow_bind_group,
};
use crate::RenderContext;
use pymol_settings::Settings;
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
        settings: &Settings,
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

/// Shared state and logic for shadow-producing pipelines (Skripkin, Full).
pub struct ShadowPipelineBase {
    pub shadow_atlas: Option<MultishadowAtlas>,
    pub shadow_pipelines: Option<ShadowPipelines>,
    pub shadow_dirty: bool,
    pub scene_bounds: Option<([f32; 3], f32)>,
    pub shadow_matrices: Vec<[[f32; 4]; 4]>,
}

impl Default for ShadowPipelineBase {
    fn default() -> Self {
        Self::new()
    }
}

impl ShadowPipelineBase {
    pub fn new() -> Self {
        Self {
            shadow_atlas: None,
            shadow_pipelines: None,
            shadow_dirty: true,
            scene_bounds: None,
            shadow_matrices: Vec::new(),
        }
    }

    pub fn ensure_pipelines(&mut self, device: &wgpu::Device) {
        if self.shadow_pipelines.is_none() {
            self.shadow_pipelines = Some(ShadowPipelines::new(device));
        }
    }

    /// Check dirty state, ensure pipelines, and resize atlas if needed.
    /// Returns `true` if the caller should compute matrices and call `finish_prepare`.
    pub fn begin_prepare(
        &mut self,
        context: &mut RenderContext,
        shadow_count: u32,
        tile_size: u32,
    ) -> bool {
        if self.shadow_atlas.is_none() {
            self.shadow_dirty = true;
        }
        if !self.shadow_dirty {
            return false;
        }

        let device = context.device();
        self.ensure_pipelines(device);

        let needs_new_atlas = self.shadow_atlas.as_ref().is_none_or(|a| {
            a.direction_count != shadow_count || a.tile_size != tile_size
        });
        if needs_new_atlas {
            self.shadow_atlas = Some(MultishadowAtlas::new(device, shadow_count, tile_size));
        }

        true
    }

    /// Upload shadow params and matrices, create bind group.
    pub fn finish_prepare(&self, context: &mut RenderContext, params: ShadowParams) {
        context
            .shadow_sampling()
            .update(context.queue(), &params, &self.shadow_matrices);

        let atlas = self.shadow_atlas.as_ref().unwrap();
        let active_bind_group = context
            .shadow_sampling()
            .create_bind_group(context.device(), &atlas.sample_view);
        context.set_shadow_bind_group(active_bind_group);
    }

    pub fn deactivate(&mut self, context: &mut RenderContext) {
        if self.shadow_atlas.is_some() {
            context.shadow_sampling().update(
                context.queue(),
                &ShadowParams {
                    shadow_count: 0,
                    grid_size: 0,
                    bias: 0.0,
                    intensity: 0.0,
                    mode: 0,
                    pcf_samples: 0,
                    _pad: [0; 2],
                },
                &[],
            );
            let disabled = create_disabled_shadow_bind_group(
                context.device(),
                context.shadow_sampling(),
            );
            context.set_shadow_bind_group(disabled);
            self.shadow_atlas = None;
        }
    }

    pub fn needs_shadow_update(&self) -> bool {
        self.shadow_dirty
    }

    pub fn invalidate_shadows(&mut self) {
        self.shadow_dirty = true;
    }

    pub fn set_scene_bounds(&mut self, center: [f32; 3], radius: f32) {
        self.scene_bounds = Some((center, radius));
    }

    pub fn mark_shadow_done(&mut self) {
        self.shadow_dirty = false;
    }

    pub fn shadow_pass_state(&self) -> Option<ShadowPassState<'_>> {
        Some(ShadowPassState {
            atlas: self.shadow_atlas.as_ref()?,
            pipelines: self.shadow_pipelines.as_ref()?,
            matrices: &self.shadow_matrices,
        })
    }
}

/// Owns all shading mode structs and routes calls to the active one.
pub struct ShadingManager {
    pub classic: ClassicPipeline,
    pub skripkin: SkripkinPipeline,
    pub full: FullPipeline,
    pub active_mode: ShadingMode,
}

impl Default for ShadingManager {
    fn default() -> Self {
        Self::new()
    }
}

impl ShadingManager {
    pub fn new() -> Self {
        Self {
            classic: ClassicPipeline,
            skripkin: SkripkinPipeline::new(),
            full: FullPipeline::new(),
            active_mode: ShadingMode::Classic,
        }
    }

    /// Switch to a new mode, calling `deactivate()` on the old one if needed.
    pub fn set_mode(&mut self, mode: ShadingMode, context: &mut RenderContext) {
        if mode != self.active_mode {
            match self.active_mode {
                ShadingMode::Classic => self.classic.deactivate(context),
                ShadingMode::Skripkin => self.skripkin.deactivate(context),
                ShadingMode::Full => self.full.deactivate(context),
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
        settings: &Settings,
    ) -> bool {
        match self.active_mode {
            ShadingMode::Classic => self.classic.prepare(context, settings),
            ShadingMode::Skripkin => self.skripkin.prepare(context, settings),
            ShadingMode::Full => self.full.prepare(context, settings),
        }
    }

    /// Provide the camera view matrix for Full mode shadow computation.
    pub fn set_camera_view(&mut self, view: [[f32; 4]; 4]) {
        self.full.set_camera_view(view);
    }

    /// Invalidate shadows in all modes (only relevant for Skripkin).
    pub fn invalidate_shadows(&mut self) {
        self.skripkin.invalidate_shadows();
        self.full.invalidate_shadows();
    }

    /// Forward scene bounds to the active pipeline.
    pub fn set_scene_bounds(&mut self, center: [f32; 3], radius: f32) {
        match self.active_mode {
            ShadingMode::Classic => self.classic.set_scene_bounds(center, radius),
            ShadingMode::Skripkin => self.skripkin.set_scene_bounds(center, radius),
            ShadingMode::Full => self.full.set_scene_bounds(center, radius),
        }
    }

    /// Notify the active pipeline that shadow passes have been rendered.
    pub fn finish_shadow_passes(&mut self) {
        match self.active_mode {
            ShadingMode::Classic => self.classic.mark_shadow_done(),
            ShadingMode::Skripkin => self.skripkin.mark_shadow_done(),
            ShadingMode::Full => self.full.mark_shadow_done(),
        }
    }

    /// Return shadow pass resources for the active pipeline, if any.
    pub fn shadow_pass_state(&self) -> Option<ShadowPassState<'_>> {
        match self.active_mode {
            ShadingMode::Classic => self.classic.shadow_pass_state(),
            ShadingMode::Skripkin => self.skripkin.shadow_pass_state(),
            ShadingMode::Full => self.full.shadow_pass_state(),
        }
    }
}
