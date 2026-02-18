//! Skripkin shading pipeline — multi-directional shadow map ambient occlusion.

use crate::RenderContext;
use crate::multishadow::{
    MultishadowAtlas, ShadowPipelines, ShadowParams,
    compute_shadow_matrix, create_disabled_shadow_bind_group,
    fibonacci_sphere_directions,
};
use pymol_settings::GlobalSettings;

use super::{ShadowPassState, ShadingPipeline};

/// All state needed for Skripkin AO shadow rendering.
pub struct SkripkinPipeline {
    pub shadow_atlas: Option<MultishadowAtlas>,
    pub shadow_pipelines: Option<ShadowPipelines>,
    pub shadow_dirty: bool,
    /// Scene bounding sphere (center, radius) for shadow matrix computation.
    scene_bounds: Option<([f32; 3], f32)>,
    /// View-projection matrices for the current set of shadow directions.
    pub shadow_matrices: Vec<[[f32; 4]; 4]>,
}

impl SkripkinPipeline {
    pub fn new() -> Self {
        Self {
            shadow_atlas: None,
            shadow_pipelines: None,
            shadow_dirty: true,
            scene_bounds: None,
            shadow_matrices: Vec::new(),
        }
    }

    /// Ensure shadow pipelines are created (idempotent).
    pub fn ensure_pipelines(&mut self, device: &wgpu::Device) {
        if self.shadow_pipelines.is_none() {
            self.shadow_pipelines = Some(ShadowPipelines::new(device));
        }
    }
}

impl ShadingPipeline for SkripkinPipeline {
    /// Allocate/rebuild the shadow atlas and upload sampling params.
    ///
    /// Does NOT create render passes — the caller must render shadow passes
    /// using `shadow_atlas`, `shadow_pipelines`, and `shadow_matrices` when
    /// this returns `true`.
    fn prepare(
        &mut self,
        context: &mut RenderContext,
        settings: &GlobalSettings,
    ) -> bool {
        // Force rebuild when atlas has been dropped (e.g., after deactivate).
        if self.shadow_atlas.is_none() {
            self.shadow_dirty = true;
        }
        if !self.shadow_dirty {
            return false;
        }

        let n_directions = settings.get_int(pymol_settings::id::skripkin_directions).max(1) as u32;
        let tile_size = settings.get_int(pymol_settings::id::skripkin_map_size).max(32) as u32;
        let bias = settings.get_float(pymol_settings::id::skripkin_bias);
        let intensity = settings.get_float(pymol_settings::id::skripkin_intensity);
        let device = context.device();

        self.ensure_pipelines(device);

        let needs_new_atlas = self.shadow_atlas.as_ref().map_or(true, |a| {
            a.direction_count != n_directions || a.tile_size != tile_size
        });
        if needs_new_atlas {
            self.shadow_atlas = Some(MultishadowAtlas::new(device, n_directions, tile_size));
        }

        let (scene_center, scene_radius) = self.scene_bounds.unwrap_or(([0.0, 0.0, 0.0], 10.0));
        let directions = fibonacci_sphere_directions(n_directions as usize);
        self.shadow_matrices = directions
            .iter()
            .map(|dir| compute_shadow_matrix(*dir, scene_center, scene_radius))
            .collect();

        // Upload sampling params and matrices so the main pass can read them.
        let atlas = self.shadow_atlas.as_ref().unwrap();
        let params = ShadowParams {
            shadow_count: n_directions,
            grid_size: atlas.grid_size,
            bias,
            intensity,
        };
        context.shadow_sampling().update(context.queue(), &params, &self.shadow_matrices);

        let active_bind_group = context.shadow_sampling().create_bind_group(
            context.device(),
            &atlas.sample_view,
        );
        context.set_shadow_bind_group(active_bind_group);

        // Return true — caller must now render shadow passes.
        // `shadow_dirty` is cleared by the caller after rendering.
        true
    }

    fn deactivate(&mut self, context: &mut RenderContext) {
        if self.shadow_atlas.is_some() {
            context.shadow_sampling().update(
                context.queue(),
                &ShadowParams {
                    shadow_count: 0,
                    grid_size: 0,
                    bias: 0.0,
                    intensity: 0.0,
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

    fn needs_shadow_update(&self) -> bool { self.shadow_dirty }

    fn invalidate_shadows(&mut self) { self.shadow_dirty = true; }

    fn set_scene_bounds(&mut self, center: [f32; 3], radius: f32) {
        self.scene_bounds = Some((center, radius));
    }

    fn mark_shadow_done(&mut self) { self.shadow_dirty = false; }

    fn shadow_pass_state(&self) -> Option<ShadowPassState<'_>> {
        Some(ShadowPassState {
            atlas: self.shadow_atlas.as_ref()?,
            pipelines: self.shadow_pipelines.as_ref()?,
            matrices: &self.shadow_matrices,
        })
    }
}
