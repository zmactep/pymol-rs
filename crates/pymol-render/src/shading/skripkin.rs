//! Skripkin shading pipeline — multi-directional shadow map ambient occlusion.

use crate::RenderContext;
use crate::multishadow::{
    ShadowParams, compute_shadow_matrix, fibonacci_sphere_directions,
};
use pymol_settings::GlobalSettings;

use super::{ShadowPassState, ShadowPipelineBase, ShadingPipeline};

pub struct SkripkinPipeline {
    pub base: ShadowPipelineBase,
}

impl Default for SkripkinPipeline {
    fn default() -> Self {
        Self::new()
    }
}

impl SkripkinPipeline {
    pub fn new() -> Self {
        Self {
            base: ShadowPipelineBase::new(),
        }
    }
}

impl ShadingPipeline for SkripkinPipeline {
    fn prepare(
        &mut self,
        context: &mut RenderContext,
        settings: &GlobalSettings,
    ) -> bool {
        let n_directions = settings.get_int(pymol_settings::id::skripkin_directions).max(1) as u32;
        let tile_size = settings.get_int(pymol_settings::id::skripkin_map_size).max(32) as u32;

        if !self.base.begin_prepare(context, n_directions, tile_size) {
            return false;
        }

        let (scene_center, scene_radius) = self.base.scene_bounds.unwrap_or(([0.0, 0.0, 0.0], 10.0));
        let directions = fibonacci_sphere_directions(n_directions as usize);
        self.base.shadow_matrices = directions
            .iter()
            .map(|dir| compute_shadow_matrix(*dir, scene_center, scene_radius))
            .collect();

        let bias = settings.get_float(pymol_settings::id::skripkin_bias);
        let intensity = settings.get_float(pymol_settings::id::skripkin_intensity);
        let atlas = self.base.shadow_atlas.as_ref().unwrap();
        let params = ShadowParams {
            shadow_count: n_directions,
            grid_size: atlas.grid_size,
            bias,
            intensity,
            mode: 1,
            pcf_samples: 1,
            _pad: [0; 2],
        };
        self.base.finish_prepare(context, params);
        true
    }

    fn deactivate(&mut self, context: &mut RenderContext) {
        self.base.deactivate(context);
    }

    fn needs_shadow_update(&self) -> bool { self.base.needs_shadow_update() }

    fn invalidate_shadows(&mut self) { self.base.invalidate_shadows(); }

    fn set_scene_bounds(&mut self, center: [f32; 3], radius: f32) {
        self.base.set_scene_bounds(center, radius);
    }

    fn mark_shadow_done(&mut self) { self.base.mark_shadow_done(); }

    fn shadow_pass_state(&self) -> Option<ShadowPassState<'_>> {
        self.base.shadow_pass_state()
    }
}
