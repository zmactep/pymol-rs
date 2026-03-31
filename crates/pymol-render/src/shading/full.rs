//! Full shading pipeline — Classic multi-light lighting + per-light directional shadow maps.
//!
//! Shadow matrices are built in view space so that the shadow direction is
//! always consistent with the view-space lighting model. The camera's view
//! matrix is folded into the shadow VP matrix so that the shader can transform
//! world-space positions directly to shadow clip space.

use crate::RenderContext;
use crate::multishadow::{ShadowParams, compute_shadow_matrix, mat4_mul};
use pymol_settings::Settings;

use super::{ShadowPassState, ShadowPipelineBase, ShadingPipeline};

pub struct FullPipeline {
    pub base: ShadowPipelineBase,
    camera_view: Option<[[f32; 4]; 4]>,
    prev_camera_view: Option<[[f32; 4]; 4]>,
}

impl Default for FullPipeline {
    fn default() -> Self {
        Self::new()
    }
}

impl FullPipeline {
    pub fn new() -> Self {
        Self {
            base: ShadowPipelineBase::new(),
            camera_view: None,
            prev_camera_view: None,
        }
    }

    /// Provide the camera view matrix. Auto-invalidates shadows when it changes.
    pub fn set_camera_view(&mut self, view: [[f32; 4]; 4]) {
        if self.prev_camera_view.as_ref() != Some(&view) {
            self.base.shadow_dirty = true;
        }
        self.camera_view = Some(view);
    }
}

/// Transform a point by a column-major 4×4 matrix (w=1 implied).
fn transform_point(m: &[[f32; 4]; 4], p: [f32; 3]) -> [f32; 3] {
    [
        m[0][0] * p[0] + m[1][0] * p[1] + m[2][0] * p[2] + m[3][0],
        m[0][1] * p[0] + m[1][1] * p[1] + m[2][1] * p[2] + m[3][1],
        m[0][2] * p[0] + m[1][2] * p[1] + m[2][2] * p[2] + m[3][2],
    ]
}

impl ShadingPipeline for FullPipeline {
    fn prepare(
        &mut self,
        context: &mut RenderContext,
        settings: &Settings,
    ) -> bool {
        let classic = &settings.shading.classic;
        let full = &settings.shading.full;
        let light_count = classic.light_count;
        let shadow_count = (light_count - 1).max(0) as u32;

        if shadow_count == 0 {
            return false;
        }

        let tile_size = full.shadow_map_size.max(64) as u32;

        if !self.base.begin_prepare(context, shadow_count, tile_size) {
            return false;
        }

        let (scene_center, scene_radius) = self.base.scene_bounds.unwrap_or(([0.0, 0.0, 0.0], 10.0));

        let light_dirs: [[f32; 3]; 9] = [
            classic.light, classic.light2, classic.light3,
            classic.light4, classic.light5, classic.light6,
            classic.light7, classic.light8, classic.light9,
        ];

        let view = self.camera_view.unwrap_or([
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ]);
        let view_center = transform_point(&view, scene_center);

        self.base.shadow_matrices = (0..shadow_count as usize)
            .map(|i| {
                let dir = light_dirs[i];
                let neg_dir = [-dir[0], -dir[1], -dir[2]];
                let shadow_vp = compute_shadow_matrix(neg_dir, view_center, scene_radius);
                mat4_mul(&shadow_vp, &view)
            })
            .collect();

        self.prev_camera_view = self.camera_view;

        let bias = full.shadow_bias;
        let intensity = full.shadow_intensity;
        let pcf_samples = full.shadow_pcf.max(1) as u32;
        let atlas = self.base.shadow_atlas.as_ref().unwrap();
        let params = ShadowParams {
            shadow_count,
            grid_size: atlas.grid_size,
            bias,
            intensity,
            mode: 2,
            pcf_samples,
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
