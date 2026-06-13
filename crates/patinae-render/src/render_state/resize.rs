use super::state::*;
use crate::frame::FrameTargets;
use crate::picking::PickingMode;

impl RenderState {
    pub fn resize(&mut self, viewport: (u32, u32)) {
        let with_picking = self.picking.picking_mode != PickingMode::Disabled;

        self.targets = FrameTargets::new_with_picking(
            &self.ctx.device,
            viewport.0,
            viewport.1,
            with_picking,
            self.ctx.color_format,
        );

        self.uniforms.set_viewport(viewport.0, viewport.1);

        self.ctx.upload_frame(&self.uniforms);

        self.screen.composite_bind_group = self.screen.composite.make_bind_group(
            &self.ctx.device,
            &self.targets.accum,
            &self.targets.reveal,
        );

        let needs_selection_overlay =
            self.screen.selection_overlay_enabled && self.scene.has_any_marker;
        let needs_overlay_id = self.screen.silhouette_params.is_some() || needs_selection_overlay;
        if needs_overlay_id {
            self.targets.ensure_overlay_id_targets(&self.ctx.device);
        }
        if needs_selection_overlay {
            self.targets.ensure_marking_targets(
                &self.ctx.device,
                self.ctx.color_format,
                self.screen.fxaa_enabled,
            );
        } else {
            self.targets.clear_marking_targets();
        }
        self.screen.marking_bindings_dirty = true;
        self.screen.marking_params_dirty = true;
        self.refresh_overlay_bind_groups();

        self.scene.scene_dirty = true;
        self.invalidate_cull_cache();

        // Fresh `FrameTargets` — picking texture starts uninitialized, force
        // a re-record on the next pick/overlay regardless of camera state.
        self.invalidate_picking_cache();
        self.invalidate_overlay_id_cache();

        if let (Some(reproject), (pw, ph)) = (
            self.picking.picking_reproject.as_mut(),
            self.targets.picking_dims(),
        ) {
            reproject.ensure_best_depth_capacity(&self.ctx.device, pw, ph);
        }

        // SSAO bind groups reference the freshly allocated depth +

        // ssao textures, so they must be rebuilt on resize.

        self.screen.ssao_bind_group = self.screen.ssao_compute.make_bind_group(
            &self.ctx.device,
            &self.ctx.frame.buffer,
            &self.screen.ssao_resources.ssao_params_buffer,
            &self.targets.depth,
            &self.targets.ssao_view,
        );

        self.screen.ssao_blur_h_bind_group = self.screen.ssao_blur.make_bind_group(
            &self.ctx.device,
            &self.screen.ssao_resources.blur_h_params_buffer,
            &self.targets.depth,
            &self.targets.ssao_view,
            &self.targets.ssao_blurred_view,
        );

        self.screen.ssao_blur_v_bind_group = self.screen.ssao_blur.make_bind_group(
            &self.ctx.device,
            &self.screen.ssao_resources.blur_v_params_buffer,
            &self.targets.depth,
            &self.targets.ssao_blurred_view,
            &self.targets.ssao_view,
        );

        self.screen.ssao_compose_bind_group = self
            .screen
            .ssao_compose
            .make_bind_group(&self.ctx.device, &self.targets.ssao_view);

        // FXAA scratch view changes on resize.

        self.screen.fxaa_bind_group = self
            .screen
            .fxaa_pass
            .make_bind_group(&self.ctx.device, &self.targets.color_scratch_view);

        // Per-rep cull bind groups reference per-rep storage buffers

        // (raw/compacted/raw_count/indirect) — those are NOT recreated by

        // resize, so the cull bind groups stay valid. The frustum params

        // uniform is refreshed every frame inside `render()`.
    }
}
