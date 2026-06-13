use super::state::*;
#[cfg(feature = "stats")]
use crate::stats::Pass as StatsPass;

impl RenderState {
    pub(super) fn refresh_marking_resources(&mut self) {
        let (Some(marking), Some(overlay_id), Some(marker_lut), Some(mask)) = (
            self.screen.marking.as_ref(),
            self.targets.overlay_id.as_ref(),
            self.scene.scene_store.marker_lut_buffer(),
            self.targets.marking_mask.as_ref(),
        ) else {
            self.screen.marking_bind_groups = None;
            self.screen.marking_bindings_dirty = true;
            self.screen.marking_params_dirty = true;
            self.screen.marking_offsets_dirty = true;

            return;
        };

        if self.screen.marking_offsets_dirty {
            marking
                .upload_object_offsets(&self.ctx.queue, self.scene.scene_store.iter_atom_offsets());
            self.screen.marking_offsets_dirty = false;
        }

        if self.screen.marking_params_dirty {
            marking.upload_params(
                &self.ctx.queue,
                (self.targets.width, self.targets.height),
                self.screen.marking_width,
            );
            self.screen.marking_params_dirty = false;
        }

        if self.screen.marking_bind_groups.is_none() || self.screen.marking_bindings_dirty {
            self.screen.marking_bind_groups = Some(marking.make_bind_groups(
                &self.ctx.device,
                overlay_id,
                marker_lut,
                &self.targets.color_scratch_view,
                mask,
            ));
            self.screen.marking_bindings_dirty = false;
        }
    }

    pub(super) fn refresh_overlay_bind_groups(&mut self) {
        self.screen.silhouette_bind_group = match (
            self.screen.silhouette.as_ref(),
            self.targets.overlay_id.as_ref(),
        ) {
            (Some(s), Some(overlay_id)) => Some(s.make_bind_group(&self.ctx.device, overlay_id)),
            _ => None,
        };

        self.screen.fxaa_overlay_bind_group = match self.targets.overlay_color_scratch.as_ref() {
            Some(view) => Some(
                self.screen
                    .fxaa_pass
                    .make_bind_group(&self.ctx.device, view),
            ),
            None => None,
        };

        self.refresh_marking_resources();
    }

    pub(super) fn record_ssao_chain(
        &mut self,
        encoder: &mut wgpu::CommandEncoder,
        target: &wgpu::TextureView,
    ) {
        if !self.screen.ssao_enabled {
            return;
        }
        self.screen.ssao_frame_counter = self.screen.ssao_frame_counter.wrapping_add(1);

        #[cfg(feature = "stats")]
        let ssao_ts = self.stats.compute_pass_timestamp_writes(StatsPass::Ssao);
        #[cfg(not(feature = "stats"))]
        let ssao_ts = None;
        self.screen.ssao_compute.dispatch(
            encoder,
            &self.screen.ssao_bind_group,
            &self.targets,
            ssao_ts,
        );

        #[cfg(feature = "stats")]
        let blur_ts = self
            .stats
            .compute_pass_timestamp_writes(StatsPass::SsaoBlur);
        #[cfg(not(feature = "stats"))]
        let blur_ts = None;
        self.screen.ssao_blur.dispatch(
            encoder,
            &self.screen.ssao_blur_h_bind_group,
            &self.targets,
            blur_ts,
        );
        self.screen.ssao_blur.dispatch(
            encoder,
            &self.screen.ssao_blur_v_bind_group,
            &self.targets,
            None,
        );

        #[cfg(feature = "stats")]
        let compose_ts = self
            .stats
            .render_pass_timestamp_writes(StatsPass::SsaoCompose);
        #[cfg(not(feature = "stats"))]
        let compose_ts = None;
        self.screen.ssao_compose.record(
            target,
            encoder,
            &self.screen.ssao_compose_bind_group,
            compose_ts,
        );
    }

    pub(super) fn record_marking_overlay(
        &self,
        encoder: &mut wgpu::CommandEncoder,
        target: &wgpu::TextureView,
        has_marker: bool,
    ) {
        if !has_marker {
            return;
        }
        let (Some(marking), Some(marking_bg), Some(mask)) = (
            self.screen.marking.as_ref(),
            self.screen.marking_bind_groups.as_ref(),
            self.targets.marking_mask.as_ref(),
        ) else {
            return;
        };

        #[cfg(feature = "stats")]
        let marking_ts = self.stats.render_pass_timestamp_writes(StatsPass::Marking);
        #[cfg(not(feature = "stats"))]
        let marking_ts = None;
        marking.record_mask(encoder, mask, &marking_bg.mask, marking_ts);
        marking.record_composite(encoder, target, &marking_bg.composite);
    }

    pub(super) fn record_silhouette_overlay(
        &self,
        encoder: &mut wgpu::CommandEncoder,
        target: &wgpu::TextureView,
    ) {
        let (Some(silhouette), Some(silhouette_bg), Some(params)) = (
            self.screen.silhouette.as_ref(),
            self.screen.silhouette_bind_group.as_ref(),
            self.screen.silhouette_params,
        ) else {
            return;
        };
        silhouette.upload_params(&self.ctx.queue, params);

        #[cfg(feature = "stats")]
        let silhouette_ts = self
            .stats
            .render_pass_timestamp_writes(StatsPass::Silhouette);
        #[cfg(not(feature = "stats"))]
        let silhouette_ts = None;
        silhouette.record(encoder, target, silhouette_bg, silhouette_ts);
    }

    pub(super) fn record_fxaa(
        &self,
        encoder: &mut wgpu::CommandEncoder,
        target: &wgpu::TextureView,
        use_overlay_source: bool,
    ) {
        if !self.screen.fxaa_enabled {
            return;
        }
        let bg = if use_overlay_source {
            self.screen
                .fxaa_overlay_bind_group
                .as_ref()
                .unwrap_or(&self.screen.fxaa_bind_group)
        } else {
            &self.screen.fxaa_bind_group
        };
        #[cfg(feature = "stats")]
        let fxaa_ts = self.stats.render_pass_timestamp_writes(StatsPass::Fxaa);
        #[cfg(not(feature = "stats"))]
        let fxaa_ts = None;
        self.screen.fxaa_pass.record(target, encoder, bg, fxaa_ts);
    }
}
