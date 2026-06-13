use super::state::*;
use crate::picking::ObjectId;
use crate::representations::BuildCtx;

impl RenderState {
    pub(super) fn dispatch_pending_compute_builds(&mut self) -> bool {
        let mut encoder = self
            .ctx
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("patinae.sync.compute_build"),
            });

        let dispatched = self.record_pending_compute_builds(&mut encoder);
        if dispatched {
            self.ctx.queue.submit(std::iter::once(encoder.finish()));
        }
        dispatched
    }

    pub(super) fn record_pending_compute_builds(
        &mut self,
        encoder: &mut wgpu::CommandEncoder,
    ) -> bool {
        let scene_bg = match self.scene.scene_store.bind_group() {
            Some(bg) => bg.clone(),
            None => return false,
        };
        let object_coords_scene_bg = self.scene.scene_store.object_coords_bind_group().cloned();

        let mut any_dispatched = false;

        for (key, entry) in self.scene.reps.iter_mut() {
            let dyn_offset = match self.scene.scene_store.slot(ObjectId(key.0)) {
                Some(s) => s.dynamic_offset(),
                None => continue,
            };
            let mut ctx = BuildCtx {
                encoder: &mut *encoder,
                scene_bg: &scene_bg,
                object_coords_scene_bg: object_coords_scene_bg.as_ref(),
                obj_dynamic_offset: dyn_offset,
                queue: &self.ctx.queue,
                device: &self.ctx.device,
                pipelines: &self.geometry,
            };
            if entry.rep.record_compute_build(&mut ctx) {
                any_dispatched = true;
            }
        }
        any_dispatched
    }
}
