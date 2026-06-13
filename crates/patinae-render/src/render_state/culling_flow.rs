use crate::compute::cull::frustum_planes_from_view_proj;
use crate::picking::ObjectId;
use crate::representations::catalog;
use crate::representations::{CullPlanCtx, ViewportLodCtx};
#[cfg(feature = "stats")]
use crate::stats::Pass as StatsPass;

use super::RenderState;

impl RenderState {
    pub(super) fn invalidate_cull_cache(&mut self) {
        self.scene.cull_pass_initialized = false;
        self.scene.last_cull_view_proj_hash = 0;
    }

    pub(super) fn poll_viewport_lod(&mut self) {
        let queue = self.ctx.queue.clone();
        let mut changed = false;
        for entry in self.scene.reps.values_mut() {
            if entry.rep.poll_viewport_lod(&queue) {
                entry.draw_phase = entry.rep.draw_phase();
                changed = true;
            }
        }

        let _ = self.ctx.device.poll(wgpu::PollType::Poll);

        for entry in self.scene.reps.values_mut() {
            if entry.rep.poll_viewport_lod(&queue) {
                entry.draw_phase = entry.rep.draw_phase();
                changed = true;
            }
        }

        if changed {
            self.scene.scene_dirty = true;
            self.invalidate_cull_cache();
            self.invalidate_picking_cache();
            self.invalidate_overlay_id_cache();
        }
    }

    /// Dispatch the per-rep cull kernel when the cached compacted buffers
    /// no longer match the scene or camera.
    pub(super) fn dispatch_cull(
        &mut self,
        encoder: &mut wgpu::CommandEncoder,
        view_proj_hash: u64,
        compute_rebuilt: bool,
    ) {
        if !should_dispatch_cull(
            self.scene.cull_pass_initialized,
            self.scene.scene_dirty,
            compute_rebuilt,
            self.scene.last_cull_view_proj_hash,
            view_proj_hash,
        ) {
            return;
        }

        let view_proj = self.uniforms.view_proj;
        let frustum_planes = frustum_planes_from_view_proj(&view_proj);
        let queue = self.ctx.queue.clone();
        let plan_ctx = CullPlanCtx {
            queue: &queue,
            view_proj,
            frustum_planes,
        };
        let mut dispatched: u32 = 0;

        for entry in self.scene.reps.values_mut() {
            let kind = entry.rep.kind();
            let Some(plan) = entry.rep.plan_cull(&plan_ctx) else {
                continue;
            };
            let Some(rep_catalog) = catalog::entry(kind).filter(|rep| rep.cullable) else {
                continue;
            };
            let Some(pipeline) = rep_catalog.cull_pipeline(&self.geometry.cull_pipeline) else {
                continue;
            };
            let label = rep_catalog.cull_label();
            let ts: Option<wgpu::ComputePassTimestampWrites<'_>> = {
                #[cfg(feature = "stats")]
                {
                    if dispatched == 0 {
                        self.stats.compute_pass_timestamp_writes(StatsPass::Cull)
                    } else {
                        None
                    }
                }
                #[cfg(not(feature = "stats"))]
                {
                    let _ = dispatched;
                    None
                }
            };
            self.geometry.cull_pipeline.dispatch_kind(
                encoder,
                pipeline,
                plan.bind_group,
                plan.upper,
                label,
                ts,
            );
            dispatched += 1;
        }

        self.scene.cull_pass_initialized = true;
        self.scene.last_cull_view_proj_hash = view_proj_hash;
    }

    pub(super) fn record_viewport_lod_readbacks(&mut self, encoder: &mut wgpu::CommandEncoder) {
        let Some(scene_bg) = self.scene.scene_store.bind_group().cloned() else {
            return;
        };
        for (key, entry) in self.scene.reps.iter_mut() {
            let Some(slot) = self.scene.scene_store.slot(ObjectId(key.0)) else {
                continue;
            };
            let mut ctx = ViewportLodCtx {
                encoder,
                scene_bg: &scene_bg,
                obj_dynamic_offset: slot.dynamic_offset(),
                pipelines: &self.geometry,
            };
            entry.rep.record_viewport_lod_readback(&mut ctx);
        }
    }
}

fn should_dispatch_cull(
    cull_pass_initialized: bool,
    scene_dirty: bool,
    compute_rebuilt: bool,
    last_view_proj_hash: u64,
    view_proj_hash: u64,
) -> bool {
    !cull_pass_initialized
        || scene_dirty
        || compute_rebuilt
        || last_view_proj_hash != view_proj_hash
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn frustum_plane_extraction_returns_normalized_planes() {
        let planes = frustum_planes_from_view_proj(&[
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ]);
        for plane in planes {
            let len = (plane[0] * plane[0] + plane[1] * plane[1] + plane[2] * plane[2]).sqrt();
            assert!((len - 1.0).abs() < 1e-5);
            assert!(plane.iter().all(|v| v.is_finite()));
        }
    }

    #[test]
    fn cull_cache_dispatches_first_frame() {
        assert!(should_dispatch_cull(false, false, false, 7, 7));
    }

    #[test]
    fn cull_cache_reuses_static_frame() {
        assert!(!should_dispatch_cull(true, false, false, 7, 7));
    }

    #[test]
    fn cull_cache_dispatches_when_camera_changes() {
        assert!(should_dispatch_cull(true, false, false, 7, 8));
    }

    #[test]
    fn cull_cache_dispatches_when_scene_dirty() {
        assert!(should_dispatch_cull(true, true, false, 7, 7));
    }

    #[test]
    fn cull_cache_dispatches_when_compute_rebuilt() {
        assert!(should_dispatch_cull(true, false, true, 7, 7));
    }

    #[test]
    fn cull_cache_reuses_marker_only_dirty() {
        assert!(!should_dispatch_cull(true, false, false, 7, 7));
    }
}
