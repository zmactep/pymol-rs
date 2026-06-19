use super::state::RenderState;
use crate::memory::{buffer_usage, GpuMemoryCategory, GpuMemoryLedger, GpuMemorySnapshot};

impl RenderState {
    /// Returns a deterministic estimated GPU-memory snapshot.
    pub fn memory_snapshot(&self) -> GpuMemorySnapshot {
        let mut ledger = GpuMemoryLedger::new();

        self.targets.record_memory(&mut ledger);
        ledger.add_usage(
            GpuMemoryCategory::FrameTargets,
            buffer_usage(&self.ctx.frame.buffer),
        );
        ledger.add_usage(
            GpuMemoryCategory::SceneStore,
            self.scene.scene_store.memory_usage(),
        );

        for entry in self.scene.reps.values() {
            ledger.add_usage(GpuMemoryCategory::Representation, entry.rep.memory_usage());
            if let Some(picking) = entry.picking.as_ref() {
                ledger.add_usage(GpuMemoryCategory::Picking, buffer_usage(&picking._buffer));
            }
        }
        for map in self.scene.maps.values() {
            ledger.add_usage(GpuMemoryCategory::Representation, map.memory_usage());
        }

        ledger.add_usage(
            GpuMemoryCategory::RepresentationScratch,
            self.geometry.dot_params_layout.memory_usage(),
        );
        ledger.add_usage(
            GpuMemoryCategory::SurfaceScratch,
            buffer_usage(&self.geometry.surface_mc_compute.tables_buf),
        );

        if let Some(readback) = self.picking.readback.as_ref() {
            ledger.add_usage(GpuMemoryCategory::Readback, readback.memory_usage());
        }
        if let Some(readback) = self.picking.hover_readback.as_ref() {
            ledger.add_usage(GpuMemoryCategory::Readback, readback.memory_usage());
        }
        if let Some(reproject) = self.picking.picking_reproject.as_ref() {
            ledger.add_usage(GpuMemoryCategory::Picking, reproject.memory_usage());
        }

        if let Some(marking) = self.screen.marking.as_ref() {
            ledger.add_usage(
                GpuMemoryCategory::Overlay,
                buffer_usage(&marking.params_buffer),
            );
            ledger.add_usage(
                GpuMemoryCategory::Overlay,
                buffer_usage(&marking.object_offsets_buffer),
            );
        }
        if let Some(silhouette) = self.screen.silhouette.as_ref() {
            ledger.add_usage(
                GpuMemoryCategory::Overlay,
                buffer_usage(&silhouette.uniform_buffer),
            );
        }

        ledger.add_usage(
            GpuMemoryCategory::Postprocess,
            buffer_usage(&self.screen.ssao_resources.ssao_params_buffer),
        );
        ledger.add_usage(
            GpuMemoryCategory::Postprocess,
            buffer_usage(&self.screen.ssao_resources.blur_h_params_buffer),
        );
        ledger.add_usage(
            GpuMemoryCategory::Postprocess,
            buffer_usage(&self.screen.ssao_resources.blur_v_params_buffer),
        );
        ledger.add_usage(
            GpuMemoryCategory::Postprocess,
            buffer_usage(&self.screen.ssao_compose.params_buffer),
        );
        ledger.add_usage(
            GpuMemoryCategory::Postprocess,
            buffer_usage(&self.screen.fxaa_pass.params_buffer),
        );

        ledger.add_usage(GpuMemoryCategory::Shadow, self.ctx.lighting.memory_usage());
        ledger.add_usage(
            GpuMemoryCategory::Shadow,
            self.lighting.directional_shadow.memory_usage(),
        );
        ledger.add_usage(
            GpuMemoryCategory::Shadow,
            self.lighting.atlas_ao.memory_usage(),
        );

        ledger.snapshot()
    }
}
