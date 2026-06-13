//! WebViewer — the main WASM-exported type.
//!
//! Mirrors patinae's viewport: owns a `Session` + `CommandExecutor` and a
//! `patinae_render::RenderState` (held by `GpuState`).

use serde::Serialize;
use wasm_bindgen::prelude::*;

use patinae_cmd::{CommandAction, CommandExecutor, MessageKind};
use patinae_render::picking::readback::{PendingPick, PickReadbackTarget};
use patinae_render::{FrameStatsHistory, PickingMode, RenderConfig};
use patinae_scene::bridge::{resolve_pick, CachedRenderScene};
use patinae_scene::{
    expand_pick_to_selection, pick_expression_for_hit, CameraDelta, InputState, MoleculeObject,
    Object, PickHit, Session, SessionAdapter,
};
use patinae_select::{build_sele_command, select};

use crate::picking::PickHitInfo;

use crate::event;
use crate::gpu::GpuState;
use crate::render_loop;

/// Browser WebGPU readbacks can disturb frame pacing when every mousemove
/// submits a hover pick. Only large Auto-LOD scenes use this cap; small scenes
/// and cartoon-only hover keep the immediate path.
const LARGE_LOD_HOVER_PICK_INTERVAL_MS: f64 = 33.0;

// ---------------------------------------------------------------------------
// Serializable types returned to JS
// ---------------------------------------------------------------------------

#[derive(Serialize)]
struct CmdOutput {
    messages: Vec<OutputMsg>,
}

#[derive(Serialize)]
struct OutputMsg {
    level: &'static str,
    text: String,
}

#[derive(Serialize)]
struct ObjectInfo {
    name: String,
    object_type: &'static str,
    atom_count: usize,
    enabled: bool,
}

#[derive(Serialize)]
struct SequenceChain {
    object_name: String,
    chain_id: String,
    residues: Vec<SequenceResidue>,
}

#[derive(Serialize)]
struct SequenceResidue {
    resn: String,
    resv: i32,
    one_letter: String,
}

#[derive(Serialize)]
struct SelectionInfo {
    name: String,
    expression: String,
    visible: bool,
}

#[derive(Serialize)]
struct LabelInfo {
    x: f32,
    y: f32,
    text: String,
    kind: &'static str,
}

#[derive(Serialize)]
struct WebPerformanceSnapshot {
    render_count: u64,
    avg_render_ms: f32,
    median_render_ms: f32,
    p95_render_ms: f32,
    last_render_ms: f32,
    last_poll_picks_ms: f32,
    last_uniforms_ms: f32,
    last_prepare_ms: f32,
    last_sync_ms: f32,
    last_sync_scene_store_object_ms: f32,
    last_sync_scene_store_flush_ms: f32,
    last_sync_marking_resources_ms: f32,
    last_sync_rep_ms: f32,
    last_sync_map_ms: f32,
    last_sync_order_bounds_ms: f32,
    last_sync_compute_dispatch_ms: f32,
    last_sync_marker_lut_upload_bytes: u64,
    last_sync_marker_lut_upload_ranges: u32,
    last_sync_marker_lut_reallocated: bool,
    last_sync_scene_store_live_atoms: u64,
    last_sync_scene_store_allocated_atoms: u64,
    last_sync_scene_store_orphaned_atoms: u64,
    last_sync_scene_store_live_bonds: u64,
    last_sync_scene_store_allocated_bonds: u64,
    last_sync_scene_store_orphaned_bonds: u64,
    last_sync_scene_store_live_table_slots: u64,
    last_sync_scene_store_allocated_table_slots: u64,
    last_sync_scene_store_orphaned_table_slots: u64,
    sphere_lod_active: bool,
    sphere_lod_sample_shift: u32,
    sphere_lod_sample_stride: u32,
    sphere_lod_base_sample_shift: u32,
    sphere_lod_source_atom_count: u64,
    sphere_lod_instance_upper_bound: u64,
    sphere_lod_cull_upper_bound: u64,
    sphere_lod_viewport_visible_count: u64,
    sphere_lod_viewport_full_detail: bool,
    stick_lod_active: bool,
    stick_lod_sample_shift: u32,
    stick_lod_sample_stride: u32,
    stick_lod_base_sample_shift: u32,
    stick_lod_source_bond_count: u64,
    stick_lod_sampled_bond_upper_bound: u64,
    stick_lod_cull_upper_bound: u64,
    stick_lod_viewport_visible_count: u64,
    stick_lod_viewport_full_detail: bool,
    overlay_id_marked_only: bool,
    last_settings_ms: f32,
    last_acquire_ms: f32,
    last_encode_ms: f32,
    last_submit_present_ms: f32,
    hover_submitted: u64,
    hover_completed: u64,
    hover_stale: u64,
    hover_queued: u64,
    hover_deferred: u64,
    hover_throttle_active: bool,
    hover_cancelled: u64,
    click_submitted: u64,
    click_completed: u64,
    hover_pending: bool,
    click_pending: bool,
}

/// True when the hover state should be re-applied (the hit identity changed).
fn hover_changed(prev: &Option<PickHit>, next: &Option<PickHit>) -> bool {
    match (prev, next) {
        (None, None) => false,
        (Some(_), None) | (None, Some(_)) => true,
        (Some(a), Some(b)) => a.object_name != b.object_name || a.atom_index != b.atom_index,
    }
}

fn performance_now_ms() -> f64 {
    web_sys::window()
        .and_then(|window| window.performance())
        .map(|performance| performance.now())
        .unwrap_or(0.0)
}

// ---------------------------------------------------------------------------
// WebViewer
// ---------------------------------------------------------------------------

/// The main web viewer — owns scene state, command executor, and GPU resources.
#[wasm_bindgen]
pub struct WebViewer {
    session: Session,
    executor: CommandExecutor,
    gpu: Option<GpuState>,
    render_scene: CachedRenderScene,
    input: InputState,
    needs_redraw: bool,
    width: u32,
    height: u32,
    picking_enabled: bool,
    selection_overlay_enabled: bool,
    hover_hit: Option<PickHit>,
    /// In-flight GPU hover pick. New hover events while this is `Some` replace
    /// `queued_hover`, so the GPU readback pipeline never has more than one
    /// hover request in flight.
    pending_hover: Option<PendingPick>,
    /// Monotonic token for hover cancellation. Pending readbacks with older
    /// epochs are discarded after cursor leave/out-of-bounds clears hover.
    hover_epoch: u64,
    pending_hover_epoch: Option<u64>,
    /// Latest hover position received while a readback is in flight.
    queued_hover: Option<(u32, u32)>,
    /// In-flight GPU click pick. Same drop-while-busy rule.
    pending_click: Option<PendingPick>,
    /// Click hit collected by `poll_pending_picks` but not yet pulled by JS.
    /// Drained by `take_completed_pick`.
    last_completed_click: Option<Option<PickHitInfo>>,
    perf_history: FrameStatsHistory,
    last_render_ms: f32,
    last_poll_picks_ms: f32,
    last_render_timings: render_loop::WebRenderTimings,
    render_count: u64,
    hover_submitted: u64,
    hover_completed: u64,
    hover_stale: u64,
    hover_queued: u64,
    hover_deferred: u64,
    hover_cancelled: u64,
    click_submitted: u64,
    click_completed: u64,
    last_hover_submit_ms: f64,
}

#[wasm_bindgen]
impl WebViewer {
    /// Create a new WebViewer bound to a `<canvas>` element.
    ///
    /// `picking_enabled = false` (default) drops hit-test picking readbacks
    /// and the half-res picking target. Selection overlay is controlled
    /// separately; silhouettes remain command/settings-driven.
    #[wasm_bindgen]
    pub async fn create(
        canvas_id: &str,
        picking_enabled: bool,
        selection_overlay_enabled: bool,
    ) -> Result<WebViewer, JsValue> {
        let render_config = RenderConfig {
            picking: if picking_enabled {
                PickingMode::Reprojected
            } else {
                PickingMode::Disabled
            },
            selection_overlay: selection_overlay_enabled,
        };
        let gpu = GpuState::from_canvas(canvas_id, render_config)
            .await
            .map_err(|e| JsValue::from_str(&e))?;

        let width = gpu.surface_config.width;
        let height = gpu.surface_config.height;

        let mut session = Session::new();
        session.camera.set_aspect(width as f32 / height as f32);

        Ok(WebViewer {
            session,
            executor: CommandExecutor::new(),
            gpu: Some(gpu),
            render_scene: CachedRenderScene::default(),
            input: InputState::new(),
            needs_redraw: true,
            width,
            height,
            picking_enabled,
            selection_overlay_enabled,
            hover_hit: None,
            pending_hover: None,
            hover_epoch: 0,
            pending_hover_epoch: None,
            queued_hover: None,
            pending_click: None,
            last_completed_click: None,
            perf_history: FrameStatsHistory::with_default_capacity(),
            last_render_ms: 0.0,
            last_poll_picks_ms: 0.0,
            last_render_timings: render_loop::WebRenderTimings::default(),
            render_count: 0,
            hover_submitted: 0,
            hover_completed: 0,
            hover_stale: 0,
            hover_queued: 0,
            hover_deferred: 0,
            hover_cancelled: 0,
            click_submitted: 0,
            click_completed: 0,
            last_hover_submit_ms: f64::NEG_INFINITY,
        })
    }

    // =======================================================================
    // Rendering
    // =======================================================================

    /// Render one frame to the canvas.
    #[wasm_bindgen]
    pub fn render_frame(&mut self) {
        let t0 = performance_now_ms();
        // Drain any GPU pick readbacks that completed since the last frame.
        let poll_t0 = performance_now_ms();
        self.poll_pending_picks();
        self.last_poll_picks_ms = (performance_now_ms() - poll_t0).max(0.0) as f32;

        let result = match &mut self.gpu {
            Some(gpu) => render_loop::render_frame(gpu, &mut self.session, &mut self.render_scene),
            None => return,
        };

        match result {
            Ok(timings) => {
                self.last_render_timings = timings;
                self.needs_redraw = false;
            }
            Err(e) => {
                log::warn!("Render error: {:?}", e);
                match e {
                    wgpu::SurfaceError::Lost | wgpu::SurfaceError::Outdated => {
                        if let Some(gpu) = &mut self.gpu {
                            gpu.resize(self.width, self.height);
                        }
                        self.needs_redraw = true;
                    }
                    wgpu::SurfaceError::Timeout | wgpu::SurfaceError::Other => {
                        self.needs_redraw = true;
                    }
                    wgpu::SurfaceError::OutOfMemory => {
                        self.needs_redraw = false;
                    }
                }
            }
        }
        let render_ms = (performance_now_ms() - t0).max(0.0) as f32;
        self.last_render_ms = render_ms;
        self.perf_history.push(render_ms);
        self.render_count = self.render_count.saturating_add(1);
    }

    /// Returns true when the scene has changed and needs a re-render.
    #[wasm_bindgen]
    pub fn needs_redraw(&self) -> bool {
        self.needs_redraw
    }

    /// Advance movie playback, rock animation, and camera interpolation.
    #[wasm_bindgen]
    pub fn update_animations(&mut self, dt: f32) {
        if self.session.update_animations(dt).needs_redraw {
            self.needs_redraw = true;
        }
    }

    /// Handle canvas resize.
    #[wasm_bindgen]
    pub fn resize(&mut self, width: u32, height: u32) {
        self.width = width.max(1);
        self.height = height.max(1);
        if let Some(gpu) = &mut self.gpu {
            gpu.resize(self.width, self.height);
        }
        self.session
            .camera
            .set_aspect(self.width as f32 / self.height as f32);
        self.needs_redraw = true;
    }

    // =======================================================================
    // Input events
    // =======================================================================

    #[wasm_bindgen]
    pub fn on_mouse_down(&mut self, x: f32, y: f32, button: u32, modifiers: u32) {
        event::handle_mouse_down(&mut self.input, x, y, button, modifiers);
        self.needs_redraw = true;
    }

    #[wasm_bindgen]
    pub fn on_mouse_move(&mut self, x: f32, y: f32, modifiers: u32) {
        event::handle_mouse_move(&mut self.input, x, y, modifiers);
        if self.input.any_button_pressed() {
            self.needs_redraw = true;
        }
    }

    #[wasm_bindgen]
    pub fn on_mouse_up(&mut self, x: f32, y: f32, button: u32) {
        event::handle_mouse_up(&mut self.input, x, y, button);
        self.needs_redraw = true;
    }

    #[wasm_bindgen]
    pub fn on_wheel(&mut self, delta_y: f32, modifiers: u32) {
        event::handle_wheel(&mut self.input, delta_y, modifiers);
        self.needs_redraw = true;
    }

    // =======================================================================
    // Picking
    // =======================================================================

    /// Enable or disable cursor-based atom picking (default: disabled).
    #[wasm_bindgen]
    pub fn set_picking_enabled(&mut self, enabled: bool) {
        self.picking_enabled = enabled;
    }

    /// Enable or disable the visible selection / hover overlay. Hit-test
    /// picking remains controlled by `set_picking_enabled`.
    #[wasm_bindgen]
    pub fn set_selection_overlay_enabled(&mut self, enabled: bool) {
        self.selection_overlay_enabled = enabled;
        if let Some(gpu) = &mut self.gpu {
            gpu.state.set_selection_overlay_enabled(enabled);
        }
        self.needs_redraw = true;
    }

    /// Update hover indicators by submitting a GPU pick at physical-pixel
    /// coordinates.
    #[wasm_bindgen]
    pub fn process_hover(&mut self, screen_x: f32, screen_y: f32) {
        if !self.picking_enabled || self.input.any_button_pressed() {
            return;
        }
        if screen_x < 0.0 || screen_y < 0.0 {
            self.cancel_hover_request();
            return;
        }
        let x = screen_x as u32;
        let y = screen_y as u32;
        if x >= self.width || y >= self.height {
            self.cancel_hover_request();
            return;
        }
        if self.gpu.is_none() {
            return;
        }

        self.queue_or_submit_hover(x, y);
    }

    fn cancel_hover_request(&mut self) {
        self.hover_epoch = self.hover_epoch.wrapping_add(1);
        self.hover_cancelled = self.hover_cancelled.saturating_add(1);
        self.queued_hover = None;
        self.apply_hover_hit(None);
    }

    /// Apply the hover-target / sele-overlap effects of a hover pick result.
    fn apply_hover_hit(&mut self, new_hit: Option<PickHit>) {
        let mode = self.session.settings.ui.mouse_selection_mode as i32;

        if let Some(hit) = new_hit.as_ref() {
            if let Some(mol_obj) = self.session.registry.get_molecule(&hit.object_name) {
                let mol = mol_obj.molecule();
                let sel = expand_pick_to_selection(hit, mode, mol);
                self.session.set_hover(patinae_scene::HoverTarget {
                    object: hit.object_name.clone(),
                    selection: sel,
                });
            } else {
                self.session.clear_hover();
            }
        } else {
            self.session.clear_hover();
        }

        self.hover_hit = new_hit;
        self.needs_redraw = true;
    }

    /// Submit a GPU click pick at physical-pixel canvas coordinates.
    /// Returns `null` immediately — the actual hit lands asynchronously
    /// via `take_completed_pick()`.
    #[wasm_bindgen]
    pub fn pick_at_screen(&mut self, screen_x: f32, screen_y: f32) -> JsValue {
        if !self.picking_enabled {
            return JsValue::NULL;
        }
        if self.pending_click.is_none() {
            let x = screen_x.max(0.0) as u32;
            let y = screen_y.max(0.0) as u32;
            if x < self.width && y < self.height {
                self.submit_gpu_pick(x, y, /*for_click=*/ true);
            }
        }
        JsValue::NULL
    }

    /// Drain the most recent GPU click pick result.
    #[wasm_bindgen]
    pub fn take_completed_pick(&mut self) -> JsValue {
        match self.last_completed_click.take() {
            Some(Some(info)) => serde_wasm_bindgen::to_value(&info).unwrap_or(JsValue::NULL),
            Some(None) => JsValue::NULL,
            None => JsValue::UNDEFINED,
        }
    }

    /// Return debug performance counters for browser-side perf harnesses.
    #[wasm_bindgen]
    pub fn get_performance_snapshot(&self) -> JsValue {
        let sphere_lod = self
            .gpu
            .as_ref()
            .map(|gpu| gpu.state.sphere_lod_diagnostics())
            .unwrap_or_default();
        let stick_lod = self
            .gpu
            .as_ref()
            .map(|gpu| gpu.state.stick_lod_diagnostics())
            .unwrap_or_default();
        let scene_fragmentation = self
            .last_render_timings
            .sync_detail
            .scene_store_fragmentation;
        let snapshot = WebPerformanceSnapshot {
            render_count: self.render_count,
            avg_render_ms: self.perf_history.avg_ms(),
            median_render_ms: self.perf_history.percentile_ms(0.5),
            p95_render_ms: self.perf_history.percentile_ms(0.95),
            last_render_ms: self.last_render_ms,
            last_poll_picks_ms: self.last_poll_picks_ms,
            last_uniforms_ms: self.last_render_timings.uniforms_ms,
            last_prepare_ms: self.last_render_timings.prepare_ms,
            last_sync_ms: self.last_render_timings.sync_ms,
            last_sync_scene_store_object_ms: self
                .last_render_timings
                .sync_detail
                .scene_store_object_sync_ms,
            last_sync_scene_store_flush_ms: self
                .last_render_timings
                .sync_detail
                .scene_store_flush_ms,
            last_sync_marking_resources_ms: self
                .last_render_timings
                .sync_detail
                .marking_resources_ms,
            last_sync_rep_ms: self.last_render_timings.sync_detail.rep_sync_ms,
            last_sync_map_ms: self.last_render_timings.sync_detail.map_sync_ms,
            last_sync_order_bounds_ms: self.last_render_timings.sync_detail.order_bounds_ms,
            last_sync_compute_dispatch_ms: self.last_render_timings.sync_detail.compute_dispatch_ms,
            last_sync_marker_lut_upload_bytes: self
                .last_render_timings
                .sync_detail
                .marker_lut_upload_bytes,
            last_sync_marker_lut_upload_ranges: self
                .last_render_timings
                .sync_detail
                .marker_lut_upload_ranges,
            last_sync_marker_lut_reallocated: self
                .last_render_timings
                .sync_detail
                .marker_lut_reallocated,
            last_sync_scene_store_live_atoms: scene_fragmentation.live_atoms,
            last_sync_scene_store_allocated_atoms: scene_fragmentation.allocated_atoms,
            last_sync_scene_store_orphaned_atoms: scene_fragmentation.orphaned_atoms,
            last_sync_scene_store_live_bonds: scene_fragmentation.live_bonds,
            last_sync_scene_store_allocated_bonds: scene_fragmentation.allocated_bonds,
            last_sync_scene_store_orphaned_bonds: scene_fragmentation.orphaned_bonds,
            last_sync_scene_store_live_table_slots: scene_fragmentation.live_table_slots,
            last_sync_scene_store_allocated_table_slots: scene_fragmentation.allocated_table_slots,
            last_sync_scene_store_orphaned_table_slots: scene_fragmentation.orphaned_table_slots,
            sphere_lod_active: sphere_lod.active,
            sphere_lod_sample_shift: sphere_lod.sample_shift,
            sphere_lod_sample_stride: sphere_lod.sample_stride,
            sphere_lod_base_sample_shift: sphere_lod.base_sample_shift,
            sphere_lod_source_atom_count: sphere_lod.source_atom_count,
            sphere_lod_instance_upper_bound: sphere_lod.instance_upper_bound,
            sphere_lod_cull_upper_bound: sphere_lod.cull_upper_bound,
            sphere_lod_viewport_visible_count: sphere_lod.viewport_visible_count,
            sphere_lod_viewport_full_detail: sphere_lod.viewport_full_detail,
            stick_lod_active: stick_lod.active,
            stick_lod_sample_shift: stick_lod.sample_shift,
            stick_lod_sample_stride: stick_lod.sample_stride,
            stick_lod_base_sample_shift: stick_lod.base_sample_shift,
            stick_lod_source_bond_count: stick_lod.source_bond_count,
            stick_lod_sampled_bond_upper_bound: stick_lod.sampled_bond_upper_bound,
            stick_lod_cull_upper_bound: stick_lod.cull_upper_bound,
            stick_lod_viewport_visible_count: stick_lod.viewport_visible_count,
            stick_lod_viewport_full_detail: stick_lod.viewport_full_detail,
            overlay_id_marked_only: false,
            last_settings_ms: self.last_render_timings.settings_ms,
            last_acquire_ms: self.last_render_timings.acquire_ms,
            last_encode_ms: self.last_render_timings.encode_ms,
            last_submit_present_ms: self.last_render_timings.submit_present_ms,
            hover_submitted: self.hover_submitted,
            hover_completed: self.hover_completed,
            hover_stale: self.hover_stale,
            hover_queued: self.hover_queued,
            hover_deferred: self.hover_deferred,
            hover_throttle_active: self.hover_throttle_active(),
            hover_cancelled: self.hover_cancelled,
            click_submitted: self.click_submitted,
            click_completed: self.click_completed,
            hover_pending: self.pending_hover.is_some(),
            click_pending: self.pending_click.is_some(),
        };
        serde_wasm_bindgen::to_value(&snapshot).unwrap_or(JsValue::NULL)
    }

    /// Clear debug performance counters for the next harness scenario.
    #[wasm_bindgen]
    pub fn reset_performance_stats(&mut self) {
        self.perf_history.clear();
        self.last_render_ms = 0.0;
        self.last_poll_picks_ms = 0.0;
        self.last_render_timings = render_loop::WebRenderTimings::default();
        self.render_count = 0;
        self.hover_submitted = 0;
        self.hover_completed = 0;
        self.hover_stale = 0;
        self.hover_queued = 0;
        self.hover_deferred = 0;
        self.hover_cancelled = 0;
        self.click_submitted = 0;
        self.click_completed = 0;
    }

    /// Run the same session-side effects that the CPU `pick_at_screen` runs
    /// for a click hit (selection command + indicator hit info).
    fn apply_click_hit(&mut self, hit: Option<PickHit>) -> Option<PickHitInfo> {
        let (cmd, hit_info) = if let Some(ref hit) = hit {
            if let Some(mol_obj) = self.session.registry.get_molecule(&hit.object_name) {
                let mol = mol_obj.molecule();
                let mode = self.session.settings.ui.mouse_selection_mode as i32;

                let cmd = pick_expression_for_hit(hit, mode, mol).and_then(|expr| {
                    let overlaps_sele = self.session.selections.get("sele").is_some_and(|entry| {
                        let sel = expand_pick_to_selection(hit, mode, mol);
                        select(mol, &entry.expression)
                            .map(|sele| sel.intersection(&sele).any())
                            .unwrap_or(false)
                    });
                    let has_sele = self.session.selections.contains("sele");
                    build_sele_command(&expr, overlaps_sele, has_sele)
                });

                let info = PickHitInfo::from_hit(hit, mode, mol);
                (cmd, Some(info))
            } else {
                (None, None)
            }
        } else {
            let cmd = if self.session.selections.contains("sele") {
                Some("deselect sele".to_string())
            } else {
                None
            };
            (cmd, None)
        };

        if let Some(ref cmd) = cmd {
            // Web viewer doesn't currently implement `CaptureRenderer`, so
            // commands that need a render context (`png`, `movie render`)
            // will fail silently. Selection / view / colour commands work
            // without one. A web-side `CaptureRenderer` impl is future work.
            let mut adapter = SessionAdapter {
                session: &mut self.session,
                render_context: None,
                default_size: (self.width, self.height),
                needs_redraw: &mut self.needs_redraw,
                async_fetch_fn: None,
            };
            let _ = self.executor.do_with_options(&mut adapter, cmd, true);
        }

        self.needs_redraw = true;
        hit_info
    }

    /// Submit a GPU pick pass via patinae-render's async API. The resulting
    /// `PendingPick` is stored in either `pending_hover` or `pending_click`;
    /// `poll_pending_picks` drains it once the GPU readback lands.
    fn submit_gpu_pick(&mut self, x: u32, y: u32, for_click: bool) {
        let gpu = match self.gpu.as_mut() {
            Some(g) => g,
            None => return,
        };
        let target = if for_click {
            PickReadbackTarget::Click
        } else {
            PickReadbackTarget::Hover
        };
        let Some(pending) = gpu.state.submit_pick_with_target(x, y, target) else {
            return;
        };
        if for_click {
            self.click_submitted = self.click_submitted.saturating_add(1);
            self.pending_click = Some(pending);
        } else {
            self.last_hover_submit_ms = performance_now_ms();
            self.hover_submitted = self.hover_submitted.saturating_add(1);
            self.pending_hover_epoch = Some(self.hover_epoch);
            self.pending_hover = Some(pending);
        }
    }

    fn hover_throttle_active(&self) -> bool {
        self.gpu
            .as_ref()
            .is_some_and(|gpu| {
                gpu.state.sphere_lod_diagnostics().active
                    || gpu.state.stick_lod_diagnostics().active
            })
    }

    fn hover_pick_due(&self, now_ms: f64) -> bool {
        !self.hover_throttle_active()
            || now_ms - self.last_hover_submit_ms >= LARGE_LOD_HOVER_PICK_INTERVAL_MS
    }

    fn queue_or_submit_hover(&mut self, x: u32, y: u32) {
        if self.pending_hover.is_some() {
            self.hover_queued = self.hover_queued.saturating_add(1);
            self.queued_hover = Some((x, y));
            return;
        }

        if !self.hover_pick_due(performance_now_ms()) {
            self.hover_deferred = self.hover_deferred.saturating_add(1);
            self.queued_hover = Some((x, y));
            return;
        }

        self.queued_hover = None;
        self.submit_gpu_pick(x, y, /*for_click=*/ false);
    }

    fn submit_queued_hover_if_due(&mut self) {
        if self.pending_hover.is_some() {
            return;
        }
        let Some((x, y)) = self.queued_hover else {
            return;
        };
        if !self.hover_pick_due(performance_now_ms()) {
            return;
        }
        self.queued_hover = None;
        self.submit_gpu_pick(x, y, /*for_click=*/ false);
    }

    /// Try to drain any in-flight GPU picks. JS calls this every rAF so
    /// readbacks complete even when no visible redraw is pending.
    #[wasm_bindgen]
    pub fn poll_pending_picks(&mut self) {
        // Hover.
        let hover_raw = match (self.pending_hover.as_ref(), self.gpu.as_ref()) {
            (Some(pending), Some(gpu)) => gpu.state.try_collect_pick(pending),
            _ => None,
        };
        if let Some(raw) = hover_raw {
            let stale = self.pending_hover_epoch != Some(self.hover_epoch);
            self.pending_hover = None;
            self.pending_hover_epoch = None;

            if !stale {
                self.hover_completed = self.hover_completed.saturating_add(1);
                let hit = raw
                    .and_then(|h| resolve_pick(h, self.render_scene.object_names(), &self.session));
                if hover_changed(&self.hover_hit, &hit) {
                    self.apply_hover_hit(hit);
                }
            } else {
                self.hover_stale = self.hover_stale.saturating_add(1);
            }

            self.submit_queued_hover_if_due();
        }
        self.submit_queued_hover_if_due();

        // Click.
        if let Some(pending) = self.pending_click.as_ref() {
            if let Some(gpu) = self.gpu.as_ref() {
                if let Some(raw) = gpu.state.try_collect_pick(pending) {
                    self.pending_click = None;
                    let hit = raw.and_then(|h| {
                        resolve_pick(h, self.render_scene.object_names(), &self.session)
                    });
                    let info = self.apply_click_hit(hit);
                    self.last_completed_click = Some(info);
                    self.click_completed = self.click_completed.saturating_add(1);
                }
            }
        }
    }

    /// Process accumulated input deltas and update the camera.
    #[wasm_bindgen]
    pub fn process_input(&mut self) {
        let deltas = self.input.take_camera_deltas();
        let screen_vertex_scale = self.session.camera.screen_vertex_scale(self.height as f32);
        for delta in deltas {
            match delta {
                CameraDelta::Rotate { x, y } => {
                    self.session.camera.rotate_x(x);
                    self.session.camera.rotate_y(y);
                }
                CameraDelta::Translate(v) => {
                    let scaled = lin_alg::f32::Vec3::new(
                        v.x * screen_vertex_scale * self.input.pan_sensitivity,
                        v.y * screen_vertex_scale * self.input.pan_sensitivity,
                        v.z * screen_vertex_scale * self.input.pan_sensitivity,
                    );
                    self.session.camera.translate(scaled);
                }
                CameraDelta::Zoom(z) => {
                    self.session.camera.zoom(z);
                }
                CameraDelta::Clip { front, back } => {
                    let view = self.session.camera.view_mut();
                    view.clip_front = (view.clip_front + front).max(0.01);
                    view.clip_back = (view.clip_back + back).max(view.clip_front + 0.01);
                }
                CameraDelta::SlabScale(raw_delta) => {
                    let mws = self.session.settings.ui.mouse_wheel_scale;
                    let scale = (1.0 + 0.04 * mws * raw_delta).clamp(0.5, 2.0);

                    let view = self.session.camera.view_mut();
                    let distance = view.position.z;
                    let half = ((view.clip_back - view.clip_front) * 0.5).max(0.1);
                    let new_half = (half * scale).max(2.0);

                    view.clip_front = (distance - new_half).max(0.01);
                    view.clip_back = (distance + new_half).max(view.clip_front + 0.1);
                }
            }
            self.needs_redraw = true;
        }
    }

    // =======================================================================
    // Commands
    // =======================================================================

    /// Execute a command string. Returns JSON with output messages.
    #[wasm_bindgen]
    pub fn execute(&mut self, command: &str) -> JsValue {
        // Web doesn't implement `CaptureRenderer` yet — commands that
        // need GPU access (png / movie render) will fail. Future work.
        let mut adapter = SessionAdapter {
            session: &mut self.session,
            render_context: None,
            default_size: (self.width, self.height),
            needs_redraw: &mut self.needs_redraw,
            async_fetch_fn: None,
        };

        let result = self.executor.do_with_options(&mut adapter, command, false);

        let messages = match result {
            Ok(output) => {
                let mut messages: Vec<OutputMsg> = output
                    .messages
                    .iter()
                    .map(|m| OutputMsg {
                        level: match m.kind {
                            MessageKind::Info => "info",
                            MessageKind::Warning => "warning",
                            MessageKind::Error => "error",
                        },
                        text: m.text.clone(),
                    })
                    .collect();

                for action in &output.actions {
                    if matches!(action, CommandAction::ClearOutput) {
                        messages.push(OutputMsg {
                            level: "clear",
                            text: String::new(),
                        });
                    }
                }

                messages
            }
            Err(e) => vec![OutputMsg {
                level: "error",
                text: e.to_string(),
            }],
        };

        let output = CmdOutput { messages };
        serde_wasm_bindgen::to_value(&output).unwrap_or(JsValue::NULL)
    }

    /// Load molecular or map data from bytes.
    #[wasm_bindgen]
    pub fn load_data(&mut self, data: &[u8], name: &str, format: &str) -> Result<(), JsValue> {
        let decompressed;
        let data = if data.len() >= 2 && data[0] == 0x1f && data[1] == 0x8b {
            use std::io::Read;
            let mut decoder = patinae_io::compress::gzip_reader(data);
            let mut buf = Vec::new();
            decoder
                .read_to_end(&mut buf)
                .map_err(|e| JsValue::from_str(&format!("Gzip decompression failed: {}", e)))?;
            decompressed = buf;
            decompressed.as_slice()
        } else {
            data
        };

        let fmt = format.to_lowercase();

        if fmt == "prs" {
            let session: Session = rmp_serde::from_slice(data)
                .map_err(|e| JsValue::from_str(&format!("PRS parse error: {}", e)))?;
            self.session = session;
            self.session.registry.mark_all_dirty();
            self.needs_redraw = true;
            return Ok(());
        }

        if matches!(fmt.as_str(), "ccp4" | "map" | "mrc") {
            let ccp4 = patinae_io::ccp4::read_ccp4_from(std::io::Cursor::new(data))
                .map_err(|e| JsValue::from_str(&format!("CCP4 parse error: {}", e)))?;
            let grid = patinae_algos::surface::Grid3D::from_dims(
                ccp4.origin,
                ccp4.spacing,
                ccp4.dims,
                ccp4.values,
            );
            let map_data = patinae_scene::MapData::new(grid);
            let map_obj = patinae_scene::MapObject::from_map_data(name, map_data);
            self.session.registry.add(map_obj);
            if let Some((min, max)) = self.session.registry.extent() {
                self.session.camera.zoom_to(min, max, 0.0);
            }
            self.needs_redraw = true;
            return Ok(());
        }

        let mol = match fmt.as_str() {
            "bcif" => patinae_io::bcif::read_bcif_bytes(data),
            _ => {
                let data_str = std::str::from_utf8(data)
                    .map_err(|_| JsValue::from_str("Data is not valid UTF-8"))?;
                match fmt.as_str() {
                    "pdb" => patinae_io::pdb::read_pdb_str(data_str),
                    "xyz" => patinae_io::xyz::read_xyz_str(data_str),
                    "cif" | "mmcif" => patinae_io::cif::read_cif_str(data_str),
                    _ => {
                        return Err(JsValue::from_str(&format!(
                            "Direct loading not yet supported for: {}. Use execute() instead.",
                            fmt
                        )))
                    }
                }
            }
        };

        match mol {
            Ok(mut molecule) => {
                let behavior = &self.session.settings.behavior;
                if behavior.auto_dss {
                    let assigner = patinae_mol::dss::assigner_for(behavior.dss_algorithm);
                    patinae_mol::dss::assign_secondary_structure(
                        &mut molecule,
                        0,
                        assigner.as_ref(),
                    );
                }

                let mol_obj = MoleculeObject::with_name(molecule, name);
                self.session.registry.add(mol_obj);
                self.session.refresh_movie_state_count();
                if let Some((min, max)) = self.session.registry.extent() {
                    self.session.camera.zoom_to(min, max, 0.0);
                }
                self.needs_redraw = true;
                Ok(())
            }
            Err(e) => Err(JsValue::from_str(&format!("Parse error: {}", e))),
        }
    }

    // =======================================================================
    // State queries
    // =======================================================================

    #[wasm_bindgen]
    pub fn get_object_names(&self) -> JsValue {
        let names: Vec<String> = self
            .session
            .registry
            .names()
            .map(|s| s.to_string())
            .collect();
        serde_wasm_bindgen::to_value(&names).unwrap_or(JsValue::NULL)
    }

    #[wasm_bindgen]
    pub fn get_object_info(&self, name: &str) -> JsValue {
        if let Some(mol_obj) = self.session.registry.get_molecule(name) {
            let info = ObjectInfo {
                name: name.to_string(),
                object_type: "molecule",
                atom_count: mol_obj.molecule().atom_count(),
                enabled: mol_obj.is_enabled(),
            };
            return serde_wasm_bindgen::to_value(&info).unwrap_or(JsValue::NULL);
        }
        if let Some(map_obj) = self.session.registry.get_map(name) {
            let info = ObjectInfo {
                name: name.to_string(),
                object_type: "map",
                atom_count: 0,
                enabled: map_obj.is_enabled(),
            };
            return serde_wasm_bindgen::to_value(&info).unwrap_or(JsValue::NULL);
        }
        JsValue::NULL
    }

    #[wasm_bindgen]
    pub fn get_sequence_data(&self) -> JsValue {
        let mut chains: Vec<SequenceChain> = Vec::new();

        for name in self.session.registry.names() {
            if let Some(mol_obj) = self.session.registry.get_molecule(name) {
                let mol = mol_obj.molecule();
                let mut chain_map: std::collections::BTreeMap<String, Vec<SequenceResidue>> =
                    std::collections::BTreeMap::new();

                let mut last_resv: Option<(String, i32)> = None;
                for atom in mol.atoms() {
                    let chain = atom.residue.key.chain.clone();
                    let resv = atom.residue.key.resv;
                    if last_resv.as_ref() == Some(&(chain.clone(), resv)) {
                        continue;
                    }
                    last_resv = Some((chain.clone(), resv));

                    if patinae_mol::is_water(&atom.residue.resn)
                        || patinae_mol::is_ion(&atom.residue.resn)
                    {
                        continue;
                    }

                    let one_letter = patinae_mol::residue_to_char(&atom.residue.resn).to_string();

                    chain_map
                        .entry(chain.clone())
                        .or_default()
                        .push(SequenceResidue {
                            resn: atom.residue.resn.clone(),
                            resv,
                            one_letter,
                        });
                }

                for (chain_id, residues) in chain_map {
                    chains.push(SequenceChain {
                        object_name: name.to_string(),
                        chain_id,
                        residues,
                    });
                }
            }
        }

        serde_wasm_bindgen::to_value(&chains).unwrap_or(JsValue::NULL)
    }

    #[wasm_bindgen]
    pub fn get_movie_state(&self) -> JsValue {
        serde_wasm_bindgen::to_value(&self.session.movie_state_snapshot()).unwrap_or(JsValue::NULL)
    }

    #[wasm_bindgen]
    pub fn get_selection_list(&self) -> JsValue {
        let mut list: Vec<SelectionInfo> = self
            .session
            .selections
            .iter()
            .map(|(name, entry)| SelectionInfo {
                name: name.clone(),
                expression: entry.expression.clone(),
                visible: entry.visible,
            })
            .collect();
        list.sort_by(|a, b| a.name.cmp(&b.name));
        serde_wasm_bindgen::to_value(&list).unwrap_or(JsValue::NULL)
    }

    #[wasm_bindgen]
    pub fn count_atoms(&self, selection: &str) -> Result<usize, JsValue> {
        let mut total = 0;
        for name in self.session.registry.names() {
            if let Some(mol_obj) = self.session.registry.get_molecule(name) {
                let mol = mol_obj.molecule();
                match select(mol, selection) {
                    Ok(mask) => total += mask.count(),
                    Err(e) => return Err(JsValue::from_str(&format!("Selection error: {}", e))),
                }
            }
        }
        Ok(total)
    }

    #[wasm_bindgen]
    pub fn get_labels(&self) -> JsValue {
        let viewport = (0.0, 0.0, self.width as f32, self.height as f32);
        let mut labels = Vec::new();

        for name in self.session.registry.names() {
            if let Some(mol_obj) = self.session.registry.get_molecule(name) {
                if !mol_obj.is_enabled() {
                    continue;
                }
                for (pos, text) in mol_obj.collect_labels() {
                    if let Some((sx, sy)) = self.session.camera.project_to_screen(pos, viewport) {
                        if sx >= 0.0
                            && sx <= self.width as f32
                            && sy >= 0.0
                            && sy <= self.height as f32
                        {
                            labels.push(LabelInfo {
                                x: sx,
                                y: sy,
                                text: text.to_string(),
                                kind: "atom",
                            });
                        }
                    }
                }
            }
        }

        for name in self.session.registry.names() {
            if let Some(meas_obj) = self.session.registry.get_measurement(name) {
                if !meas_obj.is_enabled() {
                    continue;
                }
                for (pos, text) in meas_obj.collect_labels() {
                    if let Some((sx, sy)) = self.session.camera.project_to_screen(pos, viewport) {
                        if sx >= 0.0
                            && sx <= self.width as f32
                            && sy >= 0.0
                            && sy <= self.height as f32
                        {
                            labels.push(LabelInfo {
                                x: sx,
                                y: sy,
                                text: text.to_string(),
                                kind: "measurement",
                            });
                        }
                    }
                }
            }
        }

        serde_wasm_bindgen::to_value(&labels).unwrap_or(JsValue::NULL)
    }
}
