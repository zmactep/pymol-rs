//! egui UI Construction
//!
//! Main UI loop, message dispatch, and selection synchronization.
//! Uses the Layout engine to render components into configurable panels.

use pymol_scene::Object;
use pymol_select::SelectionResult;

use pymol_framework::component::SharedContext;
use pymol_framework::message::AppMessage;
use pymol_framework::model::ResidueRef;
use crate::ui::{DragDropOverlay, NotificationOverlay};
use crate::layout::render_layout;

use super::App;

impl App {
    /// Run the egui UI and return the output for rendering later
    pub(crate) fn run_egui_ui(&mut self) -> Option<(Vec<egui::ClippedPrimitive>, egui::TexturesDelta)> {
        // Dispatch messages queued between frames (e.g. from resumed(), IPC, async tasks)
        self.dispatch_messages();

        // Start message bus frame
        self.bus.begin_frame();

        // Take egui input first
        let raw_input = {
            let (egui_state, window) = match (&mut self.view.egui.state, &self.view.window) {
                (Some(state), Some(window)) => (state, window),
                _ => return None,
            };
            egui_state.take_egui_input(window)
        };

        // Get scale factor for converting to physical pixels
        let scale_factor = self.view.scale_factor();

        let mut viewport_rect_logical = egui::Rect::NOTHING;

        // Check for pending async tasks and plugin notifications
        let mut pending_messages = self.task_runner.pending_messages();
        pending_messages.extend(self.plugin_manager.notification_messages().iter().cloned());

        // Build SharedContext for components
        let all_command_names: Vec<String> = self.executor.registry().all_names().map(|s| s.to_string()).collect();
        let setting_names = pymol_settings::setting_names();
        let mut setting_names_refs: Vec<&str> = setting_names.to_vec();
        setting_names_refs.extend(self.executor.dynamic_settings().names().iter().map(String::as_str));
        let cmd_registry = self.executor.registry();

        // Update image overlay texture if there's a viewport image
        let image_overlay_info = if let Some(vp_img) = &self.state.viewport_image {
            let needs_update = self.image_dirty || self.view.image_overlay.size != Some((vp_img.width, vp_img.height));
            if needs_update {
                self.view.update_image_overlay(&vp_img.data, vp_img.width, vp_img.height);
                self.image_dirty = false;
            }
            self.view.image_overlay.texture_id.map(|id| (id, vp_img.width, vp_img.height))
        } else {
            if self.view.image_overlay.texture_id.is_some() {
                self.view.clear_image_overlay();
            }
            None
        };

        let (gpu_device, gpu_queue) = self.view.gpu.render_context.as_ref()
            .map(|c| (c.device(), c.queue()))
            .unzip();

        let shared = SharedContext {
            registry: &self.state.registry,
            camera: &self.state.camera,
            selections: &self.state.selections,
            named_colors: &self.state.named_colors,
            movie: &self.state.movie,
            settings: &self.state.settings,
            clear_color: self.state.clear_color,
            gpu_device,
            gpu_queue,
            viewport_image: self.state.viewport_image.as_ref(),
            command_names: &all_command_names,
            command_registry: cmd_registry,
            setting_names: &setting_names_refs,
            dynamic_settings: Some(self.executor.dynamic_settings()),
            scene_generation: self.scene_generation,
        };

        // Run egui with the layout engine
        let full_output = {
            let layout = &self.layout;
            let components = &mut self.components;
            let bus = &mut self.bus;

            self.view.egui.ctx.run(raw_input, |ctx| {
                // Layout renders all panels + returns viewport rect
                let transparent_panels = self.state.settings.ui.transparent_panels;
                viewport_rect_logical = render_layout(layout, ctx, components, &shared, bus, transparent_panels);

                // Render image overlay in viewport if present
                if let Some((texture_id, img_width, img_height)) = image_overlay_info {
                    Self::render_image_overlay(ctx, viewport_rect_logical, texture_id, img_width, img_height);
                }

                // Paint atom labels and measurement labels as 2D overlay
                Self::render_labels_overlay(ctx, viewport_rect_logical, &shared);

                // Notification overlay for async tasks
                if !pending_messages.is_empty() {
                    NotificationOverlay::show(ctx, &pending_messages);
                }

                // Drag-and-drop file hint
                if let Some(ref path) = self.drag_hover_path {
                    DragDropOverlay::show(ctx, path);
                }

                // Fetch PDB dialog
                let rt_handle = self.task_runner.handle();
                self.fetch_dialog.show(ctx, &rt_handle);
            })
        };

        // Execute fetch dialog command if user clicked "Fetch"
        if let Some(cmd) = self.fetch_dialog.take_command() {
            let _ = self.execute_command(&cmd, false);
        }

        // Convert viewport rect from logical to physical pixels
        let viewport_rect_physical = egui::Rect::from_min_max(
            egui::pos2(
                viewport_rect_logical.min.x * scale_factor,
                viewport_rect_logical.min.y * scale_factor,
            ),
            egui::pos2(
                viewport_rect_logical.max.x * scale_factor,
                viewport_rect_logical.max.y * scale_factor,
            ),
        );
        self.view.viewport_rect = Some(viewport_rect_physical);

        // Handle platform output
        if let (Some(egui_state), Some(window)) = (&mut self.view.egui.state, &self.view.window) {
            egui_state.handle_platform_output(window, full_output.platform_output);
        }

        // Message Dispatch
        let has_messages = self.bus.has_pending();
        self.dispatch_messages();

        if has_messages {
            self.sync_selection_to_sequence();
        }

        // Sync sequence hover state to viewport (replaces HoverResidue message)
        if let Some(seq) = self.components.get::<crate::components::SequenceComponent>() {
            let new_hover = seq.ui_state.current_hover.clone();
            if new_hover != self.viewport.sequence_hover {
                self.viewport.sequence_hover = new_hover;
                self.scene_dirty = true;
            }
        }

        let egui_needs_repaint = full_output
            .viewport_output
            .values()
            .any(|v| v.repaint_delay.is_zero());

        if has_messages || egui_needs_repaint {
            self.view.request_redraw();
        }

        let clipped_primitives = self.view.egui.ctx.tessellate(full_output.shapes, full_output.pixels_per_point);

        Some((clipped_primitives, full_output.textures_delta))
    }

    /// Render the image overlay in the viewport area.
    fn render_image_overlay(
        ctx: &egui::Context,
        viewport_rect: egui::Rect,
        texture_id: egui::TextureId,
        img_width: u32,
        img_height: u32,
    ) {
        let available = viewport_rect.size();
        let img_aspect = img_width as f32 / img_height as f32;
        let vp_aspect = available.x / available.y;

        let (display_width, display_height) = if img_aspect > vp_aspect {
            (available.x, available.x / img_aspect)
        } else {
            (available.y * img_aspect, available.y)
        };

        let offset_x = (available.x - display_width) / 2.0;
        let offset_y = (available.y - display_height) / 2.0;

        let image_rect = egui::Rect::from_min_size(
            egui::pos2(viewport_rect.min.x + offset_x, viewport_rect.min.y + offset_y),
            egui::vec2(display_width, display_height),
        );

        let painter = ctx.layer_painter(egui::LayerId::new(
            egui::Order::Middle,
            egui::Id::new("image_overlay"),
        ));
        painter.image(
            texture_id,
            image_rect,
            egui::Rect::from_min_max(egui::pos2(0.0, 0.0), egui::pos2(1.0, 1.0)),
            egui::Color32::WHITE,
        );
    }

    /// Paint atom labels and measurement labels as 2D overlay on the viewport.
    fn render_labels_overlay(
        ctx: &egui::Context,
        viewport_rect: egui::Rect,
        shared: &SharedContext,
    ) {
        let vp = viewport_rect;
        let vp_tuple = (vp.min.x, vp.min.y, vp.width(), vp.height());
        let painter = ctx.layer_painter(egui::LayerId::new(
            egui::Order::Foreground,
            egui::Id::new("labels_overlay"),
        ));
        let font = egui::FontId::new(14.0, egui::FontFamily::Proportional);
        let label_color = egui::Color32::WHITE;

        for name in shared.registry.names() {
            let mol_obj = match shared.registry.get_molecule(name) {
                Some(m) if m.is_enabled() => m,
                _ => continue,
            };

            for (pos, text) in mol_obj.collect_labels() {
                if let Some((sx, sy)) = shared.camera.project_to_screen(pos, vp_tuple) {
                    if vp.contains(egui::pos2(sx, sy)) {
                        painter.text(
                            egui::pos2(sx, sy),
                            egui::Align2::LEFT_BOTTOM,
                            text,
                            font.clone(),
                            label_color,
                        );
                    }
                }
            }
        }

        let meas_font = egui::FontId::new(13.0, egui::FontFamily::Proportional);
        let meas_color = egui::Color32::from_rgb(255, 255, 0);
        for name in shared.registry.names() {
            if let Some(meas_obj) = shared.registry.get_measurement(name) {
                if !meas_obj.is_enabled() {
                    continue;
                }
                for (pos, text) in meas_obj.collect_labels() {
                    if let Some((sx, sy)) = shared.camera.project_to_screen(pos, vp_tuple) {
                        if vp.contains(egui::pos2(sx, sy)) {
                            painter.text(
                                egui::pos2(sx, sy),
                                egui::Align2::CENTER_CENTER,
                                text,
                                meas_font.clone(),
                                meas_color,
                            );
                        }
                    }
                }
            }
        }
    }

    /// Central message dispatcher.
    ///
    /// Loops until no new messages are generated, so that cascading messages
    /// (e.g., `ExecuteCommand` → `PrintInfo`) are processed in the same frame.
    fn dispatch_messages(&mut self) {
        for _ in 0..16 {
            let messages = self.bus.drain_outbox();
            if messages.is_empty() {
                break;
            }
            for msg in &messages {
                self.dispatch(msg);
            }
            for msg in &messages {
                self.components.broadcast(msg);
            }
            for msg in &messages {
                self.plugin_manager.broadcast(msg, &mut self.bus);
            }
        }
    }

    /// Dispatch a single message, mutating application state.
    ///
    /// Domain operations (objects, representations, coloring, camera,
    /// selections, movie) flow through `ExecuteCommand` — the command system
    /// handles them uniformly. The remaining arms cover GUI-level concerns.
    fn dispatch(&mut self, msg: &AppMessage) {
        match msg {
            // Command execution (covers all domain operations)
            AppMessage::ExecuteCommand { command, silent } => {
                let _ = self.execute_command(command, *silent);
            }

            // Layout
            AppMessage::TogglePanel(_)
            | AppMessage::ShowPanel(_)
            | AppMessage::HidePanel(_)
            | AppMessage::FloatPanel(_)
            | AppMessage::DockPanel(_)
            | AppMessage::ActivateTab(_) => self.dispatch_layout(msg),

            // Handled via broadcast to components (no app-level mutation)
            AppMessage::PrintInfo(_) | AppMessage::PrintWarning(_)
            | AppMessage::PrintError(_) | AppMessage::PrintCommand(_)
            | AppMessage::FocusPanel(_) => {}

            // Lifecycle
            AppMessage::RequestRedraw => self.scene_dirty = true,
            AppMessage::Quit => self.frame.quit_requested = true,
            AppMessage::ShowWindow => {
                self.show_window();
                self.headless = false;
                self.scene_dirty = true;
            }
            AppMessage::HideWindow => {
                self.hide_window();
                self.headless = true;
            }
            // Viewport image overlay
            AppMessage::SetViewportImage { data, width, height } => {
                self.set_viewport_image(Some(pymol_scene::ViewportImage {
                    data: data.clone(),
                    width: *width,
                    height: *height,
                }));
            }
            AppMessage::ClearViewportImage => {
                self.clear_viewport_image();
            }

            AppMessage::Custom { topic, payload } => {
                log::debug!("Custom event: {} ({} bytes)", topic, payload.len());
            }
        }
    }

    // =========================================================================
    // Grouped dispatch handlers
    // =========================================================================

    fn dispatch_layout(&mut self, msg: &AppMessage) {
        match msg {
            AppMessage::TogglePanel(id) => self.layout.toggle_expanded(id),
            AppMessage::ShowPanel(id) => self.layout.show_panel(id),
            AppMessage::HidePanel(id) => self.layout.hide_panel(id),
            AppMessage::FloatPanel(id) => self.layout.float_panel(id),
            AppMessage::DockPanel(id) => self.layout.dock_panel(id),
            AppMessage::ActivateTab(id) => self.layout.activate_tab(id),
            _ => unreachable!(),
        }
        self.scene_dirty = true;
    }

    /// Evaluate all visible selections for all molecules.
    pub(crate) fn evaluate_visible_selections(&mut self) -> Vec<(String, SelectionResult)> {
        self.state.selections.evaluate_visible(
            &self.state.registry,
            Default::default(),
        )
    }

    /// Sync the current "sele" selection to the sequence panel highlights
    pub(crate) fn sync_selection_to_sequence(&mut self) {
        use std::collections::HashSet;

        let mut highlighted = HashSet::new();

        if let Some(entry) = self.state.selections.get("sele") {
            let expression = entry.expression.clone();

            for name in self.state.registry.names() {
                if let Some(mol_obj) = self.state.registry.get_molecule(name) {
                    let mol = mol_obj.molecule();

                    let result = match pymol_select::select(mol, &expression) {
                        Ok(r) => r,
                        Err(_) => continue,
                    };

                    if !result.any() {
                        continue;
                    }

                    for residue in mol.residues() {
                        let any_selected = residue.atom_range.clone().any(|idx| {
                            result.contains_index(idx)
                        });
                        if any_selected {
                            highlighted.insert(ResidueRef {
                                object_name: name.to_string(),
                                chain_id: residue.chain().to_string(),
                                resv: residue.resv(),
                            });
                        }
                    }
                }
            }
        }

        if let Some(seq) = self
            .components
            .get_mut::<crate::components::SequenceComponent>()
        {
            seq.model.update_highlights(highlighted);
        }
    }
}
