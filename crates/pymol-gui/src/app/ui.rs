//! egui UI Construction
//!
//! Main UI loop, message dispatch, and selection synchronization.
//! Uses the Layout engine to render components into configurable panels.

use pymol_mol::RepMask;
use pymol_scene::Object;
use pymol_select::SelectionResult;

use crate::component::SharedContext;
use crate::message::AppMessage;
use crate::model::ResidueRef;
use crate::ui::{DragDropOverlay, NotificationOverlay};

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

        // Check for pending async tasks and get their messages before entering the closure
        let pending_messages = self.task_runner.pending_messages();

        // Build SharedContext for components
        let builtin_names: Vec<String> = self.executor.registry().all_names().map(|s| s.to_string()).collect();
        let external_names: Vec<String> = self.ipc.external_commands.names().map(|s| s.to_string()).collect();
        let all_command_names: Vec<String> = builtin_names.into_iter().chain(external_names).collect();
        let setting_names = pymol_settings::setting_names();
        let setting_names_refs: Vec<&str> = setting_names.iter().copied().collect();
        let cmd_registry = self.executor.registry();

        // Update raytraced overlay texture if there's a new image
        let ray_overlay_info = if let Some(ray_img) = &self.state.raytraced_image {
            let needs_update = self.view.ray_overlay.size != Some((ray_img.width, ray_img.height));
            if needs_update {
                self.view.update_ray_overlay(&ray_img.data, ray_img.width, ray_img.height);
            }
            self.view.ray_overlay.texture_id.map(|id| (id, ray_img.width, ray_img.height))
        } else {
            if self.view.ray_overlay.texture_id.is_some() {
                self.view.clear_ray_overlay();
            }
            None
        };

        let shared = SharedContext {
            registry: &self.state.registry,
            camera: &self.state.camera,
            selections: &self.state.selections,
            named_colors: &self.state.named_colors,
            movie: &self.state.movie,
            command_names: &all_command_names,
            command_registry: &cmd_registry,
            setting_names: &setting_names_refs,
        };

        // Run egui with the layout engine
        let full_output = {
            let layout = &self.layout;
            let components = &mut self.components;
            let bus = &mut self.bus;

            self.view.egui.ctx.run(raw_input, |ctx| {
                // Layout renders all panels + returns viewport rect
                viewport_rect_logical = layout.show(ctx, components, &shared, bus);

                // Render ray overlay in viewport if present
                if let Some((texture_id, img_width, img_height)) = ray_overlay_info {
                    Self::render_ray_overlay(ctx, viewport_rect_logical, texture_id, img_width, img_height);
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
            })
        };

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

    /// Render the raytraced image overlay in the viewport area.
    fn render_ray_overlay(
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
            egui::Id::new("ray_overlay"),
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
        }
    }

    /// Dispatch a single message, mutating application state.
    ///
    /// The exhaustive match is kept for compiler safety; each arm delegates
    /// to a focused handler method.
    fn dispatch(&mut self, msg: &AppMessage) {
        match msg {
            // Command execution
            AppMessage::ExecuteCommand { command, silent } => {
                let _ = self.execute_command(command, *silent);
            }

            // Object visibility / deletion
            AppMessage::ToggleObject(_)
            | AppMessage::EnableAllObjects
            | AppMessage::DisableAllObjects
            | AppMessage::DeleteObject(_) => self.dispatch_object(msg),

            // Representation toggling
            AppMessage::ShowRepresentation { .. }
            | AppMessage::HideRepresentation { .. }
            | AppMessage::ShowAllRepresentations(_)
            | AppMessage::HideAllRepresentations(_)
            | AppMessage::SetColor { .. } => self.dispatch_representation(msg),

            // Camera navigation
            AppMessage::ZoomTo(name) => {
                let _ = self.execute_command(&format!("zoom {}", name), false);
            }
            AppMessage::CenterOn(name) => {
                let _ = self.execute_command(&format!("center {}", name), false);
            }

            // Selections
            AppMessage::DeleteSelection(_)
            | AppMessage::ToggleSelectionVisibility(_)
            | AppMessage::SelectionChanged { .. } => self.dispatch_selection(msg),

            // Movie
            AppMessage::MoviePlay
            | AppMessage::MoviePause
            | AppMessage::MovieStop
            | AppMessage::MovieGotoFrame(_)
            | AppMessage::MovieSetFps(_)
            | AppMessage::MovieSetLoopMode(_) => self.dispatch_movie(msg),

            // Layout
            AppMessage::TogglePanel(_)
            | AppMessage::FloatPanel(_)
            | AppMessage::DockPanel(_)
            | AppMessage::ActivateTab(_) => self.dispatch_layout(msg),

            // Viewport
            AppMessage::HoverResidue(info) => {
                self.viewport.sequence_hover = info.clone();
                self.scene_dirty = true;
            }

            // Handled via broadcast to components (no app-level mutation)
            AppMessage::PrintInfo(_) | AppMessage::PrintWarning(_)
            | AppMessage::PrintError(_) | AppMessage::PrintCommand(_)
            | AppMessage::UpdateHighlights(_) | AppMessage::FocusPanel(_) => {}

            // Lifecycle
            AppMessage::RequestRedraw => self.scene_dirty = true,
            AppMessage::Quit => self.frame.quit_requested = true,
            AppMessage::Custom { topic, payload } => {
                log::debug!("Custom event: {} ({} bytes)", topic, payload.len());
            }
        }
    }

    // =========================================================================
    // Grouped dispatch handlers
    // =========================================================================

    fn dispatch_object(&mut self, msg: &AppMessage) {
        match msg {
            AppMessage::ToggleObject(name) => {
                if let Some(obj) = self.state.registry.get_mut(name) {
                    if obj.is_enabled() { obj.disable(); } else { obj.enable(); }
                }
            }
            AppMessage::EnableAllObjects => {
                let names: Vec<String> = self.state.registry.names().map(|s| s.to_string()).collect();
                for name in names {
                    if let Some(obj) = self.state.registry.get_mut(&name) { obj.enable(); }
                }
            }
            AppMessage::DisableAllObjects => {
                let names: Vec<String> = self.state.registry.names().map(|s| s.to_string()).collect();
                for name in names {
                    if let Some(obj) = self.state.registry.get_mut(&name) { obj.disable(); }
                }
            }
            AppMessage::DeleteObject(name) => {
                self.state.registry.remove(name);
            }
            _ => unreachable!(),
        }
        self.scene_dirty = true;
    }

    fn dispatch_representation(&mut self, msg: &AppMessage) {
        match msg {
            AppMessage::ShowRepresentation { object, rep } => {
                if let Some(mol) = self.state.registry.get_molecule_mut(object) {
                    mol.show(*rep);
                }
            }
            AppMessage::HideRepresentation { object, rep } => {
                if let Some(mol) = self.state.registry.get_molecule_mut(object) {
                    mol.hide(*rep);
                }
            }
            AppMessage::ShowAllRepresentations(name) => {
                if let Some(mol) = self.state.registry.get_molecule_mut(name) {
                    mol.show(RepMask::ALL);
                }
            }
            AppMessage::HideAllRepresentations(name) => {
                if let Some(mol) = self.state.registry.get_molecule_mut(name) {
                    mol.hide_all();
                }
            }
            AppMessage::SetColor { object, color } => {
                if let Some(color_idx) = self.state.named_colors.get_by_name(color).map(|(idx, _)| idx) {
                    if let Some(mol_obj) = self.state.registry.get_molecule_mut(object) {
                        for atom in mol_obj.molecule_mut().atoms_mut() {
                            atom.repr.colors.base = color_idx as i32;
                        }
                    }
                }
            }
            _ => unreachable!(),
        }
        self.scene_dirty = true;
    }

    fn dispatch_selection(&mut self, msg: &AppMessage) {
        match msg {
            AppMessage::DeleteSelection(name) => {
                self.state.selections.remove(name);
                self.bus.print_info(format!("Deleted selection \"{}\"", name));
                self.scene_dirty = true;
            }
            AppMessage::ToggleSelectionVisibility(name) => {
                if let Some(entry) = self.state.selections.get_mut(name) {
                    entry.visible = !entry.visible;
                    log::debug!("Selection '{}' indicators toggled", name);
                    self.scene_dirty = true;
                }
            }
            AppMessage::SelectionChanged { .. } => {
                self.sync_selection_to_sequence();
            }
            _ => unreachable!(),
        }
    }

    fn dispatch_movie(&mut self, msg: &AppMessage) {
        match msg {
            AppMessage::MoviePlay => self.state.movie.play(),
            AppMessage::MoviePause => self.state.movie.pause(),
            AppMessage::MovieStop => self.state.movie.stop(),
            AppMessage::MovieGotoFrame(frame) => {
                self.state.movie.goto_frame(*frame);
                self.scene_dirty = true;
            }
            AppMessage::MovieSetFps(fps) => self.state.movie.set_fps(*fps),
            AppMessage::MovieSetLoopMode(mode) => self.state.movie.set_loop_mode(*mode),
            _ => unreachable!(),
        }
    }

    fn dispatch_layout(&mut self, msg: &AppMessage) {
        match msg {
            AppMessage::TogglePanel(id) => self.layout.toggle_expanded(id),
            AppMessage::FloatPanel(id) => self.layout.float_panel(id),
            AppMessage::DockPanel(id) => self.layout.dock_panel(id),
            AppMessage::ActivateTab(id) => self.layout.activate_tab(id),
            _ => unreachable!(),
        }
        self.scene_dirty = true;
    }

    /// Evaluate all visible selections for all molecules.
    pub(crate) fn evaluate_visible_selections(&self) -> Vec<(String, SelectionResult)> {
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
