//! egui UI Construction
//!
//! Main UI loop, object action handling, and selection synchronization.

use pymol_mol::RepMask;
use pymol_scene::Object;
use pymol_select::SelectionResult;

use crate::ui::command::CommandAction;
use crate::ui::objects::ObjectAction;
use crate::ui::sequence::{ResidueRef, SequenceAction};
use crate::ui::{DragDropOverlay, NotificationOverlay, OutputPanel};

use super::App;

impl App {
    /// Run the egui UI and return the output for rendering later
    pub(crate) fn run_egui_ui(&mut self) -> Option<(Vec<egui::ClippedPrimitive>, egui::TexturesDelta)> {
        // Take egui input first
        let raw_input = {
            let (egui_state, window) = match (&mut self.view.egui_state, &self.view.window) {
                (Some(state), Some(window)) => (state, window),
                _ => return None,
            };
            egui_state.take_egui_input(window)
        };

        // Get scale factor for converting to physical pixels
        let scale_factor = self.view.scale_factor();

        // Collect data needed for UI without borrowing self in the closure
        let mut commands_to_execute = Vec::new();
        let mut object_actions = Vec::new();
        let mut sequence_actions = Vec::new();
        let mut viewport_rect_logical = egui::Rect::NOTHING;

        // Check for pending async tasks and get their messages before entering the closure
        let pending_messages = self.task_runner.pending_messages();

        // Get command names for autocomplete (combine built-in + external)
        let builtin_names: Vec<&str> = self.executor.registry().all_names().collect();
        let external_names: Vec<&str> = self.external_commands.names().collect();
        let all_command_names: Vec<&str> = builtin_names.into_iter().chain(external_names).collect();

        // Build completion context
        let setting_names = pymol_settings::setting_names();
        let color_names: Vec<String> = self.state.named_colors.names().iter().map(|s| s.to_string()).collect();
        let object_names: Vec<String> = self.state.registry.names().map(|s| s.to_string()).collect();
        let selection_names: Vec<String> = self.state.selections.names();

        // Update raytraced overlay texture if there's a new image
        // This must happen outside the egui closure to avoid borrow conflicts
        let ray_overlay_info = if let Some(ray_img) = &self.state.raytraced_image {
            // Check if we need to update the texture (new image or size changed)
            let needs_update = self.view.ray_overlay_size != Some((ray_img.width, ray_img.height));
            if needs_update {
                self.view.update_ray_overlay(&ray_img.data, ray_img.width, ray_img.height);
            }
            self.view.ray_overlay_texture_id.map(|id| (id, ray_img.width, ray_img.height))
        } else {
            // Clear overlay if raytraced image was removed
            if self.view.ray_overlay_texture_id.is_some() {
                self.view.clear_ray_overlay();
            }
            None
        };

        // Run egui with a closure that captures only what it needs
        let full_output = {
            let output = &mut self.output;
            let command_line = &mut self.command_line;
            let ui_config = &mut self.view.ui_config;
            let registry = &self.state.registry;
            let camera = &self.state.camera;
            let selections = &self.state.selections;
            let named_colors = &self.state.named_colors;
            let object_list_panel = &mut self.object_list_panel;
            let sequence_panel = &mut self.sequence_panel;
            let cmd_registry = self.executor.registry();
            let setting_names_refs: Vec<&str> = setting_names.iter().copied().collect();

            let completion_ctx = crate::ui::completion::CompletionContext {
                command_names: &all_command_names,
                registry: &cmd_registry,
                setting_names: &setting_names_refs,
                color_names: &color_names,
                object_names: &object_names,
                selection_names: &selection_names,
            };

            self.view.egui_ctx.run(raw_input, |ctx| {
                // Top panel - output and command line
                egui::TopBottomPanel::top("top_panel").show(ctx, |ui| {
                    if ui_config.show_output_panel {
                        OutputPanel::show(ui, output, ui_config.output_panel_height);

                        // Draggable resize handle between output and command line
                        let resize_id = ui.id().with("output_resize_handle");
                        let sep_response = ui.separator();
                        let handle_rect = sep_response.rect.expand2(egui::vec2(0.0, 2.0));
                        let response =
                            ui.interact(handle_rect, resize_id, egui::Sense::drag());
                        if response.dragged() {
                            ui_config.output_panel_height =
                                (ui_config.output_panel_height + response.drag_delta().y)
                                    .clamp(30.0, 600.0);
                        }
                        if response.hovered() || response.dragged() {
                            ui.ctx().set_cursor_icon(egui::CursorIcon::ResizeVertical);
                        }
                    }

                    match crate::ui::CommandLinePanel::show(ui, command_line, &completion_ctx) {
                        CommandAction::Execute(cmd) => {
                            commands_to_execute.push(cmd);
                        }
                        CommandAction::None => {}
                    }
                });

                // Bottom panel - sequence viewer (docked or floating)
                if ui_config.show_sequence_panel {
                    let has_sele = selections.contains("sele");

                    if ui_config.sequence_panel_floating {
                        // Floating mode: egui::Window
                        let mut open = true;
                        egui::Window::new("Sequence")
                            .id(egui::Id::new("sequence_panel_window"))
                            .open(&mut open)
                            .default_size(ui_config.sequence_window_size)
                            .min_size([200.0, 60.0])
                            .resizable(true)
                            .collapsible(true)
                            .show(ctx, |ui| {
                                ui.horizontal(|ui| {
                                    ui.with_layout(
                                        egui::Layout::right_to_left(egui::Align::Center),
                                        |ui| {
                                            if ui
                                                .small_button("\u{2B07}")
                                                .on_hover_text("Dock panel")
                                                .clicked()
                                            {
                                                ui_config.sequence_panel_floating = false;
                                                sequence_panel.clear_drag();
                                            }
                                        },
                                    );
                                });
                                ui.separator();
                                sequence_actions =
                                    sequence_panel.show(ui, registry, named_colors, has_sele);
                            });
                        if !open {
                            ui_config.show_sequence_panel = false;
                        }
                    } else {
                        // Docked mode: TopBottomPanel at bottom
                        let panel_height = sequence_panel.desired_height(registry);
                        egui::TopBottomPanel::bottom("sequence_panel")
                            .resizable(true)
                            .default_height(panel_height)
                            .height_range(30.0..=200.0)
                            .show(ctx, |ui| {
                                // Float button in top-right corner
                                let rect = ui.max_rect();
                                let btn_size = egui::vec2(20.0, 16.0);
                                let btn_pos = egui::pos2(
                                    rect.right() - btn_size.x - 4.0,
                                    rect.top() + 2.0,
                                );
                                let btn_rect = egui::Rect::from_min_size(btn_pos, btn_size);
                                if ui
                                    .put(
                                        btn_rect,
                                        egui::Button::new(
                                            egui::RichText::new("\u{2197}").size(12.0),
                                        )
                                        .frame(false),
                                    )
                                    .on_hover_text("Float panel")
                                    .clicked()
                                {
                                    ui_config.sequence_panel_floating = true;
                                    sequence_panel.clear_drag();
                                }

                                sequence_actions =
                                    sequence_panel.show(ui, registry, named_colors, has_sele);
                            });
                    }
                }

                // Right panel - object list only
                if ui_config.show_control_panel {
                    egui::SidePanel::right("right_panel")
                        .default_width(ui_config.right_panel_width)
                        .show(ctx, |ui| {
                            let panel_x_range = ui.max_rect().x_range();

                            // Bottom toolbar (pinned below the scrollable list)
                            egui::TopBottomPanel::bottom("right_panel_toolbar")
                                .show_separator_line(false)
                                .show_inside(ui, |ui| {
                                    // Full-width separator using the outer panel's x range
                                    let y = ui.max_rect().top();
                                    let stroke = ui.visuals().widgets.noninteractive.bg_stroke;
                                    ui.painter().hline(panel_x_range, y, stroke);
                                    ui.add_space(4.0);

                                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                                        let size = ui.spacing().interact_size.y;
                                        let btn = egui::Button::new("S")
                                            .min_size(egui::vec2(size, size))
                                            .selected(ui_config.show_sequence_panel);
                                        if ui.add(btn).on_hover_text("Sequence viewer").clicked() {
                                            ui_config.show_sequence_panel = !ui_config.show_sequence_panel;
                                        }
                                    });
                                });

                            // Scrollable object list fills remaining space
                            egui::ScrollArea::vertical().show(ui, |ui| {
                                object_actions = object_list_panel.show(ui, registry, selections);
                            });
                        });
                }

                // Central panel to capture the viewport rect (area for 3D rendering)
                let central_response = egui::CentralPanel::default()
                    .frame(egui::Frame::NONE)
                    .show(ctx, |ui| {
                        let available = ui.available_size();

                        // If we have a raytraced overlay, render it scaled to fit the viewport
                        if let Some((texture_id, img_width, img_height)) = ray_overlay_info {
                            // Calculate scaling to fit the image in the viewport while maintaining aspect ratio
                            let img_aspect = img_width as f32 / img_height as f32;
                            let vp_aspect = available.x / available.y;

                            let (display_width, display_height) = if img_aspect > vp_aspect {
                                // Image is wider than viewport - fit to width
                                (available.x, available.x / img_aspect)
                            } else {
                                // Image is taller than viewport - fit to height
                                (available.y * img_aspect, available.y)
                            };

                            // Center the image in the viewport
                            let offset_x = (available.x - display_width) / 2.0;
                            let offset_y = (available.y - display_height) / 2.0;

                            // Allocate space for interaction
                            let (rect, response) = ui.allocate_exact_size(available, egui::Sense::hover());

                            let image_rect = egui::Rect::from_min_size(
                                egui::pos2(rect.min.x + offset_x, rect.min.y + offset_y),
                                egui::vec2(display_width, display_height),
                            );

                            // Render the raytraced image
                            ui.painter().image(
                                texture_id,
                                image_rect,
                                egui::Rect::from_min_max(egui::pos2(0.0, 0.0), egui::pos2(1.0, 1.0)),
                                egui::Color32::WHITE,
                            );

                            response
                        } else {
                            // No overlay - just allocate the space for 3D rendering
                            ui.allocate_response(available, egui::Sense::hover())
                        }
                    });
                viewport_rect_logical = central_response.response.rect;

                // Paint atom labels as 2D overlay on the 3D viewport
                {
                    let vp = viewport_rect_logical;
                    let vp_tuple = (vp.min.x, vp.min.y, vp.width(), vp.height());
                    let painter = ctx.layer_painter(egui::LayerId::new(
                        egui::Order::Foreground,
                        egui::Id::new("labels_overlay"),
                    ));
                    let font = egui::FontId::new(14.0, egui::FontFamily::Proportional);
                    let label_color = egui::Color32::WHITE;

                    for name in registry.names() {
                        let mol_obj = match registry.get_molecule(name) {
                            Some(m) if m.is_enabled() => m,
                            _ => continue,
                        };

                        for (pos, text) in mol_obj.collect_labels() {
                            if let Some((sx, sy)) = camera.project_to_screen(pos, vp_tuple) {
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

                    // Paint measurement labels
                    let meas_font = egui::FontId::new(13.0, egui::FontFamily::Proportional);
                    let meas_color = egui::Color32::from_rgb(255, 255, 0);
                    for name in registry.names() {
                        if let Some(meas_obj) = registry.get_measurement(name) {
                            if !meas_obj.is_enabled() {
                                continue;
                            }
                            for (pos, text) in meas_obj.collect_labels() {
                                if let Some((sx, sy)) = camera.project_to_screen(pos, vp_tuple) {
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

                // Show notification overlay when async tasks are in progress
                if !pending_messages.is_empty() {
                    NotificationOverlay::show(ctx, &pending_messages);
                }

                // Show drag-and-drop hint when a file is being dragged over the window
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
        if let (Some(egui_state), Some(window)) = (&mut self.view.egui_state, &self.view.window) {
            egui_state.handle_platform_output(window, full_output.platform_output);
        }

        // Process commands and actions
        let has_actions = !commands_to_execute.is_empty()
            || !object_actions.is_empty()
            || !sequence_actions.is_empty();
        for cmd in commands_to_execute {
            // Errors are displayed in the GUI output, no need to propagate
            let _ = self.execute_command(&cmd, false);
        }
        for action in object_actions {
            self.handle_object_action(action);
        }
        for action in sequence_actions {
            match action {
                SequenceAction::Execute(cmd) => {
                    let _ = self.execute_command(&cmd, false);
                }
                SequenceAction::Notify(msg) => {
                    self.output.print_info(msg);
                }
                SequenceAction::Hover(info) => {
                    self.sequence_hover = info;
                    self.needs_redraw = true;
                }
            }
        }

        // Sync 3D selection highlights to sequence panel
        if has_actions {
            self.sync_selection_to_sequence();
        }

        // Check if egui needs a repaint (e.g., for popup menus, animations, etc.)
        // repaint_delay of Duration::ZERO means immediate repaint needed
        let egui_needs_repaint = full_output
            .viewport_output
            .values()
            .any(|v| v.repaint_delay.is_zero());

        // Request window redraw if we processed any actions OR if egui needs repaint
        // This is needed because the event loop is in Wait mode and won't redraw
        // unless explicitly requested
        if has_actions || egui_needs_repaint {
            self.view.request_redraw();
        }

        let clipped_primitives = self.view.egui_ctx.tessellate(full_output.shapes, full_output.pixels_per_point);

        Some((clipped_primitives, full_output.textures_delta))
    }

    /// Evaluate all visible selections for all molecules.
    ///
    /// Delegates to `SelectionManager::evaluate_visible` which builds full
    /// evaluation contexts (implicit object-name selections, pre-evaluated
    /// named selections) for each molecule.
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

        // Look for the "sele" selection (the default selection name)
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

                    // Map selected atoms to residues
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

        self.sequence_panel.update_highlights(highlighted);
    }

    /// Handle object list actions
    pub(crate) fn handle_object_action(&mut self, action: ObjectAction) {
        match action {
            ObjectAction::None => {}
            ObjectAction::ToggleEnabled(name) => {
                if let Some(obj) = self.state.registry.get_mut(&name) {
                    let enabled = obj.is_enabled();
                    if enabled {
                        obj.disable();
                    } else {
                        obj.enable();
                    }
                    self.needs_redraw = true;
                }
            }
            ObjectAction::EnableAll => {
                let names: Vec<String> = self.state.registry.names().map(|s| s.to_string()).collect();
                for name in names {
                    if let Some(obj) = self.state.registry.get_mut(&name) {
                        obj.enable();
                    }
                }
                self.needs_redraw = true;
            }
            ObjectAction::DisableAll => {
                let names: Vec<String> = self.state.registry.names().map(|s| s.to_string()).collect();
                for name in names {
                    if let Some(obj) = self.state.registry.get_mut(&name) {
                        obj.disable();
                    }
                }
                self.needs_redraw = true;
            }
            ObjectAction::ShowAll(name) => {
                if let Some(mol) = self.state.registry.get_molecule_mut(&name) {
                    mol.show(RepMask::ALL);
                    self.needs_redraw = true;
                }
            }
            ObjectAction::HideAll(name) => {
                if let Some(mol) = self.state.registry.get_molecule_mut(&name) {
                    mol.hide_all();
                    self.needs_redraw = true;
                }
            }
            ObjectAction::ShowRep(name, rep) => {
                if let Some(mol) = self.state.registry.get_molecule_mut(&name) {
                    mol.show(rep);
                    self.needs_redraw = true;
                }
            }
            ObjectAction::HideRep(name, rep) => {
                if let Some(mol) = self.state.registry.get_molecule_mut(&name) {
                    mol.hide(rep);
                    self.needs_redraw = true;
                }
            }
            ObjectAction::SetColor(name, color) => {
                if let Some(color_idx) = self.state.named_colors.get_by_name(&color).map(|(idx, _)| idx) {
                    if let Some(mol_obj) = self.state.registry.get_molecule_mut(&name) {
                        for atom in mol_obj.molecule_mut().atoms_mut() {
                            atom.repr.colors.base = color_idx as i32;
                        }
                        self.needs_redraw = true;
                    }
                }
            }
            ObjectAction::Delete(name) => {
                self.state.registry.remove(&name);
                self.needs_redraw = true;
            }
            ObjectAction::ZoomTo(name) => {
                let _ = self.execute_command(&format!("zoom {}", name), false);
            }
            ObjectAction::CenterOn(name) => {
                let _ = self.execute_command(&format!("center {}", name), false);
            }
            ObjectAction::DeleteSelection(name) => {
                self.state.selections.remove(&name);
                self.output.print_info(format!("Deleted selection \"{}\"", name));
                self.needs_redraw = true;
            }
            ObjectAction::ToggleSelectionEnabled(name) => {
                if let Some(entry) = self.state.selections.get_mut(&name) {
                    entry.visible = !entry.visible;
                    let state = if entry.visible { "shown" } else { "hidden" };
                    log::debug!("Selection '{}' indicators {}", name, state);
                    self.needs_redraw = true;
                }
            }
        }
    }
}
