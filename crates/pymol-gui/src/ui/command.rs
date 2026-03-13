//! Command Line Input Panel
//!
//! Provides a PyMOL-style command input with history navigation and tab completion.

use egui::{Key, TextEdit, Ui, Color32, RichText, ScrollArea, Sense, Id};

use pymol_framework::message::{AppMessage, MessageBus};
use pymol_framework::completion::{generate_completions, CompletionContext};
use pymol_framework::model::CommandLineModel;

// Re-export CommandLineState as CommandLineUiState for backwards compatibility
pub use pymol_framework::completion::CommandLineState as CommandLineUiState;

/// ID for the command line text edit widget
const COMMAND_INPUT_ID: &str = "pymol_command_input";

// =============================================================================
// View
// =============================================================================

/// Command line input panel
pub struct CommandLinePanel;

impl CommandLinePanel {
    /// Draw the command line panel, sending messages directly to the bus.
    pub fn show(
        ui: &mut Ui,
        model: &mut CommandLineModel,
        ui_state: &mut CommandLineUiState,
        ctx: &CompletionContext,
        bus: &mut MessageBus,
    ) {
        let mut text_edit_rect = egui::Rect::NOTHING;
        let text_edit_id = Id::new(COMMAND_INPUT_ID);

        // Consume Tab key before widgets to prevent focus navigation
        let tab_pressed = if ui_state.has_focus || ui_state.completion.visible {
            ui.ctx().input_mut(|i| i.consume_key(egui::Modifiers::NONE, Key::Tab))
        } else {
            false
        };

        // Track if we need to restore focus and cursor after completion
        let mut restore_focus = false;

        // Handle completion navigation keys when popup is visible
        if ui_state.completion.visible {
            if ui.ctx().input_mut(|i| i.consume_key(egui::Modifiers::NONE, Key::ArrowDown)) {
                ui_state.completion.select_next();
            }
            if ui.ctx().input_mut(|i| i.consume_key(egui::Modifiers::NONE, Key::ArrowUp)) {
                ui_state.completion.select_previous();
            }
            if ui.ctx().input_mut(|i| i.consume_key(egui::Modifiers::NONE, Key::Escape)) {
                ui_state.completion.reset();
            }
            if ui.ctx().input_mut(|i| i.consume_key(egui::Modifiers::NONE, Key::Enter)) {
                ui_state.completion.apply_to_input(&mut model.input);
                restore_focus = true;
            }
        }

        ui.horizontal(|ui| {
            ui.label(
                RichText::new("PyMOL>")
                    .family(egui::FontFamily::Monospace)
                    .color(Color32::from_rgb(150, 255, 150)),
            );

            let response = ui.add(
                TextEdit::singleline(&mut model.input)
                    .id(text_edit_id)
                    .font(egui::TextStyle::Monospace)
                    .desired_width(ui.available_width() - 10.0)
                    .hint_text("Enter command... (Tab for completion)")
                    .lock_focus(true),
            );

            text_edit_rect = response.rect;
            ui_state.has_focus = response.has_focus();

            // Handle Tab key for completion
            if tab_pressed {
                if ui_state.completion.visible {
                    ui_state.completion.apply_to_input(&mut model.input);
                } else {
                    Self::trigger_completion(model, ui_state, ctx);
                }
                restore_focus = true;
            }

            // Handle Enter key for command execution (only when completion not visible)
            if !ui_state.completion.visible && response.lost_focus() && ui.input(|i| i.key_pressed(Key::Enter)) {
                let command = model.take_command();
                if !command.is_empty() {
                    model.add_to_history(command.clone());
                    bus.send(AppMessage::ExecuteCommand {
                        command,
                        silent: false,
                    });
                }
                restore_focus = true;
            }

            // Handle history navigation only when focused and completion not visible
            if response.has_focus() && !ui_state.completion.visible {
                if ui.input(|i| i.key_pressed(Key::ArrowUp)) {
                    model.history_previous();
                } else if ui.input(|i| i.key_pressed(Key::ArrowDown)) {
                    model.history_next();
                }
            }

            // Hide completion when input changes (user is typing)
            if response.changed() {
                ui_state.completion.reset();
            }

            // Restore focus and cursor position after completion
            if restore_focus {
                ui.memory_mut(|mem| mem.request_focus(text_edit_id));
                Self::set_cursor_to_end(ui, text_edit_id, &model.input);
            }

            // Startup focus
            if ui_state.wants_focus {
                response.request_focus();
                ui_state.wants_focus = false;
            }
        });

        // Show completion popup
        if ui_state.completion.visible && !ui_state.completion.suggestions.is_empty() {
            Self::show_completion_popup(ui, model, ui_state, text_edit_rect);
        }
    }

    /// Generate completions for current input
    fn trigger_completion(model: &mut CommandLineModel, ui_state: &mut CommandLineUiState, ctx: &CompletionContext) {
        let cursor_pos = model.input.len();
        let result = generate_completions(&model.input, cursor_pos, ctx);

        if result.suggestions.len() == 1 {
            ui_state.completion.show(result.start_pos, result.suggestions);
            ui_state.completion.apply_to_input(&mut model.input);
        } else if !result.suggestions.is_empty() {
            ui_state.completion.show(result.start_pos, result.suggestions);
        }
    }

    /// Set the TextEdit cursor to the end of the text
    fn set_cursor_to_end(ui: &mut Ui, id: Id, text: &str) {
        if let Some(mut te_state) = egui::TextEdit::load_state(ui.ctx(), id) {
            let ccursor = egui::text::CCursor::new(text.len());
            te_state.cursor.set_char_range(Some(egui::text::CCursorRange::one(ccursor)));
            te_state.store(ui.ctx(), id);
        }
    }

    /// Show the completion popup
    fn show_completion_popup(ui: &mut Ui, model: &mut CommandLineModel, ui_state: &mut CommandLineUiState, text_edit_rect: egui::Rect) {
        let popup_pos = egui::pos2(text_edit_rect.left(), text_edit_rect.bottom() + 2.0);
        let mut clicked_idx: Option<usize> = None;

        egui::Area::new(ui.make_persistent_id("completion_popup"))
            .fixed_pos(popup_pos)
            .order(egui::Order::Foreground)
            .show(ui.ctx(), |ui| {
                egui::Frame::popup(ui.style())
                    .fill(Color32::from_rgb(30, 30, 40))
                    .stroke(egui::Stroke::new(1.0, Color32::from_rgb(80, 80, 100)))
                    .corner_radius(4.0)
                    .inner_margin(egui::Margin::same(4))
                    .show(ui, |ui| {
                        let max_len = ui_state.completion.suggestions.iter().map(|s| s.len()).max().unwrap_or(10);
                        let popup_width = (max_len as f32 * 8.0 + 24.0).max(200.0).min(text_edit_rect.width());
                        ui.set_min_width(popup_width);

                        let row_height = 20.0;
                        let visible_count = ui_state.completion.suggestions.len().min(10);
                        let scroll_height = visible_count as f32 * row_height;

                        ScrollArea::vertical()
                            .max_height(scroll_height)
                            .auto_shrink([false, true])
                            .show(ui, |ui| {
                                ui.set_min_width(popup_width - 16.0);

                                for (idx, suggestion) in ui_state.completion.suggestions.iter().enumerate() {
                                    let is_selected = idx == ui_state.completion.selected;
                                    let bg = if is_selected { Color32::from_rgb(60, 100, 150) } else { Color32::TRANSPARENT };
                                    let fg = if is_selected { Color32::WHITE } else { Color32::from_rgb(200, 200, 200) };

                                    let resp = ui.allocate_ui_with_layout(
                                        egui::vec2(ui.available_width(), row_height),
                                        egui::Layout::left_to_right(egui::Align::Center),
                                        |ui| {
                                            ui.painter().rect_filled(ui.max_rect(), 2.0, bg);
                                            ui.add_space(8.0);
                                            ui.label(RichText::new(suggestion).family(egui::FontFamily::Monospace).color(fg).size(13.0));
                                        },
                                    );

                                    if is_selected && ui_state.completion.scroll_to_selected {
                                        resp.response.scroll_to_me(Some(egui::Align::Center));
                                    }
                                    if resp.response.interact(Sense::click()).clicked() {
                                        clicked_idx = Some(idx);
                                    }
                                }
                            });

                        ui.add_space(4.0);
                        ui.label(RichText::new("Tab/Enter: select | Up/Down: navigate | Esc: close").size(10.0).color(Color32::from_rgb(120, 120, 140)));
                    });
            });

        // Clear scroll flag after rendering
        ui_state.completion.scroll_to_selected = false;

        if let Some(idx) = clicked_idx {
            ui_state.completion.selected = idx;
            ui_state.completion.apply_to_input(&mut model.input);
        }
    }
}
