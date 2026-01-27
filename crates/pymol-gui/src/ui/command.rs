//! Command Line Input Panel
//!
//! Provides a PyMOL-style command input with history navigation.

use egui::{Key, TextEdit, Ui};

use crate::state::GuiState;

/// Command to execute after UI processing
pub enum CommandAction {
    /// No action needed
    None,
    /// Execute the given command
    Execute(String),
}

/// Command line input panel
pub struct CommandLinePanel;

impl CommandLinePanel {
    /// Draw the command line panel and return any action to take
    pub fn show(ui: &mut Ui, state: &mut GuiState) -> CommandAction {
        let mut action = CommandAction::None;

        ui.horizontal(|ui| {
            ui.label(
                egui::RichText::new("PyMOL>")
                    .family(egui::FontFamily::Monospace)
                    .color(egui::Color32::from_rgb(150, 255, 150)),
            );

            let response = ui.add(
                TextEdit::singleline(&mut state.command_input)
                    .font(egui::TextStyle::Monospace)
                    .desired_width(ui.available_width() - 10.0)
                    .hint_text("Enter command..."),
            );

            // Track whether the command input has focus
            state.command_has_focus = response.has_focus();

            // Handle Enter key - check if text field lost focus due to Enter
            // This is the reliable way to detect Enter submission in egui
            if response.lost_focus() && ui.input(|i| i.key_pressed(Key::Enter)) {
                let command = state.take_command();
                if !command.is_empty() {
                    state.add_to_history(command.clone());
                    state.print_command(format!("PyMOL> {}", command));
                    action = CommandAction::Execute(command);
                }
                // Request re-focus after executing command
                state.command_wants_focus = true;
            }

            // Handle history navigation only when focused
            if response.has_focus() {
                if ui.input(|i| i.key_pressed(Key::ArrowUp)) {
                    state.history_previous();
                } else if ui.input(|i| i.key_pressed(Key::ArrowDown)) {
                    state.history_next();
                }
            }

            // Only request focus when explicitly requested (startup or after command)
            // This prevents stealing focus every frame
            if state.command_wants_focus {
                response.request_focus();
                state.command_wants_focus = false;
            }
        });

        action
    }
}
