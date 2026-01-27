//! Output/Log Panel
//!
//! Displays command output, info messages, warnings, and errors in a scrollable panel.

use egui::{Color32, RichText, ScrollArea, Ui};

use crate::state::{GuiState, OutputKind};

/// Output panel that displays log messages
pub struct OutputPanel;

impl OutputPanel {
    /// Draw the output panel
    pub fn show(ui: &mut Ui, state: &mut GuiState) {
        let available_height = state.output_panel_height;

        ScrollArea::vertical()
            .max_height(available_height)
            .auto_shrink([false, false])
            .stick_to_bottom(state.output_auto_scroll)
            .show(ui, |ui| {
                ui.set_min_width(ui.available_width());

                for message in &state.output_buffer {
                    let color = match message.kind {
                        OutputKind::Normal => Color32::LIGHT_GRAY,
                        OutputKind::Info => Color32::from_rgb(100, 180, 255),
                        OutputKind::Warning => Color32::YELLOW,
                        OutputKind::Error => Color32::from_rgb(255, 100, 100),
                        OutputKind::Command => Color32::from_rgb(150, 255, 150),
                    };

                    let text = RichText::new(&message.text)
                        .color(color)
                        .family(egui::FontFamily::Monospace)
                        .size(12.0);

                    ui.label(text);
                }
            });
    }
}
