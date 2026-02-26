//! Output/Log Panel
//!
//! Displays command output, info messages, warnings, and errors in a scrollable panel.

use egui::{Color32, RichText, ScrollArea, Ui};

use crate::state::{OutputBufferState, OutputKind};

/// Output panel that displays log messages
pub struct OutputPanel;

impl OutputPanel {
    /// Draw the output panel
    ///
    /// # Arguments
    /// * `ui` - The egui UI context
    /// * `state` - The output buffer state
    pub fn show(ui: &mut Ui, state: &OutputBufferState, max_height: f32) {
        ScrollArea::vertical()
            .max_height(max_height)
            .auto_shrink([false, false])
            .stick_to_bottom(state.auto_scroll)
            .show(ui, |ui| {
                ui.set_min_width(ui.available_width());

                for message in &state.buffer {
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
