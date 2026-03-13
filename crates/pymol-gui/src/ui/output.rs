//! Output/Log Panel
//!
//! Displays command output, info messages, warnings, and errors in a scrollable panel.

use egui::{Color32, RichText, ScrollArea, Ui};

use pymol_framework::model::{OutputModel, OutputKind};

/// Output panel that displays log messages
pub struct OutputPanel;

impl OutputPanel {
    /// Draw the output panel.
    ///
    /// The caller is responsible for constraining the available space
    /// (e.g., via `allocate_ui_with_layout`). The scroll area fills
    /// whatever space the parent `Ui` provides.
    pub fn show(ui: &mut Ui, model: &OutputModel) {
        ScrollArea::vertical()
            .auto_shrink([false, false])
            .stick_to_bottom(model.auto_scroll)
            .show(ui, |ui| {
                ui.set_min_width(ui.available_width());

                for message in &model.buffer {
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
