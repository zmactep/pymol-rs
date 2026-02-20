//! Drag-and-drop overlay
//!
//! Shown when a file is dragged over the application window.

use std::path::Path;

use egui::{Align2, Color32, CornerRadius, Frame, Id, Margin, Stroke};

/// Overlay that displays a "Drop to Open" hint when a file is dragged over the window.
pub struct DragDropOverlay;

impl DragDropOverlay {
    /// Show a centered "Drop to open" overlay.
    ///
    /// * `ctx`  — egui context
    /// * `path` — the file currently being hovered
    pub fn show(ctx: &egui::Context, path: &Path) {
        // Dim the whole screen slightly
        let screen_rect = ctx.input(|i| i.content_rect());
        let painter = ctx.layer_painter(egui::LayerId::new(
            egui::Order::Foreground,
            Id::new("drag_drop_dim"),
        ));
        painter.rect_filled(
            screen_rect,
            0.0,
            Color32::from_rgba_unmultiplied(0, 0, 0, 100),
        );

        // Centered card
        egui::Area::new(Id::new("drag_drop_overlay"))
            .anchor(Align2::CENTER_CENTER, egui::vec2(0.0, 0.0))
            .interactable(false)
            .order(egui::Order::Foreground)
            .show(ctx, |ui| {
                Frame::NONE
                    .fill(Color32::from_rgba_unmultiplied(30, 30, 30, 230))
                    .stroke(Stroke::new(2.0, Color32::from_rgb(100, 180, 255)))
                    .corner_radius(CornerRadius::same(12))
                    .inner_margin(Margin::symmetric(40, 30))
                    .show(ui, |ui| {
                        ui.vertical_centered(|ui| {
                            ui.label(
                                egui::RichText::new("Drop to Open")
                                    .color(Color32::WHITE)
                                    .size(28.0)
                                    .strong(),
                            );

                            ui.add_space(8.0);

                            let filename = path
                                .file_name()
                                .map(|f| f.to_string_lossy().into_owned())
                                .unwrap_or_else(|| path.to_string_lossy().into_owned());

                            ui.label(
                                egui::RichText::new(&filename)
                                    .color(Color32::from_gray(180))
                                    .size(14.0),
                            );
                        });
                    });
            });

        ctx.request_repaint();
    }
}
