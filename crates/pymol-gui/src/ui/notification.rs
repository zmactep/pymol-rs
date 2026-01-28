//! Notification Overlay Component
//!
//! Displays a notification overlay in the bottom-left corner of the screen
//! to indicate that async operations are in progress.

use egui::{Align2, Color32, CornerRadius, Frame, Id, Margin, Stroke};

/// Notification overlay that displays in the bottom-left corner
pub struct NotificationOverlay;

impl NotificationOverlay {
    /// Show the notification overlay with messages from pending tasks
    ///
    /// Renders a rounded rectangle notification in the bottom-left corner
    /// with a spinner and the provided messages. If multiple messages are
    /// provided, they are displayed on separate lines.
    ///
    /// # Arguments
    ///
    /// * `ctx` - The egui context
    /// * `messages` - Slice of messages to display (one per pending task)
    pub fn show(ctx: &egui::Context, messages: &[String]) {
        if messages.is_empty() {
            return;
        }

        // Use an Area for absolute positioning
        egui::Area::new(Id::new("async_notification_overlay"))
            .anchor(Align2::LEFT_BOTTOM, egui::vec2(10.0, -10.0))
            .interactable(false) // Don't block mouse events
            .order(egui::Order::Foreground)
            .show(ctx, |ui| {
                // Create a frame with rounded corners and semi-transparent background
                Frame::NONE
                    .fill(Color32::from_rgba_unmultiplied(40, 40, 40, 220))
                    .stroke(Stroke::new(1.0, Color32::from_gray(80)))
                    .corner_radius(CornerRadius::same(6))
                    .inner_margin(Margin::symmetric(12, 8))
                    .show(ui, |ui| {
                        ui.horizontal(|ui| {
                            // Add a spinner
                            ui.spinner();
                            ui.add_space(8.0);

                            if messages.len() == 1 {
                                // Single message - show inline
                                ui.label(
                                    egui::RichText::new(&messages[0])
                                        .color(Color32::from_gray(220))
                                        .size(13.0),
                                );
                            } else {
                                // Multiple messages - show vertically
                                ui.vertical(|ui| {
                                    for message in messages {
                                        ui.label(
                                            egui::RichText::new(message)
                                                .color(Color32::from_gray(220))
                                                .size(13.0),
                                        );
                                    }
                                });
                            }
                        });
                    });
            });

        // Request continuous repaint while notification is shown (for spinner animation)
        ctx.request_repaint();
    }
}
