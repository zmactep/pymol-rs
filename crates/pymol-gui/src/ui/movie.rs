//! Movie Panel Component
//!
//! Provides playback controls for movie/animation sequences.
//! Sends `AppMessage` directly to the `MessageBus`.

use pymol_scene::Movie;

use crate::message::{AppMessage, MessageBus};

/// Movie panel with playback controls
pub struct MoviePanel;

impl MoviePanel {
    /// Show the movie panel, sending messages directly to the bus.
    pub fn show(ui: &mut egui::Ui, movie: &Movie, bus: &mut MessageBus) {
        ui.horizontal(|ui| {
            // Play/Pause toggle button
            let play_text = if movie.is_playing() { "⏸" } else { "▶" };
            let play_tooltip = if movie.is_playing() { "Pause" } else { "Play" };
            if ui.button(play_text).on_hover_text(play_tooltip).clicked() {
                if movie.is_playing() {
                    bus.send(AppMessage::MoviePause);
                } else {
                    bus.send(AppMessage::MoviePlay);
                }
            }

            // Stop button
            if ui.button("⏹").on_hover_text("Stop").clicked() {
                bus.send(AppMessage::MovieStop);
            }

            // Rewind button
            if ui.button("⏮").on_hover_text("First frame").clicked() {
                bus.send(AppMessage::MovieGotoFrame(0));
            }

            // Previous frame button
            if ui.button("⏪").on_hover_text("Previous frame").clicked() {
                let current = movie.current_frame();
                if current > 0 {
                    bus.send(AppMessage::MovieGotoFrame(current - 1));
                }
            }

            // Next frame button
            if ui.button("⏩").on_hover_text("Next frame").clicked() {
                let current = movie.current_frame();
                let count = movie.frame_count();
                if current + 1 < count {
                    bus.send(AppMessage::MovieGotoFrame(current + 1));
                }
            }

            // Fast forward to end button
            if ui.button("⏭").on_hover_text("Last frame").clicked() {
                let count = movie.frame_count();
                if count > 0 {
                    bus.send(AppMessage::MovieGotoFrame(count - 1));
                }
            }

            ui.separator();

            // Frame counter display
            let count = movie.frame_count();
            let current = movie.current_frame() + 1; // 1-indexed for display
            ui.label(format!("{} / {}", current, count));
        });

        // Frame slider (only if movie has frames)
        let count = movie.frame_count();
        if count > 0 {
            let mut frame = movie.current_frame();
            let slider = egui::Slider::new(&mut frame, 0..=count.saturating_sub(1))
                .show_value(false)
                .text("Frame");
            if ui.add(slider).changed() {
                bus.send(AppMessage::MovieGotoFrame(frame));
            }
        }
    }

    /// Show a compact horizontal movie control bar
    pub fn show_compact(ui: &mut egui::Ui, movie: &Movie, bus: &mut MessageBus) {
        ui.horizontal(|ui| {
            // Play/Pause toggle
            let play_text = if movie.is_playing() { "⏸" } else { "▶" };
            if ui.small_button(play_text).clicked() {
                if movie.is_playing() {
                    bus.send(AppMessage::MoviePause);
                } else {
                    bus.send(AppMessage::MoviePlay);
                }
            }

            // Stop
            if ui.small_button("⏹").clicked() {
                bus.send(AppMessage::MovieStop);
            }

            // Frame counter
            let count = movie.frame_count();
            let current = movie.current_frame() + 1;
            ui.label(format!("{}/{}", current, count));
        });
    }

    /// Show the movie panel as a collapsing header
    pub fn show_collapsible(ui: &mut egui::Ui, movie: &Movie, bus: &mut MessageBus) {
        let header = if movie.is_playing() {
            "▶ Movie (playing)"
        } else {
            "Movie"
        };

        egui::CollapsingHeader::new(header)
            .default_open(false)
            .show(ui, |ui| {
                Self::show(ui, movie, bus);
            });
    }
}
