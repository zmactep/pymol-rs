//! Movie Panel Component
//!
//! Provides playback controls for movie/animation sequences.
//! Sends commands via `MessageBus::execute_command()`.

use pymol_scene::Movie;

use crate::message::MessageBus;

/// Movie panel with playback controls
pub struct MoviePanel;

impl MoviePanel {
    /// Show the movie panel, sending commands directly to the bus.
    pub fn show(ui: &mut egui::Ui, movie: &Movie, bus: &mut MessageBus) {
        ui.horizontal(|ui| {
            // Play/Pause toggle button
            let play_text = if movie.is_playing() { "⏸" } else { "▶" };
            let play_tooltip = if movie.is_playing() { "Pause" } else { "Play" };
            if ui.button(play_text).on_hover_text(play_tooltip).clicked() {
                if movie.is_playing() {
                    bus.execute_command("mpause");
                } else {
                    bus.execute_command("mplay");
                }
            }

            // Stop button
            if ui.button("⏹").on_hover_text("Stop").clicked() {
                bus.execute_command("mstop");
            }

            // Rewind button
            if ui.button("⏮").on_hover_text("First frame").clicked() {
                bus.execute_command("rewind");
            }

            // Previous frame button
            if ui.button("⏪").on_hover_text("Previous frame").clicked() {
                bus.execute_command("backward");
            }

            // Next frame button
            if ui.button("⏩").on_hover_text("Next frame").clicked() {
                bus.execute_command("forward");
            }

            // Fast forward to end button
            if ui.button("⏭").on_hover_text("Last frame").clicked() {
                bus.execute_command("ending");
            }

            ui.separator();

            // Frame counter display
            let count = movie.effective_frame_count();
            let current = movie.current_frame() + 1; // 1-indexed for display
            ui.label(format!("{} / {}", current, count));
        });

        // Frame slider (only if movie has frames)
        let count = movie.effective_frame_count();
        if count > 0 {
            let mut frame = movie.current_frame();
            let slider = egui::Slider::new(&mut frame, 0..=count.saturating_sub(1))
                .show_value(false)
                .text("Frame");
            if ui.add(slider).changed() {
                // `frame` command uses 1-based indexing; silent to avoid flooding output
                bus.execute_command_silent(format!("frame {}", frame + 1));
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
                    bus.execute_command("mpause");
                } else {
                    bus.execute_command("mplay");
                }
            }

            // Stop
            if ui.small_button("⏹").clicked() {
                bus.execute_command("mstop");
            }

            // Frame counter
            let count = movie.effective_frame_count();
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
