//! Movie Panel Component
//!
//! Provides playback controls for movie/animation sequences.

use pymol_scene::Movie;

/// Action returned by the movie panel when user interacts with controls
#[derive(Debug, Clone, PartialEq)]
pub enum MovieAction {
    /// No action taken
    None,
    /// Play button clicked
    Play,
    /// Pause button clicked
    Pause,
    /// Stop button clicked
    Stop,
    /// Go to specific frame (0-indexed)
    GotoFrame(usize),
    /// Set playback FPS
    SetFps(f32),
    /// Set loop mode (0=once, 1=loop, 2=swing)
    SetLoopMode(u8),
}

/// Movie panel with playback controls
pub struct MoviePanel;

impl MoviePanel {
    /// Show the movie panel and return any action taken by the user
    ///
    /// This panel displays playback controls (play/pause/stop), a frame counter,
    /// and a frame slider when the movie has frames.
    ///
    /// # Arguments
    ///
    /// * `ui` - The egui UI context to render into
    /// * `movie` - Reference to the current movie state
    ///
    /// # Returns
    ///
    /// The action taken by the user, or `MovieAction::None` if no interaction
    pub fn show(ui: &mut egui::Ui, movie: &Movie) -> MovieAction {
        let mut action = MovieAction::None;

        ui.horizontal(|ui| {
            // Play/Pause toggle button
            let play_text = if movie.is_playing() { "⏸" } else { "▶" };
            let play_tooltip = if movie.is_playing() { "Pause" } else { "Play" };
            if ui.button(play_text).on_hover_text(play_tooltip).clicked() {
                action = if movie.is_playing() {
                    MovieAction::Pause
                } else {
                    MovieAction::Play
                };
            }

            // Stop button
            if ui.button("⏹").on_hover_text("Stop").clicked() {
                action = MovieAction::Stop;
            }

            // Rewind button
            if ui.button("⏮").on_hover_text("First frame").clicked() {
                action = MovieAction::GotoFrame(0);
            }

            // Previous frame button
            if ui.button("⏪").on_hover_text("Previous frame").clicked() {
                let current = movie.current_frame();
                if current > 0 {
                    action = MovieAction::GotoFrame(current - 1);
                }
            }

            // Next frame button
            if ui.button("⏩").on_hover_text("Next frame").clicked() {
                let current = movie.current_frame();
                let count = movie.frame_count();
                if current + 1 < count {
                    action = MovieAction::GotoFrame(current + 1);
                }
            }

            // Fast forward to end button
            if ui.button("⏭").on_hover_text("Last frame").clicked() {
                let count = movie.frame_count();
                if count > 0 {
                    action = MovieAction::GotoFrame(count - 1);
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
                action = MovieAction::GotoFrame(frame);
            }
        }

        action
    }

    /// Show a compact horizontal movie control bar
    ///
    /// This is a more compact version suitable for embedding in a toolbar.
    pub fn show_compact(ui: &mut egui::Ui, movie: &Movie) -> MovieAction {
        let mut action = MovieAction::None;

        ui.horizontal(|ui| {
            // Play/Pause toggle
            let play_text = if movie.is_playing() { "⏸" } else { "▶" };
            if ui.small_button(play_text).clicked() {
                action = if movie.is_playing() {
                    MovieAction::Pause
                } else {
                    MovieAction::Play
                };
            }

            // Stop
            if ui.small_button("⏹").clicked() {
                action = MovieAction::Stop;
            }

            // Frame counter
            let count = movie.frame_count();
            let current = movie.current_frame() + 1;
            ui.label(format!("{}/{}", current, count));
        });

        action
    }

    /// Show the movie panel as a collapsing header
    ///
    /// Wraps the full panel in a collapsible section.
    pub fn show_collapsible(ui: &mut egui::Ui, movie: &Movie) -> MovieAction {
        let mut action = MovieAction::None;

        let header = if movie.is_playing() {
            "▶ Movie (playing)"
        } else {
            "Movie"
        };

        egui::CollapsingHeader::new(header)
            .default_open(false)
            .show(ui, |ui| {
                action = Self::show(ui, movie);
            });

        action
    }
}
