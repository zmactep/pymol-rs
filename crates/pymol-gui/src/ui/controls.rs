//! Control Buttons Panel
//!
//! Provides view controls, selection tools, and playback buttons.

use egui::{Button, RichText, Ui};

/// Control action to perform
#[derive(Debug, Clone)]
pub enum ControlAction {
    /// No action
    None,
    /// Reset view to default
    Reset,
    /// Zoom to fit all objects
    Zoom,
    /// Orient to principal axes
    Orient,
    /// Toggle orthographic projection
    ToggleOrtho,
    /// Take screenshot
    Screenshot,
    /// Ray trace image
    Ray,
    /// Clear all picks
    Unpick,
    /// Deselect all
    Deselect,
    /// Toggle rocking animation
    Rock,
    /// Get current view matrix
    GetView,
    /// Go to first state
    FirstState,
    /// Go to previous state
    PreviousState,
    /// Stop playback
    Stop,
    /// Start playback
    Play,
    /// Go to next state
    NextState,
    /// Go to last state
    LastState,
    /// Clear movie
    MovieClear,
    /// Open builder panel
    Builder,
    /// Open properties panel
    Properties,
    /// Rebuild representations
    Rebuild,
}

/// Control buttons panel
pub struct ControlsPanel;

impl ControlsPanel {
    /// Draw the controls panel and return any action to take
    pub fn show(ui: &mut Ui, is_playing: bool, is_rocking: bool) -> Vec<ControlAction> {
        let mut actions = Vec::new();

        ui.vertical(|ui| {
            // View controls row 1
            ui.horizontal(|ui| {
                if ui.add(Button::new("Reset")).on_hover_text("Reset view (R)").clicked() {
                    actions.push(ControlAction::Reset);
                }
                if ui.add(Button::new("Zoom")).on_hover_text("Zoom to fit all").clicked() {
                    actions.push(ControlAction::Zoom);
                }
                if ui.add(Button::new("Orient")).on_hover_text("Orient to principal axes").clicked() {
                    actions.push(ControlAction::Orient);
                }
                ui.menu_button("Draw/Ray", |ui| {
                    if ui.button("Draw").clicked() {
                        ui.close_menu();
                    }
                    if ui.button("Ray").clicked() {
                        actions.push(ControlAction::Ray);
                        ui.close_menu();
                    }
                    if ui.button("Screenshot").clicked() {
                        actions.push(ControlAction::Screenshot);
                        ui.close_menu();
                    }
                });
            });

            // View controls row 2
            ui.horizontal(|ui| {
                if ui.add(Button::new("Unpick")).on_hover_text("Clear all picks").clicked() {
                    actions.push(ControlAction::Unpick);
                }
                if ui.add(Button::new("Deselect")).on_hover_text("Clear selection").clicked() {
                    actions.push(ControlAction::Deselect);
                }
                let rock_text = if is_rocking {
                    RichText::new("Rock").color(egui::Color32::GREEN)
                } else {
                    RichText::new("Rock")
                };
                if ui.add(Button::new(rock_text)).on_hover_text("Toggle rocking").clicked() {
                    actions.push(ControlAction::Rock);
                }
                if ui.add(Button::new("Get View")).on_hover_text("Copy view matrix").clicked() {
                    actions.push(ControlAction::GetView);
                }
            });

            ui.separator();

            // Playback controls
            ui.horizontal(|ui| {
                if ui.add(Button::new("|<")).on_hover_text("First state").clicked() {
                    actions.push(ControlAction::FirstState);
                }
                if ui.add(Button::new("<")).on_hover_text("Previous state").clicked() {
                    actions.push(ControlAction::PreviousState);
                }
                if ui.add(Button::new("Stop")).on_hover_text("Stop playback").clicked() {
                    actions.push(ControlAction::Stop);
                }
                let play_text = if is_playing {
                    RichText::new("Pause").color(egui::Color32::YELLOW)
                } else {
                    RichText::new("Play").color(egui::Color32::GREEN)
                };
                if ui.add(Button::new(play_text)).on_hover_text("Play/Pause").clicked() {
                    actions.push(ControlAction::Play);
                }
                if ui.add(Button::new(">")).on_hover_text("Next state").clicked() {
                    actions.push(ControlAction::NextState);
                }
                if ui.add(Button::new(">|")).on_hover_text("Last state").clicked() {
                    actions.push(ControlAction::LastState);
                }
                if ui.add(Button::new("MClear")).on_hover_text("Clear movie").clicked() {
                    actions.push(ControlAction::MovieClear);
                }
            });

            ui.separator();

            // Tool buttons
            ui.horizontal(|ui| {
                if ui.add(Button::new("Builder")).on_hover_text("Open builder").clicked() {
                    actions.push(ControlAction::Builder);
                }
                if ui.add(Button::new("Properties")).on_hover_text("Object properties").clicked() {
                    actions.push(ControlAction::Properties);
                }
                if ui.add(Button::new("Rebuild")).on_hover_text("Rebuild representations").clicked() {
                    actions.push(ControlAction::Rebuild);
                }
            });
        });

        actions
    }
}
