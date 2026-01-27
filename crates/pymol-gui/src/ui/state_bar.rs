//! State/Playback Controls Bar
//!
//! Bottom bar with playback controls and state information.

use egui::{Button, Color32, RichText, Ui};

use crate::state::GuiState;

/// State bar action
#[derive(Debug, Clone)]
pub enum StateBarAction {
    /// No action
    None,
    /// Go to first state
    FirstState,
    /// Go to previous state
    PreviousState,
    /// Stop playback
    Stop,
    /// Toggle play/pause
    TogglePlay,
    /// Go to next state
    NextState,
    /// Go to last state
    LastState,
    /// Toggle scene cycling
    SceneCycle,
    /// Toggle fullscreen
    Fullscreen,
}

/// State bar at the bottom of the window
pub struct StateBar;

impl StateBar {
    /// Draw the state bar and return any action to take
    pub fn show(ui: &mut Ui, state: &GuiState) -> Vec<StateBarAction> {
        let mut actions = Vec::new();

        ui.horizontal(|ui| {
            // Playback controls
            if ui.add(Button::new("|<").small()).on_hover_text("First state").clicked() {
                actions.push(StateBarAction::FirstState);
            }
            if ui.add(Button::new("◀").small()).on_hover_text("Previous state").clicked() {
                actions.push(StateBarAction::PreviousState);
            }
            if ui.add(Button::new("■").small()).on_hover_text("Stop").clicked() {
                actions.push(StateBarAction::Stop);
            }

            let play_btn = if state.is_playing {
                Button::new(RichText::new("▐▐").color(Color32::YELLOW)).small()
            } else {
                Button::new(RichText::new("▶").color(Color32::GREEN)).small()
            };
            if ui.add(play_btn).on_hover_text("Play/Pause").clicked() {
                actions.push(StateBarAction::TogglePlay);
            }

            if ui.add(Button::new("▶").small()).on_hover_text("Next state").clicked() {
                actions.push(StateBarAction::NextState);
            }
            if ui.add(Button::new(">|").small()).on_hover_text("Last state").clicked() {
                actions.push(StateBarAction::LastState);
            }

            // Spacer
            ui.add_space(ui.available_width() - 150.0);

            // State indicator
            ui.label(
                RichText::new(format!("State: {}/{}", state.current_state, state.total_states))
                    .color(Color32::LIGHT_GRAY)
                    .family(egui::FontFamily::Monospace),
            );

            // Scene/Fullscreen buttons
            if ui.add(Button::new("S").small()).on_hover_text("Scene cycle").clicked() {
                actions.push(StateBarAction::SceneCycle);
            }
            if ui.add(Button::new("F").small()).on_hover_text("Fullscreen").clicked() {
                actions.push(StateBarAction::Fullscreen);
            }
        });

        actions
    }
}
