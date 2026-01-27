//! Mouse Mode Info Panel
//!
//! Displays current mouse mode and key bindings reference.

use egui::{Color32, RichText, Ui};

use crate::state::GuiState;

/// Mouse info panel
pub struct MouseInfoPanel;

impl MouseInfoPanel {
    /// Draw the mouse info panel
    pub fn show(ui: &mut Ui, state: &mut GuiState) {
        ui.vertical(|ui| {
            // Mouse mode header
            ui.horizontal(|ui| {
                ui.label(
                    RichText::new("Mouse Mode")
                        .color(Color32::from_rgb(100, 180, 255))
                        .strong(),
                );
                ui.label(
                    RichText::new(state.mouse_mode.display_name())
                        .color(Color32::WHITE),
                );
            });

            // Button mapping
            ui.horizontal(|ui| {
                ui.label(RichText::new("Buttons").color(Color32::from_rgb(100, 180, 255)));
                ui.label(RichText::new("L").color(Color32::WHITE));
                ui.label(RichText::new("M").color(Color32::WHITE));
                ui.label(RichText::new("R").color(Color32::WHITE));
                ui.label(RichText::new("Wheel").color(Color32::WHITE));
            });

            // Key combinations
            Self::show_key_binding(ui, "& Keys", "Rota", "Move", "MovZ", "Slab");
            Self::show_key_binding(ui, "Shft", "+Box", "-Box", "Clip", "MovS");
            Self::show_key_binding(ui, "Ctrl", "Move", "PkAt", "Pk1", "MvSZ");
            Self::show_key_binding(ui, "CtSh", "Sele", "Orig", "Clip", "MovZ");
            Self::show_key_binding(ui, "SnglClk", "+/-", "", "Cent", "Menu");
            Self::show_key_binding(ui, "DblClk", "Menu", "", "-", "PkAt");

            ui.separator();

            // Selection mode
            ui.horizontal(|ui| {
                ui.label(
                    RichText::new("Selecting")
                        .color(Color32::from_rgb(100, 180, 255)),
                );

                let mode_text = RichText::new(state.selecting_mode.display_name())
                    .color(Color32::from_rgb(255, 200, 100));

                if ui.add(egui::Label::new(mode_text).sense(egui::Sense::click())).clicked() {
                    state.selecting_mode = state.selecting_mode.next();
                }
            });

            // State info
            ui.horizontal(|ui| {
                ui.label(RichText::new("State").color(Color32::from_rgb(100, 180, 255)));
                ui.label(
                    RichText::new(format!("{}/", state.current_state))
                        .color(Color32::WHITE),
                );
                ui.label(
                    RichText::new(format!("{}", state.total_states))
                        .color(Color32::WHITE),
                );
            });
        });
    }

    fn show_key_binding(ui: &mut Ui, modifier: &str, l: &str, m: &str, r: &str, wheel: &str) {
        ui.horizontal(|ui| {
            ui.label(
                RichText::new(format!("{:8}", modifier))
                    .family(egui::FontFamily::Monospace)
                    .color(Color32::from_rgb(200, 200, 100))
                    .size(10.0),
            );
            for action in [l, m, r, wheel] {
                ui.label(
                    RichText::new(format!("{:5}", action))
                        .family(egui::FontFamily::Monospace)
                        .color(Color32::LIGHT_GRAY)
                        .size(10.0),
                );
            }
        });
    }
}
