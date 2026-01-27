//! Object List Panel
//!
//! Displays loaded objects with action buttons (A/S/H/L/C) for each.

use std::collections::HashMap;

use egui::{Color32, RichText, Ui};
use pymol_mol::RepMask;
use pymol_scene::ObjectRegistry;

/// Action to perform on an object
pub enum ObjectAction {
    /// No action
    None,
    /// Toggle object visibility
    ToggleEnabled(String),
    /// Show all representations
    ShowAll(String),
    /// Hide all representations
    HideAll(String),
    /// Show specific representation
    ShowRep(String, u32),
    /// Hide specific representation  
    HideRep(String, u32),
    /// Set color
    SetColor(String, String),
    /// Delete object
    Delete(String),
    /// Zoom to object
    ZoomTo(String),
    /// Center on object
    CenterOn(String),
    /// Delete a named selection
    DeleteSelection(String),
    /// Toggle selection enabled state
    ToggleSelectionEnabled(String),
}

/// Object list panel
pub struct ObjectListPanel;

/// Color for selections (pink/magenta like PyMOL)
const SELECTION_COLOR: Color32 = Color32::from_rgb(255, 85, 255);

impl ObjectListPanel {
    /// Draw the object list panel
    /// 
    /// # Arguments
    /// * `ui` - The egui UI context
    /// * `registry` - Object registry containing molecules and other objects
    /// * `selections` - Named selections (name -> expression)
    pub fn show(ui: &mut Ui, registry: &ObjectRegistry, selections: &HashMap<String, String>) -> Vec<ObjectAction> {
        let mut actions = Vec::new();

        ui.vertical(|ui| {
            // Header
            ui.horizontal(|ui| {
                ui.label(
                    RichText::new("Objects")
                        .strong()
                        .color(Color32::WHITE),
                );
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    ui.label(
                        RichText::new("A  S  H  L  C")
                            .small()
                            .color(Color32::GRAY),
                    );
                });
            });

            ui.separator();

            // "all" row
            ui.horizontal(|ui| {
                ui.label(
                    RichText::new("all")
                        .family(egui::FontFamily::Monospace)
                        .color(Color32::WHITE),
                );

                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    // C - Color (all objects)
                    if ui
                        .add(egui::Button::new(RichText::new("C").color(Color32::from_rgb(255, 200, 100))).small())
                        .on_hover_text("Color all objects")
                        .clicked()
                    {
                        // Color picker would go here
                    }

                    // L - Label
                    if ui
                        .add(egui::Button::new(RichText::new("L").color(Color32::from_rgb(200, 200, 255))).small())
                        .on_hover_text("Toggle labels")
                        .clicked()
                    {
                        // Toggle labels for all
                    }

                    // H - Hide
                    if ui
                        .add(egui::Button::new(RichText::new("H").color(Color32::from_rgb(200, 100, 100))).small())
                        .on_hover_text("Hide all representations")
                        .clicked()
                    {
                        for name in registry.names() {
                            actions.push(ObjectAction::HideAll(name.to_string()));
                        }
                    }

                    // S - Show
                    if ui
                        .add(egui::Button::new(RichText::new("S").color(Color32::from_rgb(100, 200, 100))).small())
                        .on_hover_text("Show representations")
                        .clicked()
                    {
                        for name in registry.names() {
                            actions.push(ObjectAction::ShowRep(name.to_string(), RepMask::LINES));
                        }
                    }

                    // A - Actions
                    if ui
                        .add(egui::Button::new(RichText::new("A").color(Color32::from_rgb(150, 150, 255))).small())
                        .on_hover_text("Actions menu")
                        .clicked()
                    {
                        // Actions menu would go here
                    }
                });
            });

            ui.separator();

            // Individual objects (molecules, etc.)
            let names: Vec<_> = registry.names().map(|s| s.to_string()).collect();
            for name in names {
                let is_enabled = registry
                    .get(&name)
                    .map(|o| o.is_enabled())
                    .unwrap_or(false);

                ui.horizontal(|ui| {
                    // Object name with enabled/disabled indicator
                    let name_color = if is_enabled {
                        Color32::from_rgb(100, 255, 100)
                    } else {
                        Color32::GRAY
                    };

                    let name_response = ui.add(
                        egui::Label::new(
                            RichText::new(&name)
                                .family(egui::FontFamily::Monospace)
                                .color(name_color),
                        )
                        .sense(egui::Sense::click()),
                    );

                    if name_response.clicked() {
                        actions.push(ObjectAction::ToggleEnabled(name.clone()));
                    }
                    if name_response.secondary_clicked() {
                        actions.push(ObjectAction::ZoomTo(name.clone()));
                    }

                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        // C - Color
                        ui.menu_button(
                            RichText::new("C").color(Color32::from_rgb(255, 200, 100)).small(),
                            |ui| {
                                Self::color_menu(ui, &name, &mut actions);
                            },
                        ).response.on_hover_text("Color");

                        // L - Label
                        if ui
                            .add(egui::Button::new(RichText::new("L").color(Color32::from_rgb(200, 200, 255))).small())
                            .on_hover_text("Labels")
                            .clicked()
                        {
                            // Toggle labels
                        }

                        // H - Hide menu
                        ui.menu_button(
                            RichText::new("H").color(Color32::from_rgb(200, 100, 100)).small(),
                            |ui| {
                                Self::hide_menu(ui, &name, &mut actions);
                            },
                        ).response.on_hover_text("Hide representations");

                        // S - Show menu
                        ui.menu_button(
                            RichText::new("S").color(Color32::from_rgb(100, 200, 100)).small(),
                            |ui| {
                                Self::show_menu(ui, &name, &mut actions);
                            },
                        ).response.on_hover_text("Show representations");

                        // A - Actions menu
                        ui.menu_button(
                            RichText::new("A").color(Color32::from_rgb(150, 150, 255)).small(),
                            |ui| {
                                Self::actions_menu(ui, &name, &mut actions);
                            },
                        ).response.on_hover_text("Actions");
                    });
                });
            }

            // Named selections (shown after objects, with parentheses like PyMOL)
            if !selections.is_empty() {
                ui.separator();
                
                // Sort selection names for consistent display
                let mut sel_names: Vec<_> = selections.keys().collect();
                sel_names.sort();
                
                for sel_name in sel_names {
                    ui.horizontal(|ui| {
                        // Display selection with parentheses: (name)
                        let display_name = format!("({})", sel_name);
                        
                        let name_response = ui.add(
                            egui::Label::new(
                                RichText::new(&display_name)
                                    .family(egui::FontFamily::Monospace)
                                    .color(SELECTION_COLOR),
                            )
                            .sense(egui::Sense::click()),
                        );

                        // Check clicked before using response for hover (to avoid move)
                        let clicked = name_response.clicked();

                        // Show expression on hover
                        if let Some(expr) = selections.get(sel_name.as_str()) {
                            name_response.on_hover_text(format!("Selection: {}", expr));
                        }

                        if clicked {
                            actions.push(ObjectAction::ToggleSelectionEnabled(sel_name.to_string()));
                        }

                        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                            // C - Color (apply color to selection)
                            ui.menu_button(
                                RichText::new("C").color(Color32::from_rgb(255, 200, 100)).small(),
                                |ui| {
                                    Self::color_menu(ui, sel_name, &mut actions);
                                },
                            ).response.on_hover_text("Color selection");

                            // L - Label (placeholder)
                            if ui
                                .add(egui::Button::new(RichText::new("L").color(Color32::from_rgb(200, 200, 255))).small())
                                .on_hover_text("Labels")
                                .clicked()
                            {
                                // Toggle labels for selection
                            }

                            // H - Hide
                            ui.menu_button(
                                RichText::new("H").color(Color32::from_rgb(200, 100, 100)).small(),
                                |ui| {
                                    Self::hide_menu(ui, sel_name, &mut actions);
                                },
                            ).response.on_hover_text("Hide representations");

                            // S - Show
                            ui.menu_button(
                                RichText::new("S").color(Color32::from_rgb(100, 200, 100)).small(),
                                |ui| {
                                    Self::show_menu(ui, sel_name, &mut actions);
                                },
                            ).response.on_hover_text("Show representations");

                            // A - Actions menu for selections
                            ui.menu_button(
                                RichText::new("A").color(Color32::from_rgb(150, 150, 255)).small(),
                                |ui| {
                                    Self::selection_actions_menu(ui, sel_name, &mut actions);
                                },
                            ).response.on_hover_text("Actions");
                        });
                    });
                }
            }
        });

        actions
    }

    /// Actions menu for selections
    fn selection_actions_menu(ui: &mut Ui, name: &str, actions: &mut Vec<ObjectAction>) {
        if ui.button("delete").clicked() {
            actions.push(ObjectAction::DeleteSelection(name.to_string()));
            ui.close_menu();
        }
    }

    fn show_menu(ui: &mut Ui, name: &str, actions: &mut Vec<ObjectAction>) {
        if ui.button("lines").clicked() {
            actions.push(ObjectAction::ShowRep(name.to_string(), RepMask::LINES));
            ui.close_menu();
        }
        if ui.button("sticks").clicked() {
            actions.push(ObjectAction::ShowRep(name.to_string(), RepMask::STICKS));
            ui.close_menu();
        }
        if ui.button("spheres").clicked() {
            actions.push(ObjectAction::ShowRep(name.to_string(), RepMask::SPHERES));
            ui.close_menu();
        }
        if ui.button("cartoon").clicked() {
            actions.push(ObjectAction::ShowRep(name.to_string(), RepMask::CARTOON));
            ui.close_menu();
        }
        if ui.button("surface").clicked() {
            actions.push(ObjectAction::ShowRep(name.to_string(), RepMask::SURFACE));
            ui.close_menu();
        }
        if ui.button("mesh").clicked() {
            actions.push(ObjectAction::ShowRep(name.to_string(), RepMask::MESH));
            ui.close_menu();
        }
        if ui.button("dots").clicked() {
            actions.push(ObjectAction::ShowRep(name.to_string(), RepMask::DOTS));
            ui.close_menu();
        }
        if ui.button("ribbon").clicked() {
            actions.push(ObjectAction::ShowRep(name.to_string(), RepMask::RIBBON));
            ui.close_menu();
        }
        ui.separator();
        if ui.button("all").clicked() {
            actions.push(ObjectAction::ShowAll(name.to_string()));
            ui.close_menu();
        }
    }

    fn hide_menu(ui: &mut Ui, name: &str, actions: &mut Vec<ObjectAction>) {
        if ui.button("lines").clicked() {
            actions.push(ObjectAction::HideRep(name.to_string(), RepMask::LINES));
            ui.close_menu();
        }
        if ui.button("sticks").clicked() {
            actions.push(ObjectAction::HideRep(name.to_string(), RepMask::STICKS));
            ui.close_menu();
        }
        if ui.button("spheres").clicked() {
            actions.push(ObjectAction::HideRep(name.to_string(), RepMask::SPHERES));
            ui.close_menu();
        }
        if ui.button("cartoon").clicked() {
            actions.push(ObjectAction::HideRep(name.to_string(), RepMask::CARTOON));
            ui.close_menu();
        }
        if ui.button("surface").clicked() {
            actions.push(ObjectAction::HideRep(name.to_string(), RepMask::SURFACE));
            ui.close_menu();
        }
        if ui.button("mesh").clicked() {
            actions.push(ObjectAction::HideRep(name.to_string(), RepMask::MESH));
            ui.close_menu();
        }
        if ui.button("dots").clicked() {
            actions.push(ObjectAction::HideRep(name.to_string(), RepMask::DOTS));
            ui.close_menu();
        }
        if ui.button("ribbon").clicked() {
            actions.push(ObjectAction::HideRep(name.to_string(), RepMask::RIBBON));
            ui.close_menu();
        }
        ui.separator();
        if ui.button("everything").clicked() {
            actions.push(ObjectAction::HideAll(name.to_string()));
            ui.close_menu();
        }
    }

    fn actions_menu(ui: &mut Ui, name: &str, actions: &mut Vec<ObjectAction>) {
        if ui.button("zoom").clicked() {
            actions.push(ObjectAction::ZoomTo(name.to_string()));
            ui.close_menu();
        }
        if ui.button("center").clicked() {
            actions.push(ObjectAction::CenterOn(name.to_string()));
            ui.close_menu();
        }
        ui.separator();
        if ui.button("enable").clicked() {
            actions.push(ObjectAction::ToggleEnabled(name.to_string()));
            ui.close_menu();
        }
        if ui.button("disable").clicked() {
            actions.push(ObjectAction::ToggleEnabled(name.to_string()));
            ui.close_menu();
        }
        ui.separator();
        if ui.button("delete").clicked() {
            actions.push(ObjectAction::Delete(name.to_string()));
            ui.close_menu();
        }
    }

    fn color_menu(ui: &mut Ui, name: &str, actions: &mut Vec<ObjectAction>) {
        let colors = [
            ("red", Color32::RED),
            ("green", Color32::GREEN),
            ("blue", Color32::BLUE),
            ("yellow", Color32::YELLOW),
            ("cyan", Color32::from_rgb(0, 255, 255)),
            ("magenta", Color32::from_rgb(255, 0, 255)),
            ("orange", Color32::from_rgb(255, 165, 0)),
            ("white", Color32::WHITE),
            ("gray", Color32::GRAY),
        ];

        for (color_name, color) in colors {
            if ui
                .add(egui::Button::new(RichText::new(color_name).color(color)))
                .clicked()
            {
                actions.push(ObjectAction::SetColor(name.to_string(), color_name.to_string()));
                ui.close_menu();
            }
        }

        ui.separator();

        let by_colors = ["by element", "by chain", "by ss", "by residue"];
        for by_color in by_colors {
            if ui.button(by_color).clicked() {
                actions.push(ObjectAction::SetColor(
                    name.to_string(),
                    by_color.replace(' ', "_"),
                ));
                ui.close_menu();
            }
        }
    }
}
