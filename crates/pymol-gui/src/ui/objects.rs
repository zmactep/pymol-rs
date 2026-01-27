//! Object List Panel
//!
//! Displays loaded objects with action buttons (A/S/H/L/C) for each.
//! Uses a single global popup to ensure only one menu is visible at a time.

use std::collections::HashMap;

use egui::{Color32, RichText, Ui, Vec2, Id};
use pymol_mol::RepMask;
use pymol_scene::{ObjectRegistry, SelectionEntry};

/// Minimum size for action buttons (A, S, H, L, C) to ensure consistent sizing
const BUTTON_MIN_SIZE: Vec2 = Vec2::splat(22.0);

/// Button colors for the action buttons
mod button_colors {
    use egui::Color32;
    pub const ACTION: Color32 = Color32::from_rgb(150, 150, 255);   // A - blue
    pub const SHOW: Color32 = Color32::from_rgb(100, 200, 100);     // S - green
    pub const HIDE: Color32 = Color32::from_rgb(200, 100, 100);     // H - red
    pub const LABEL: Color32 = Color32::from_rgb(200, 200, 255);    // L - light blue
    pub const COLOR: Color32 = Color32::from_rgb(255, 200, 100);    // C - orange/yellow
}

/// Create a consistent action button with fixed size
fn action_button(text: &str, color: Color32) -> egui::Button<'static> {
    egui::Button::new(RichText::new(text.to_string()).color(color))
        .min_size(BUTTON_MIN_SIZE)
}

/// Key for the single global popup ID
const GLOBAL_POPUP_KEY: &str = "object_list_global_popup";

/// Get the currently open popup ID
fn get_open_popup(ui: &Ui) -> Option<Id> {
    ui.ctx().data_mut(|d| d.get_persisted::<Id>(Id::new(GLOBAL_POPUP_KEY)))
}

/// Set the currently open popup ID (closing any other)
fn set_open_popup(ui: &Ui, id: Option<Id>) {
    ui.ctx().data_mut(|d| {
        if let Some(id) = id {
            d.insert_persisted(Id::new(GLOBAL_POPUP_KEY), id);
        } else {
            d.remove::<Id>(Id::new(GLOBAL_POPUP_KEY));
        }
    });
}

/// Menu button with popup - ensures only ONE popup is open at a time globally
fn menu_button_with_popup<R>(
    ui: &mut Ui,
    text: &str,
    color: Color32,
    hover_text: &str,
    add_contents: impl FnOnce(&mut Ui) -> R,
) -> Option<R> {
    let popup_id = ui.make_persistent_id(text);
    
    // Create and render button
    let button = egui::Button::new(RichText::new(text).color(color))
        .min_size(BUTTON_MIN_SIZE);
    let response = ui.add(button);
    
    // Check if THIS popup is currently open
    let is_this_popup_open = get_open_popup(ui) == Some(popup_id);
    
    // Toggle popup on click - this is the ONLY place we change popup state
    if response.clicked() {
        if is_this_popup_open {
            // Close this popup
            set_open_popup(ui, None);
        } else {
            // Open this popup (automatically closes any other)
            set_open_popup(ui, Some(popup_id));
        }
    }
    
    // Re-check after potential state change
    let is_this_popup_open = get_open_popup(ui) == Some(popup_id);
    
    // Show popup if this one is open
    let mut result = None;
    if is_this_popup_open {
        let popup_pos = response.rect.left_bottom();
        
        egui::Area::new(popup_id)
            .order(egui::Order::Foreground)
            .fixed_pos(popup_pos)
            .constrain(true)
            .show(ui.ctx(), |ui| {
                egui::Frame::popup(ui.style())
                    .show(ui, |ui| {
                        ui.set_min_width(100.0);
                        result = Some(add_contents(ui));
                    });
            });
        
        // Close on click outside (check if clicked anywhere but not on this button)
        let clicked_elsewhere = ui.input(|i| i.pointer.any_click()) && !response.clicked();
        if clicked_elsewhere {
            set_open_popup(ui, None);
        }
    }
    
    // Show tooltip when no popup is open
    if get_open_popup(ui).is_none() {
        response.on_hover_text(hover_text);
    }
    
    result
}

/// Action to perform on an object
pub enum ObjectAction {
    /// No action
    None,
    /// Toggle object visibility
    ToggleEnabled(String),
    /// Enable all objects
    EnableAll,
    /// Disable all objects
    DisableAll,
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
/// Color for hidden selections (dimmed)
const SELECTION_HIDDEN_COLOR: Color32 = Color32::from_rgb(128, 43, 128);

impl ObjectListPanel {
    /// Draw the object list panel
    /// 
    /// # Arguments
    /// * `ui` - The egui UI context
    /// * `registry` - Object registry containing molecules and other objects
    /// * `selections` - Named selections (name -> expression)
    pub fn show(ui: &mut Ui, registry: &ObjectRegistry, selections: &HashMap<String, SelectionEntry>) -> Vec<ObjectAction> {
        let mut actions = Vec::new();

        ui.vertical(|ui| {
            // Header
            ui.label(
                RichText::new("Objects")
                    .strong()
                    .color(Color32::WHITE),
            );

            ui.separator();

            // ========== "all" row ==========
            let all_enabled = registry.names().all(|name| {
                registry.get(name).map(|o| o.is_enabled()).unwrap_or(false)
            });
            let has_objects = registry.names().next().is_some();
            
            ui.horizontal(|ui| {
                let all_color = if !has_objects {
                    Color32::GRAY
                } else if all_enabled {
                    Color32::from_rgb(100, 255, 100)
                } else {
                    Color32::WHITE
                };
                
                let all_response = ui.add(
                    egui::Label::new(
                        RichText::new("all")
                            .family(egui::FontFamily::Monospace)
                            .color(all_color),
                    )
                    .sense(egui::Sense::click()),
                );
                
                if all_response.clicked() && has_objects {
                    if all_enabled {
                        actions.push(ObjectAction::DisableAll);
                    } else {
                        actions.push(ObjectAction::EnableAll);
                    }
                }
                
                all_response.on_hover_text(if all_enabled { "Click to disable all" } else { "Click to enable all" });

                // "all" row buttons - simple action buttons, no popups
                if ui
                    .add(action_button("A", button_colors::ACTION))
                    .on_hover_text("Actions menu")
                    .clicked()
                {
                    // Actions menu would go here
                }

                if ui
                    .add(action_button("S", button_colors::SHOW))
                    .on_hover_text("Show representations")
                    .clicked()
                {
                    for name in registry.names() {
                        actions.push(ObjectAction::ShowRep(name.to_string(), RepMask::LINES));
                    }
                }

                if ui
                    .add(action_button("H", button_colors::HIDE))
                    .on_hover_text("Hide all representations")
                    .clicked()
                {
                    for name in registry.names() {
                        actions.push(ObjectAction::HideAll(name.to_string()));
                    }
                }

                if ui
                    .add(action_button("L", button_colors::LABEL))
                    .on_hover_text("Toggle labels")
                    .clicked()
                {
                    // Toggle labels for all
                }

                if ui
                    .add(action_button("C", button_colors::COLOR))
                    .on_hover_text("Color all objects")
                    .clicked()
                {
                    // Color picker would go here
                }
            });

            ui.separator();

            // ========== Object rows ==========
            let names: Vec<_> = registry.names().map(|s| s.to_string()).collect();
            for name in names {
                let is_enabled = registry
                    .get(&name)
                    .map(|o| o.is_enabled())
                    .unwrap_or(false);

                ui.push_id(&name, |ui| {
                    ui.horizontal(|ui| {
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

                        // Menu buttons with popups
                        menu_button_with_popup(ui, "A", button_colors::ACTION, "Actions", |ui| {
                            Self::actions_menu(ui, &name, &mut actions);
                        });

                        menu_button_with_popup(ui, "S", button_colors::SHOW, "Show representations", |ui| {
                            Self::show_menu(ui, &name, &mut actions);
                        });

                        menu_button_with_popup(ui, "H", button_colors::HIDE, "Hide representations", |ui| {
                            Self::hide_menu(ui, &name, &mut actions);
                        });

                        if ui
                            .add(action_button("L", button_colors::LABEL))
                            .on_hover_text("Labels")
                            .clicked()
                        {
                            // Toggle labels
                        }

                        menu_button_with_popup(ui, "C", button_colors::COLOR, "Color", |ui| {
                            Self::color_menu(ui, &name, &mut actions);
                        });
                    });
                });
            }

            // ========== Selection rows ==========
            if !selections.is_empty() {
                ui.separator();
                
                let mut sel_entries: Vec<(&String, &SelectionEntry)> = selections.iter().collect();
                sel_entries.sort_by(|a, b| a.0.cmp(b.0));
                
                for (sel_name, entry) in sel_entries {
                    ui.push_id(format!("sel_{}", sel_name), |ui| {
                        ui.horizontal(|ui| {
                            let display_name = format!("({})", sel_name);
                            
                            // Use different color based on visibility
                            let color = if entry.visible {
                                SELECTION_COLOR  // Bright pink when visible
                            } else {
                                SELECTION_HIDDEN_COLOR  // Dimmed when hidden
                            };
                            
                            let name_response = ui.add(
                                egui::Label::new(
                                    RichText::new(&display_name)
                                        .family(egui::FontFamily::Monospace)
                                        .color(color),
                                )
                                .sense(egui::Sense::click()),
                            );

                            let clicked = name_response.clicked();

                            let visibility_hint = if entry.visible { "visible" } else { "hidden" };
                            name_response.on_hover_text(format!("Selection: {} ({})\nClick to toggle visibility", entry.expression, visibility_hint));

                            if clicked {
                                actions.push(ObjectAction::ToggleSelectionEnabled(sel_name.to_string()));
                            }

                            // Selection menu buttons
                            menu_button_with_popup(ui, "A", button_colors::ACTION, "Actions", |ui| {
                                Self::selection_actions_menu(ui, sel_name, &mut actions);
                            });

                            menu_button_with_popup(ui, "S", button_colors::SHOW, "Show representations", |ui| {
                                Self::show_menu(ui, sel_name, &mut actions);
                            });

                            menu_button_with_popup(ui, "H", button_colors::HIDE, "Hide representations", |ui| {
                                Self::hide_menu(ui, sel_name, &mut actions);
                            });

                            if ui
                                .add(action_button("L", button_colors::LABEL))
                                .on_hover_text("Labels")
                                .clicked()
                            {
                                // Toggle labels for selection
                            }

                            menu_button_with_popup(ui, "C", button_colors::COLOR, "Color selection", |ui| {
                                Self::color_menu(ui, sel_name, &mut actions);
                            });
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
            set_open_popup(ui, None);
        }
    }

    fn show_menu(ui: &mut Ui, name: &str, actions: &mut Vec<ObjectAction>) {
        if ui.button("lines").clicked() {
            actions.push(ObjectAction::ShowRep(name.to_string(), RepMask::LINES));
            set_open_popup(ui, None);
        }
        if ui.button("sticks").clicked() {
            actions.push(ObjectAction::ShowRep(name.to_string(), RepMask::STICKS));
            set_open_popup(ui, None);
        }
        if ui.button("spheres").clicked() {
            actions.push(ObjectAction::ShowRep(name.to_string(), RepMask::SPHERES));
            set_open_popup(ui, None);
        }
        if ui.button("cartoon").clicked() {
            actions.push(ObjectAction::ShowRep(name.to_string(), RepMask::CARTOON));
            set_open_popup(ui, None);
        }
        if ui.button("surface").clicked() {
            actions.push(ObjectAction::ShowRep(name.to_string(), RepMask::SURFACE));
            set_open_popup(ui, None);
        }
        if ui.button("mesh").clicked() {
            actions.push(ObjectAction::ShowRep(name.to_string(), RepMask::MESH));
            set_open_popup(ui, None);
        }
        if ui.button("dots").clicked() {
            actions.push(ObjectAction::ShowRep(name.to_string(), RepMask::DOTS));
            set_open_popup(ui, None);
        }
        if ui.button("ribbon").clicked() {
            actions.push(ObjectAction::ShowRep(name.to_string(), RepMask::RIBBON));
            set_open_popup(ui, None);
        }
        ui.separator();
        if ui.button("all").clicked() {
            actions.push(ObjectAction::ShowAll(name.to_string()));
            set_open_popup(ui, None);
        }
    }

    fn hide_menu(ui: &mut Ui, name: &str, actions: &mut Vec<ObjectAction>) {
        if ui.button("lines").clicked() {
            actions.push(ObjectAction::HideRep(name.to_string(), RepMask::LINES));
            set_open_popup(ui, None);
        }
        if ui.button("sticks").clicked() {
            actions.push(ObjectAction::HideRep(name.to_string(), RepMask::STICKS));
            set_open_popup(ui, None);
        }
        if ui.button("spheres").clicked() {
            actions.push(ObjectAction::HideRep(name.to_string(), RepMask::SPHERES));
            set_open_popup(ui, None);
        }
        if ui.button("cartoon").clicked() {
            actions.push(ObjectAction::HideRep(name.to_string(), RepMask::CARTOON));
            set_open_popup(ui, None);
        }
        if ui.button("surface").clicked() {
            actions.push(ObjectAction::HideRep(name.to_string(), RepMask::SURFACE));
            set_open_popup(ui, None);
        }
        if ui.button("mesh").clicked() {
            actions.push(ObjectAction::HideRep(name.to_string(), RepMask::MESH));
            set_open_popup(ui, None);
        }
        if ui.button("dots").clicked() {
            actions.push(ObjectAction::HideRep(name.to_string(), RepMask::DOTS));
            set_open_popup(ui, None);
        }
        if ui.button("ribbon").clicked() {
            actions.push(ObjectAction::HideRep(name.to_string(), RepMask::RIBBON));
            set_open_popup(ui, None);
        }
        ui.separator();
        if ui.button("everything").clicked() {
            actions.push(ObjectAction::HideAll(name.to_string()));
            set_open_popup(ui, None);
        }
    }

    fn actions_menu(ui: &mut Ui, name: &str, actions: &mut Vec<ObjectAction>) {
        if ui.button("zoom").clicked() {
            actions.push(ObjectAction::ZoomTo(name.to_string()));
            set_open_popup(ui, None);
        }
        if ui.button("center").clicked() {
            actions.push(ObjectAction::CenterOn(name.to_string()));
            set_open_popup(ui, None);
        }
        ui.separator();
        if ui.button("enable").clicked() {
            actions.push(ObjectAction::ToggleEnabled(name.to_string()));
            set_open_popup(ui, None);
        }
        if ui.button("disable").clicked() {
            actions.push(ObjectAction::ToggleEnabled(name.to_string()));
            set_open_popup(ui, None);
        }
        ui.separator();
        if ui.button("delete").clicked() {
            actions.push(ObjectAction::Delete(name.to_string()));
            set_open_popup(ui, None);
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
                set_open_popup(ui, None);
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
                set_open_popup(ui, None);
            }
        }
    }
}
