//! Object List Panel
//!
//! Displays loaded objects with action buttons (A/S/H/L/C) for each.
//! Uses a stateful menu system where buttons only update state,
//! and a single menu is rendered at the end to avoid overlapping popups.

use std::collections::HashMap;

use egui::{Color32, RichText, Ui, Vec2, Pos2, Id};
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

/// What menu is currently showing
#[derive(Clone, PartialEq, Debug)]
pub enum ActiveMenu {
    /// No menu open
    None,
    /// Actions menu for an object (zoom, center, enable, disable, delete)
    Actions { object: String },
    /// Show representations menu for an object
    Show { object: String },
    /// Hide representations menu for an object
    Hide { object: String },
    /// Color menu for an object
    Color { object: String },
    /// Actions menu for a selection (just delete)
    SelectionActions { selection: String },
}

impl Default for ActiveMenu {
    fn default() -> Self {
        Self::None
    }
}

/// Stores which menu is active and where to position it
#[derive(Default)]
pub struct MenuState {
    /// Which menu is currently active
    pub active: ActiveMenu,
    /// Position to anchor the menu popup
    pub anchor_pos: Pos2,
}

/// Create a consistent action button with fixed size
fn action_button(text: &str, color: Color32) -> egui::Button<'static> {
    egui::Button::new(RichText::new(text.to_string()).color(color))
        .min_size(BUTTON_MIN_SIZE)
}

/// Menu button that only updates state - NEVER renders popups directly
fn menu_button(
    ui: &mut Ui,
    text: &str,
    color: Color32,
    hover_text: &str,
    menu_state: &mut MenuState,
    target_menu: ActiveMenu,
) -> egui::Response {
    let button = egui::Button::new(RichText::new(text).color(color))
        .min_size(BUTTON_MIN_SIZE);
    let response = ui.add(button);
    
    if response.clicked() {
        if menu_state.active == target_menu {
            // Clicking the same button closes the menu
            menu_state.active = ActiveMenu::None;
        } else {
            // Open this menu (automatically closes any other)
            menu_state.active = target_menu;
            menu_state.anchor_pos = response.rect.left_bottom();
        }
    }
    
    // Show tooltip only when no menu is open
    if menu_state.active == ActiveMenu::None {
        response.clone().on_hover_text(hover_text);
    }
    
    response
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

/// Object list panel - stateful to track menu state
#[derive(Default)]
pub struct ObjectListPanel {
    /// Current menu state
    menu_state: MenuState,
}

/// Color for selections (pink/magenta like PyMOL)
const SELECTION_COLOR: Color32 = Color32::from_rgb(255, 85, 255);
/// Color for hidden selections (dimmed)
const SELECTION_HIDDEN_COLOR: Color32 = Color32::from_rgb(128, 43, 128);

impl ObjectListPanel {
    /// Create a new ObjectListPanel
    pub fn new() -> Self {
        Self::default()
    }

    /// Draw the object list panel
    /// 
    /// # Arguments
    /// * `ui` - The egui UI context
    /// * `registry` - Object registry containing molecules and other objects
    /// * `selections` - Named selections (name -> expression)
    pub fn show(&mut self, ui: &mut Ui, registry: &ObjectRegistry, selections: &HashMap<String, SelectionEntry>) -> Vec<ObjectAction> {
        let mut actions = Vec::new();
        
        // Track if any menu button was clicked this frame
        let mut menu_button_clicked = false;

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

                // "all" row buttons - simple action buttons, no menus
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

                        // Menu buttons - they only update state, NEVER render popups
                        if menu_button(
                            ui, "A", button_colors::ACTION, "Actions",
                            &mut self.menu_state,
                            ActiveMenu::Actions { object: name.clone() }
                        ).clicked() {
                            menu_button_clicked = true;
                        }

                        if menu_button(
                            ui, "S", button_colors::SHOW, "Show representations",
                            &mut self.menu_state,
                            ActiveMenu::Show { object: name.clone() }
                        ).clicked() {
                            menu_button_clicked = true;
                        }

                        if menu_button(
                            ui, "H", button_colors::HIDE, "Hide representations",
                            &mut self.menu_state,
                            ActiveMenu::Hide { object: name.clone() }
                        ).clicked() {
                            menu_button_clicked = true;
                        }

                        if ui
                            .add(action_button("L", button_colors::LABEL))
                            .on_hover_text("Labels")
                            .clicked()
                        {
                            // Toggle labels
                        }

                        if menu_button(
                            ui, "C", button_colors::COLOR, "Color",
                            &mut self.menu_state,
                            ActiveMenu::Color { object: name.clone() }
                        ).clicked() {
                            menu_button_clicked = true;
                        }
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
                            if menu_button(
                                ui, "A", button_colors::ACTION, "Actions",
                                &mut self.menu_state,
                                ActiveMenu::SelectionActions { selection: sel_name.to_string() }
                            ).clicked() {
                                menu_button_clicked = true;
                            }

                            if menu_button(
                                ui, "S", button_colors::SHOW, "Show representations",
                                &mut self.menu_state,
                                ActiveMenu::Show { object: sel_name.to_string() }
                            ).clicked() {
                                menu_button_clicked = true;
                            }

                            if menu_button(
                                ui, "H", button_colors::HIDE, "Hide representations",
                                &mut self.menu_state,
                                ActiveMenu::Hide { object: sel_name.to_string() }
                            ).clicked() {
                                menu_button_clicked = true;
                            }

                            if ui
                                .add(action_button("L", button_colors::LABEL))
                                .on_hover_text("Labels")
                                .clicked()
                            {
                                // Toggle labels for selection
                            }

                            if menu_button(
                                ui, "C", button_colors::COLOR, "Color selection",
                                &mut self.menu_state,
                                ActiveMenu::Color { object: sel_name.to_string() }
                            ).clicked() {
                                menu_button_clicked = true;
                            }
                        });
                    });
                }
            }
        });

        // ========== RENDER THE SINGLE ACTIVE MENU ==========
        // This happens AFTER all buttons are processed, so state is final
        if self.menu_state.active != ActiveMenu::None {
            let menu_clicked = self.render_active_menu(ui, &mut actions);
            
            // Close menu on click outside (but not if we clicked a menu button or menu item)
            if !menu_button_clicked && !menu_clicked {
                if ui.input(|i| i.pointer.any_click()) {
                    self.menu_state.active = ActiveMenu::None;
                }
            }
        }

        actions
    }

    /// Render the currently active menu - called ONCE at the end
    /// Returns true if a menu item was clicked
    fn render_active_menu(&mut self, ui: &mut Ui, actions: &mut Vec<ObjectAction>) -> bool {
        let mut item_clicked = false;
        
        egui::Area::new(Id::new("object_list_menu_panel"))
            .order(egui::Order::Foreground)
            .fixed_pos(self.menu_state.anchor_pos)
            .constrain(true)
            .show(ui.ctx(), |ui| {
                egui::Frame::popup(ui.style())
                    .show(ui, |ui| {
                        ui.set_min_width(100.0);
                        
                        match &self.menu_state.active.clone() {
                            ActiveMenu::None => {}
                            ActiveMenu::Actions { object } => {
                                item_clicked = Self::render_actions_menu(ui, object, actions);
                            }
                            ActiveMenu::Show { object } => {
                                item_clicked = Self::render_show_menu(ui, object, actions);
                            }
                            ActiveMenu::Hide { object } => {
                                item_clicked = Self::render_hide_menu(ui, object, actions);
                            }
                            ActiveMenu::Color { object } => {
                                item_clicked = Self::render_color_menu(ui, object, actions);
                            }
                            ActiveMenu::SelectionActions { selection } => {
                                item_clicked = Self::render_selection_actions_menu(ui, selection, actions);
                            }
                        }
                    });
            });
        
        // Close menu if an item was clicked
        if item_clicked {
            self.menu_state.active = ActiveMenu::None;
        }
        
        item_clicked
    }

    /// Render actions menu for selections - returns true if item clicked
    fn render_selection_actions_menu(ui: &mut Ui, name: &str, actions: &mut Vec<ObjectAction>) -> bool {
        if ui.button("delete").clicked() {
            actions.push(ObjectAction::DeleteSelection(name.to_string()));
            return true;
        }
        false
    }

    /// Render show menu - returns true if item clicked
    fn render_show_menu(ui: &mut Ui, name: &str, actions: &mut Vec<ObjectAction>) -> bool {
        if ui.button("lines").clicked() {
            actions.push(ObjectAction::ShowRep(name.to_string(), RepMask::LINES));
            return true;
        }
        if ui.button("sticks").clicked() {
            actions.push(ObjectAction::ShowRep(name.to_string(), RepMask::STICKS));
            return true;
        }
        if ui.button("spheres").clicked() {
            actions.push(ObjectAction::ShowRep(name.to_string(), RepMask::SPHERES));
            return true;
        }
        if ui.button("cartoon").clicked() {
            actions.push(ObjectAction::ShowRep(name.to_string(), RepMask::CARTOON));
            return true;
        }
        if ui.button("surface").clicked() {
            actions.push(ObjectAction::ShowRep(name.to_string(), RepMask::SURFACE));
            return true;
        }
        if ui.button("mesh").clicked() {
            actions.push(ObjectAction::ShowRep(name.to_string(), RepMask::MESH));
            return true;
        }
        if ui.button("dots").clicked() {
            actions.push(ObjectAction::ShowRep(name.to_string(), RepMask::DOTS));
            return true;
        }
        if ui.button("ribbon").clicked() {
            actions.push(ObjectAction::ShowRep(name.to_string(), RepMask::RIBBON));
            return true;
        }
        ui.separator();
        if ui.button("all").clicked() {
            actions.push(ObjectAction::ShowAll(name.to_string()));
            return true;
        }
        false
    }

    /// Render hide menu - returns true if item clicked
    fn render_hide_menu(ui: &mut Ui, name: &str, actions: &mut Vec<ObjectAction>) -> bool {
        if ui.button("lines").clicked() {
            actions.push(ObjectAction::HideRep(name.to_string(), RepMask::LINES));
            return true;
        }
        if ui.button("sticks").clicked() {
            actions.push(ObjectAction::HideRep(name.to_string(), RepMask::STICKS));
            return true;
        }
        if ui.button("spheres").clicked() {
            actions.push(ObjectAction::HideRep(name.to_string(), RepMask::SPHERES));
            return true;
        }
        if ui.button("cartoon").clicked() {
            actions.push(ObjectAction::HideRep(name.to_string(), RepMask::CARTOON));
            return true;
        }
        if ui.button("surface").clicked() {
            actions.push(ObjectAction::HideRep(name.to_string(), RepMask::SURFACE));
            return true;
        }
        if ui.button("mesh").clicked() {
            actions.push(ObjectAction::HideRep(name.to_string(), RepMask::MESH));
            return true;
        }
        if ui.button("dots").clicked() {
            actions.push(ObjectAction::HideRep(name.to_string(), RepMask::DOTS));
            return true;
        }
        if ui.button("ribbon").clicked() {
            actions.push(ObjectAction::HideRep(name.to_string(), RepMask::RIBBON));
            return true;
        }
        ui.separator();
        if ui.button("everything").clicked() {
            actions.push(ObjectAction::HideAll(name.to_string()));
            return true;
        }
        false
    }

    /// Render actions menu - returns true if item clicked
    fn render_actions_menu(ui: &mut Ui, name: &str, actions: &mut Vec<ObjectAction>) -> bool {
        if ui.button("zoom").clicked() {
            actions.push(ObjectAction::ZoomTo(name.to_string()));
            return true;
        }
        if ui.button("center").clicked() {
            actions.push(ObjectAction::CenterOn(name.to_string()));
            return true;
        }
        ui.separator();
        if ui.button("enable").clicked() {
            actions.push(ObjectAction::ToggleEnabled(name.to_string()));
            return true;
        }
        if ui.button("disable").clicked() {
            actions.push(ObjectAction::ToggleEnabled(name.to_string()));
            return true;
        }
        ui.separator();
        if ui.button("delete").clicked() {
            actions.push(ObjectAction::Delete(name.to_string()));
            return true;
        }
        false
    }

    /// Render color menu - returns true if item clicked
    fn render_color_menu(ui: &mut Ui, name: &str, actions: &mut Vec<ObjectAction>) -> bool {
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
                return true;
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
                return true;
            }
        }
        
        false
    }
}
