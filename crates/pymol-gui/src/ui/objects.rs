//! Object List Panel
//!
//! Displays loaded objects with action buttons (A/S/H/L/C) for each.
//! Uses a stateful menu system where buttons only update state,
//! and a single menu is rendered at the end to avoid overlapping popups.

use egui::{Color32, RichText, Ui, Vec2, Pos2, Id, Align2};
use pymol_mol::RepMask;
use pymol_scene::{ObjectRegistry, SelectionEntry, SelectionManager};

/// Minimum size for action buttons (A, S, H, L, C) to ensure consistent sizing
const BUTTON_MIN_SIZE: Vec2 = Vec2::splat(22.0);

/// All color constants used in the object list panel
mod colors {
    use egui::Color32;

    // Button colors (A, S, H, L, C)
    pub const ACTION: Color32 = Color32::from_rgb(150, 150, 255);   // A - blue
    pub const SHOW: Color32 = Color32::from_rgb(100, 200, 100);     // S - green
    pub const HIDE: Color32 = Color32::from_rgb(200, 100, 100);     // H - red
    pub const LABEL: Color32 = Color32::from_rgb(200, 200, 255);    // L - light blue
    pub const COLOR: Color32 = Color32::from_rgb(255, 200, 100);    // C - orange/yellow

    // Object state colors
    pub const ENABLED: Color32 = Color32::from_rgb(100, 255, 100);  // Green for enabled
    pub const DISABLED: Color32 = Color32::GRAY;                     // Gray for disabled

    // Selection colors
    pub const SELECTION: Color32 = Color32::from_rgb(255, 85, 255);      // Bright pink
    pub const SELECTION_HIDDEN: Color32 = Color32::from_rgb(128, 43, 128); // Dimmed pink
}

/// Available molecular representations for show/hide menus
const REPRESENTATIONS: &[(&str, RepMask)] = &[
    ("lines", RepMask::LINES),
    ("sticks", RepMask::STICKS),
    ("spheres", RepMask::SPHERES),
    ("cartoon", RepMask::CARTOON),
    ("surface", RepMask::SURFACE),
    ("mesh", RepMask::MESH),
    ("dots", RepMask::DOTS),
    ("ribbon", RepMask::RIBBON),
];

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

    // Check if this button's menu is currently active
    let is_this_menu_active = menu_state.active == target_menu;

    if response.clicked() {
        if is_this_menu_active {
            // Clicking the same button closes the menu
            menu_state.active = ActiveMenu::None;
        } else {
            // Open this menu (automatically closes any other)
            menu_state.active = target_menu;
            // Use right_bottom so popup opens to the left (better for right-side panel)
            menu_state.anchor_pos = response.rect.right_bottom();
        }
        // Request immediate repaint to show menu without delay
        ui.ctx().request_repaint();
    } else if is_this_menu_active {
        // Keep anchor_pos updated every frame for the active menu's button
        // This prevents "jumping" when scroll position or layout changes
        menu_state.anchor_pos = response.rect.right_bottom();
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
    ShowRep(String, RepMask),
    /// Hide specific representation
    HideRep(String, RepMask),
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


impl ObjectListPanel {
    /// Create a new ObjectListPanel
    pub fn new() -> Self {
        Self::default()
    }

    /// Render the right-aligned button bar (A, S, H, L, C) for an object or selection
    ///
    /// # Arguments
    /// * `ui` - The egui UI context
    /// * `target_name` - Name of the object or selection
    /// * `is_selection` - true if this is a selection, false if it's an object
    /// * `menu_button_clicked` - Mutable flag to track if any menu button was clicked
    fn render_button_bar(
        &mut self,
        ui: &mut Ui,
        target_name: &str,
        is_selection: bool,
        menu_button_clicked: &mut bool,
    ) {
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            // Buttons are added in reverse order (right-to-left): C, L, H, S, A

            // C - Color
            let color_tooltip = if is_selection { "Color selection" } else { "Color" };
            if menu_button(
                ui, "C", colors::COLOR, color_tooltip,
                &mut self.menu_state,
                ActiveMenu::Color { object: target_name.to_string() }
            ).clicked() {
                *menu_button_clicked = true;
            }

            // L - Labels (simple button, no menu)
            if ui
                .add(action_button("L", colors::LABEL))
                .on_hover_text("Labels")
                .clicked()
            {
                // Toggle labels (TODO: implement)
            }

            // H - Hide
            if menu_button(
                ui, "H", colors::HIDE, "Hide representations",
                &mut self.menu_state,
                ActiveMenu::Hide { object: target_name.to_string() }
            ).clicked() {
                *menu_button_clicked = true;
            }

            // S - Show
            if menu_button(
                ui, "S", colors::SHOW, "Show representations",
                &mut self.menu_state,
                ActiveMenu::Show { object: target_name.to_string() }
            ).clicked() {
                *menu_button_clicked = true;
            }

            // A - Actions (different menu for selections vs objects)
            let action_menu = if is_selection {
                ActiveMenu::SelectionActions { selection: target_name.to_string() }
            } else {
                ActiveMenu::Actions { object: target_name.to_string() }
            };
            if menu_button(
                ui, "A", colors::ACTION, "Actions",
                &mut self.menu_state,
                action_menu
            ).clicked() {
                *menu_button_clicked = true;
            }
        });
    }

    /// Render the panel header
    fn render_header(_ui: &mut Ui) {
    }

    /// Render the "all" row with aggregate actions for all objects
    fn render_all_row(
        ui: &mut Ui,
        registry: &ObjectRegistry,
        actions: &mut Vec<ObjectAction>,
    ) {
        let all_enabled = registry.names().all(|name| {
            registry.get(name).map(|o| o.is_enabled()).unwrap_or(false)
        });
        let has_objects = registry.names().next().is_some();

        ui.horizontal(|ui| {
            let all_color = if !has_objects {
                colors::DISABLED
            } else if all_enabled {
                colors::ENABLED
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

            // Right-align the buttons (simple action buttons, no menus)
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                // Buttons are added in reverse order (right-to-left)
                if ui
                    .add(action_button("C", colors::COLOR))
                    .on_hover_text("Color all objects")
                    .clicked()
                {
                    // Color picker would go here
                }

                if ui
                    .add(action_button("L", colors::LABEL))
                    .on_hover_text("Toggle labels")
                    .clicked()
                {
                    // Toggle labels for all
                }

                if ui
                    .add(action_button("H", colors::HIDE))
                    .on_hover_text("Hide all representations")
                    .clicked()
                {
                    for name in registry.names() {
                        actions.push(ObjectAction::HideAll(name.to_string()));
                    }
                }

                if ui
                    .add(action_button("S", colors::SHOW))
                    .on_hover_text("Show representations")
                    .clicked()
                {
                    for name in registry.names() {
                        actions.push(ObjectAction::ShowRep(name.to_string(), RepMask::LINES));
                    }
                }

                if ui
                    .add(action_button("A", colors::ACTION))
                    .on_hover_text("Actions menu")
                    .clicked()
                {
                    // Actions menu would go here
                }
            });
        });
    }

    /// Render a single object row
    fn render_object_row(
        &mut self,
        ui: &mut Ui,
        name: &str,
        is_enabled: bool,
        actions: &mut Vec<ObjectAction>,
        menu_button_clicked: &mut bool,
    ) {
        ui.push_id(name, |ui| {
            ui.horizontal(|ui| {
                let name_color = if is_enabled {
                    colors::ENABLED
                } else {
                    colors::DISABLED
                };

                let name_response = ui.add(
                    egui::Label::new(
                        RichText::new(name)
                            .family(egui::FontFamily::Monospace)
                            .color(name_color),
                    )
                    .sense(egui::Sense::click()),
                );

                if name_response.clicked() {
                    actions.push(ObjectAction::ToggleEnabled(name.to_string()));
                }
                if name_response.secondary_clicked() {
                    actions.push(ObjectAction::ZoomTo(name.to_string()));
                }

                self.render_button_bar(ui, name, false, menu_button_clicked);
            });
        });
    }

    /// Render a single selection row
    fn render_selection_row(
        &mut self,
        ui: &mut Ui,
        name: &str,
        entry: &SelectionEntry,
        actions: &mut Vec<ObjectAction>,
        menu_button_clicked: &mut bool,
    ) {
        ui.push_id(format!("sel_{}", name), |ui| {
            ui.horizontal(|ui| {
                let display_name = format!("({})", name);

                let color = if entry.visible {
                    colors::SELECTION
                } else {
                    colors::SELECTION_HIDDEN
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
                name_response.on_hover_text(format!(
                    "Selection: {} ({})\nClick to toggle visibility",
                    entry.expression, visibility_hint
                ));

                if clicked {
                    actions.push(ObjectAction::ToggleSelectionEnabled(name.to_string()));
                }

                self.render_button_bar(ui, name, true, menu_button_clicked);
            });
        });
    }

    /// Render the active menu and handle click-outside-to-close behavior
    fn render_menu_and_handle_close(
        &mut self,
        ui: &mut Ui,
        actions: &mut Vec<ObjectAction>,
        menu_button_clicked: bool,
    ) {
        if self.menu_state.active == ActiveMenu::None {
            return;
        }

        let menu_clicked = self.render_active_menu(ui, actions);

        // Close menu on click outside (but not if we clicked a menu button or menu item)
        if !menu_button_clicked && !menu_clicked {
            if ui.input(|i| i.pointer.any_click()) {
                self.menu_state.active = ActiveMenu::None;
            }
        }
    }

    /// Draw the object list panel
    ///
    /// # Arguments
    /// * `ui` - The egui UI context
    /// * `registry` - Object registry containing molecules and other objects
    /// * `selections` - Named selections (name -> expression)
    pub fn show(
        &mut self,
        ui: &mut Ui,
        registry: &ObjectRegistry,
        selections: &SelectionManager,
    ) -> Vec<ObjectAction> {
        let mut actions = Vec::new();
        let mut menu_button_clicked = false;

        ui.vertical(|ui| {
            Self::render_header(ui);
            Self::render_all_row(ui, registry, &mut actions);

            // Object rows
            for name in registry.names() {
                let is_enabled = registry
                    .get(name)
                    .map(|o| o.is_enabled())
                    .unwrap_or(false);
                self.render_object_row(ui, name, is_enabled, &mut actions, &mut menu_button_clicked);
            }

            // Selection rows
            if !selections.is_empty() {
                ui.separator();
                let mut sel_entries: Vec<_> = selections.iter().collect();
                sel_entries.sort_by(|a, b| a.0.cmp(b.0));

                for (sel_name, entry) in sel_entries {
                    self.render_selection_row(ui, sel_name, entry, &mut actions, &mut menu_button_clicked);
                }
            }
        });

        self.render_menu_and_handle_close(ui, &mut actions, menu_button_clicked);
        actions
    }

    /// Render the currently active menu - called ONCE at the end
    /// Returns true if a menu item was clicked
    fn render_active_menu(&mut self, ui: &mut Ui, actions: &mut Vec<ObjectAction>) -> bool {
        let mut item_clicked = false;

        // Request repaint to ensure the popup is fully rendered
        // This is critical for the first frame after opening
        ui.ctx().request_repaint();

        // Use a stable Area configuration:
        // - fixed_pos: anchor to button's right-bottom corner
        // - pivot RIGHT_TOP: popup's right-top aligns with anchor, so it expands LEFT
        // - movable(false): prevent dragging
        // - constrain(false): don't auto-adjust position (causes jumping)
        // - Order::Tooltip: highest rendering priority, ensures immediate visibility
        egui::Area::new(Id::new("object_list_menu_panel"))
            .order(egui::Order::Tooltip)
            .fixed_pos(self.menu_state.anchor_pos)
            .pivot(Align2::RIGHT_TOP)
            .movable(false)
            .constrain(false)
            .show(ui.ctx(), |ui| {
                // Use explicit solid colors to avoid semi-transparency on first frame
                let popup_frame = egui::Frame::popup(ui.style())
                    .fill(Color32::from_rgb(40, 40, 45))  // Solid dark background
                    .stroke(egui::Stroke::new(1.0, Color32::from_rgb(80, 80, 85)));

                popup_frame.show(ui, |ui| {
                        ui.set_min_width(100.0);

                        match &self.menu_state.active.clone() {
                            ActiveMenu::None => {}
                            ActiveMenu::Actions { object } => {
                                item_clicked = Self::render_actions_menu(ui, object, actions);
                            }
                            ActiveMenu::Show { object } => {
                                item_clicked = Self::render_representation_menu(ui, object, actions, true);
                            }
                            ActiveMenu::Hide { object } => {
                                item_clicked = Self::render_representation_menu(ui, object, actions, false);
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

    /// Render representation menu (show or hide) - returns true if item clicked
    ///
    /// # Arguments
    /// * `is_show` - true for show menu, false for hide menu
    fn render_representation_menu(
        ui: &mut Ui,
        name: &str,
        actions: &mut Vec<ObjectAction>,
        is_show: bool,
    ) -> bool {
        // Render buttons for each representation
        for (rep_name, rep_mask) in REPRESENTATIONS {
            if ui.button(*rep_name).clicked() {
                let action = if is_show {
                    ObjectAction::ShowRep(name.to_string(), *rep_mask)
                } else {
                    ObjectAction::HideRep(name.to_string(), *rep_mask)
                };
                actions.push(action);
                return true;
            }
        }

        // Separator and "all/everything" option
        ui.separator();
        let all_label = if is_show { "all" } else { "everything" };
        if ui.button(all_label).clicked() {
            let action = if is_show {
                ObjectAction::ShowAll(name.to_string())
            } else {
                ObjectAction::HideAll(name.to_string())
            };
            actions.push(action);
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
