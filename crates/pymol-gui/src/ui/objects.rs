//! Object List Panel
//!
//! Displays loaded objects with action buttons (A/S/H/L/C) for each.
//! Sends `AppMessage` directly to the `MessageBus` — no intermediate action enums.

use egui::{Color32, RichText, Ui, Vec2, Pos2, Id, Align2};
use pymol_mol::RepMask;
use pymol_scene::{ObjectRegistry, SelectionEntry, SelectionManager};

use crate::message::{AppMessage, MessageBus};

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

// =============================================================================
// UI State
// =============================================================================

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

/// Object list UI state (egui-specific: menu popups and anchor positions)
#[derive(Default)]
pub struct ObjectListUiState {
    /// Which menu is currently active
    pub active_menu: ActiveMenu,
    /// Position to anchor the menu popup
    pub anchor_pos: Pos2,
}

// =============================================================================
// View
// =============================================================================

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
    ui_state: &mut ObjectListUiState,
    target_menu: ActiveMenu,
) -> egui::Response {
    let button = egui::Button::new(RichText::new(text).color(color))
        .min_size(BUTTON_MIN_SIZE);
    let response = ui.add(button);

    let is_this_menu_active = ui_state.active_menu == target_menu;

    if response.clicked() {
        if is_this_menu_active {
            ui_state.active_menu = ActiveMenu::None;
        } else {
            ui_state.active_menu = target_menu;
            ui_state.anchor_pos = response.rect.right_bottom();
        }
        ui.ctx().request_repaint();
    } else if is_this_menu_active {
        ui_state.anchor_pos = response.rect.right_bottom();
    }

    if ui_state.active_menu == ActiveMenu::None {
        response.clone().on_hover_text(hover_text);
    }

    response
}

/// Object list panel
pub struct ObjectListPanel;

impl ObjectListPanel {
    /// Draw the object list panel, sending messages directly to the bus.
    pub fn show(
        ui: &mut Ui,
        ui_state: &mut ObjectListUiState,
        registry: &ObjectRegistry,
        selections: &SelectionManager,
        bus: &mut MessageBus,
    ) {
        let mut menu_button_clicked = false;

        ui.vertical(|ui| {
            Self::render_all_row(ui, registry, bus);

            // Object rows
            for name in registry.names() {
                let is_enabled = registry
                    .get(name)
                    .map(|o| o.is_enabled())
                    .unwrap_or(false);
                Self::render_object_row(ui, ui_state, name, is_enabled, bus, &mut menu_button_clicked);
            }

            // Selection rows
            if !selections.is_empty() {
                ui.separator();
                let mut sel_entries: Vec<_> = selections.iter().collect();
                sel_entries.sort_by(|a, b| a.0.cmp(b.0));

                for (sel_name, entry) in sel_entries {
                    Self::render_selection_row(ui, ui_state, sel_name, entry, bus, &mut menu_button_clicked);
                }
            }
        });

        Self::render_menu_and_handle_close(ui, ui_state, bus, menu_button_clicked);
    }

    /// Render the "all" row with aggregate actions for all objects
    fn render_all_row(
        ui: &mut Ui,
        registry: &ObjectRegistry,
        bus: &mut MessageBus,
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
                    bus.send(AppMessage::DisableAllObjects);
                } else {
                    bus.send(AppMessage::EnableAllObjects);
                }
            }

            all_response.on_hover_text(if all_enabled { "Click to disable all" } else { "Click to enable all" });

            // Right-align the buttons (simple action buttons, no menus)
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
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
                        bus.send(AppMessage::HideAllRepresentations(name.to_string()));
                    }
                }

                if ui
                    .add(action_button("S", colors::SHOW))
                    .on_hover_text("Show representations")
                    .clicked()
                {
                    for name in registry.names() {
                        bus.send(AppMessage::ShowRepresentation {
                            object: name.to_string(),
                            rep: RepMask::LINES,
                        });
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

    /// Render the button bar (A, S, H, L, C) for an object or selection
    fn render_button_bar(
        ui: &mut Ui,
        ui_state: &mut ObjectListUiState,
        target_name: &str,
        is_selection: bool,
        menu_button_clicked: &mut bool,
    ) {
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            // C - Color
            let color_tooltip = if is_selection { "Color selection" } else { "Color" };
            if menu_button(
                ui, "C", colors::COLOR, color_tooltip,
                ui_state,
                ActiveMenu::Color { object: target_name.to_string() }
            ).clicked() {
                *menu_button_clicked = true;
            }

            // L - Labels
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
                ui_state,
                ActiveMenu::Hide { object: target_name.to_string() }
            ).clicked() {
                *menu_button_clicked = true;
            }

            // S - Show
            if menu_button(
                ui, "S", colors::SHOW, "Show representations",
                ui_state,
                ActiveMenu::Show { object: target_name.to_string() }
            ).clicked() {
                *menu_button_clicked = true;
            }

            // A - Actions
            let action_menu = if is_selection {
                ActiveMenu::SelectionActions { selection: target_name.to_string() }
            } else {
                ActiveMenu::Actions { object: target_name.to_string() }
            };
            if menu_button(
                ui, "A", colors::ACTION, "Actions",
                ui_state,
                action_menu
            ).clicked() {
                *menu_button_clicked = true;
            }
        });
    }

    /// Render a single object row
    fn render_object_row(
        ui: &mut Ui,
        ui_state: &mut ObjectListUiState,
        name: &str,
        is_enabled: bool,
        bus: &mut MessageBus,
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
                    bus.send(AppMessage::ToggleObject(name.to_string()));
                }
                if name_response.secondary_clicked() {
                    bus.send(AppMessage::ZoomTo(name.to_string()));
                }

                Self::render_button_bar(ui, ui_state, name, false, menu_button_clicked);
            });
        });
    }

    /// Render a single selection row
    fn render_selection_row(
        ui: &mut Ui,
        ui_state: &mut ObjectListUiState,
        name: &str,
        entry: &SelectionEntry,
        bus: &mut MessageBus,
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
                    bus.send(AppMessage::ToggleSelectionVisibility(name.to_string()));
                }

                Self::render_button_bar(ui, ui_state, name, true, menu_button_clicked);
            });
        });
    }

    /// Render the active menu and handle click-outside-to-close behavior
    fn render_menu_and_handle_close(
        ui: &mut Ui,
        ui_state: &mut ObjectListUiState,
        bus: &mut MessageBus,
        menu_button_clicked: bool,
    ) {
        if ui_state.active_menu == ActiveMenu::None {
            return;
        }

        let menu_clicked = Self::render_active_menu(ui, ui_state, bus);

        if !menu_button_clicked && !menu_clicked {
            if ui.input(|i| i.pointer.any_click()) {
                ui_state.active_menu = ActiveMenu::None;
            }
        }
    }

    /// Render the currently active menu - called ONCE at the end
    /// Returns true if a menu item was clicked
    fn render_active_menu(ui: &mut Ui, ui_state: &mut ObjectListUiState, bus: &mut MessageBus) -> bool {
        let mut item_clicked = false;

        ui.ctx().request_repaint();

        egui::Area::new(Id::new("object_list_menu_panel"))
            .order(egui::Order::Tooltip)
            .fixed_pos(ui_state.anchor_pos)
            .pivot(Align2::RIGHT_TOP)
            .movable(false)
            .constrain(false)
            .show(ui.ctx(), |ui| {
                let popup_frame = egui::Frame::popup(ui.style())
                    .fill(Color32::from_rgb(40, 40, 45))
                    .stroke(egui::Stroke::new(1.0, Color32::from_rgb(80, 80, 85)));

                popup_frame.show(ui, |ui| {
                        ui.set_min_width(100.0);

                        match &ui_state.active_menu.clone() {
                            ActiveMenu::None => {}
                            ActiveMenu::Actions { object } => {
                                item_clicked = Self::render_actions_menu(ui, object, bus);
                            }
                            ActiveMenu::Show { object } => {
                                item_clicked = Self::render_representation_menu(ui, object, bus, true);
                            }
                            ActiveMenu::Hide { object } => {
                                item_clicked = Self::render_representation_menu(ui, object, bus, false);
                            }
                            ActiveMenu::Color { object } => {
                                item_clicked = Self::render_color_menu(ui, object, bus);
                            }
                            ActiveMenu::SelectionActions { selection } => {
                                item_clicked = Self::render_selection_actions_menu(ui, selection, bus);
                            }
                        }
                    });
            });

        if item_clicked {
            ui_state.active_menu = ActiveMenu::None;
        }

        item_clicked
    }

    /// Render actions menu for selections
    fn render_selection_actions_menu(ui: &mut Ui, name: &str, bus: &mut MessageBus) -> bool {
        if ui.button("delete").clicked() {
            bus.send(AppMessage::DeleteSelection(name.to_string()));
            return true;
        }
        false
    }

    /// Render representation menu (show or hide)
    fn render_representation_menu(
        ui: &mut Ui,
        name: &str,
        bus: &mut MessageBus,
        is_show: bool,
    ) -> bool {
        for (rep_name, rep_mask) in REPRESENTATIONS {
            if ui.button(*rep_name).clicked() {
                if is_show {
                    bus.send(AppMessage::ShowRepresentation {
                        object: name.to_string(),
                        rep: *rep_mask,
                    });
                } else {
                    bus.send(AppMessage::HideRepresentation {
                        object: name.to_string(),
                        rep: *rep_mask,
                    });
                }
                return true;
            }
        }

        ui.separator();
        let all_label = if is_show { "all" } else { "everything" };
        if ui.button(all_label).clicked() {
            if is_show {
                bus.send(AppMessage::ShowAllRepresentations(name.to_string()));
            } else {
                bus.send(AppMessage::HideAllRepresentations(name.to_string()));
            }
            return true;
        }

        false
    }

    /// Render actions menu
    fn render_actions_menu(ui: &mut Ui, name: &str, bus: &mut MessageBus) -> bool {
        if ui.button("zoom").clicked() {
            bus.send(AppMessage::ZoomTo(name.to_string()));
            return true;
        }
        if ui.button("center").clicked() {
            bus.send(AppMessage::CenterOn(name.to_string()));
            return true;
        }
        ui.separator();
        if ui.button("enable").clicked() {
            bus.send(AppMessage::ToggleObject(name.to_string()));
            return true;
        }
        if ui.button("disable").clicked() {
            bus.send(AppMessage::ToggleObject(name.to_string()));
            return true;
        }
        ui.separator();
        if ui.button("delete").clicked() {
            bus.send(AppMessage::DeleteObject(name.to_string()));
            return true;
        }
        false
    }

    /// Render color menu
    fn render_color_menu(ui: &mut Ui, name: &str, bus: &mut MessageBus) -> bool {
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
                bus.send(AppMessage::SetColor {
                    object: name.to_string(),
                    color: color_name.to_string(),
                });
                return true;
            }
        }

        ui.separator();

        let by_colors = ["by element", "by chain", "by ss", "by residue"];
        for by_color in by_colors {
            if ui.button(by_color).clicked() {
                bus.send(AppMessage::SetColor {
                    object: name.to_string(),
                    color: by_color.replace(' ', "_"),
                });
                return true;
            }
        }

        false
    }
}
