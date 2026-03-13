//! Object List Panel
//!
//! Displays loaded objects with action buttons (A/S/H/L/C) for each.
//! Sends commands via `MessageBus::execute_command()`.

use egui::{Align2, Color32, Id, Pos2, Rect, RichText, Ui, Vec2};
use pymol_scene::{ObjectRegistry, SelectionEntry, SelectionManager};

use pymol_framework::message::MessageBus;

/// Minimum size for action buttons (A, S, H, L, C) to ensure consistent sizing
const BUTTON_MIN_SIZE: Vec2 = Vec2::splat(16.0);

/// All color constants used in the object list panel
mod colors {
    use egui::Color32;

    // Button colors (A, S, H, L, C)
    pub const ACTION: Color32 = Color32::from_rgb(150, 150, 255); // A - blue
    pub const SHOW: Color32 = Color32::from_rgb(100, 200, 100); // S - green
    pub const HIDE: Color32 = Color32::from_rgb(200, 100, 100); // H - red
    pub const LABEL: Color32 = Color32::from_rgb(200, 200, 255); // L - light blue
    pub const COLOR: Color32 = Color32::from_rgb(255, 200, 100); // C - orange/yellow

    // Object state colors
    pub const ENABLED: Color32 = Color32::from_rgb(100, 255, 100); // Green for enabled
    pub const DISABLED: Color32 = Color32::GRAY; // Gray for disabled

    // Selection colors
    pub const SELECTION: Color32 = Color32::from_rgb(255, 85, 255); // Bright pink
    pub const SELECTION_HIDDEN: Color32 = Color32::from_rgb(128, 43, 128); // Dimmed pink
}

/// Named colors for color menus, with display swatches
const NAMED_COLORS: &[(&str, Color32)] = &[
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

/// Representation groups for Show/Hide menus (separated by group)
const REP_GROUPS: &[&[&str]] = &[
    &["lines", "sticks"],
    &["cartoon", "ribbon"],
    &["surface", "mesh"],
    &["spheres", "dots"],
];

/// Representations available for "by rep" color submenu
const COLOR_REPS: &[(&str, &str)] = &[
    ("cartoon", "cartoon_color"),
    ("ribbon", "ribbon_color"),
    ("sticks", "stick_color"),
    ("lines", "line_color"),
    ("surface", "surface_color"),
    ("mesh", "mesh_color"),
    ("spheres", "sphere_color"),
    ("dots", "dot_color"),
];

// =============================================================================
// Command formatting helpers
// =============================================================================

/// Format a command with an optional positional argument.
/// `cmd_arg("zoom", "")` → `"zoom"`, `cmd_arg("zoom", "obj")` → `"zoom obj"`
fn cmd_arg(cmd: &str, arg: &str) -> String {
    if arg.is_empty() {
        cmd.to_string()
    } else {
        format!("{} {}", cmd, arg)
    }
}

/// Format a command with an optional comma-separated selection.
/// `cmd_sel("show lines", "")` → `"show lines"`, `cmd_sel("show lines", "obj")` → `"show lines, obj"`
fn cmd_sel(cmd: &str, sel: &str) -> String {
    if sel.is_empty() {
        cmd.to_string()
    } else {
        format!("{}, {}", cmd, sel)
    }
}

/// CA atom selection, optionally scoped to an object.
/// `ca_sel("")` → `"name CA"`, `ca_sel("obj")` → `"name CA and obj"`
fn ca_sel(name: &str) -> String {
    if name.is_empty() {
        "name CA".to_string()
    } else {
        format!("name CA and {}", name)
    }
}

/// Atom selection: the object name, or "all" if empty.
fn atom_sel(name: &str) -> &str {
    if name.is_empty() {
        "all"
    } else {
        name
    }
}

// =============================================================================
// UI State
// =============================================================================

/// What menu is currently showing
#[derive(Clone, PartialEq, Debug, Default)]
pub enum ActiveMenu {
    /// No menu open
    #[default]
    None,
    /// Actions menu for an object (zoom, orient, center, origin)
    Actions { object: String, is_group: bool },
    /// Show representations menu for an object
    Show { object: String },
    /// Hide representations menu for an object
    Hide { object: String },
    /// Label menu for an object
    Label { object: String },
    /// Color menu for an object
    Color { object: String },
    /// Actions menu for a selection (just delete)
    SelectionActions { selection: String },
}

/// Submenu state for nested menus (Show > as, Color > by rep)
#[derive(Clone, PartialEq, Debug, Default)]
pub enum SubMenuState {
    #[default]
    None,
    /// Show > as (show_as submenu)
    ShowAs,
    /// Color > by rep (list of representations)
    ColorByRep,
    /// Color > by rep > specific rep (color picker for that rep)
    ColorByRepColors { rep: String, setting: String },
}


/// Object list UI state (egui-specific: menu popups and anchor positions)
pub struct ObjectListUiState {
    /// Which menu is currently active
    pub active_menu: ActiveMenu,
    /// Position to anchor the menu popup
    pub anchor_pos: Pos2,
    /// Current submenu state
    pub submenu: SubMenuState,
    /// Position to anchor submenu
    pub submenu_anchor: Pos2,
    /// Position to anchor sub-submenu (3rd level)
    pub sub_submenu_anchor: Pos2,
    /// Rects of all rendered menu areas this frame (for click-outside detection)
    menu_rects: Vec<Rect>,
}

impl Default for ObjectListUiState {
    fn default() -> Self {
        Self {
            active_menu: ActiveMenu::None,
            anchor_pos: Pos2::ZERO,
            submenu: SubMenuState::None,
            submenu_anchor: Pos2::ZERO,
            sub_submenu_anchor: Pos2::ZERO,
            menu_rects: Vec::new(),
        }
    }
}

// =============================================================================
// View
// =============================================================================

/// Menu button that only updates state - NEVER renders popups directly
fn menu_button(
    ui: &mut Ui,
    text: &str,
    color: Color32,
    hover_text: &str,
    ui_state: &mut ObjectListUiState,
    target_menu: ActiveMenu,
) -> egui::Response {
    let button = egui::Button::new(RichText::new(text).color(color)).min_size(BUTTON_MIN_SIZE);
    let response = ui.add(button);

    let is_this_menu_active = ui_state.active_menu == target_menu;

    if response.clicked() {
        if is_this_menu_active {
            ui_state.active_menu = ActiveMenu::None;
        } else {
            ui_state.active_menu = target_menu;
            ui_state.anchor_pos = response.rect.right_bottom();
            ui_state.submenu = SubMenuState::None;
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

/// Popup frame styling shared by all menu panels
fn popup_frame(ui: &Ui) -> egui::Frame {
    egui::Frame::popup(ui.style())
        .fill(Color32::from_rgb(40, 40, 45))
        .stroke(egui::Stroke::new(1.0, Color32::from_rgb(80, 80, 85)))
}

/// A button with a ">" arrow indicating a submenu
fn submenu_button(ui: &mut Ui, label: &str) -> egui::Response {
    ui.add(
        egui::Button::new(RichText::new(format!("{label}  >")))
            .min_size(Vec2::new(100.0, 0.0)),
    )
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
            Self::render_all_row(ui, ui_state, registry, bus, &mut menu_button_clicked);

            // Object rows (hierarchical: only top-level, groups render children)
            for name in registry.top_level_names() {
                if registry.get_group(name).is_some() {
                    Self::render_group_row(
                        ui,
                        ui_state,
                        registry,
                        name,
                        0,
                        bus,
                        &mut menu_button_clicked,
                    );
                } else {
                    let is_enabled = registry
                        .get(name)
                        .map(|o| o.is_enabled())
                        .unwrap_or(false);
                    Self::render_object_row(
                        ui,
                        ui_state,
                        name,
                        is_enabled,
                        bus,
                        &mut menu_button_clicked,
                    );
                }
            }

            // Selection rows
            if !selections.is_empty() {
                ui.separator();
                let mut sel_entries: Vec<_> = selections.iter().collect();
                sel_entries.sort_by(|a, b| a.0.cmp(b.0));

                for (sel_name, entry) in sel_entries {
                    Self::render_selection_row(
                        ui,
                        ui_state,
                        sel_name,
                        entry,
                        bus,
                        &mut menu_button_clicked,
                    );
                }
            }
        });

        Self::render_menu_and_handle_close(ui, ui_state, bus, menu_button_clicked);
    }

    /// Render the "all" row with aggregate actions for all objects
    fn render_all_row(
        ui: &mut Ui,
        ui_state: &mut ObjectListUiState,
        registry: &ObjectRegistry,
        bus: &mut MessageBus,
        menu_button_clicked: &mut bool,
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
                    bus.execute_command("disable");
                } else {
                    bus.execute_command("enable");
                }
            }

            all_response.on_hover_text(if all_enabled {
                "Click to disable all"
            } else {
                "Click to enable all"
            });

            // "all" uses the same button bar with empty name = commands apply to everything
            Self::render_button_bar(ui, ui_state, "", false, false, menu_button_clicked);
        });
    }

    /// Render the button bar (A, S, H, L, C) for an object or selection.
    /// Pass empty `target_name` for the "all" row (commands without object scope).
    fn render_button_bar(
        ui: &mut Ui,
        ui_state: &mut ObjectListUiState,
        target_name: &str,
        is_selection: bool,
        is_group: bool,
        menu_button_clicked: &mut bool,
    ) {
        ui.with_layout(
            egui::Layout::right_to_left(egui::Align::Center),
            |ui| {
                ui.spacing_mut().item_spacing.x = 2.0;
                // C - Color
                let color_tooltip = if is_selection {
                    "Color selection"
                } else {
                    "Color"
                };
                if menu_button(
                    ui,
                    "C",
                    colors::COLOR,
                    color_tooltip,
                    ui_state,
                    ActiveMenu::Color {
                        object: target_name.to_string(),
                    },
                )
                .clicked()
                {
                    *menu_button_clicked = true;
                }

                // L - Labels
                if menu_button(
                    ui,
                    "L",
                    colors::LABEL,
                    "Labels",
                    ui_state,
                    ActiveMenu::Label {
                        object: target_name.to_string(),
                    },
                )
                .clicked()
                {
                    *menu_button_clicked = true;
                }

                // H - Hide
                if menu_button(
                    ui,
                    "H",
                    colors::HIDE,
                    "Hide representations",
                    ui_state,
                    ActiveMenu::Hide {
                        object: target_name.to_string(),
                    },
                )
                .clicked()
                {
                    *menu_button_clicked = true;
                }

                // S - Show
                if menu_button(
                    ui,
                    "S",
                    colors::SHOW,
                    "Show representations",
                    ui_state,
                    ActiveMenu::Show {
                        object: target_name.to_string(),
                    },
                )
                .clicked()
                {
                    *menu_button_clicked = true;
                }

                // A - Actions
                let action_menu = if is_selection {
                    ActiveMenu::SelectionActions {
                        selection: target_name.to_string(),
                    }
                } else {
                    ActiveMenu::Actions {
                        object: target_name.to_string(),
                        is_group,
                    }
                };
                if menu_button(ui, "A", colors::ACTION, "Actions", ui_state, action_menu)
                    .clicked()
                {
                    *menu_button_clicked = true;
                }
            },
        );
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
                    bus.execute_command(format!("toggle {}", name));
                }
                if name_response.secondary_clicked() {
                    bus.execute_command(format!("zoom {}", name));
                }

                Self::render_button_bar(ui, ui_state, name, false, false, menu_button_clicked);
            });
        });
    }

    /// Render a group row with expand/collapse and recursive children
    fn render_group_row(
        ui: &mut Ui,
        ui_state: &mut ObjectListUiState,
        registry: &ObjectRegistry,
        name: &str,
        indent_level: usize,
        bus: &mut MessageBus,
        menu_button_clicked: &mut bool,
    ) {
        let is_open = registry
            .get_group(name)
            .map(|g| g.is_open())
            .unwrap_or(false);
        let is_enabled = registry
            .get(name)
            .map(|o| o.is_enabled())
            .unwrap_or(false);
        let children: Vec<String> = registry
            .get_group(name)
            .map(|g| g.children().to_vec())
            .unwrap_or_default();

        ui.push_id(name, |ui| {
            ui.horizontal(|ui| {
                // Indentation for nested groups
                if indent_level > 0 {
                    ui.add_space(indent_level as f32 * 16.0);
                }

                // Group name
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
                    bus.execute_command(format!("toggle {}", name));
                }
                if name_response.secondary_clicked() {
                    bus.execute_command(format!("zoom {}", name));
                }

                // Expand/collapse triangle (after name)
                let (rect, response) = ui.allocate_exact_size(
                    Vec2::splat(ui.spacing().icon_width),
                    egui::Sense::click(),
                );
                if response.clicked() {
                    bus.execute_command(format!("group {}, action=toggle", name));
                }
                {
                    let color = if response.hovered() {
                        Color32::WHITE
                    } else {
                        Color32::GRAY
                    };
                    let center = rect.center();
                    let half = rect.width() * 0.25;
                    let points = if is_open {
                        // Down-pointing triangle
                        vec![
                            Pos2::new(center.x - half, center.y - half * 0.5),
                            Pos2::new(center.x + half, center.y - half * 0.5),
                            Pos2::new(center.x, center.y + half * 0.5),
                        ]
                    } else {
                        // Right-pointing triangle
                        vec![
                            Pos2::new(center.x - half * 0.5, center.y - half),
                            Pos2::new(center.x + half * 0.5, center.y),
                            Pos2::new(center.x - half * 0.5, center.y + half),
                        ]
                    };
                    ui.painter()
                        .add(egui::Shape::convex_polygon(points, color, egui::Stroke::NONE));
                }

                Self::render_button_bar(ui, ui_state, name, false, true, menu_button_clicked);
            });

            // Render children if open
            if is_open {
                for child_name in &children {
                    if registry.get_group(child_name).is_some() {
                        Self::render_group_row(
                            ui,
                            ui_state,
                            registry,
                            child_name,
                            indent_level + 1,
                            bus,
                            menu_button_clicked,
                        );
                    } else {
                        let child_enabled = registry
                            .get(child_name)
                            .map(|o| o.is_enabled())
                            .unwrap_or(false);
                        ui.horizontal(|ui| {
                            ui.add_space((indent_level + 1) as f32 * 16.0);
                            let name_color = if child_enabled {
                                colors::ENABLED
                            } else {
                                colors::DISABLED
                            };

                            let name_response = ui.add(
                                egui::Label::new(
                                    RichText::new(child_name.as_str())
                                        .family(egui::FontFamily::Monospace)
                                        .color(name_color),
                                )
                                .sense(egui::Sense::click()),
                            );

                            if name_response.clicked() {
                                bus.execute_command(format!("toggle {}", child_name));
                            }
                            if name_response.secondary_clicked() {
                                bus.execute_command(format!("zoom {}", child_name));
                            }

                            Self::render_button_bar(
                                ui,
                                ui_state,
                                child_name,
                                false,
                                false,
                                menu_button_clicked,
                            );
                        });
                    }
                }
            }
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
                    bus.execute_command(format!("toggle {}", name));
                }

                Self::render_button_bar(ui, ui_state, name, true, false, menu_button_clicked);
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

        // Clear rects from previous frame, will be populated by render_active_menu
        ui_state.menu_rects.clear();

        let item_clicked = Self::render_active_menu(ui, ui_state, bus);

        if item_clicked {
            ui_state.active_menu = ActiveMenu::None;
            ui_state.submenu = SubMenuState::None;
            return;
        }

        // Close on click outside all menu areas
        if !menu_button_clicked
            && ui.input(|i| i.pointer.any_click())
        {
            let pointer_in_menu = ui
                .input(|i| i.pointer.interact_pos())
                .is_some_and(|pos| {
                    ui_state.menu_rects.iter().any(|r| r.contains(pos))
                });
            if !pointer_in_menu {
                ui_state.active_menu = ActiveMenu::None;
                ui_state.submenu = SubMenuState::None;
            }
        }
    }

    /// Render the currently active menu - called ONCE at the end.
    /// Returns true if a final action item was clicked.
    fn render_active_menu(
        ui: &mut Ui,
        ui_state: &mut ObjectListUiState,
        bus: &mut MessageBus,
    ) -> bool {
        let mut item_clicked = false;

        ui.ctx().request_repaint();

        // -- Main menu panel --
        let area_resp = egui::Area::new(Id::new("object_list_menu_panel"))
            .order(egui::Order::Tooltip)
            .fixed_pos(ui_state.anchor_pos)
            .pivot(Align2::RIGHT_TOP)
            .movable(false)
            .constrain(false)
            .show(ui.ctx(), |ui| {
                popup_frame(ui).show(ui, |ui| {
                    ui.set_min_width(100.0);

                    match &ui_state.active_menu.clone() {
                        ActiveMenu::None => {}
                        ActiveMenu::Actions { object, is_group } => {
                            item_clicked = Self::render_actions_menu(ui, object, *is_group, bus);
                        }
                        ActiveMenu::Show { object } => {
                            item_clicked =
                                Self::render_show_menu(ui, ui_state, object, bus);
                        }
                        ActiveMenu::Hide { object } => {
                            item_clicked = Self::render_hide_menu(ui, object, bus);
                        }
                        ActiveMenu::Label { object } => {
                            item_clicked = Self::render_label_menu(ui, object, bus);
                        }
                        ActiveMenu::Color { object } => {
                            item_clicked =
                                Self::render_color_menu(ui, ui_state, object, bus);
                        }
                        ActiveMenu::SelectionActions { selection } => {
                            item_clicked =
                                Self::render_selection_actions_menu(ui, selection, bus);
                        }
                    }
                });
            });
        ui_state.menu_rects.push(area_resp.response.rect);

        // -- Submenus (rendered as separate Areas) --
        if !item_clicked {
            item_clicked = Self::render_submenus(ui, ui_state, bus);
        }

        item_clicked
    }

    // =========================================================================
    // Menu renderers
    // =========================================================================

    /// Actions menu: zoom, orient, center, origin (+ group actions for groups)
    fn render_actions_menu(ui: &mut Ui, name: &str, is_group: bool, bus: &mut MessageBus) -> bool {
        if ui.button("zoom").clicked() {
            bus.execute_command(cmd_arg("zoom", name));
            return true;
        }
        if ui.button("orient").clicked() {
            bus.execute_command(cmd_arg("orient", name));
            return true;
        }
        if ui.button("center").clicked() {
            bus.execute_command(cmd_arg("center", name));
            return true;
        }
        if ui.button("origin").clicked() {
            bus.execute_command(cmd_arg("origin", name));
            return true;
        }

        if is_group {
            ui.separator();
            if ui.button("open").clicked() {
                bus.execute_command(format!("group {}, action=open", name));
                return true;
            }
            if ui.button("close").clicked() {
                bus.execute_command(format!("group {}, action=close", name));
                return true;
            }
            if ui.button("ungroup").clicked() {
                bus.execute_command(format!("ungroup {}", name));
                return true;
            }
        }

        false
    }

    /// Show menu: "as" submenu trigger, then grouped representations
    fn render_show_menu(
        ui: &mut Ui,
        ui_state: &mut ObjectListUiState,
        name: &str,
        bus: &mut MessageBus,
    ) -> bool {
        // "as" submenu trigger
        let as_resp = submenu_button(ui, "as");
        if as_resp.clicked() {
            if ui_state.submenu == SubMenuState::ShowAs {
                ui_state.submenu = SubMenuState::None;
            } else {
                ui_state.submenu = SubMenuState::ShowAs;
            }
        }
        if matches!(ui_state.submenu, SubMenuState::ShowAs) {
            ui_state.submenu_anchor = as_resp.rect.left_center();
        }

        ui.separator();

        Self::render_rep_groups(ui, name, bus, "show")
    }

    /// Hide menu: "everything" first, then grouped representations
    fn render_hide_menu(ui: &mut Ui, name: &str, bus: &mut MessageBus) -> bool {
        if ui.button("everything").clicked() {
            bus.execute_command(cmd_sel("hide everything", name));
            return true;
        }
        ui.separator();

        Self::render_rep_groups(ui, name, bus, "hide")
    }

    /// Label menu: hide, separator, label options
    fn render_label_menu(ui: &mut Ui, name: &str, bus: &mut MessageBus) -> bool {
        if ui.button("hide").clicked() {
            bus.execute_command(cmd_sel("hide labels", name));
            return true;
        }
        ui.separator();

        let ca = ca_sel(name);
        let sel = atom_sel(name);

        if ui.button("residue (3)").clicked() {
            bus.execute_command(format!("label {}, resn", ca));
            return true;
        }
        if ui.button("residue (1)").clicked() {
            bus.execute_command(format!("label {}, oneletter", ca));
            return true;
        }
        if ui.button("position").clicked() {
            bus.execute_command(format!("label {}, resi", ca));
            return true;
        }
        ui.separator();

        if ui.button("atom name").clicked() {
            bus.execute_command(format!("label {}, name", sel));
            return true;
        }
        if ui.button("element").clicked() {
            bus.execute_command(format!("label {}, elem", sel));
            return true;
        }
        ui.separator();

        if ui.button("b-factor").clicked() {
            bus.execute_command(format!("label {}, b", sel));
            return true;
        }
        if ui.button("vdw radius").clicked() {
            bus.execute_command(format!("label {}, vdw", sel));
            return true;
        }

        false
    }

    /// Color menu: "by rep" submenu trigger, named colors, color schemes
    fn render_color_menu(
        ui: &mut Ui,
        ui_state: &mut ObjectListUiState,
        name: &str,
        bus: &mut MessageBus,
    ) -> bool {
        // "by rep" submenu trigger
        let by_rep_resp = submenu_button(ui, "by rep");
        if by_rep_resp.clicked() {
            if matches!(
                ui_state.submenu,
                SubMenuState::ColorByRep | SubMenuState::ColorByRepColors { .. }
            ) {
                ui_state.submenu = SubMenuState::None;
            } else {
                ui_state.submenu = SubMenuState::ColorByRep;
            }
        }
        if matches!(
            ui_state.submenu,
            SubMenuState::ColorByRep | SubMenuState::ColorByRepColors { .. }
        ) {
            ui_state.submenu_anchor = by_rep_resp.rect.left_center();
        }

        ui.separator();

        // Named colors
        for &(color_name, color) in NAMED_COLORS {
            if ui
                .add(egui::Button::new(RichText::new(color_name).color(color)))
                .clicked()
            {
                bus.execute_command(cmd_sel(&format!("color {}", color_name), name));
                return true;
            }
        }

        ui.separator();

        // Color schemes
        let schemes = [
            ("by element", "atomic"),
            ("by chain", "chain"),
            ("by ss", "ss"),
            ("by b-factor", "b"),
            ("by residue", "residue"),
        ];
        for (label, scheme) in schemes {
            if ui.button(label).clicked() {
                bus.execute_command(cmd_sel(&format!("color {}", scheme), name));
                return true;
            }
        }

        false
    }

    /// Actions menu for selections
    fn render_selection_actions_menu(ui: &mut Ui, name: &str, bus: &mut MessageBus) -> bool {
        if ui.button("delete").clicked() {
            bus.execute_command(format!("delete {}", name));
            return true;
        }
        false
    }

    // =========================================================================
    // Helpers
    // =========================================================================

    /// Render representation groups with separators between groups
    fn render_rep_groups(
        ui: &mut Ui,
        name: &str,
        bus: &mut MessageBus,
        action: &str,
    ) -> bool {
        for (i, group) in REP_GROUPS.iter().enumerate() {
            if i > 0 {
                ui.separator();
            }
            for rep in *group {
                if ui.button(*rep).clicked() {
                    bus.execute_command(cmd_sel(&format!("{} {}", action, rep), name));
                    return true;
                }
            }
        }
        false
    }

    // =========================================================================
    // Submenus
    // =========================================================================

    /// Render active submenus. Returns true if a final action item was clicked.
    fn render_submenus(
        ui: &mut Ui,
        ui_state: &mut ObjectListUiState,
        bus: &mut MessageBus,
    ) -> bool {
        match ui_state.submenu.clone() {
            SubMenuState::None => false,
            SubMenuState::ShowAs => {
                let object = match &ui_state.active_menu {
                    ActiveMenu::Show { object } => object.clone(),
                    _ => return false,
                };
                Self::render_submenu_panel(
                    ui,
                    ui_state,
                    "submenu_show_as",
                    ui_state.submenu_anchor,
                    |ui| {
                        for (i, group) in REP_GROUPS.iter().enumerate() {
                            if i > 0 {
                                ui.separator();
                            }
                            for rep in *group {
                                if ui.button(*rep).clicked() {
                                    bus.execute_command(cmd_sel(
                                        &format!("show_as {}", rep),
                                        &object,
                                    ));
                                    return true;
                                }
                            }
                        }
                        false
                    },
                )
            }
            SubMenuState::ColorByRep | SubMenuState::ColorByRepColors { .. } => {
                let object = match &ui_state.active_menu {
                    ActiveMenu::Color { object } => object.clone(),
                    _ => return false,
                };
                let mut sub_sub_anchor = ui_state.sub_submenu_anchor;
                let mut new_sub = ui_state.submenu.clone();

                let clicked = Self::render_submenu_panel(
                    ui,
                    ui_state,
                    "submenu_color_by_rep",
                    ui_state.submenu_anchor,
                    |ui| {
                        for &(rep_name, setting_name) in COLOR_REPS {
                            let resp = submenu_button(ui, rep_name);
                            if resp.clicked() {
                                let target = SubMenuState::ColorByRepColors {
                                    rep: rep_name.to_string(),
                                    setting: setting_name.to_string(),
                                };
                                if new_sub == target {
                                    new_sub = SubMenuState::ColorByRep;
                                } else {
                                    new_sub = target;
                                    sub_sub_anchor = resp.rect.left_center();
                                }
                            }
                            // Keep updating anchor for the active rep
                            if matches!(&new_sub, SubMenuState::ColorByRepColors { rep, .. } if rep == rep_name)
                            {
                                sub_sub_anchor = resp.rect.left_center();
                            }
                        }
                        false
                    },
                );

                // Update sub-submenu state
                if new_sub != ui_state.submenu {
                    ui_state.submenu = new_sub;
                    ui_state.sub_submenu_anchor = sub_sub_anchor;
                }

                if clicked {
                    return true;
                }

                // Render 3rd-level color picker if active
                if let SubMenuState::ColorByRepColors { setting, .. } = &ui_state.submenu {
                    let setting = setting.clone();
                    Self::render_submenu_panel(
                        ui,
                        ui_state,
                        "submenu_color_rep_colors",
                        ui_state.sub_submenu_anchor,
                        |ui| {
                            for &(color_name, color) in NAMED_COLORS {
                                if ui
                                    .add(egui::Button::new(
                                        RichText::new(color_name).color(color),
                                    ))
                                    .clicked()
                                {
                                    bus.execute_command(cmd_sel(
                                        &format!("set {}, {}", setting, color_name),
                                        &object,
                                    ));
                                    return true;
                                }
                            }
                            false
                        },
                    )
                } else {
                    false
                }
            }
        }
    }

    /// Render a submenu panel at the given anchor position.
    /// Tracks the area rect for click-outside detection.
    fn render_submenu_panel(
        ui: &mut Ui,
        ui_state: &mut ObjectListUiState,
        id: &str,
        anchor: Pos2,
        content: impl FnOnce(&mut Ui) -> bool,
    ) -> bool {
        let mut clicked = false;
        let area_resp = egui::Area::new(Id::new(id))
            .order(egui::Order::Tooltip)
            .fixed_pos(anchor)
            .pivot(Align2::RIGHT_CENTER)
            .movable(false)
            .constrain(false)
            .show(ui.ctx(), |ui| {
                popup_frame(ui).show(ui, |ui| {
                    ui.set_min_width(100.0);
                    clicked = content(ui);
                });
            });
        ui_state.menu_rects.push(area_resp.response.rect);
        clicked
    }
}
