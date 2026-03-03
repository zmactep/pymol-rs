//! Layout Engine
//!
//! Configurable panel placement around the central 3D viewport.
//! Panels are rendered in a fixed order (top -> bottom -> left -> right -> floating)
//! so that egui's panel system correctly allocates space.

use serde::{Deserialize, Serialize};

use crate::component::SharedContext;
use crate::component_store::ComponentStore;
use crate::message::MessageBus;

// ---------------------------------------------------------------------------
// Slot types
// ---------------------------------------------------------------------------

/// Which side of the window a docked panel belongs to.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum SlotSide {
    Top,
    Bottom,
    Left,
    Right,
}

/// Where a panel is positioned relative to the 3D viewport.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum Slot {
    /// Fixed at top of window. Panels stack top-to-bottom in declaration order.
    Top { height: f32 },
    /// Fixed at bottom of window.
    Bottom { height: f32 },
    /// Fixed at left side.
    Left { width: f32 },
    /// Fixed at right side.
    Right { width: f32 },
    /// Floating window at an arbitrary position.
    Floating { pos: [f32; 2], size: [f32; 2] },
}

impl Slot {
    /// Which side this slot belongs to, or `None` for floating panels.
    fn side(&self) -> Option<SlotSide> {
        match self {
            Slot::Top { .. } => Some(SlotSide::Top),
            Slot::Bottom { .. } => Some(SlotSide::Bottom),
            Slot::Left { .. } => Some(SlotSide::Left),
            Slot::Right { .. } => Some(SlotSide::Right),
            Slot::Floating { .. } => None,
        }
    }

    /// Convert to an orientation for the panel shell renderer.
    fn orientation(&self) -> Option<Orientation> {
        match *self {
            Slot::Top { height } => Some(Orientation::Horizontal { is_top: true, size: height }),
            Slot::Bottom { height } => Some(Orientation::Horizontal { is_top: false, size: height }),
            Slot::Left { width } => Some(Orientation::Vertical { is_left: true, size: width }),
            Slot::Right { width } => Some(Orientation::Vertical { is_left: false, size: width }),
            Slot::Floating { .. } => None,
        }
    }
}

/// Axis + direction, abstracting over the four docked slot types.
#[derive(Debug, Clone, Copy)]
enum Orientation {
    /// Top or Bottom panel — uses `TopBottomPanel`, sized by height.
    Horizontal { is_top: bool, size: f32 },
    /// Left or Right panel — uses `SidePanel`, sized by width.
    Vertical { is_left: bool, size: f32 },
}

// ---------------------------------------------------------------------------
// Panel configuration
// ---------------------------------------------------------------------------

/// One panel placement in the layout.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PanelSlot {
    /// Component ID — must match [`Component::id()`].
    pub component_id: String,
    /// Where the panel is placed.
    pub slot: Slot,
    /// Whether the panel is currently expanded (false = collapsed tab).
    pub expanded: bool,
    /// Whether the user can drag-resize the panel edge.
    pub resizable: bool,
    /// Whether the panel can be detached into a floating window.
    pub floatable: bool,
    /// Original docked slot, saved when the panel is floated so it can dock back.
    #[serde(skip)]
    pub docked_slot: Option<Slot>,
}

/// The full window layout configuration.
///
/// Panels are rendered in declaration order within each slot type.
/// The order of slot types is: Top, Bottom, Left, Right, then Floating.
/// After all panels are rendered, `CentralPanel` claims the remaining space
/// as the 3D viewport.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Layout {
    pub panels: Vec<PanelSlot>,
}

impl Default for Layout {
    fn default() -> Self {
        Self::pymol_classic()
    }
}

impl Layout {
    /// Empty layout with no panels.
    pub fn empty() -> Self {
        Self {
            panels: Vec::new(),
        }
    }

    /// Default PyMOL-like layout: REPL on top, objects on right.
    pub fn pymol_classic() -> Self {
        Self {
            panels: vec![
                PanelSlot {
                    component_id: "repl".into(),
                    slot: Slot::Top { height: 150.0 },
                    expanded: true,
                    resizable: true,
                    floatable: true,
                    docked_slot: None,
                },
                PanelSlot {
                    component_id: "movie".into(),
                    slot: Slot::Bottom { height: 40.0 },
                    expanded: false,
                    resizable: false,
                    floatable: true,
                    docked_slot: None,
                },
                PanelSlot {
                    component_id: "sequence".into(),
                    slot: Slot::Bottom { height: 80.0 },
                    expanded: false,
                    resizable: true,
                    floatable: true,
                    docked_slot: None,
                },
                PanelSlot {
                    component_id: "object_list".into(),
                    slot: Slot::Right { width: 200.0 },
                    expanded: true,
                    resizable: true,
                    floatable: true,
                    docked_slot: None,
                },
            ],
        }
    }

    /// Toggle expanded/collapsed state of a panel by component ID.
    ///
    /// Group-aware: in a multi-panel group (same slot side), expanding
    /// one panel collapses its siblings (tab switching). Collapsing the active
    /// tab collapses the entire group.
    pub fn toggle_expanded(&mut self, id: &str) {
        let side = self
            .panels
            .iter()
            .find(|p| p.component_id == id)
            .and_then(|p| p.slot.side());

        let is_currently_expanded = self
            .panels
            .iter()
            .find(|p| p.component_id == id)
            .map_or(false, |p| p.expanded);

        if let Some(side) = side {
            if is_currently_expanded {
                // Collapse the active tab -> collapses the group
                if let Some(panel) = self.panels.iter_mut().find(|p| p.component_id == id) {
                    panel.expanded = false;
                }
            } else {
                // Expand this tab -> collapse siblings in the same group
                for panel in &mut self.panels {
                    if panel.slot.side() == Some(side) {
                        panel.expanded = panel.component_id == id;
                    }
                }
            }
        } else {
            // Floating panel -- simple toggle
            if let Some(panel) = self.panels.iter_mut().find(|p| p.component_id == id) {
                panel.expanded = !panel.expanded;
            }
        }
    }

    /// Set expanded state of a panel by component ID.
    pub fn set_expanded(&mut self, id: &str, expanded: bool) {
        if let Some(panel) = self.panels.iter_mut().find(|p| p.component_id == id) {
            panel.expanded = expanded;
        }
    }

    /// Check if a panel is expanded.
    pub fn is_expanded(&self, id: &str) -> bool {
        self.panels
            .iter()
            .find(|p| p.component_id == id)
            .map_or(false, |p| p.expanded)
    }

    /// Detach a docked panel into a floating window.
    pub fn float_panel(&mut self, id: &str) {
        if let Some(panel) = self.panels.iter_mut().find(|p| p.component_id == id) {
            if matches!(panel.slot, Slot::Floating { .. }) {
                return;
            }
            let original = std::mem::replace(
                &mut panel.slot,
                Slot::Floating {
                    pos: [100.0, 100.0],
                    size: [400.0, 300.0],
                },
            );
            panel.docked_slot = Some(original);
            panel.expanded = true;
        }
    }

    /// Dock a floating panel back to its original slot.
    pub fn dock_panel(&mut self, id: &str) {
        if let Some(panel) = self.panels.iter_mut().find(|p| p.component_id == id) {
            if let Some(original) = panel.docked_slot.take() {
                panel.slot = original;
            }
        }
    }

    /// Make a panel the active tab in its group (expand it, collapse siblings).
    pub fn activate_tab(&mut self, id: &str) {
        let side = self
            .panels
            .iter()
            .find(|p| p.component_id == id)
            .and_then(|p| p.slot.side());
        if let Some(side) = side {
            for panel in &mut self.panels {
                if panel.slot.side() == Some(side) {
                    panel.expanded = panel.component_id == id;
                }
            }
        }
    }

    /// Render all panels around the central viewport and return the viewport rect.
    ///
    /// Panels are rendered in a specific order to satisfy egui's layout rules:
    /// Top, Bottom, Left, Right, then Floating. CentralPanel claims the rest.
    ///
    /// Panels sharing the same slot type are rendered as tabs within a single
    /// egui panel. A solo panel renders with expand/collapse as before.
    pub fn show(
        &self,
        egui_ctx: &egui::Context,
        store: &mut ComponentStore,
        shared: &SharedContext,
        bus: &mut MessageBus,
    ) -> egui::Rect {
        // Collect panel descriptors so closures don't borrow self.
        let descriptors: Vec<_> = self
            .panels
            .iter()
            .map(|p| PanelDescriptor {
                id: p.component_id.clone(),
                slot: p.slot.clone(),
                resizable: p.resizable,
                expanded: p.expanded,
                floatable: p.floatable,
                is_floating: p.docked_slot.is_some(),
            })
            .collect();

        // Group docked panels by side, preserving declaration order.
        let mut groups: [Vec<&PanelDescriptor>; 4] = Default::default();
        let mut floating: Vec<&PanelDescriptor> = Vec::new();

        for desc in &descriptors {
            match desc.slot.side() {
                Some(SlotSide::Top) => groups[0].push(desc),
                Some(SlotSide::Bottom) => groups[1].push(desc),
                Some(SlotSide::Left) => groups[2].push(desc),
                Some(SlotSide::Right) => groups[3].push(desc),
                None => floating.push(desc),
            }
        }

        // Render docked groups in egui's required order: Top, Bottom, Left, Right
        for group in &groups {
            render_group(egui_ctx, store, group, shared, bus);
        }

        // Floating panels
        for desc in &floating {
            render_floating_panel(egui_ctx, store, desc, shared, bus);
        }

        // CentralPanel claims remaining space = 3D viewport
        let response = egui::CentralPanel::default()
            .frame(egui::Frame::NONE)
            .show(egui_ctx, |ui| {
                ui.allocate_response(ui.available_size(), egui::Sense::hover())
            });

        response.inner.rect
    }
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// Snapshot of a panel's properties for rendering (avoids borrowing Layout).
struct PanelDescriptor {
    id: String,
    slot: Slot,
    resizable: bool,
    expanded: bool,
    floatable: bool,
    /// True when the panel was docked and is now floating (can dock back).
    is_floating: bool,
}

/// Height of a collapsed horizontal tab (Top/Bottom panels).
const COLLAPSED_TAB_HEIGHT: f32 = 22.0;

/// Width of a collapsed vertical tab (Left/Right panels).
const COLLAPSED_TAB_WIDTH: f32 = 24.0;

/// Icon for the "detach to float" button.
const FLOAT_ICON: &str = "\u{2197}"; // arrow upper-right

/// Icon for the "dock back" button.
const DOCK_ICON: &str = "\u{2199}"; // arrow lower-left

// ---------------------------------------------------------------------------
// Panel shell — the core abstraction that eliminates the 4-way duplication
// ---------------------------------------------------------------------------

/// Show an egui panel shell (TopBottomPanel or SidePanel) based on orientation.
///
/// When `expanded`, the panel uses its configured size and is optionally resizable.
/// When collapsed, it shrinks to a minimal tab strip.
fn show_panel_shell(
    egui_ctx: &egui::Context,
    orientation: Orientation,
    panel_id: egui::Id,
    resizable: bool,
    expanded: bool,
    render_inner: impl FnOnce(&mut egui::Ui),
) {
    match orientation {
        Orientation::Horizontal { is_top, size } => {
            let panel = if is_top {
                egui::TopBottomPanel::top(panel_id)
            } else {
                egui::TopBottomPanel::bottom(panel_id)
            };
            if expanded {
                panel
                    .resizable(resizable)
                    .default_height(size)
                    .min_height(size)
                    .show(egui_ctx, render_inner);
            } else {
                panel
                    .resizable(false)
                    .exact_height(COLLAPSED_TAB_HEIGHT)
                    .show(egui_ctx, render_inner);
            }
        }
        Orientation::Vertical { is_left, size } => {
            let panel = if is_left {
                egui::SidePanel::left(panel_id)
            } else {
                egui::SidePanel::right(panel_id)
            };
            if expanded {
                panel
                    .resizable(resizable)
                    .default_width(size)
                    .show(egui_ctx, render_inner);
            } else {
                panel
                    .resizable(false)
                    .exact_width(COLLAPSED_TAB_WIDTH)
                    .show(egui_ctx, render_inner);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Group rendering
// ---------------------------------------------------------------------------

/// Render a group of panels sharing the same slot type.
///
/// - Empty group: no-op.
/// - Single panel: rendered with expand/collapse as a standalone panel.
/// - Multiple panels: rendered as tabs within a single egui panel.
fn render_group(
    egui_ctx: &egui::Context,
    store: &mut ComponentStore,
    group: &[&PanelDescriptor],
    shared: &SharedContext,
    bus: &mut MessageBus,
) {
    match group.len() {
        0 => {}
        1 => render_single_panel(egui_ctx, store, group[0], shared, bus),
        _ => render_tabbed_group(egui_ctx, store, group, shared, bus),
    }
}

/// Render a single docked panel -- expanded with header or collapsed tab.
fn render_single_panel(
    egui_ctx: &egui::Context,
    store: &mut ComponentStore,
    desc: &PanelDescriptor,
    shared: &SharedContext,
    bus: &mut MessageBus,
) {
    let orientation = match desc.slot.orientation() {
        Some(o) => o,
        None => return, // Floating handled separately
    };
    let panel_id = egui::Id::new(desc.id.as_str());
    let expanded = desc.expanded;
    let floatable = desc.floatable;
    let resizable = desc.resizable;

    store.show_in_panel(
        &desc.id,
        |component, shared, bus| {
            let title = component.title().to_owned();
            let comp_id = component.id();

            show_panel_shell(egui_ctx, orientation, panel_id, resizable, expanded, |ui| {
                if expanded {
                    show_panel_header(ui, &title, comp_id, floatable, bus);
                    component.show(ui, shared, bus);
                } else {
                    show_collapsed_tab(ui, orientation, &title, comp_id, bus);
                }
            });
        },
        shared,
        bus,
    );
}

/// Render multiple panels sharing a slot type as tabs within a single egui panel.
fn render_tabbed_group(
    egui_ctx: &egui::Context,
    store: &mut ComponentStore,
    group: &[&PanelDescriptor],
    shared: &SharedContext,
    bus: &mut MessageBus,
) {
    let orientation = match group[0].slot.orientation() {
        Some(o) => o,
        None => return,
    };
    let side = group[0].slot.side().expect("tabbed group must be docked");
    let group_id = egui::Id::new(format!("tab_group_{side:?}"));
    let any_expanded = group.iter().any(|d| d.expanded);
    let active_id = group
        .iter()
        .find(|d| d.expanded)
        .map(|d| d.id.clone())
        .unwrap_or_else(|| group[0].id.clone());

    // Collect IDs and floatable flags for the tab bar.
    let tabs: Vec<_> = group
        .iter()
        .map(|d| (d.id.clone(), d.floatable))
        .collect();

    // Use active descriptor's sizing and resizable setting.
    let active_desc = group.iter().find(|d| d.id == active_id).unwrap_or(&group[0]);
    let resizable = active_desc.resizable;

    // Override orientation size from the active descriptor.
    let orientation = match active_desc.slot.orientation() {
        Some(o) => o,
        None => orientation,
    };

    let render_inner = &mut |ui: &mut egui::Ui| {
        if any_expanded {
            let active_floatable = tabs
                .iter()
                .find(|(id, _)| *id == active_id)
                .map_or(false, |(_, f)| *f);

            show_tab_bar(ui, &tabs, &active_id, active_floatable, store, bus);
            ui.separator();

            store.show_in_panel(
                &active_id,
                |component, shared, bus| {
                    component.show(ui, shared, bus);
                },
                shared,
                bus,
            );
        } else {
            show_collapsed_tab_bar(ui, &tabs, store, bus);
        }
    };

    show_panel_shell(egui_ctx, orientation, group_id, resizable, any_expanded, render_inner);
}

/// Render a floating panel as an egui::Window.
fn render_floating_panel(
    egui_ctx: &egui::Context,
    store: &mut ComponentStore,
    desc: &PanelDescriptor,
    shared: &SharedContext,
    bus: &mut MessageBus,
) {
    if !desc.expanded {
        return;
    }
    let (pos, size) = match desc.slot {
        Slot::Floating { pos, size } => (pos, size),
        _ => return,
    };
    let panel_id = egui::Id::new(desc.id.as_str());
    let can_dock = desc.is_floating;

    store.show_in_panel(
        &desc.id,
        |component, shared, bus| {
            egui::Window::new(component.title())
                .id(panel_id)
                .default_pos(egui::pos2(pos[0], pos[1]))
                .default_size(egui::vec2(size[0], size[1]))
                .resizable(desc.resizable)
                .collapsible(true)
                .show(egui_ctx, |ui| {
                    if can_dock {
                        ui.horizontal(|ui| {
                            ui.with_layout(
                                egui::Layout::right_to_left(egui::Align::Center),
                                |ui| {
                                    if ui
                                        .add(egui::Button::new(DOCK_ICON).frame(false))
                                        .on_hover_text("Dock panel back")
                                        .on_hover_cursor(egui::CursorIcon::PointingHand)
                                        .clicked()
                                    {
                                        bus.send(crate::message::AppMessage::DockPanel(
                                            component.id().to_string(),
                                        ));
                                    }
                                },
                            );
                        });
                        ui.separator();
                    }
                    component.show(ui, shared, bus);
                });
        },
        shared,
        bus,
    );
}

// ---------------------------------------------------------------------------
// Tab bar rendering
// ---------------------------------------------------------------------------

/// Render the tab bar for an expanded tabbed group.
fn show_tab_bar(
    ui: &mut egui::Ui,
    tabs: &[(String, bool)],
    active_id: &str,
    active_floatable: bool,
    store: &mut ComponentStore,
    bus: &mut MessageBus,
) {
    ui.horizontal(|ui| {
        for (id, _floatable) in tabs {
            let title = store.get_title(id).unwrap_or_else(|| id.clone());
            let is_active = id == active_id;

            if ui.selectable_label(is_active, &title).clicked() {
                if is_active {
                    bus.send(crate::message::AppMessage::TogglePanel(id.clone()));
                } else {
                    bus.send(crate::message::AppMessage::ActivateTab(id.clone()));
                }
            }
        }

        if active_floatable {
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                if ui
                    .add(egui::Button::new(FLOAT_ICON).frame(false))
                    .on_hover_text("Detach to floating window")
                    .on_hover_cursor(egui::CursorIcon::PointingHand)
                    .clicked()
                {
                    bus.send(crate::message::AppMessage::FloatPanel(
                        active_id.to_string(),
                    ));
                }
            });
        }
    });
}

/// Render a collapsed tab bar showing all tab titles in a horizontal strip.
fn show_collapsed_tab_bar(
    ui: &mut egui::Ui,
    tabs: &[(String, bool)],
    store: &mut ComponentStore,
    bus: &mut MessageBus,
) {
    ui.horizontal_centered(|ui| {
        for (id, _) in tabs {
            let title = store.get_title(id).unwrap_or_else(|| id.clone());
            if ui
                .add(egui::Button::new(&title).frame(false))
                .on_hover_cursor(egui::CursorIcon::PointingHand)
                .clicked()
            {
                bus.send(crate::message::AppMessage::ActivateTab(id.clone()));
            }
        }
    });
}

// ---------------------------------------------------------------------------
// Panel header & collapsed tab rendering
// ---------------------------------------------------------------------------

/// Render a clickable title bar for an expanded docked panel.
///
/// Clicking the title collapses the panel. If `floatable`, a float button
/// detaches the panel into a floating window.
fn show_panel_header(
    ui: &mut egui::Ui,
    title: &str,
    component_id: &str,
    floatable: bool,
    bus: &mut MessageBus,
) {
    ui.horizontal(|ui| {
        if ui
            .add(egui::Button::new(egui::RichText::new(title).strong()).frame(false))
            .on_hover_cursor(egui::CursorIcon::PointingHand)
            .clicked()
        {
            bus.send(crate::message::AppMessage::TogglePanel(
                component_id.to_string(),
            ));
        }

        if floatable {
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                if ui
                    .add(egui::Button::new(FLOAT_ICON).frame(false))
                    .on_hover_text("Detach to floating window")
                    .on_hover_cursor(egui::CursorIcon::PointingHand)
                    .clicked()
                {
                    bus.send(crate::message::AppMessage::FloatPanel(
                        component_id.to_string(),
                    ));
                }
            });
        }
    });
    ui.separator();
}

/// Render a collapsed tab for a docked panel.
///
/// Horizontal panels (Top/Bottom) show a clickable text label.
/// Vertical panels (Left/Right) show rotated text.
fn show_collapsed_tab(
    ui: &mut egui::Ui,
    orientation: Orientation,
    title: &str,
    component_id: &str,
    bus: &mut MessageBus,
) {
    match orientation {
        Orientation::Horizontal { .. } => {
            ui.horizontal_centered(|ui| {
                if ui
                    .add(egui::Button::new(title).frame(false))
                    .on_hover_cursor(egui::CursorIcon::PointingHand)
                    .clicked()
                {
                    bus.send(crate::message::AppMessage::TogglePanel(
                        component_id.to_string(),
                    ));
                }
            });
        }
        Orientation::Vertical { .. } => {
            let rect = ui.available_rect_before_wrap();
            let response = ui.allocate_rect(rect, egui::Sense::click());

            if response.clicked() {
                bus.send(crate::message::AppMessage::TogglePanel(
                    component_id.to_string(),
                ));
            }

            if response.hovered() {
                ui.painter().rect_filled(
                    rect,
                    0.0,
                    ui.visuals().widgets.hovered.bg_fill,
                );
                ui.ctx().set_cursor_icon(egui::CursorIcon::PointingHand);
            }

            // Draw rotated title text (bottom-to-top)
            let painter = ui.painter();
            let font_id = egui::TextStyle::Button.resolve(ui.style());
            let color = ui.visuals().text_color();
            let galley = painter.layout_no_wrap(title.to_owned(), font_id, color);
            let center = rect.center();

            let mut text_shape = egui::epaint::TextShape::new(center, galley, color);
            text_shape.angle = -std::f32::consts::FRAC_PI_2;
            // Pivot around center so the text is visually centered in the strip.
            let half = text_shape.galley.rect.size() * 0.5;
            text_shape.pos -= egui::vec2(half.y, -half.x);

            painter.add(text_shape);
        }
    }
}
