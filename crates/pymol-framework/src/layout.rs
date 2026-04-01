//! Layout Data Types and State
//!
//! Configurable panel placement around the central 3D viewport.
//! This module contains the layout **data model** and **state mutation** methods.
//! Rendering (egui-specific) stays in `pymol-gui`.

use serde::{Deserialize, Serialize};

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
    pub fn side(&self) -> Option<SlotSide> {
        match self {
            Slot::Top { .. } => Some(SlotSide::Top),
            Slot::Bottom { .. } => Some(SlotSide::Bottom),
            Slot::Left { .. } => Some(SlotSide::Left),
            Slot::Right { .. } => Some(SlotSide::Right),
            Slot::Floating { .. } => None,
        }
    }
}

// ---------------------------------------------------------------------------
// Panel configuration
// ---------------------------------------------------------------------------

/// One panel placement in the layout.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PanelSlot {
    /// Component ID — must match [`Component::id()`](crate::component::Component::id).
    pub component_id: String,
    /// Where the panel is placed.
    pub slot: Slot,
    /// Whether the panel is currently expanded (false = collapsed tab).
    pub expanded: bool,
    /// Whether the user can drag-resize the panel edge.
    pub resizable: bool,
    /// Whether the panel can be detached into a floating window.
    pub floatable: bool,
    /// Whether the panel is visible at all (when false, not rendered even as a collapsed tab).
    #[serde(default = "default_visible")]
    pub visible: bool,
    /// Original docked slot, saved when the panel is floated so it can dock back.
    #[serde(skip)]
    pub docked_slot: Option<Slot>,
}

/// Convenience type for plugin component registration.
///
/// Like [`PanelSlot`] but without `component_id` or runtime state (`docked_slot`).
/// Use [`to_panel_slot`](PanelConfig::to_panel_slot) to convert.
pub struct PanelConfig {
    pub slot: Slot,
    pub expanded: bool,
    pub resizable: bool,
    pub floatable: bool,
    pub visible: bool,
}

impl PanelConfig {
    /// Panel at the top of the window.
    pub fn top(height: f32) -> Self {
        Self {
            slot: Slot::Top { height },
            expanded: true,
            resizable: true,
            floatable: true,
            visible: true,
        }
    }

    /// Panel at the bottom of the window.
    pub fn bottom(height: f32) -> Self {
        Self {
            slot: Slot::Bottom { height },
            expanded: true,
            resizable: true,
            floatable: true,
            visible: true,
        }
    }

    /// Panel at the left side of the window.
    pub fn left(width: f32) -> Self {
        Self {
            slot: Slot::Left { width },
            expanded: true,
            resizable: true,
            floatable: true,
            visible: true,
        }
    }

    /// Panel at the right side of the window.
    pub fn right(width: f32) -> Self {
        Self {
            slot: Slot::Right { width },
            expanded: true,
            resizable: true,
            floatable: true,
            visible: true,
        }
    }

    /// Floating window at default position.
    pub fn floating(width: f32, height: f32) -> Self {
        Self {
            slot: Slot::Floating {
                pos: [100.0, 100.0],
                size: [width, height],
            },
            expanded: true,
            resizable: true,
            floatable: false,
            visible: true,
        }
    }

    /// Convert into a [`PanelSlot`] with the given component ID.
    pub fn to_panel_slot(self, component_id: String) -> PanelSlot {
        PanelSlot {
            component_id,
            slot: self.slot,
            expanded: self.expanded,
            resizable: self.resizable,
            floatable: self.floatable,
            visible: self.visible,
            docked_slot: None,
        }
    }
}

// ---------------------------------------------------------------------------
// Layout
// ---------------------------------------------------------------------------

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

fn default_visible() -> bool {
    true
}

impl Layout {
    /// Empty layout with no panels.
    pub fn empty() -> Self {
        Self {
            panels: Vec::new(),
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
            .is_some_and(|p| p.expanded);

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
            .is_some_and(|p| p.expanded)
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

    /// Show a panel (make visible and expanded).
    pub fn show_panel(&mut self, id: &str) {
        if let Some(panel) = self.panels.iter_mut().find(|p| p.component_id == id) {
            panel.visible = true;
            panel.expanded = true;
            // Collapse siblings on the same side
            let side = panel.slot.side();
            let id_owned = id.to_string();
            if let Some(side) = side {
                for panel in &mut self.panels {
                    if panel.slot.side() == Some(side) && panel.component_id != id_owned {
                        panel.expanded = false;
                    }
                }
            }
        }
    }

    /// Hide a panel (make invisible — not even rendered as a collapsed tab).
    pub fn hide_panel(&mut self, id: &str) {
        if let Some(panel) = self.panels.iter_mut().find(|p| p.component_id == id) {
            panel.visible = false;
            panel.expanded = false;
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
}
