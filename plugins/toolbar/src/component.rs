//! Toolbar Component
//!
//! Wraps [`ToolbarPanel`] into a self-contained component with lazy-loaded
//! icon textures. Groups and buttons are configured at construction time.

use std::collections::HashMap;

use pymol_plugin::prelude::{Component, EguiComponent, SharedContext, MessageBus};
use crate::model::{ToolbarGroup, default_toolbar_groups};
use crate::panel::ToolbarPanel;

/// Self-contained toolbar component with icon texture cache.
pub struct ToolbarComponent {
    groups: Vec<ToolbarGroup>,
    textures: HashMap<String, egui::TextureHandle>,
}

impl ToolbarComponent {
    /// Create with default toolbar groups.
    pub fn new() -> Self {
        Self {
            groups: default_toolbar_groups(),
            textures: HashMap::new(),
        }
    }
}

impl Default for ToolbarComponent {
    fn default() -> Self {
        Self::new()
    }
}

impl Component for ToolbarComponent {
    fn id(&self) -> &'static str {
        "toolbar"
    }

    fn title(&self) -> &str {
        "Toolbar"
    }
}

impl EguiComponent for ToolbarComponent {
    fn show(&mut self, ui: &mut egui::Ui, _ctx: &SharedContext, bus: &mut MessageBus) {
        ToolbarPanel::show(ui, &self.groups, &mut self.textures, bus);
    }
}
