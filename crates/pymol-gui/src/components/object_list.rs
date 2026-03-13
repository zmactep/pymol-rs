//! Object List Component
//!
//! Bundles [`ObjectListUiState`] + view into a single self-contained component.
//! This component has no domain model — it reads objects from [`SharedContext`].

use pymol_framework::component::{Component, EguiComponent, SharedContext};
use pymol_framework::message::MessageBus;
use crate::ui::objects::ObjectListUiState;
use crate::ui::ObjectListPanel;

/// Self-contained object list component.
pub struct ObjectListComponent {
    pub ui_state: ObjectListUiState,
}

impl ObjectListComponent {
    pub fn new() -> Self {
        Self {
            ui_state: ObjectListUiState::default(),
        }
    }
}

impl Default for ObjectListComponent {
    fn default() -> Self {
        Self::new()
    }
}

impl Component for ObjectListComponent {
    fn id(&self) -> &'static str {
        "object_list"
    }

    fn title(&self) -> &str {
        "Objects"
    }
}

impl EguiComponent for ObjectListComponent {
    fn show(&mut self, ui: &mut egui::Ui, ctx: &SharedContext, bus: &mut MessageBus) {
        egui::ScrollArea::vertical().show(ui, |ui| {
            ObjectListPanel::show(ui, &mut self.ui_state, ctx.registry, ctx.selections, bus);
        });
    }
}
