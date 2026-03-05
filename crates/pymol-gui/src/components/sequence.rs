//! Sequence Viewer Component
//!
//! Bundles [`SequenceModel`] + [`SequenceUiState`] + view into a single
//! self-contained component.

use pymol_framework::component::{Component, SharedContext};
use pymol_framework::message::MessageBus;
use crate::model::SequenceModel;
use crate::ui::sequence::SequenceUiState;
use crate::ui::SequencePanel;

/// Self-contained sequence viewer component.
pub struct SequenceComponent {
    pub model: SequenceModel,
    pub ui_state: SequenceUiState,
}

impl SequenceComponent {
    pub fn new() -> Self {
        Self {
            model: SequenceModel::new(),
            ui_state: SequenceUiState::new(),
        }
    }
}

impl Default for SequenceComponent {
    fn default() -> Self {
        Self::new()
    }
}

impl Component for SequenceComponent {
    fn id(&self) -> &'static str {
        "sequence"
    }

    fn title(&self) -> &str {
        "Sequence"
    }

    fn show(&mut self, ui: &mut egui::Ui, ctx: &SharedContext, bus: &mut MessageBus) {
        let has_sele = ctx.selections.contains("sele");
        SequencePanel::show(
            ui,
            &mut self.model,
            &mut self.ui_state,
            ctx.registry,
            ctx.named_colors,
            has_sele,
            bus,
        );
    }
}
