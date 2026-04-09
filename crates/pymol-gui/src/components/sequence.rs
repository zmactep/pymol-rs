//! Sequence Viewer Component
//!
//! Bundles [`SequenceModel`] + [`SequenceUiState`] + view into a single
//! self-contained component.

use pymol_color::ElementColors;
use pymol_framework::component::{Component, EguiComponent, SharedContext};
use pymol_framework::message::MessageBus;
use pymol_framework::model::{SequenceColorContext, SequenceModel, SequenceUiState};
use crate::ui::SequencePanel;

/// Self-contained sequence viewer component.
pub struct SequenceComponent {
    pub model: SequenceModel,
    pub ui_state: SequenceUiState,
    /// Element color table (created once, reused across rebuilds)
    element_colors: ElementColors,
}

impl SequenceComponent {
    pub fn new() -> Self {
        Self {
            model: SequenceModel::new(),
            ui_state: SequenceUiState::new(),
            element_colors: ElementColors::default(),
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
}

impl EguiComponent for SequenceComponent {
    fn show(&mut self, ui: &mut egui::Ui, ctx: &SharedContext, bus: &mut MessageBus) {
        if self.model.needs_rebuild(ctx.registry, self.ui_state.viewer_colors) {
            let color_ctx = SequenceColorContext {
                named_colors: ctx.named_colors,
                element_colors: &self.element_colors,
            };
            self.model.rebuild_cache(ctx.registry, &color_ctx);
        }
        let has_sele = ctx.selections.contains("sele");
        SequencePanel::show(ui, &self.model, &mut self.ui_state, has_sele, bus);
    }
}
