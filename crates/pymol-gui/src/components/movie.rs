//! Movie Controls Component
//!
//! Wraps [`MoviePanel`] view into a self-contained component.
//! This component has no model — it reads movie state from [`SharedContext`].

use pymol_framework::component::{Component, SharedContext};
use pymol_framework::message::MessageBus;
use crate::ui::MoviePanel;

/// Self-contained movie controls component.
pub struct MovieComponent;

impl MovieComponent {
    pub fn new() -> Self {
        Self
    }
}

impl Default for MovieComponent {
    fn default() -> Self {
        Self::new()
    }
}

impl Component for MovieComponent {
    fn id(&self) -> &'static str {
        "movie"
    }

    fn title(&self) -> &str {
        "Movie"
    }

    fn show(&mut self, ui: &mut egui::Ui, ctx: &SharedContext, bus: &mut MessageBus) {
        MoviePanel::show(ui, ctx.movie, bus);
    }
}
