//! Component System
//!
//! Defines the [`Component`] trait for self-contained UI panels and
//! [`SharedContext`] for read-only access to application state.

use pymol_cmd::CommandRegistry;
use pymol_color::NamedColors;
use pymol_scene::{Camera, Movie, ObjectRegistry, SelectionManager, ViewportImage};

use crate::message::{AppMessage, MessageBus};

/// Read-only snapshot of application state, shared with all components.
///
/// Components read this context during rendering but never mutate it directly.
/// All mutations go through [`MessageBus`].
pub struct SharedContext<'a> {
    // Scene state
    pub registry: &'a ObjectRegistry,
    pub camera: &'a Camera,
    pub selections: &'a SelectionManager,
    pub named_colors: &'a NamedColors,
    pub movie: &'a Movie,

    // Viewport image overlay (e.g. from `ray` command or plugins)
    pub viewport_image: Option<&'a ViewportImage>,

    // Command system (for autocomplete)
    pub command_names: &'a [String],
    pub command_registry: &'a CommandRegistry,
    pub setting_names: &'a [&'a str],
}

/// A self-contained UI component (panel).
///
/// Each component owns its domain model and egui-specific UI state.
/// It reads shared application state via [`SharedContext`] and writes
/// actions via [`MessageBus`].
///
/// # Design
///
/// - `show()` is the main render method, called by the layout engine
/// - `on_message()` lets components react to dispatched messages
/// - The 3D viewport is intentionally **not** a Component (it needs
///   wgpu rendering and CentralPanel placement)
pub trait Component: std::any::Any {
    /// Unique string identifier (e.g., `"command_line"`, `"sequence"`).
    fn id(&self) -> &'static str;

    /// Display name for panel headers and tabs.
    fn title(&self) -> &str;

    /// Render the component's UI into the given region.
    fn show(&mut self, ui: &mut egui::Ui, ctx: &SharedContext, bus: &mut MessageBus);

    /// React to a dispatched message. Called once per message after dispatch.
    ///
    /// Default implementation: no-op.
    fn on_message(&mut self, _msg: &AppMessage) {}

    // =========================================================================
    // Lifecycle hooks (default no-ops)
    // =========================================================================

    /// Called once after the component is registered.
    fn on_init(&mut self) {}

    /// Called when the component's panel gains focus (e.g. tab selected).
    fn on_focus(&mut self) {}

    /// Called when the component's panel loses focus (e.g. tab deselected).
    fn on_blur(&mut self) {}
}
