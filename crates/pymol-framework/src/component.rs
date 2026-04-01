//! Component System
//!
//! Defines the [`Component`] trait for self-contained UI panels and
//! [`SharedContext`] for read-only access to application state.
//!
//! The base [`Component`] trait is UI-framework agnostic. Rendering is handled
//! by backend-specific extension traits (e.g. [`EguiComponent`] behind the
//! `egui` feature flag).

use pymol_cmd::{CommandRegistry, DynamicSettingRegistry};
use pymol_color::NamedColors;
use pymol_scene::{Camera, Movie, ObjectRegistry, SelectionManager, ViewportImage};
use pymol_settings::Settings;

use crate::message::AppMessage;
#[cfg(feature = "egui")]
use crate::message::MessageBus;

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
    pub settings: &'a Settings,
    pub clear_color: [f32; 3],

    // GPU (optional — available when a wgpu surface is initialised)
    pub gpu_device: Option<&'a wgpu::Device>,
    pub gpu_queue: Option<&'a wgpu::Queue>,

    // Change detection
    /// Monotonically increasing counter bumped on every scene change
    /// (camera, objects, commands, etc.). Compare against a cached value
    /// to detect viewport changes.
    pub scene_generation: u64,

    // Viewport image overlay (e.g. from `ray` command or plugins)
    pub viewport_image: Option<&'a ViewportImage>,

    // Command system (for autocomplete)
    pub command_names: &'a [String],
    pub command_registry: &'a CommandRegistry,
    pub setting_names: &'a [&'a str],
    pub dynamic_settings: Option<&'a DynamicSettingRegistry>,
}

/// A self-contained UI component (panel).
///
/// Each component owns its domain model and UI state.
/// It reads shared application state via [`SharedContext`] and writes
/// actions via [`MessageBus`].
///
/// # Design
///
/// - Rendering is handled by backend-specific traits (e.g. [`EguiComponent`])
/// - `on_message()` lets components react to dispatched messages
/// - The 3D viewport is intentionally **not** a Component (it needs
///   wgpu rendering and CentralPanel placement)
pub trait Component: std::any::Any {
    /// Unique string identifier (e.g., `"command_line"`, `"sequence"`).
    fn id(&self) -> &'static str;

    /// Display name for panel headers and tabs.
    fn title(&self) -> &str;

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

/// Extension trait for components that render via egui.
///
/// Enabled with the `egui` feature flag. Desktop applications and plugins
/// implement this trait; web frontends use their own rendering layer.
#[cfg(feature = "egui")]
pub trait EguiComponent: Component {
    /// Render the component's UI into the given egui region.
    fn show(&mut self, ui: &mut egui::Ui, ctx: &SharedContext, bus: &mut MessageBus);
}
