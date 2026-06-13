//! Component System
//!
//! Defines the [`Component`] trait for self-contained UI panels and
//! [`SharedContext`] for read-only access to application state.

use patinae_cmd::{CommandRegistry, DynamicSettingRegistry};
use patinae_color::NamedPalette;
use patinae_scene::{Camera, Movie, ObjectRegistry, SelectionManager, ViewportImage};
use patinae_settings::Settings;

use crate::atom_stream::{AtomStreamPlan, AtomStreamRequest, AtomStreamRows};
use crate::message::AppMessage;

/// Read-only snapshot of application state, shared with all components.
///
/// Components read this context during rendering but never mutate it directly.
/// All mutations go through [`crate::message::MessageBus`].
pub struct SharedContext<'a> {
    // Scene state
    pub registry: &'a ObjectRegistry,
    pub camera: &'a Camera,
    pub selections: &'a SelectionManager,
    pub named_palette: &'a NamedPalette,
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

impl<'a> SharedContext<'a> {
    /// Opens a borrowed atom iterator.
    ///
    /// # Errors
    /// Returns an error when the request references missing objects or invalid
    /// selection expressions.
    pub fn iter_atoms(&self, request: &AtomStreamRequest) -> Result<AtomStreamRows<'a>, String> {
        let plan = AtomStreamPlan::open(self, request)?;
        Ok(plan.iter(self.registry))
    }

    /// Opens a reusable atom stream plan.
    ///
    /// # Errors
    /// Returns an error when the request cannot be resolved against this scene.
    pub fn atom_stream_plan(&self, request: &AtomStreamRequest) -> Result<AtomStreamPlan, String> {
        AtomStreamPlan::open(self, request)
    }
}

/// A self-contained UI component (panel).
///
/// Each component owns its domain model and UI state.
/// It reads shared application state via [`SharedContext`] and writes
/// actions via [`crate::message::MessageBus`].
///
/// # Design
///
/// - `on_message()` lets components react to dispatched messages
/// - The 3D viewport is intentionally **not** a Component
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
