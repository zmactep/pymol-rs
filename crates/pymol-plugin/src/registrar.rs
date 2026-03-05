//! Plugin Registrar
//!
//! The host creates a [`PluginRegistrar`] and passes it to the plugin's
//! registration function. The plugin populates it with commands, components,
//! and message handlers. After registration, the host drains the registrar
//! and integrates everything into the running application.

use pymol_cmd::Command;
use pymol_framework::component::Component;
use pymol_framework::layout::PanelConfig;
use pymol_framework::message::{AppMessage, MessageBus};

/// Plugin metadata — name, version, description.
pub struct PluginMetadata {
    pub name: &'static str,
    pub version: &'static str,
    pub description: &'static str,
}

/// Trait for headless plugins that need to react to messages without a GUI.
pub trait MessageHandler: Send {
    /// Called for every dispatched message after components have been notified.
    fn on_message(&mut self, msg: &AppMessage, bus: &mut MessageBus);
}

/// Accumulator filled by the plugin's registration function.
///
/// The host creates this, passes it to the plugin via FFI, then drains
/// the collected commands, components, and handlers.
pub struct PluginRegistrar {
    pub(crate) metadata: Option<PluginMetadata>,
    pub(crate) commands: Vec<Box<dyn Command>>,
    pub(crate) components: Vec<(Box<dyn Component>, PanelConfig)>,
    pub(crate) message_handler: Option<Box<dyn MessageHandler>>,
}

impl PluginRegistrar {
    /// Create an empty registrar.
    pub fn new() -> Self {
        Self {
            metadata: None,
            commands: Vec::new(),
            components: Vec::new(),
            message_handler: None,
        }
    }

    /// Set the plugin's metadata.
    pub fn set_metadata(&mut self, metadata: PluginMetadata) {
        self.metadata = Some(metadata);
    }

    /// Register a command implementation.
    pub fn register_command(&mut self, cmd: impl Command + 'static) {
        self.commands.push(Box::new(cmd));
    }

    /// Register a GUI component with its panel configuration.
    pub fn register_component(&mut self, comp: impl Component + 'static, config: PanelConfig) {
        self.components.push((Box::new(comp), config));
    }

    /// Set a message handler for headless message processing.
    pub fn set_message_handler(&mut self, handler: impl MessageHandler + 'static) {
        self.message_handler = Some(Box::new(handler));
    }

    // =================================================================
    // Host-side drain methods
    // =================================================================

    /// Take the plugin metadata (returns `None` if already taken or not set).
    pub fn take_metadata(&mut self) -> Option<PluginMetadata> {
        self.metadata.take()
    }

    /// Drain all registered commands.
    pub fn drain_commands(&mut self) -> Vec<Box<dyn Command>> {
        std::mem::take(&mut self.commands)
    }

    /// Drain all registered components with their panel configurations.
    pub fn drain_components(&mut self) -> Vec<(Box<dyn Component>, PanelConfig)> {
        std::mem::take(&mut self.components)
    }

    /// Take the message handler (returns `None` if not set or already taken).
    pub fn take_message_handler(&mut self) -> Option<Box<dyn MessageHandler>> {
        self.message_handler.take()
    }
}

impl Default for PluginRegistrar {
    fn default() -> Self {
        Self::new()
    }
}
