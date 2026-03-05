//! Plugin Registrar
//!
//! The host creates a [`PluginRegistrar`] and passes it to the plugin's
//! registration function. The plugin populates it with commands, components,
//! and message handlers. After registration, the host drains the registrar
//! and integrates everything into the running application.

use std::sync::Arc;

use pymol_cmd::Command;
pub use pymol_cmd::FileHandler;
pub use pymol_cmd::DynamicCommandInvocation;
use pymol_framework::component::{Component, SharedContext};
use pymol_framework::layout::PanelConfig;
use pymol_framework::message::{AppMessage, MessageBus};

// =============================================================================
// Polling types
// =============================================================================

/// Result of a command execution requested via [`PollContext::execute_command`].
///
/// Delivered to the plugin in the next `poll()` call.
pub struct CommandResult {
    /// Correlation ID (set by the plugin when requesting execution).
    pub id: u64,
    /// `Ok(())` on success, `Err(message)` on failure.
    pub result: Result<(), String>,
}

/// Queued command execution request (internal).
pub struct CommandExecRequest {
    pub id: u64,
    pub command: String,
    pub silent: bool,
}

/// Queued dynamic command registration (internal).
pub struct DynCmdRegistration {
    pub name: String,
    pub help: String,
}

/// Context provided to plugins during periodic polling.
///
/// Plugins that return `needs_poll() == true` receive this each frame.
/// It provides read-only state access, message bus, deferred command
/// execution, and dynamic command registration.
pub struct PollContext<'a> {
    /// Read-only application state.
    pub shared: &'a SharedContext<'a>,
    /// Message bus for sending messages (print, quit, execute, etc.).
    pub bus: &'a mut MessageBus,
    /// Results from command executions requested in previous poll cycles.
    pub command_results: &'a [CommandResult],
    /// Invocations of dynamic commands since last poll.
    pub dynamic_invocations: &'a [DynamicCommandInvocation],
    // Internal queues — accessed via methods
    pub(crate) exec_queue: &'a mut Vec<CommandExecRequest>,
    pub(crate) reg_queue: &'a mut Vec<DynCmdRegistration>,
    pub(crate) unreg_queue: &'a mut Vec<String>,
}

impl<'a> PollContext<'a> {
    /// Create a new poll context.
    ///
    /// This is used by the host to build the context before calling `poll()`.
    pub fn new(
        shared: &'a SharedContext<'a>,
        bus: &'a mut MessageBus,
        command_results: &'a [CommandResult],
        dynamic_invocations: &'a [DynamicCommandInvocation],
        exec_queue: &'a mut Vec<CommandExecRequest>,
        reg_queue: &'a mut Vec<DynCmdRegistration>,
        unreg_queue: &'a mut Vec<String>,
    ) -> Self {
        Self {
            shared,
            bus,
            command_results,
            dynamic_invocations,
            exec_queue,
            reg_queue,
            unreg_queue,
        }
    }

    /// Queue a command for execution.
    ///
    /// The command is executed by the host after `poll()` returns.
    /// The result (success or error message) is delivered in the next
    /// `poll()` call via [`PollContext::command_results`].
    pub fn execute_command(&mut self, id: u64, command: &str, silent: bool) {
        self.exec_queue.push(CommandExecRequest {
            id,
            command: command.to_string(),
            silent,
        });
    }

    /// Register a dynamic command.
    ///
    /// The command appears in autocomplete and help immediately (next frame).
    /// When a user invokes it, the invocation is delivered to the plugin
    /// via [`PollContext::dynamic_invocations`] in the next `poll()` call.
    pub fn register_dynamic_command(&mut self, name: String, help: String) {
        self.reg_queue.push(DynCmdRegistration { name, help });
    }

    /// Unregister a dynamic command.
    pub fn unregister_dynamic_command(&mut self, name: &str) {
        self.unreg_queue.push(name.to_string());
    }
}

// =============================================================================
// MessageHandler trait
// =============================================================================

/// Plugin metadata — name, version, description.
pub struct PluginMetadata {
    pub name: &'static str,
    pub version: &'static str,
    pub description: &'static str,
}

/// Trait for headless plugins that need to react to messages without a GUI.
///
/// Plugins that also need periodic polling should override [`needs_poll`] and
/// [`poll`] — the host will call `poll()` each frame with a [`PollContext`].
pub trait MessageHandler: Send {
    /// Called for every dispatched message after components have been notified.
    fn on_message(&mut self, msg: &AppMessage, bus: &mut MessageBus);

    /// Whether this handler needs periodic `poll()` calls (default: `false`).
    ///
    /// When `true`, the host will call [`poll`] each frame and set a 50 ms
    /// wake-up timer so the event loop doesn't sleep indefinitely.
    fn needs_poll(&self) -> bool {
        false
    }

    /// Called each frame when [`needs_poll`] returns `true`.
    ///
    /// Use this for non-blocking I/O (e.g., polling a socket), deferred
    /// command execution, and dynamic command management.
    fn poll(&mut self, _ctx: &mut PollContext<'_>) {}
}

// =============================================================================
// PluginRegistrar
// =============================================================================

/// Accumulator filled by the plugin's registration function.
///
/// The host creates this, passes it to the plugin via FFI, then drains
/// the collected commands, components, and handlers.
pub struct PluginRegistrar {
    pub(crate) metadata: Option<PluginMetadata>,
    pub(crate) commands: Vec<Box<dyn Command>>,
    pub(crate) components: Vec<(Box<dyn Component>, PanelConfig)>,
    pub(crate) message_handler: Option<Box<dyn MessageHandler>>,
    pub(crate) file_handlers: Vec<(String, FileHandler)>,
}

impl PluginRegistrar {
    /// Create an empty registrar.
    pub fn new() -> Self {
        Self {
            metadata: None,
            commands: Vec::new(),
            components: Vec::new(),
            message_handler: None,
            file_handlers: Vec::new(),
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

    /// Register a file handler for a specific extension.
    ///
    /// Used by the builtin `run` command to dispatch non-.pml files
    /// to the appropriate plugin (e.g., `.py` → Python plugin).
    pub fn register_file_handler(
        &mut self,
        extension: &str,
        handler: impl Fn(&str) -> Result<(), String> + Send + Sync + 'static,
    ) {
        self.file_handlers
            .push((extension.to_string(), Arc::new(handler)));
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

    /// Drain all registered file handlers.
    pub fn drain_file_handlers(&mut self) -> Vec<(String, FileHandler)> {
        std::mem::take(&mut self.file_handlers)
    }
}

impl Default for PluginRegistrar {
    fn default() -> Self {
        Self::new()
    }
}
