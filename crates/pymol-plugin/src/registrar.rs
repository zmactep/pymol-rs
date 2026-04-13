//! Plugin Registrar
//!
//! The host creates a [`PluginRegistrar`] and passes it to the plugin's
//! registration function. The plugin populates it with commands, components,
//! and message handlers. After registration, the host drains the registrar
//! and integrates everything into the running application.

use std::sync::Arc;

use pymol_cmd::Command;
pub use pymol_cmd::DynamicCommandInvocation;
pub use pymol_cmd::{FormatHandler, PluginReaderFn, PluginWriterFn, ScriptHandler, ViewerLike};
use pymol_framework::component::{EguiComponent, SharedContext};
use pymol_framework::layout::PanelConfig;
use pymol_framework::message::{AppMessage, MessageBus};
pub use pymol_scene::{parse_key_string, KeyBinding, KeyCode};
pub use pymol_settings::{DynamicSettingDescriptor, DynamicSettingStore, SharedSettingStore};

/// A boxed closure that mutates the viewer on the main thread.
///
/// Queued by plugins via [`PollContext::queue_viewer_mutation`] during `poll()`,
/// then executed by the host in phase 2 with full mutable viewer access.
pub type ViewerMutation = Box<dyn FnOnce(&mut dyn ViewerLike) + Send>;

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

// =============================================================================
// Hotkey types
// =============================================================================

/// Callback for hotkey actions — runs during the plugin poll phase.
pub type HotkeyCallback = Box<dyn FnMut(&mut PollContext<'_>) + Send>;

/// Action to perform when a plugin hotkey is triggered.
pub enum PluginKeyAction {
    /// Execute a PyMOL command string.
    Command(String),
    /// Invoke a dynamic command (delivered via `dynamic_invocations` in next poll).
    DynamicCommand { name: String, args: Vec<String> },
    /// Publish a custom message to the message bus.
    Custom { topic: String, payload: Vec<u8> },
    /// Run arbitrary code during the next poll phase with `PollContext` access.
    Callback(HotkeyCallback),
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
    /// Hotkey bindings triggered since last poll (read-only).
    pub triggered_hotkeys: &'a [KeyBinding],
    // Internal queues — accessed via methods
    pub(crate) exec_queue: &'a mut Vec<CommandExecRequest>,
    pub(crate) reg_queue: &'a mut Vec<DynCmdRegistration>,
    pub(crate) unreg_queue: &'a mut Vec<String>,
    pub(crate) notification_queue: &'a mut Vec<String>,
    pub(crate) hotkey_reg_queue: &'a mut Vec<(String, PluginKeyAction)>,
    pub(crate) hotkey_unreg_queue: &'a mut Vec<String>,
    pub(crate) mutation_queue: &'a mut Vec<ViewerMutation>,
}

impl<'a> PollContext<'a> {
    /// Create a new poll context.
    ///
    /// This is used by the host to build the context before calling `poll()`.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        shared: &'a SharedContext<'a>,
        bus: &'a mut MessageBus,
        command_results: &'a [CommandResult],
        dynamic_invocations: &'a [DynamicCommandInvocation],
        triggered_hotkeys: &'a [KeyBinding],
        exec_queue: &'a mut Vec<CommandExecRequest>,
        reg_queue: &'a mut Vec<DynCmdRegistration>,
        unreg_queue: &'a mut Vec<String>,
        notification_queue: &'a mut Vec<String>,
        hotkey_reg_queue: &'a mut Vec<(String, PluginKeyAction)>,
        hotkey_unreg_queue: &'a mut Vec<String>,
        mutation_queue: &'a mut Vec<ViewerMutation>,
    ) -> Self {
        Self {
            shared,
            bus,
            command_results,
            dynamic_invocations,
            triggered_hotkeys,
            exec_queue,
            reg_queue,
            unreg_queue,
            notification_queue,
            hotkey_reg_queue,
            hotkey_unreg_queue,
            mutation_queue,
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

    /// Show a notification message in the overlay (spinner + text).
    ///
    /// Call this during `poll()` while a background operation is in progress.
    /// The notification is cleared automatically when `poll()` returns without
    /// calling this method.
    pub fn set_notification(&mut self, msg: impl Into<String>) {
        self.notification_queue.push(msg.into());
    }

    /// Register a hotkey binding by key string (e.g. `"ctrl+s"`).
    ///
    /// The string is parsed on the host side (which has the correct `KeyCode`
    /// enum from winit). Applied after `poll()` returns.
    pub fn register_hotkey(&mut self, key: impl Into<String>, action: PluginKeyAction) {
        self.hotkey_reg_queue.push((key.into(), action));
    }

    /// Unregister a hotkey binding by key string (e.g. `"ctrl+s"`).
    ///
    /// Applied after `poll()` returns.
    pub fn unregister_hotkey(&mut self, key: impl Into<String>) {
        self.hotkey_unreg_queue.push(key.into());
    }

    /// Queue a mutation to be applied to the viewer on the main thread.
    ///
    /// The closure is executed by the host after `poll()` returns, during
    /// phase 2, with full mutable access to the viewer. Use this for
    /// operations that require modifying viewer state (e.g., atom properties).
    pub fn queue_viewer_mutation(
        &mut self,
        f: impl FnOnce(&mut dyn ViewerLike) + Send + 'static,
    ) {
        self.mutation_queue.push(Box::new(f));
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
    pub(crate) components: Vec<(Box<dyn EguiComponent>, PanelConfig)>,
    pub(crate) message_handler: Option<Box<dyn MessageHandler>>,
    pub(crate) script_handlers: Vec<(String, ScriptHandler)>,
    pub(crate) format_handlers: Vec<FormatHandler>,
    pub(crate) hotkeys: Vec<(KeyBinding, PluginKeyAction)>,
    pub(crate) settings: Vec<(Vec<DynamicSettingDescriptor>, SharedSettingStore)>,
}

impl PluginRegistrar {
    /// Create an empty registrar.
    pub fn new() -> Self {
        Self {
            metadata: None,
            commands: Vec::new(),
            components: Vec::new(),
            message_handler: None,
            script_handlers: Vec::new(),
            format_handlers: Vec::new(),
            hotkeys: Vec::new(),
            settings: Vec::new(),
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
    pub fn register_component(&mut self, comp: impl EguiComponent + 'static, config: PanelConfig) {
        self.components.push((Box::new(comp), config));
    }

    /// Set a message handler for headless message processing.
    pub fn set_message_handler(&mut self, handler: impl MessageHandler + 'static) {
        self.message_handler = Some(Box::new(handler));
    }

    /// Register a script handler for a specific extension.
    ///
    /// Used by the builtin `run` command to dispatch non-.pml files
    /// to the appropriate plugin (e.g., `.py` → Python plugin).
    pub fn register_script_handler(
        &mut self,
        extension: &str,
        handler: impl Fn(&str) -> Result<(), String> + Send + Sync + 'static,
    ) {
        self.script_handlers
            .push((extension.to_string(), Arc::new(handler)));
    }

    /// Register a file format handler for `load` and `save` commands.
    ///
    /// The handler specifies supported extensions and optional reader/writer
    /// factories. When a user runs `load file.ext` or `save file.ext` and
    /// the extension matches, the plugin's reader or writer is used.
    pub fn register_format_handler(&mut self, handler: FormatHandler) {
        self.format_handlers.push(handler);
    }

    /// Register plugin settings with their shared store.
    ///
    /// The plugin creates a [`SharedSettingStore`], populates it with defaults,
    /// and passes the descriptors + store here. The host will dispatch
    /// `set`/`get`/`unset` commands to the store when a matching setting
    /// name is used.
    ///
    /// Use the [`define_plugin_settings!`](crate::define_plugin_settings) macro
    /// to generate descriptors and store initialization from a struct definition.
    pub fn register_settings(
        &mut self,
        descriptors: Vec<DynamicSettingDescriptor>,
        store: SharedSettingStore,
    ) {
        self.settings.push((descriptors, store));
    }

    /// Register a keyboard shortcut.
    ///
    /// Plugin hotkeys are checked after the main application bindings.
    /// For `Callback` actions, the closure runs during the next poll phase
    /// with access to [`PollContext`].
    pub fn register_hotkey(&mut self, key: impl Into<KeyBinding>, action: PluginKeyAction) {
        self.hotkeys.push((key.into(), action));
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
    pub fn drain_components(&mut self) -> Vec<(Box<dyn EguiComponent>, PanelConfig)> {
        std::mem::take(&mut self.components)
    }

    /// Take the message handler (returns `None` if not set or already taken).
    pub fn take_message_handler(&mut self) -> Option<Box<dyn MessageHandler>> {
        self.message_handler.take()
    }

    /// Drain all registered script handlers.
    pub fn drain_script_handlers(&mut self) -> Vec<(String, ScriptHandler)> {
        std::mem::take(&mut self.script_handlers)
    }

    /// Drain all registered format handlers.
    pub fn drain_format_handlers(&mut self) -> Vec<FormatHandler> {
        std::mem::take(&mut self.format_handlers)
    }

    /// Drain all registered hotkey bindings.
    pub fn drain_hotkeys(&mut self) -> Vec<(KeyBinding, PluginKeyAction)> {
        std::mem::take(&mut self.hotkeys)
    }

    /// Drain all registered plugin settings (descriptors + shared stores).
    pub fn drain_settings(&mut self) -> Vec<(Vec<DynamicSettingDescriptor>, SharedSettingStore)> {
        std::mem::take(&mut self.settings)
    }
}

impl Default for PluginRegistrar {
    fn default() -> Self {
        Self::new()
    }
}
