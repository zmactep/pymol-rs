//! Plugin Manager
//!
//! Loads plugin shared libraries at runtime and integrates their commands,
//! components, and message handlers into the application.

use std::any::Any;
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::path::Path;
use std::sync::{Arc, Mutex};

use libloading::Library;

use pymol_cmd::CommandExecutor;
use pymol_plugin::ffi::{ABI_VERSION, SDK_VERSION};
use pymol_cmd::DynamicCommandInvocation;
use pymol_plugin::registrar::{
    CommandExecRequest, CommandResult, DynCmdRegistration, KeyBinding,
    MessageHandler, PluginKeyAction, PluginMetadata, PluginRegistrar, PollContext,
    ViewerMutation,
};
use pymol_scene::KeyBindings;

use pymol_framework::component::SharedContext;
use pymol_framework::component_store::ComponentStore;
use pymol_framework::message::{AppMessage, MessageBus};

use crate::layout::Layout;

/// A loaded plugin, keeping the shared library alive.
struct LoadedPlugin {
    /// Must stay alive for the duration of the plugin.
    _library: Library,
    _metadata: PluginMetadata,
    message_handler: Option<Box<dyn MessageHandler>>,
    /// Set to `true` after a panic; the plugin will be skipped for all future calls.
    faulted: bool,
    /// Keyboard shortcut bindings registered by this plugin.
    hotkeys: KeyBindings<PluginKeyAction>,
}

/// Manages dynamically loaded plugins.
pub struct PluginManager {
    plugins: Vec<LoadedPlugin>,
    // Deferred command execution
    pending_executions: Vec<CommandExecRequest>,
    command_results: Vec<CommandResult>,
    // Dynamic commands
    dynamic_invocations: Arc<Mutex<Vec<DynamicCommandInvocation>>>,
    pending_registrations: Vec<DynCmdRegistration>,
    pending_unregistrations: Vec<String>,
    // Notification messages from plugins (cleared each poll cycle)
    notification_messages: Vec<String>,
    // Hotkeys
    /// Hotkey bindings triggered since last poll (filled by handle_key, drained by poll_all).
    triggered_hotkeys: Vec<KeyBinding>,
    pending_hotkey_registrations: Vec<(KeyBinding, PluginKeyAction)>,
    pending_hotkey_unregistrations: Vec<KeyBinding>,
    // Viewer mutations queued by plugins
    pending_mutations: Vec<ViewerMutation>,
}

impl PluginManager {
    pub fn new() -> Self {
        Self {
            plugins: Vec::new(),
            pending_executions: Vec::new(),
            command_results: Vec::new(),
            dynamic_invocations: Arc::new(Mutex::new(Vec::new())),
            pending_registrations: Vec::new(),
            pending_unregistrations: Vec::new(),
            notification_messages: Vec::new(),
            triggered_hotkeys: Vec::new(),
            pending_hotkey_registrations: Vec::new(),
            pending_hotkey_unregistrations: Vec::new(),
            pending_mutations: Vec::new(),
        }
    }

    /// Load all plugin libraries from a directory.
    ///
    /// Scans for `.dylib` (macOS), `.so` (Linux), or `.dll` (Windows) files
    /// and attempts to load each one.
    pub fn load_dir(
        &mut self,
        dir: &Path,
        executor: &mut CommandExecutor,
        components: &mut ComponentStore,
        layout: &mut Layout,
    ) {
        let entries = match std::fs::read_dir(dir) {
            Ok(entries) => entries,
            Err(e) => {
                log::debug!("Plugin directory {:?}: {}", dir, e);
                return;
            }
        };

        for entry in entries.flatten() {
            let path = entry.path();
            let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("");

            let is_plugin = match std::env::consts::OS {
                "macos" => ext == "dylib",
                "linux" => ext == "so",
                "windows" => ext == "dll",
                _ => false,
            };

            if is_plugin {
                match self.load_library(&path, executor, components, layout) {
                    Ok(name) => log::info!("Loaded plugin: {}", name),
                    Err(e) => log::warn!("Failed to load plugin {:?}: {}", path, e),
                }
            }
        }
    }

    /// Load a single plugin library.
    pub fn load_library(
        &mut self,
        path: &Path,
        executor: &mut CommandExecutor,
        components: &mut ComponentStore,
        layout: &mut Layout,
    ) -> Result<String, String> {
        // Safety: loading a shared library is inherently unsafe.
        // We require same-compiler builds and check ABI + SDK versions.
        let library = unsafe {
            Library::new(path).map_err(|e| format!("Failed to load library: {}", e))?
        };

        // Look up the plugin declaration symbol
        let declaration = unsafe {
            let symbol = library
                .get::<*const pymol_plugin::ffi::PluginDeclaration>(b"PYMOL_PLUGIN_DECLARATION")
                .map_err(|e| format!("Missing PYMOL_PLUGIN_DECLARATION symbol: {}", e))?;
            &**symbol
        };

        // Check ABI version
        if declaration.abi_version != ABI_VERSION {
            return Err(format!(
                "ABI version mismatch: plugin has {}, host expects {}",
                declaration.abi_version, ABI_VERSION
            ));
        }

        // Check SDK version
        let plugin_sdk_version = unsafe {
            std::str::from_utf8(std::slice::from_raw_parts(
                declaration.sdk_version_ptr,
                declaration.sdk_version_len,
            ))
            .map_err(|_| "Invalid SDK version string")?
        };
        if plugin_sdk_version != SDK_VERSION {
            return Err(format!(
                "SDK version mismatch: plugin has {}, host expects {}",
                plugin_sdk_version, SDK_VERSION
            ));
        }

        // Forward the host's logger so log macros work inside the plugin
        let init_result = catch_unwind(AssertUnwindSafe(|| unsafe {
            (declaration.init)(log::logger(), log::max_level());
        }));
        if let Err(panic_info) = init_result {
            return Err(format!(
                "Plugin panicked during init: {}",
                panic_payload_to_string(&panic_info)
            ));
        }

        // Create registrar and call the plugin's register function
        let mut registrar = PluginRegistrar::new();
        let register_result = catch_unwind(AssertUnwindSafe(|| unsafe {
            (declaration.register)(&mut registrar as *mut PluginRegistrar);
        }));
        if let Err(panic_info) = register_result {
            return Err(format!(
                "Plugin panicked during register: {}",
                panic_payload_to_string(&panic_info)
            ));
        }

        // Extract metadata
        let metadata = registrar.take_metadata().ok_or("Plugin did not set metadata")?;
        let plugin_name = format!("{} v{}", metadata.name, metadata.version);

        // Register commands
        for cmd in registrar.drain_commands() {
            executor.registry_mut().register_boxed(cmd);
        }

        // Register script handlers
        for (ext, handler) in registrar.drain_script_handlers() {
            executor.register_script_handler(ext, handler);
        }

        // Register format handlers (for load/save)
        for handler in registrar.drain_format_handlers() {
            executor.register_format_handler(handler);
        }

        // Register components
        for (component, config) in registrar.drain_components() {
            let id = component.id().to_string();
            let panel_slot = config.to_panel_slot(id);
            components.add_boxed(component);
            layout.panels.push(panel_slot);
        }

        // Register hotkeys
        let mut hotkeys = KeyBindings::new();
        for (key, action) in registrar.drain_hotkeys() {
            hotkeys.bind(key, action);
        }

        // Store plugin
        self.plugins.push(LoadedPlugin {
            _library: library,
            _metadata: metadata,
            message_handler: registrar.take_message_handler(),
            faulted: false,
            hotkeys,
        });

        Ok(plugin_name)
    }

    /// Broadcast a message to all plugin message handlers.
    pub fn broadcast(&mut self, msg: &AppMessage, bus: &mut MessageBus) {
        for plugin in &mut self.plugins {
            if plugin.faulted {
                continue;
            }
            if let Some(handler) = &mut plugin.message_handler {
                let result = catch_unwind(AssertUnwindSafe(|| {
                    handler.on_message(msg, bus);
                }));
                if let Err(panic_info) = result {
                    log::error!(
                        "Plugin '{}' panicked during on_message: {}. Plugin disabled.",
                        plugin._metadata.name,
                        panic_payload_to_string(&panic_info),
                    );
                    plugin.faulted = true;
                }
            }
        }
    }

    /// Whether any loaded plugin needs periodic polling.
    pub fn any_needs_poll(&self) -> bool {
        self.plugins.iter().any(|p| {
            !p.faulted
                && p.message_handler
                    .as_ref()
                    .map_or(false, |h| h.needs_poll())
        })
    }

    /// Poll all handlers that need it (Phase 1 of the three-phase cycle).
    ///
    /// After this call, use `take_pending_executions()`, etc. to process
    /// the queued requests in subsequent phases.
    pub fn poll_all(&mut self, shared: &SharedContext<'_>, bus: &mut MessageBus) {
        // Drain invocations captured since last poll
        let invocations: Vec<DynamicCommandInvocation> = self
            .dynamic_invocations
            .lock()
            .map(|mut v| std::mem::take(&mut *v))
            .unwrap_or_default();

        // Swap out results so plugins can read them
        let results = std::mem::take(&mut self.command_results);

        // Drain triggered hotkeys
        let triggered = std::mem::take(&mut self.triggered_hotkeys);

        let mut exec_queue = Vec::new();
        let mut reg_queue = Vec::new();
        let mut unreg_queue = Vec::new();
        let mut notification_queue = Vec::new();
        let mut hotkey_reg_queue = Vec::new();
        let mut hotkey_unreg_queue = Vec::new();
        let mut mutation_queue = Vec::new();

        // Phase 0: Execute triggered hotkey callbacks
        if !triggered.is_empty() {
            let mut ctx = PollContext::new(
                shared,
                bus,
                &results,
                &invocations,
                &triggered,
                &mut exec_queue,
                &mut reg_queue,
                &mut unreg_queue,
                &mut notification_queue,
                &mut hotkey_reg_queue,
                &mut hotkey_unreg_queue,
                &mut mutation_queue,
            );

            for binding in &triggered {
                for plugin in &mut self.plugins {
                    if plugin.faulted {
                        continue;
                    }
                    if let Some(PluginKeyAction::Callback(cb)) = plugin.hotkeys.get_mut(binding) {
                        let result = catch_unwind(AssertUnwindSafe(|| {
                            cb(&mut ctx);
                        }));
                        if let Err(panic_info) = result {
                            log::error!(
                                "Plugin '{}' panicked in hotkey callback: {}. Plugin disabled.",
                                plugin._metadata.name,
                                panic_payload_to_string(&panic_info),
                            );
                            plugin.faulted = true;
                        }
                        break; // first match wins
                    }
                }
            }
        }

        // Phase 1: Regular plugin polling
        let mut ctx = PollContext::new(
            shared,
            bus,
            &results,
            &invocations,
            &triggered,
            &mut exec_queue,
            &mut reg_queue,
            &mut unreg_queue,
            &mut notification_queue,
            &mut hotkey_reg_queue,
            &mut hotkey_unreg_queue,
            &mut mutation_queue,
        );

        for plugin in &mut self.plugins {
            if plugin.faulted {
                continue;
            }
            if let Some(handler) = &mut plugin.message_handler {
                if handler.needs_poll() {
                    let result = catch_unwind(AssertUnwindSafe(|| {
                        handler.poll(&mut ctx);
                    }));
                    if let Err(panic_info) = result {
                        log::error!(
                            "Plugin '{}' panicked during poll: {}. Plugin disabled.",
                            plugin._metadata.name,
                            panic_payload_to_string(&panic_info),
                        );
                        plugin.faulted = true;
                    }
                }
            }
        }

        // Store queued requests for the host to process
        self.pending_executions = exec_queue;
        self.pending_registrations = reg_queue;
        self.pending_unregistrations = unreg_queue;
        self.notification_messages = notification_queue;
        self.pending_hotkey_registrations = hotkey_reg_queue;
        self.pending_hotkey_unregistrations = hotkey_unreg_queue;
        self.pending_mutations = mutation_queue;
    }

    /// Take pending command execution requests (Phase 2 input).
    pub fn take_pending_executions(&mut self) -> Vec<CommandExecRequest> {
        std::mem::take(&mut self.pending_executions)
    }

    /// Store command results for delivery in the next poll cycle.
    pub fn store_command_results(&mut self, results: Vec<CommandResult>) {
        self.command_results = results;
    }

    /// Take pending viewer mutations (Phase 2 input, after commands).
    pub fn take_pending_mutations(&mut self) -> Vec<ViewerMutation> {
        std::mem::take(&mut self.pending_mutations)
    }

    /// Take pending dynamic command registrations (Phase 3 input).
    pub fn take_pending_registrations(&mut self) -> Vec<DynCmdRegistration> {
        std::mem::take(&mut self.pending_registrations)
    }

    /// Take pending dynamic command unregistrations (Phase 3 input).
    pub fn take_pending_unregistrations(&mut self) -> Vec<String> {
        std::mem::take(&mut self.pending_unregistrations)
    }

    /// Get a handle to the shared invocations sink for `DynamicCommand`.
    pub fn invocations_handle(&self) -> Arc<Mutex<Vec<DynamicCommandInvocation>>> {
        self.dynamic_invocations.clone()
    }

    /// Get the number of loaded plugins.
    pub fn plugin_count(&self) -> usize {
        self.plugins.len()
    }

    /// Get notification messages from plugins (set during the last poll cycle).
    pub fn notification_messages(&self) -> &[String] {
        &self.notification_messages
    }

    // =========================================================================
    // Hotkeys
    // =========================================================================

    /// Look up a plugin hotkey action. First loaded plugin wins.
    pub fn lookup_hotkey(&self, binding: &KeyBinding) -> Option<&PluginKeyAction> {
        for plugin in &self.plugins {
            if plugin.faulted {
                continue;
            }
            if let Some(action) = plugin.hotkeys.get(binding) {
                return Some(action);
            }
        }
        None
    }

    /// Queue a hotkey binding as triggered (called from handle_key).
    ///
    /// The callback will be executed during the next poll phase.
    pub fn trigger_hotkey(&mut self, binding: KeyBinding) {
        self.triggered_hotkeys.push(binding);
    }

    /// Whether any hotkeys were triggered since last poll.
    pub fn has_triggered_hotkeys(&self) -> bool {
        !self.triggered_hotkeys.is_empty()
    }

    /// Apply pending hotkey registration/unregistration changes (Phase 3).
    pub fn apply_hotkey_changes(&mut self) {
        let registrations = std::mem::take(&mut self.pending_hotkey_registrations);
        let unregistrations = std::mem::take(&mut self.pending_hotkey_unregistrations);

        // Apply unregistrations — remove from all plugins
        for key in &unregistrations {
            for plugin in &mut self.plugins {
                plugin.hotkeys.unbind(*key);
            }
        }

        // Apply registrations — bind to the last loaded plugin
        // (runtime registrations go to last plugin since we don't track which plugin requested it)
        if let Some(plugin) = self.plugins.last_mut() {
            for (key, action) in registrations {
                plugin.hotkeys.bind(key, action);
            }
        }
    }
}

impl Default for PluginManager {
    fn default() -> Self {
        Self::new()
    }
}

/// Extract a human-readable message from a panic payload.
fn panic_payload_to_string(payload: &Box<dyn Any + Send>) -> String {
    if let Some(s) = payload.downcast_ref::<&str>() {
        s.to_string()
    } else if let Some(s) = payload.downcast_ref::<String>() {
        s.clone()
    } else {
        "non-string panic".to_string()
    }
}
