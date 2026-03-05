//! Plugin Manager
//!
//! Loads plugin shared libraries at runtime and integrates their commands,
//! components, and message handlers into the application.

use std::path::Path;

use libloading::Library;

use pymol_cmd::CommandRegistry;
use pymol_plugin::ffi::{ABI_VERSION, SDK_VERSION};
use pymol_plugin::registrar::{MessageHandler, PluginMetadata, PluginRegistrar};

use pymol_framework::component_store::ComponentStore;
use pymol_framework::message::{AppMessage, MessageBus};

use crate::layout::Layout;

/// A loaded plugin, keeping the shared library alive.
struct LoadedPlugin {
    /// Must stay alive for the duration of the plugin.
    _library: Library,
    _metadata: PluginMetadata,
    message_handler: Option<Box<dyn MessageHandler>>,
}

/// Manages dynamically loaded plugins.
pub struct PluginManager {
    plugins: Vec<LoadedPlugin>,
}

impl PluginManager {
    pub fn new() -> Self {
        Self {
            plugins: Vec::new(),
        }
    }

    /// Load all plugin libraries from a directory.
    ///
    /// Scans for `.dylib` (macOS), `.so` (Linux), or `.dll` (Windows) files
    /// and attempts to load each one.
    pub fn load_dir(
        &mut self,
        dir: &Path,
        registry: &mut CommandRegistry,
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
                match self.load_library(&path, registry, components, layout) {
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
        registry: &mut CommandRegistry,
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

        // Create registrar and call the plugin's register function
        let mut registrar = PluginRegistrar::new();
        unsafe {
            (declaration.register)(&mut registrar as *mut PluginRegistrar);
        }

        // Extract metadata
        let metadata = registrar.take_metadata().ok_or("Plugin did not set metadata")?;
        let plugin_name = format!("{} v{}", metadata.name, metadata.version);

        // Register commands
        for cmd in registrar.drain_commands() {
            registry.register_boxed(cmd);
        }

        // Register components
        for (component, config) in registrar.drain_components() {
            let id = component.id().to_string();
            let panel_slot = config.to_panel_slot(id);
            components.add_boxed(component);
            layout.panels.push(panel_slot);
        }

        // Store plugin
        self.plugins.push(LoadedPlugin {
            _library: library,
            _metadata: metadata,
            message_handler: registrar.take_message_handler(),
        });

        Ok(plugin_name)
    }

    /// Broadcast a message to all plugin message handlers.
    pub fn broadcast(&mut self, msg: &AppMessage, bus: &mut MessageBus) {
        for plugin in &mut self.plugins {
            if let Some(handler) = &mut plugin.message_handler {
                handler.on_message(msg, bus);
            }
        }
    }

    /// Get the number of loaded plugins.
    pub fn plugin_count(&self) -> usize {
        self.plugins.len()
    }
}

impl Default for PluginManager {
    fn default() -> Self {
        Self::new()
    }
}
