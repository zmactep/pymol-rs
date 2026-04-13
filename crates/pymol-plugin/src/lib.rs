//! PyMOL-RS Plugin SDK
//!
//! This crate provides everything a plugin author needs to extend PyMOL-RS
//! with new commands, GUI components, and message handlers.
//!
//! # Quick Start
//!
//! ```rust,ignore
//! use pymol_plugin::prelude::*;
//!
//! pymol_plugin! {
//!     name: "my-plugin",
//!     version: "0.1.0",
//!     description: "My awesome plugin",
//!     commands: [MyCommand],
//! }
//! ```
//!
//! # Architecture
//!
//! Plugins are compiled as `cdylib` shared libraries and loaded at runtime
//! by the host application. The `pymol_plugin!` macro generates a C-compatible
//! entry point that the host discovers via `libloading`.

pub mod ffi;
pub mod macros;
pub mod registrar;

// Re-export log so the pymol_plugin! macro can reference it without
// plugins needing an explicit `log` dependency.
pub use log;

/// Hidden re-exports used by `define_plugin_settings!` macro internals.
/// Not part of the public API.
#[doc(hidden)]
pub mod __private {
    pub use pymol_settings::{SettingType, SettingValue, SideEffectCategory};
}

/// Convenient re-exports for plugin authors.
///
/// A single `use pymol_plugin::prelude::*;` gives access to all types
/// needed to implement commands, components, and message handlers.
pub mod prelude {
    // Framework types
    pub use pymol_framework::component::{Component, EguiComponent, SharedContext};
    pub use pymol_framework::layout::{PanelConfig, Slot};
    pub use pymol_framework::message::{AppMessage, MessageBus};
    pub use pymol_framework::topics::{publish, subscribe};

    // Command system
    pub use pymol_cmd::{ArgHint, CmdError, CmdResult, Command, CommandContext, ParsedCommand, ViewerLike};

    // Domain types
    pub use pymol_mol::{Atom, Bond, CoordSet, ObjectMolecule};
    pub use pymol_scene::{Camera, ObjectRegistry, SelectionManager};
    pub use pymol_select::select;

    // Plugin API
    pub use crate::registrar::{
        parse_key_string, CommandResult, DynamicCommandInvocation, DynamicSettingDescriptor,
        DynamicSettingStore, FormatHandler, HotkeyCallback, KeyBinding, KeyCode,
        MessageHandler, PluginKeyAction, PluginMetadata, PluginReaderFn,
        PluginRegistrar, PluginWriterFn, PollContext, SharedSettingStore,
    };
    // Settings side effects (needed for define_plugin_settings! macro)
    pub use pymol_settings::SideEffectCategory;
    // pymol_plugin! macro is auto-exported via #[macro_export]
    // define_plugin_settings! macro is auto-exported via #[macro_export]
}
