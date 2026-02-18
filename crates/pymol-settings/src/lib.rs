//! PyMOL-RS Settings System
//!
//! This crate provides the configuration and settings management for PyMOL-RS.
//! It supports global settings, per-object settings, per-atom/bond settings,
//! and hierarchical setting inheritance.
//!
//! # Architecture
//!
//! The settings system mirrors PyMOL's original design with 798 settings organized
//! into a hierarchical inheritance model:
//!
//! ```text
//! global < object < object-state
//!                   object-state < atom < atom-state
//!                   object-state < bond < bond-state
//! ```
//!
//! # Setting Types
//!
//! - `Bool` - Boolean values
//! - `Int` - Integer values
//! - `Float` - Single-precision floats
//! - `Float3` - Three-component float vectors (positions, directions)
//! - `Color` - Color indices (resolved via pymol-color)
//! - `String` - String values
//!
//! # Example
//!
//! ```rust
//! use pymol_settings::{GlobalSettings, SettingValue, id};
//!
//! let mut settings = GlobalSettings::new();
//!
//! // Set a float setting
//! settings.set_float(id::ambient, 0.2).unwrap();
//!
//! // Get with type-safe getter
//! let ambient = settings.get_float(id::ambient);
//! assert_eq!(ambient, 0.2);
//! ```

mod definitions;
mod error;
mod setting;
pub mod shading_mode;
mod side_effects;
mod store;

// Re-export main types
pub use definitions::{
    create_default_settings, get_setting, get_setting_id, get_string_default, SETTINGS,
    SETTING_COUNT,
};
pub use error::SettingError;
pub use setting::{Setting, SettingLevel, SettingType, SettingValue};
pub use side_effects::{
    FnSideEffect, SettingSideEffect, SideEffectCategory, SideEffectRegistry,
};
pub use shading_mode::ShadingMode;
pub use store::{
    GlobalSettings, SerializedSetting, SettingResolver, SettingStore, UniqueId, UniqueSettings,
};

/// Setting ID constants
pub mod id {
    pub use crate::definitions::id::*;
}

/// Re-export commonly used types for convenience
pub mod prelude {
    pub use crate::definitions::id;
    pub use crate::{
        GlobalSettings, Setting, SettingError, SettingLevel, SettingResolver, SettingStore,
        SettingType, SettingValue, UniqueSettings,
    };
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic_usage() {
        let mut settings = GlobalSettings::new();

        // Test setting and getting values
        settings.set_float(id::ambient, 0.25).unwrap();
        assert_eq!(settings.get_float(id::ambient), 0.25);

        // Test reset
        settings.reset(id::ambient).unwrap();
        assert_eq!(settings.get_float(id::ambient), 0.14); // default
    }

    #[test]
    fn test_resolver_chain() {
        let mut global = GlobalSettings::new();
        let mut object = GlobalSettings::new();

        global.set_float(id::ambient, 0.1).unwrap();
        object.set_float(id::ambient, 0.2).unwrap();

        let resolver = SettingResolver::with_object(&global, &object);
        assert_eq!(resolver.get_float(id::ambient), 0.2);
    }

    #[test]
    fn test_serialization_roundtrip() {
        let mut settings = GlobalSettings::new();
        settings.set_float(id::ambient, 0.3).unwrap();
        settings.set_bool(id::orthoscopic, true).unwrap();

        let list = settings.to_session_list();
        assert!(!list.is_empty());

        let mut restored = GlobalSettings::new();
        restored.from_session_list(&list);

        assert_eq!(restored.get_float(id::ambient), 0.3);
        assert_eq!(restored.get_bool(id::orthoscopic), true);
    }
}
