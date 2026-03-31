//! Dynamic settings — runtime-registered settings for plugins and extensions.
//!
//! Unlike the typed settings groups (which are compile-time structs), dynamic
//! settings are registered at runtime and stored in a shared key-value store.
//! The host dispatches `set`/`get`/`unset` commands to the store via
//! [`SharedSettingStore`] (an `Arc<RwLock<DynamicSettingStore>>`).
//!
//! Values live on the plugin side. The host only holds an Arc reference
//! for command dispatch.

use std::sync::{Arc, RwLock};

use ahash::AHashMap;

use crate::setting::{SettingType, SettingValue};
use crate::side_effects::SideEffectCategory;

// =============================================================================
// Descriptor
// =============================================================================

/// Metadata for a dynamically registered setting.
///
/// Contains everything the `set`/`get` commands need to validate and display
/// a setting, but no function pointers — values are stored externally in a
/// [`DynamicSettingStore`].
#[derive(Debug, Clone)]
pub struct DynamicSettingDescriptor {
    /// Setting name (e.g., `"my_plugin_threshold"`).
    pub name: String,
    /// Value type.
    pub setting_type: SettingType,
    /// Default value (used by `get` when no value is stored, and by `unset`).
    pub default: SettingValue,
    /// Min constraint (numeric).
    pub min: Option<f32>,
    /// Max constraint (numeric).
    pub max: Option<f32>,
    /// Named value variants (owned, unlike the `&'static` slice in built-in descriptors).
    pub value_hints: Vec<(String, SettingValue)>,
    /// Side effects triggered when this setting changes.
    pub side_effects: Vec<SideEffectCategory>,
    /// Whether this setting supports per-object overrides.
    pub object_overridable: bool,
}

impl DynamicSettingDescriptor {
    /// Whether this setting has named value variants.
    pub fn has_value_hints(&self) -> bool {
        !self.value_hints.is_empty()
    }

    /// Look up a setting value from a string alias (case-insensitive).
    pub fn resolve_hint(&self, name: &str) -> Option<&SettingValue> {
        let name_lower = name.to_lowercase();
        self.value_hints
            .iter()
            .find(|(n, _)| n.to_lowercase() == name_lower)
            .map(|(_, v)| v)
    }

    /// Look up a display name from a setting value.
    pub fn hint_name(&self, value: &SettingValue) -> Option<&str> {
        self.value_hints
            .iter()
            .find(|(_, v)| v == value)
            .map(|(n, _)| n.as_str())
    }

    /// Get all hint names (for autocomplete).
    pub fn hint_names(&self) -> impl Iterator<Item = &str> + '_ {
        self.value_hints.iter().map(|(n, _)| n.as_str())
    }
}

// =============================================================================
// Store
// =============================================================================

/// Key-value store for dynamic setting values (global + per-object).
///
/// Owned by the plugin, shared with the host via [`SharedSettingStore`].
/// The host writes values when `set` commands target dynamic settings;
/// the plugin reads them during `poll()`.
#[derive(Debug, Clone, Default)]
pub struct DynamicSettingStore {
    /// Global values: setting name → value.
    global: AHashMap<String, SettingValue>,
    /// Per-object overrides: object name → (setting name → value).
    per_object: AHashMap<String, AHashMap<String, SettingValue>>,
}

impl DynamicSettingStore {
    /// Create an empty store.
    pub fn new() -> Self {
        Self::default()
    }

    // ----- Global -----

    /// Get a global setting value. Returns `None` if not explicitly set.
    pub fn get(&self, name: &str) -> Option<&SettingValue> {
        self.global.get(name)
    }

    /// Set a global setting value.
    pub fn set(&mut self, name: &str, value: SettingValue) {
        self.global.insert(name.to_string(), value);
    }

    /// Remove a global setting value (reverts to descriptor default).
    pub fn remove(&mut self, name: &str) -> Option<SettingValue> {
        self.global.remove(name)
    }

    // ----- Per-object -----

    /// Get a per-object override. Returns `None` if not set for this object.
    pub fn get_object(&self, obj_name: &str, setting_name: &str) -> Option<&SettingValue> {
        self.per_object
            .get(obj_name)
            .and_then(|m| m.get(setting_name))
    }

    /// Set a per-object override.
    pub fn set_object(&mut self, obj_name: &str, setting_name: &str, value: SettingValue) {
        self.per_object
            .entry(obj_name.to_string())
            .or_default()
            .insert(setting_name.to_string(), value);
    }

    /// Remove a per-object override.
    pub fn remove_object(&mut self, obj_name: &str, setting_name: &str) -> Option<SettingValue> {
        let removed = self
            .per_object
            .get_mut(obj_name)
            .and_then(|m| m.remove(setting_name));
        // Clean up empty sub-maps
        if let Some(m) = self.per_object.get(obj_name) {
            if m.is_empty() {
                self.per_object.remove(obj_name);
            }
        }
        removed
    }

    /// Iterate over all global values.
    pub fn iter(&self) -> impl Iterator<Item = (&str, &SettingValue)> {
        self.global.iter().map(|(k, v)| (k.as_str(), v))
    }
}

/// Thread-safe shared handle to a [`DynamicSettingStore`].
///
/// The plugin creates this and gives a clone to the host during registration.
/// Both sides share the same underlying store.
pub type SharedSettingStore = Arc<RwLock<DynamicSettingStore>>;

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    fn make_descriptor(name: &str) -> DynamicSettingDescriptor {
        DynamicSettingDescriptor {
            name: name.to_string(),
            setting_type: SettingType::Float,
            default: SettingValue::Float(0.5),
            min: Some(0.0),
            max: Some(1.0),
            value_hints: vec![],
            side_effects: vec![],
            object_overridable: false,
        }
    }

    #[test]
    fn test_store_global_get_set_remove() {
        let mut store = DynamicSettingStore::new();
        assert!(store.get("foo").is_none());

        store.set("foo", SettingValue::Float(0.7));
        assert_eq!(store.get("foo"), Some(&SettingValue::Float(0.7)));

        let removed = store.remove("foo");
        assert_eq!(removed, Some(SettingValue::Float(0.7)));
        assert!(store.get("foo").is_none());
    }

    #[test]
    fn test_store_per_object() {
        let mut store = DynamicSettingStore::new();

        store.set_object("mol1", "threshold", SettingValue::Float(0.3));
        assert_eq!(
            store.get_object("mol1", "threshold"),
            Some(&SettingValue::Float(0.3))
        );
        assert!(store.get_object("mol2", "threshold").is_none());

        store.remove_object("mol1", "threshold");
        assert!(store.get_object("mol1", "threshold").is_none());
        // Sub-map should be cleaned up
        assert!(!store.per_object.contains_key("mol1"));
    }

    #[test]
    fn test_shared_store() {
        let store: SharedSettingStore = Arc::new(RwLock::new(DynamicSettingStore::new()));
        let clone = store.clone();

        store.write().unwrap().set("x", SettingValue::Int(42));
        assert_eq!(
            clone.read().unwrap().get("x"),
            Some(&SettingValue::Int(42))
        );
    }

    #[test]
    fn test_descriptor_hints() {
        let desc = DynamicSettingDescriptor {
            name: "mode".to_string(),
            setting_type: SettingType::Int,
            default: SettingValue::Int(0),
            min: None,
            max: None,
            value_hints: vec![
                ("fast".to_string(), SettingValue::Int(0)),
                ("quality".to_string(), SettingValue::Int(1)),
            ],
            side_effects: vec![],
            object_overridable: false,
        };

        assert!(desc.has_value_hints());
        assert_eq!(desc.resolve_hint("Fast"), Some(&SettingValue::Int(0)));
        assert_eq!(desc.hint_name(&SettingValue::Int(1)), Some("quality"));
        assert!(desc.resolve_hint("unknown").is_none());
    }

    #[test]
    fn test_descriptor_no_hints() {
        let desc = make_descriptor("foo");
        assert!(!desc.has_value_hints());
        assert!(desc.resolve_hint("any").is_none());
        assert!(desc.hint_name(&SettingValue::Float(0.5)).is_none());
    }
}
