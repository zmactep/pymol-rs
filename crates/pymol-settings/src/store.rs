//! Setting storage and retrieval with hierarchical inheritance

use ahash::AHashMap;
use parking_lot::RwLock;

use crate::definitions::{get_setting, SETTINGS, SETTING_COUNT};
use crate::error::SettingError;
use crate::setting::{SettingLevel, SettingValue};

// =============================================================================
// Global Settings Store
// =============================================================================

/// Global settings store - holds session-wide settings
/// Uses a Vec indexed by setting ID for O(1) access
#[derive(Debug)]
pub struct GlobalSettings {
    /// Current values (None = use default)
    values: Vec<Option<SettingValue>>,
    /// Tracks which settings have been modified
    changed: Vec<bool>,
}

impl GlobalSettings {
    /// Create a new global settings store with defaults
    pub fn new() -> Self {
        GlobalSettings {
            values: vec![None; SETTING_COUNT],
            changed: vec![false; SETTING_COUNT],
        }
    }

    /// Get a setting value, returning the default if not explicitly set
    pub fn get(&self, id: u16) -> Option<SettingValue> {
        let idx = id as usize;
        if idx >= SETTING_COUNT {
            return None;
        }

        self.values
            .get(idx)
            .and_then(|v| v.clone())
            .or_else(|| get_setting(id).map(|s| s.default.clone()))
    }

    /// Check if a setting is explicitly defined (not using default)
    pub fn is_defined(&self, id: u16) -> bool {
        let idx = id as usize;
        idx < SETTING_COUNT && self.values.get(idx).map_or(false, |v| v.is_some())
    }

    /// Get only if explicitly defined (not default)
    pub fn get_if_defined(&self, id: u16) -> Option<&SettingValue> {
        let idx = id as usize;
        if idx >= SETTING_COUNT {
            return None;
        }
        self.values.get(idx).and_then(|v| v.as_ref())
    }

    /// Set a setting value
    pub fn set(&mut self, id: u16, value: SettingValue) -> Result<(), SettingError> {
        let idx = id as usize;
        if idx >= SETTING_COUNT {
            return Err(SettingError::NotFound(format!("id:{}", id)));
        }

        let setting = get_setting(id).ok_or_else(|| SettingError::NotFound(format!("id:{}", id)))?;

        // Type check
        if !value.is_compatible_with(setting.setting_type) {
            return Err(SettingError::TypeMismatch {
                expected: setting.setting_type.to_string().leak(),
                actual: value.setting_type().to_string().leak(),
            });
        }

        // Range check for numeric types
        if let (Some(min), Some(max)) = (setting.min, setting.max) {
            if min != max {
                if let Some(v) = value.as_float() {
                    if v < min || v > max {
                        return Err(SettingError::InvalidValue {
                            name: setting.name.to_string(),
                            reason: format!("value {} is outside range [{}, {}]", v, min, max),
                        });
                    }
                }
            }
        }

        self.values[idx] = Some(value);
        self.changed[idx] = true;
        Ok(())
    }

    /// Unset a setting (revert to default)
    pub fn unset(&mut self, id: u16) -> bool {
        let idx = id as usize;
        if idx >= SETTING_COUNT {
            return false;
        }
        let was_set = self.values[idx].is_some();
        self.values[idx] = None;
        if was_set {
            self.changed[idx] = true;
        }
        was_set
    }

    /// Reset a setting to its default value
    pub fn reset(&mut self, id: u16) -> Result<(), SettingError> {
        self.unset(id);
        Ok(())
    }

    /// Reset all settings to defaults
    pub fn reset_all(&mut self) {
        for i in 0..SETTING_COUNT {
            if self.values[i].is_some() {
                self.values[i] = None;
                self.changed[i] = true;
            }
        }
    }

    /// Check if a setting has changed since last check
    pub fn check_changed(&mut self, id: u16) -> bool {
        let idx = id as usize;
        if idx >= SETTING_COUNT {
            return false;
        }
        let changed = self.changed[idx];
        self.changed[idx] = false;
        changed
    }

    /// Get all changed setting IDs and clear changed flags
    pub fn get_changed_and_clear(&mut self) -> Vec<u16> {
        let mut changed = Vec::new();
        for (idx, c) in self.changed.iter_mut().enumerate() {
            if *c {
                changed.push(idx as u16);
                *c = false;
            }
        }
        changed
    }

    /// Iterate over all defined (non-default) settings
    pub fn iter_defined(&self) -> impl Iterator<Item = (u16, &SettingValue)> {
        self.values
            .iter()
            .enumerate()
            .filter_map(|(idx, v)| v.as_ref().map(|val| (idx as u16, val)))
    }

    // =========================================================================
    // Type-Safe Getters
    // =========================================================================

    /// Get a boolean setting
    pub fn get_bool(&self, id: u16) -> bool {
        self.get(id).and_then(|v| v.as_bool()).unwrap_or(false)
    }

    /// Get an integer setting
    pub fn get_int(&self, id: u16) -> i32 {
        self.get(id).and_then(|v| v.as_int()).unwrap_or(0)
    }

    /// Get a float setting (with automatic int-to-float coercion)
    pub fn get_float(&self, id: u16) -> f32 {
        self.get(id).and_then(|v| v.as_float()).unwrap_or(0.0)
    }

    /// Get a float3 setting
    pub fn get_float3(&self, id: u16) -> [f32; 3] {
        self.get(id)
            .and_then(|v| v.as_float3())
            .unwrap_or([0.0, 0.0, 0.0])
    }

    /// Get a color index setting
    pub fn get_color(&self, id: u16) -> i32 {
        self.get(id).and_then(|v| v.as_color()).unwrap_or(-1)
    }

    /// Get a string setting
    pub fn get_string(&self, id: u16) -> String {
        self.get(id)
            .and_then(|v| v.as_string().map(|s| s.to_string()))
            .unwrap_or_default()
    }

    // =========================================================================
    // Type-Safe Setters
    // =========================================================================

    /// Set a boolean setting
    pub fn set_bool(&mut self, id: u16, value: bool) -> Result<(), SettingError> {
        self.set(id, SettingValue::Bool(value))
    }

    /// Set an integer setting
    pub fn set_int(&mut self, id: u16, value: i32) -> Result<(), SettingError> {
        self.set(id, SettingValue::Int(value))
    }

    /// Set a float setting
    pub fn set_float(&mut self, id: u16, value: f32) -> Result<(), SettingError> {
        self.set(id, SettingValue::Float(value))
    }

    /// Set a float3 setting
    pub fn set_float3(&mut self, id: u16, value: [f32; 3]) -> Result<(), SettingError> {
        self.set(id, SettingValue::Float3(value))
    }

    /// Set a color setting
    pub fn set_color(&mut self, id: u16, value: i32) -> Result<(), SettingError> {
        self.set(id, SettingValue::Color(value))
    }

    /// Set a string setting
    pub fn set_string(&mut self, id: u16, value: impl Into<String>) -> Result<(), SettingError> {
        self.set(id, SettingValue::String(value.into()))
    }
}

impl Default for GlobalSettings {
    fn default() -> Self {
        Self::new()
    }
}

impl Clone for GlobalSettings {
    fn clone(&self) -> Self {
        GlobalSettings {
            values: self.values.clone(),
            changed: vec![false; SETTING_COUNT], // Don't clone changed state
        }
    }
}

// =============================================================================
// Serialization Support
// =============================================================================

/// Settings that should not be saved in sessions (system-dependent)
const SESSION_BLACKLIST: &[u16] = &[
    // Add system-dependent settings here that shouldn't persist
    // e.g., internal_gui_width, window positions, etc.
];

/// Check if a setting should be excluded from sessions
fn is_session_blacklisted(id: u16) -> bool {
    SESSION_BLACKLIST.contains(&id)
}

/// Serialized setting entry
#[derive(Debug, Clone)]
pub struct SerializedSetting {
    pub id: u16,
    pub value: SettingValue,
}

impl GlobalSettings {
    /// Serialize settings to a list for session saving
    /// Only includes explicitly set (non-default) settings
    pub fn to_session_list(&self) -> Vec<SerializedSetting> {
        self.values
            .iter()
            .enumerate()
            .filter_map(|(idx, v)| {
                let id = idx as u16;
                if is_session_blacklisted(id) {
                    return None;
                }
                v.as_ref().map(|value| SerializedSetting {
                    id,
                    value: value.clone(),
                })
            })
            .collect()
    }

    /// Restore settings from a session list
    pub fn from_session_list(&mut self, list: &[SerializedSetting]) {
        for entry in list {
            if !is_session_blacklisted(entry.id) {
                let idx = entry.id as usize;
                if idx < SETTING_COUNT {
                    self.values[idx] = Some(entry.value.clone());
                    self.changed[idx] = true;
                }
            }
        }
    }

    /// Create a new GlobalSettings from a session list
    pub fn new_from_session_list(list: &[SerializedSetting]) -> Self {
        let mut settings = Self::new();
        settings.from_session_list(list);
        settings
    }

    /// Export all settings (including defaults) for debugging
    pub fn export_all(&self) -> Vec<(u16, String, SettingValue)> {
        SETTINGS
            .iter()
            .map(|s| {
                let value = self.get(s.id).unwrap_or_else(|| s.default.clone());
                (s.id, s.name.to_string(), value)
            })
            .collect()
    }
}

impl UniqueSettings {
    /// Serialize unique settings to a list
    pub fn to_session_list(&self) -> Vec<(UniqueId, u16, SettingValue)> {
        let mut result = Vec::new();
        for (&unique_id, &first_offset) in &self.id_to_offset {
            let mut offset = first_offset;
            loop {
                if let Some(entry) = &self.entries[offset] {
                    if !is_session_blacklisted(entry.setting_id) {
                        result.push((unique_id, entry.setting_id, entry.value.clone()));
                    }
                }
                offset = self.next[offset];
                if offset == 0 {
                    break;
                }
            }
        }
        result
    }

    /// Restore unique settings from a session list
    pub fn from_session_list(&mut self, list: &[(UniqueId, u16, SettingValue)]) {
        for (unique_id, setting_id, value) in list {
            if !is_session_blacklisted(*setting_id) {
                // We use a simplified set that ignores errors
                let _ = self.set(*unique_id, *setting_id, value.clone());
            }
        }
    }
}

// =============================================================================
// Unique Settings (Per-atom/bond/object-state)
// =============================================================================

/// Unique ID for per-atom/bond settings
pub type UniqueId = i32;

/// Entry in the unique settings storage
#[derive(Debug, Clone)]
struct UniqueEntry {
    setting_id: u16,
    value: SettingValue,
}

/// Storage for per-atom/bond settings using a sparse representation
/// Uses a linked-list style approach similar to PyMOL for memory efficiency
#[derive(Debug, Default)]
pub struct UniqueSettings {
    /// Map from unique_id to first entry index
    id_to_offset: AHashMap<UniqueId, usize>,
    /// Linked list of entries
    entries: Vec<Option<UniqueEntry>>,
    /// Next links for entries (0 = end of chain)
    next: Vec<usize>,
    /// Free list head
    next_free: usize,
}

impl UniqueSettings {
    /// Create a new unique settings store
    pub fn new() -> Self {
        let initial_capacity = 64;
        let mut settings = UniqueSettings {
            id_to_offset: AHashMap::new(),
            entries: Vec::with_capacity(initial_capacity),
            next: Vec::with_capacity(initial_capacity),
            next_free: 0,
        };
        // Initialize with some free entries
        for i in 0..initial_capacity {
            settings.entries.push(None);
            settings.next.push(if i + 1 < initial_capacity {
                i + 1
            } else {
                0
            });
        }
        settings.next_free = 0;
        settings
    }

    /// Expand the storage pool
    fn expand(&mut self) {
        let old_size = self.entries.len();
        let new_size = old_size + old_size / 2 + 16;

        for i in old_size..new_size {
            self.entries.push(None);
            self.next.push(if i + 1 < new_size { i + 1 } else { 0 });
        }
        self.next_free = old_size;
    }

    /// Allocate a new entry
    fn alloc(&mut self) -> usize {
        if self.next_free == 0 && self.entries.iter().all(|e| e.is_some()) {
            self.expand();
        }

        // Find a free slot
        if self.next_free != 0 || self.entries[0].is_none() {
            let offset = self.next_free;
            self.next_free = self.next[offset];
            offset
        } else {
            // Expand and try again
            self.expand();
            let offset = self.next_free;
            self.next_free = self.next[offset];
            offset
        }
    }

    /// Free an entry
    fn free(&mut self, offset: usize) {
        self.entries[offset] = None;
        self.next[offset] = self.next_free;
        self.next_free = offset;
    }

    /// Get a setting value for a unique ID
    pub fn get(&self, unique_id: UniqueId, setting_id: u16) -> Option<&SettingValue> {
        let mut offset = *self.id_to_offset.get(&unique_id)?;

        while offset != 0 || self.entries[offset].is_some() {
            if let Some(entry) = &self.entries[offset] {
                if entry.setting_id == setting_id {
                    return Some(&entry.value);
                }
            }
            offset = self.next[offset];
            if offset == 0 {
                break;
            }
        }
        None
    }

    /// Check if a setting is defined for a unique ID
    pub fn is_defined(&self, unique_id: UniqueId, setting_id: u16) -> bool {
        self.get(unique_id, setting_id).is_some()
    }

    /// Set a setting value for a unique ID
    pub fn set(
        &mut self,
        unique_id: UniqueId,
        setting_id: u16,
        value: SettingValue,
    ) -> Result<(), SettingError> {
        // Validate the setting
        let setting =
            get_setting(setting_id).ok_or_else(|| SettingError::NotFound(format!("id:{}", setting_id)))?;

        if !value.is_compatible_with(setting.setting_type) {
            return Err(SettingError::TypeMismatch {
                expected: setting.setting_type.to_string().leak(),
                actual: value.setting_type().to_string().leak(),
            });
        }

        // Check if we already have this setting for this unique_id
        if let Some(&first_offset) = self.id_to_offset.get(&unique_id) {
            let mut offset = first_offset;
            loop {
                if let Some(entry) = &mut self.entries[offset] {
                    if entry.setting_id == setting_id {
                        entry.value = value;
                        return Ok(());
                    }
                }
                let next = self.next[offset];
                if next == 0 {
                    // Append to end of chain
                    let new_offset = self.alloc();
                    self.entries[new_offset] = Some(UniqueEntry { setting_id, value });
                    self.next[new_offset] = 0;
                    self.next[offset] = new_offset;
                    return Ok(());
                }
                offset = next;
            }
        } else {
            // New unique_id
            let offset = self.alloc();
            self.entries[offset] = Some(UniqueEntry { setting_id, value });
            self.next[offset] = 0;
            self.id_to_offset.insert(unique_id, offset);
            Ok(())
        }
    }

    /// Unset a setting for a unique ID
    pub fn unset(&mut self, unique_id: UniqueId, setting_id: u16) -> bool {
        let first_offset = match self.id_to_offset.get(&unique_id) {
            Some(&o) => o,
            None => return false,
        };

        let mut prev = 0;
        let mut offset = first_offset;

        loop {
            if let Some(entry) = &self.entries[offset] {
                if entry.setting_id == setting_id {
                    let next = self.next[offset];

                    if prev == 0 {
                        // First entry in chain
                        if next == 0 {
                            // Only entry - remove the whole chain
                            self.id_to_offset.remove(&unique_id);
                        } else {
                            // Update chain head
                            self.id_to_offset.insert(unique_id, next);
                        }
                    } else {
                        // Middle or end of chain
                        self.next[prev] = next;
                    }

                    self.free(offset);
                    return true;
                }
            }

            let next = self.next[offset];
            if next == 0 {
                break;
            }
            prev = offset;
            offset = next;
        }

        false
    }

    /// Remove all settings for a unique ID
    pub fn remove_all(&mut self, unique_id: UniqueId) {
        let first_offset = match self.id_to_offset.remove(&unique_id) {
            Some(o) => o,
            None => return,
        };

        let mut offset = first_offset;
        loop {
            let next = self.next[offset];
            self.free(offset);
            if next == 0 {
                break;
            }
            offset = next;
        }
    }

    /// Get all settings for a unique ID
    pub fn get_all(&self, unique_id: UniqueId) -> Vec<(u16, &SettingValue)> {
        let mut result = Vec::new();
        let first_offset = match self.id_to_offset.get(&unique_id) {
            Some(&o) => o,
            None => return result,
        };

        let mut offset = first_offset;
        loop {
            if let Some(entry) = &self.entries[offset] {
                result.push((entry.setting_id, &entry.value));
            }
            let next = self.next[offset];
            if next == 0 {
                break;
            }
            offset = next;
        }
        result
    }

    /// Clear all unique settings
    pub fn clear(&mut self) {
        self.id_to_offset.clear();
        for entry in &mut self.entries {
            *entry = None;
        }
        // Rebuild free list
        for (i, next) in self.next.iter_mut().enumerate() {
            *next = if i + 1 < self.entries.len() {
                i + 1
            } else {
                0
            };
        }
        self.next_free = 0;
    }
}

// =============================================================================
// Setting Resolver (Hierarchical Lookup)
// =============================================================================

/// A resolver that looks up settings with hierarchical inheritance
/// 
/// Inheritance order:
/// - Global < Object < ObjectState
/// - ObjectState < Atom < AtomState
/// - ObjectState < Bond < BondState
pub struct SettingResolver<'a> {
    /// Global settings (always present)
    pub global: &'a GlobalSettings,
    /// Object-level settings (optional)
    pub object: Option<&'a GlobalSettings>,
    /// Object-state-level settings (optional)
    pub state: Option<&'a GlobalSettings>,
    /// Atom/bond-level unique settings (optional)
    pub unique: Option<(&'a UniqueSettings, UniqueId)>,
}

impl<'a> SettingResolver<'a> {
    /// Create a resolver for global settings only
    pub fn global(global: &'a GlobalSettings) -> Self {
        SettingResolver {
            global,
            object: None,
            state: None,
            unique: None,
        }
    }

    /// Create a resolver with object-level settings
    pub fn with_object(global: &'a GlobalSettings, object: &'a GlobalSettings) -> Self {
        SettingResolver {
            global,
            object: Some(object),
            state: None,
            unique: None,
        }
    }

    /// Create a resolver with object and state-level settings
    pub fn with_state(
        global: &'a GlobalSettings,
        object: &'a GlobalSettings,
        state: &'a GlobalSettings,
    ) -> Self {
        SettingResolver {
            global,
            object: Some(object),
            state: Some(state),
            unique: None,
        }
    }

    /// Create a full resolver with unique settings
    pub fn with_unique(
        global: &'a GlobalSettings,
        object: &'a GlobalSettings,
        state: &'a GlobalSettings,
        unique: &'a UniqueSettings,
        unique_id: UniqueId,
    ) -> Self {
        SettingResolver {
            global,
            object: Some(object),
            state: Some(state),
            unique: Some((unique, unique_id)),
        }
    }

    /// Get a setting value with full inheritance resolution
    pub fn get(&self, id: u16) -> Option<SettingValue> {
        let setting = get_setting(id)?;

        // Check unique settings first (atom/bond level)
        if let Some((unique, unique_id)) = self.unique {
            if matches!(
                setting.level,
                SettingLevel::Atom
                    | SettingLevel::AtomState
                    | SettingLevel::Bond
                    | SettingLevel::BondState
            ) {
                if let Some(value) = unique.get(unique_id, id) {
                    return Some(value.clone());
                }
            }
        }

        // Check state settings
        if let Some(state) = self.state {
            if let Some(value) = state.get_if_defined(id) {
                return Some(value.clone());
            }
        }

        // Check object settings
        if let Some(object) = self.object {
            if let Some(value) = object.get_if_defined(id) {
                return Some(value.clone());
            }
        }

        // Fall back to global
        self.global.get(id)
    }

    /// Get a setting value, returning the default if not found at any level
    pub fn get_or_default(&self, id: u16) -> SettingValue {
        self.get(id)
            .or_else(|| get_setting(id).map(|s| s.default.clone()))
            .unwrap_or_default()
    }

    // =========================================================================
    // Type-Safe Getters
    // =========================================================================

    /// Get a boolean setting
    pub fn get_bool(&self, id: u16) -> bool {
        self.get_or_default(id).as_bool().unwrap_or(false)
    }

    /// Get an integer setting
    pub fn get_int(&self, id: u16) -> i32 {
        self.get_or_default(id).as_int().unwrap_or(0)
    }

    /// Get a float setting (with automatic int-to-float coercion)
    pub fn get_float(&self, id: u16) -> f32 {
        self.get_or_default(id).as_float().unwrap_or(0.0)
    }

    /// Get a float3 setting
    pub fn get_float3(&self, id: u16) -> [f32; 3] {
        self.get_or_default(id).as_float3().unwrap_or([0.0, 0.0, 0.0])
    }

    /// Get a color index setting
    pub fn get_color(&self, id: u16) -> i32 {
        self.get_or_default(id).as_color().unwrap_or(-1)
    }

    /// Get a string setting
    pub fn get_string(&self, id: u16) -> String {
        self.get_or_default(id)
            .as_string()
            .map(|s| s.to_string())
            .unwrap_or_default()
    }

    // =========================================================================
    // Optional Type-Safe Getters (only if defined)
    // =========================================================================

    /// Get a boolean setting if defined at any level
    pub fn get_bool_if_defined(&self, id: u16) -> Option<bool> {
        self.get(id).and_then(|v| v.as_bool())
    }

    /// Get an integer setting if defined at any level
    pub fn get_int_if_defined(&self, id: u16) -> Option<i32> {
        self.get(id).and_then(|v| v.as_int())
    }

    /// Get a float setting if defined at any level
    pub fn get_float_if_defined(&self, id: u16) -> Option<f32> {
        self.get(id).and_then(|v| v.as_float())
    }

    /// Get a float3 setting if defined at any level
    pub fn get_float3_if_defined(&self, id: u16) -> Option<[f32; 3]> {
        self.get(id).and_then(|v| v.as_float3())
    }

    /// Get a color index setting if defined at any level
    pub fn get_color_if_defined(&self, id: u16) -> Option<i32> {
        self.get(id).and_then(|v| v.as_color())
    }

    /// Get a string setting if defined at any level
    pub fn get_string_if_defined(&self, id: u16) -> Option<String> {
        self.get(id).and_then(|v| v.as_string().map(|s| s.to_string()))
    }
}

// =============================================================================
// Thread-Safe Setting Store
// =============================================================================

/// Thread-safe wrapper around GlobalSettings
#[derive(Debug)]
pub struct SettingStore {
    inner: RwLock<GlobalSettings>,
}

impl SettingStore {
    /// Create a new setting store
    pub fn new() -> Self {
        SettingStore {
            inner: RwLock::new(GlobalSettings::new()),
        }
    }

    /// Get a setting value
    pub fn get(&self, id: u16) -> Option<SettingValue> {
        self.inner.read().get(id)
    }

    /// Get a setting by name
    pub fn get_by_name(&self, name: &str) -> Result<SettingValue, SettingError> {
        let id = crate::definitions::get_setting_id(name)
            .ok_or_else(|| SettingError::NotFound(name.to_string()))?;
        self.get(id)
            .ok_or_else(|| SettingError::NotFound(name.to_string()))
    }

    /// Set a setting value
    pub fn set(&self, id: u16, value: SettingValue) -> Result<(), SettingError> {
        self.inner.write().set(id, value)
    }

    /// Set a setting by name
    pub fn set_by_name(&self, name: &str, value: SettingValue) -> Result<(), SettingError> {
        let id = crate::definitions::get_setting_id(name)
            .ok_or_else(|| SettingError::NotFound(name.to_string()))?;
        self.set(id, value)
    }

    /// Reset a setting to default
    pub fn reset(&self, id: u16) -> Result<(), SettingError> {
        self.inner.write().reset(id)
    }

    /// Reset all settings to defaults
    pub fn reset_all(&self) {
        self.inner.write().reset_all()
    }

    /// Get read access to the inner settings
    pub fn read(&self) -> parking_lot::RwLockReadGuard<'_, GlobalSettings> {
        self.inner.read()
    }

    /// Get write access to the inner settings
    pub fn write(&self) -> parking_lot::RwLockWriteGuard<'_, GlobalSettings> {
        self.inner.write()
    }

    // =========================================================================
    // Type-Safe Getters
    // =========================================================================

    /// Get a boolean setting
    pub fn get_bool(&self, id: u16) -> bool {
        self.inner.read().get_bool(id)
    }

    /// Get an integer setting
    pub fn get_int(&self, id: u16) -> i32 {
        self.inner.read().get_int(id)
    }

    /// Get a float setting
    pub fn get_float(&self, id: u16) -> f32 {
        self.inner.read().get_float(id)
    }

    /// Get a float3 setting
    pub fn get_float3(&self, id: u16) -> [f32; 3] {
        self.inner.read().get_float3(id)
    }

    /// Get a color index setting
    pub fn get_color(&self, id: u16) -> i32 {
        self.inner.read().get_color(id)
    }

    /// Get a string setting
    pub fn get_string(&self, id: u16) -> String {
        self.inner.read().get_string(id)
    }

    // =========================================================================
    // Type-Safe Setters
    // =========================================================================

    /// Set a boolean setting
    pub fn set_bool(&self, id: u16, value: bool) -> Result<(), SettingError> {
        self.inner.write().set_bool(id, value)
    }

    /// Set an integer setting
    pub fn set_int(&self, id: u16, value: i32) -> Result<(), SettingError> {
        self.inner.write().set_int(id, value)
    }

    /// Set a float setting
    pub fn set_float(&self, id: u16, value: f32) -> Result<(), SettingError> {
        self.inner.write().set_float(id, value)
    }

    /// Set a float3 setting
    pub fn set_float3(&self, id: u16, value: [f32; 3]) -> Result<(), SettingError> {
        self.inner.write().set_float3(id, value)
    }

    /// Set a color setting
    pub fn set_color(&self, id: u16, value: i32) -> Result<(), SettingError> {
        self.inner.write().set_color(id, value)
    }

    /// Set a string setting
    pub fn set_string(&self, id: u16, value: impl Into<String>) -> Result<(), SettingError> {
        self.inner.write().set_string(id, value)
    }
}

impl Default for SettingStore {
    fn default() -> Self {
        Self::new()
    }
}

impl Clone for SettingStore {
    fn clone(&self) -> Self {
        SettingStore {
            inner: RwLock::new(self.inner.read().clone()),
        }
    }
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;

    #[test]
    fn test_global_settings_default() {
        let settings = GlobalSettings::new();
        let value = settings.get(crate::definitions::id::sphere_quality);
        assert!(value.is_some());
    }

    #[test]
    fn test_global_settings_set_get() {
        let mut settings = GlobalSettings::new();
        settings
            .set(crate::definitions::id::ambient, SettingValue::Float(0.5))
            .unwrap();
        let value = settings.get(crate::definitions::id::ambient).unwrap();
        assert_eq!(value.as_float(), Some(0.5));
    }

    #[test]
    fn test_global_settings_reset() {
        let mut settings = GlobalSettings::new();
        settings
            .set(crate::definitions::id::ambient, SettingValue::Float(0.5))
            .unwrap();
        settings.reset(crate::definitions::id::ambient).unwrap();
        let value = settings.get(crate::definitions::id::ambient).unwrap();
        // Should be back to default (0.14)
        assert_eq!(value.as_float(), Some(0.14));
    }

    #[test]
    fn test_unique_settings() {
        let mut unique = UniqueSettings::new();
        unique
            .set(1, crate::definitions::id::sphere_quality, SettingValue::Int(3))
            .unwrap();

        let value = unique.get(1, crate::definitions::id::sphere_quality);
        assert!(value.is_some());
        assert_eq!(value.unwrap().as_int(), Some(3));

        // Different unique_id should not have this setting
        assert!(unique.get(2, crate::definitions::id::sphere_quality).is_none());
    }

    #[test]
    fn test_resolver_inheritance() {
        let mut global = GlobalSettings::new();
        let mut object = GlobalSettings::new();

        // Set different values at different levels
        global
            .set(crate::definitions::id::ambient, SettingValue::Float(0.1))
            .unwrap();
        object
            .set(crate::definitions::id::ambient, SettingValue::Float(0.2))
            .unwrap();

        // Global resolver
        let resolver = SettingResolver::global(&global);
        assert_eq!(
            resolver.get(crate::definitions::id::ambient).unwrap().as_float(),
            Some(0.1)
        );

        // Object resolver should override global
        let resolver = SettingResolver::with_object(&global, &object);
        assert_eq!(
            resolver.get(crate::definitions::id::ambient).unwrap().as_float(),
            Some(0.2)
        );
    }

    #[test]
    fn test_setting_store_thread_safe() {
        let store = Arc::new(SettingStore::new());
        let store2 = store.clone();

        std::thread::spawn(move || {
            store2
                .set(crate::definitions::id::ambient, SettingValue::Float(0.3))
                .unwrap();
        })
        .join()
        .unwrap();

        let value = store.get(crate::definitions::id::ambient).unwrap();
        assert_eq!(value.as_float(), Some(0.3));
    }
}
