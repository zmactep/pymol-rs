//! Side effects system for setting changes
//!
//! When settings change, other parts of the system may need to react:
//! - Re-render the scene
//! - Reload shaders
//! - Update UI
//! - Recalculate representations
//!
//! This module provides a trait-based callback system that allows other crates
//! to register handlers for setting changes.

use ahash::AHashMap;
use parking_lot::RwLock;
use std::sync::Arc;

use crate::setting::SettingValue;

// =============================================================================
// Side Effect Handler Trait
// =============================================================================

/// Trait for handling setting change side effects
pub trait SettingSideEffect: Send + Sync {
    /// Called when a setting changes
    ///
    /// # Arguments
    /// * `id` - The setting ID that changed
    /// * `value` - The new value (None if reset to default)
    /// * `selection` - Optional selection string (for per-object/atom changes)
    /// * `state` - State index (-1 for all states)
    fn on_change(
        &self,
        id: u16,
        value: Option<&SettingValue>,
        selection: Option<&str>,
        state: i32,
    );

    /// Returns a list of setting IDs this handler cares about
    /// Return None to receive all setting changes
    fn watched_settings(&self) -> Option<Vec<u16>> {
        None
    }
}

// =============================================================================
// Side Effect Categories
// =============================================================================

/// Categories of side effects for batch processing
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SideEffectCategory {
    /// Scene needs to be invalidated/redrawn
    SceneInvalidate,
    /// Scene has changed (more significant than invalidate)
    SceneChanged,
    /// Shader variables need to be reloaded
    ShaderReload,
    /// Shader computation (lighting) needs to be recalculated
    ShaderComputeLighting,
    /// Ortho/UI needs to be marked dirty
    OrthoDirty,
    /// Sequence view needs update
    SeqChanged,
    /// Stereo mode update
    StereoUpdate,
    /// Representations need to be rebuilt
    RepresentationRebuild,
    /// Full rebuild of all objects
    FullRebuild,
    /// Viewport update
    ViewportUpdate,
}

/// Pre-defined side effect mappings for common settings
/// This mirrors PyMOL's SettingGenerateSideEffects function
#[allow(non_upper_case_globals)]
pub fn get_side_effects(id: u16) -> Vec<SideEffectCategory> {
    use crate::definitions::id::*;

    match id {
        // Stereo settings
        stereo | stereo_mode | anaglyph_mode => vec![
            SideEffectCategory::StereoUpdate,
            SideEffectCategory::ShaderReload,
            SideEffectCategory::OrthoDirty,
        ],

        // Lighting settings
        light_count | spec_count | precomputed_lighting => vec![
            SideEffectCategory::ShaderComputeLighting,
            SideEffectCategory::SceneInvalidate,
        ],

        dot_lighting | mesh_lighting | field_of_view | fog_start |
        two_sided_lighting | transparency_global_sort | dot_normals | mesh_normals => vec![
            SideEffectCategory::SceneInvalidate,
        ],

        // Sequence view settings
        seq_view | seq_view_label_spacing | seq_view_label_mode |
        seq_view_label_start | seq_view_format | seq_view_color => vec![
            SideEffectCategory::SeqChanged,
        ],

        // UI settings
        show_frame_rate | group_full_member_names | group_arrow_prefix => vec![
            SideEffectCategory::OrthoDirty,
        ],

        // Build settings
        defer_builds_mode => vec![
            SideEffectCategory::FullRebuild,
        ],

        // Grid settings
        grid_mode | grid_slot | grid_max => vec![
            SideEffectCategory::SceneChanged,
        ],

        // Pick surface
        pick_surface | pickable => vec![
            SideEffectCategory::SceneChanged,
        ],

        // Representation geometry settings - need rep rebuild
        sphere_scale | transparency | surface_quality |
        stick_radius | line_width |
        cartoon_fancy_helices | cartoon_fancy_sheets |
        cartoon_oval_width | cartoon_oval_length |
        cartoon_rect_width | cartoon_rect_length |
        cartoon_loop_radius | cartoon_dumbbell_width |
        cartoon_dumbbell_length | cartoon_dumbbell_radius |
        cartoon_round_helices | cartoon_sampling |
        cartoon_smooth_loops => vec![
            SideEffectCategory::RepresentationRebuild,
        ],

        // Representation color settings - need rep rebuild
        stick_color | line_color | cartoon_color |
        surface_color | mesh_color | sphere_color |
        ribbon_color => vec![
            SideEffectCategory::RepresentationRebuild,
        ],

        // Background color settings
        bg_rgb | bg_rgb_top | bg_rgb_bottom => vec![
            SideEffectCategory::ViewportUpdate,
        ],

        _ => vec![],
    }
}

// =============================================================================
// Side Effect Registry
// =============================================================================

/// Registry for side effect handlers
pub struct SideEffectRegistry {
    /// Handlers for specific settings (id -> handlers)
    specific: RwLock<AHashMap<u16, Vec<Arc<dyn SettingSideEffect>>>>,
    /// Handlers that receive all setting changes
    global: RwLock<Vec<Arc<dyn SettingSideEffect>>>,
    /// Category handlers
    category_handlers: RwLock<AHashMap<SideEffectCategory, Vec<Arc<dyn Fn() + Send + Sync>>>>,
}

impl SideEffectRegistry {
    /// Create a new empty registry
    pub fn new() -> Self {
        SideEffectRegistry {
            specific: RwLock::new(AHashMap::new()),
            global: RwLock::new(Vec::new()),
            category_handlers: RwLock::new(AHashMap::new()),
        }
    }

    /// Register a side effect handler
    pub fn register(&self, handler: Arc<dyn SettingSideEffect>) {
        if let Some(ids) = handler.watched_settings() {
            let mut specific = self.specific.write();
            for id in ids {
                specific.entry(id).or_default().push(handler.clone());
            }
        } else {
            self.global.write().push(handler);
        }
    }

    /// Register a handler for a specific setting
    pub fn register_for_setting(&self, id: u16, handler: Arc<dyn SettingSideEffect>) {
        self.specific.write().entry(id).or_default().push(handler);
    }

    /// Register a category handler
    pub fn register_category_handler<F>(&self, category: SideEffectCategory, handler: F)
    where
        F: Fn() + Send + Sync + 'static,
    {
        self.category_handlers
            .write()
            .entry(category)
            .or_default()
            .push(Arc::new(handler));
    }

    /// Unregister all handlers (useful for cleanup)
    pub fn clear(&self) {
        self.specific.write().clear();
        self.global.write().clear();
        self.category_handlers.write().clear();
    }

    /// Trigger side effects for a setting change
    pub fn trigger(
        &self,
        id: u16,
        value: Option<&SettingValue>,
        selection: Option<&str>,
        state: i32,
    ) {
        // Call specific handlers
        if let Some(handlers) = self.specific.read().get(&id) {
            for handler in handlers {
                handler.on_change(id, value, selection, state);
            }
        }

        // Call global handlers
        for handler in self.global.read().iter() {
            handler.on_change(id, value, selection, state);
        }

        // Trigger category handlers based on the setting
        let categories = get_side_effects(id);
        let category_handlers = self.category_handlers.read();
        for category in categories {
            if let Some(handlers) = category_handlers.get(&category) {
                for handler in handlers {
                    handler();
                }
            }
        }
    }

    /// Trigger side effects for multiple setting changes (batch)
    pub fn trigger_batch(&self, changes: &[(u16, Option<&SettingValue>, Option<&str>, i32)]) {
        // Collect unique categories to avoid duplicate triggers
        let mut triggered_categories = ahash::AHashSet::new();

        for (id, value, selection, state) in changes {
            // Call specific handlers
            if let Some(handlers) = self.specific.read().get(id) {
                for handler in handlers {
                    handler.on_change(*id, *value, *selection, *state);
                }
            }

            // Call global handlers
            for handler in self.global.read().iter() {
                handler.on_change(*id, *value, *selection, *state);
            }

            // Collect categories
            for category in get_side_effects(*id) {
                triggered_categories.insert(category);
            }
        }

        // Trigger category handlers once per category
        let category_handlers = self.category_handlers.read();
        for category in triggered_categories {
            if let Some(handlers) = category_handlers.get(&category) {
                for handler in handlers {
                    handler();
                }
            }
        }
    }
}

impl Default for SideEffectRegistry {
    fn default() -> Self {
        Self::new()
    }
}

// =============================================================================
// Convenience Side Effect Handler Implementations
// =============================================================================

/// A simple closure-based side effect handler
pub struct FnSideEffect<F>
where
    F: Fn(u16, Option<&SettingValue>, Option<&str>, i32) + Send + Sync,
{
    handler: F,
    watched: Option<Vec<u16>>,
}

impl<F> FnSideEffect<F>
where
    F: Fn(u16, Option<&SettingValue>, Option<&str>, i32) + Send + Sync,
{
    /// Create a new closure-based handler
    pub fn new(handler: F) -> Self {
        FnSideEffect {
            handler,
            watched: None,
        }
    }

    /// Create a handler that only watches specific settings
    pub fn watching(mut self, settings: Vec<u16>) -> Self {
        self.watched = Some(settings);
        self
    }
}

impl<F> SettingSideEffect for FnSideEffect<F>
where
    F: Fn(u16, Option<&SettingValue>, Option<&str>, i32) + Send + Sync,
{
    fn on_change(
        &self,
        id: u16,
        value: Option<&SettingValue>,
        selection: Option<&str>,
        state: i32,
    ) {
        (self.handler)(id, value, selection, state);
    }

    fn watched_settings(&self) -> Option<Vec<u16>> {
        self.watched.clone()
    }
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicU32, Ordering};

    #[test]
    fn test_registry_basic() {
        let registry = SideEffectRegistry::new();
        let counter = Arc::new(AtomicU32::new(0));
        let counter_clone = counter.clone();

        let handler = FnSideEffect::new(move |_id, _value, _sel, _state| {
            counter_clone.fetch_add(1, Ordering::SeqCst);
        });

        registry.register(Arc::new(handler));
        registry.trigger(1, None, None, 0);

        assert_eq!(counter.load(Ordering::SeqCst), 1);
    }

    #[test]
    fn test_registry_specific_setting() {
        let registry = SideEffectRegistry::new();
        let counter = Arc::new(AtomicU32::new(0));
        let counter_clone = counter.clone();

        let handler = FnSideEffect::new(move |_id, _value, _sel, _state| {
            counter_clone.fetch_add(1, Ordering::SeqCst);
        })
        .watching(vec![42]);

        registry.register(Arc::new(handler));

        // Should not trigger for setting 1
        registry.trigger(1, None, None, 0);
        assert_eq!(counter.load(Ordering::SeqCst), 0);

        // Should trigger for setting 42
        registry.trigger(42, None, None, 0);
        assert_eq!(counter.load(Ordering::SeqCst), 1);
    }

    #[test]
    fn test_category_handlers() {
        let registry = SideEffectRegistry::new();
        let counter = Arc::new(AtomicU32::new(0));
        let counter_clone = counter.clone();

        registry.register_category_handler(SideEffectCategory::SceneInvalidate, move || {
            counter_clone.fetch_add(1, Ordering::SeqCst);
        });

        // Trigger a setting that causes SceneInvalidate
        registry.trigger(crate::definitions::id::dot_lighting, None, None, 0);

        assert_eq!(counter.load(Ordering::SeqCst), 1);
    }

    #[test]
    fn test_batch_trigger() {
        let registry = SideEffectRegistry::new();
        let counter = Arc::new(AtomicU32::new(0));
        let counter_clone = counter.clone();

        // Category handler should only be called once even for multiple settings
        registry.register_category_handler(SideEffectCategory::SceneInvalidate, move || {
            counter_clone.fetch_add(1, Ordering::SeqCst);
        });

        let changes: Vec<(u16, Option<&SettingValue>, Option<&str>, i32)> = vec![
            (crate::definitions::id::dot_lighting, None, None, 0),
            (crate::definitions::id::mesh_lighting, None, None, 0),
        ];

        registry.trigger_batch(&changes);

        // Should only be 1 because both settings trigger the same category
        assert_eq!(counter.load(Ordering::SeqCst), 1);
    }
}
