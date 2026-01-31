//! Named view management
//!
//! This module provides named view storage for storing and recalling
//! camera states, similar to PyMOL's `view` command.
//!
//! Unlike scenes, views only store the camera state (18 values):
//! rotation, position, origin, clipping planes, and FOV.

use ahash::AHashMap;

use crate::camera::{Camera, SceneView};
use crate::error::{SceneError, SceneResult};

/// View manager for storing and recalling named camera views
///
/// Views are simpler than scenes - they only store the camera state,
/// not colors, representations, or object visibility.
#[derive(Debug, Default, Clone)]
pub struct ViewManager {
    /// Stored views by name
    views: AHashMap<String, SceneView>,
    /// Ordered list of view names
    view_order: Vec<String>,
    /// Currently active view (if any)
    current: Option<String>,
}

impl ViewManager {
    /// Create a new empty view manager
    pub fn new() -> Self {
        Self::default()
    }

    /// Get the number of stored views
    pub fn len(&self) -> usize {
        self.views.len()
    }

    /// Check if there are no stored views
    pub fn is_empty(&self) -> bool {
        self.views.is_empty()
    }

    /// Store the current camera view under a name
    ///
    /// If a view with this name already exists, it will be updated.
    pub fn store(&mut self, key: &str, camera: &Camera) {
        let view = camera.current_view();

        // Update order if new view
        if !self.views.contains_key(key) {
            self.view_order.push(key.to_string());
        }

        self.views.insert(key.to_string(), view);
    }

    /// Store a view directly (without capturing from camera)
    pub fn store_view(&mut self, key: &str, view: SceneView) {
        if !self.views.contains_key(key) {
            self.view_order.push(key.to_string());
        }
        self.views.insert(key.to_string(), view);
    }

    /// Recall a view by name, applying it to the camera
    ///
    /// # Arguments
    /// * `key` - View name to recall
    /// * `camera` - Camera to apply the view to
    /// * `animate` - Animation duration in seconds (0 = instant)
    pub fn recall(
        &mut self,
        key: &str,
        camera: &mut Camera,
        animate: f32,
    ) -> SceneResult<()> {
        let view = self
            .views
            .get(key)
            .ok_or_else(|| SceneError::ViewNotFound(key.to_string()))?
            .clone();

        if animate > 0.0 {
            camera.animate_to(view, animate);
        } else {
            camera.set_view(view);
        }

        self.current = Some(key.to_string());
        Ok(())
    }

    /// Get a view by name without recalling it
    pub fn get(&self, key: &str) -> Option<&SceneView> {
        self.views.get(key)
    }

    /// Get a mutable view by name
    pub fn get_mut(&mut self, key: &str) -> Option<&mut SceneView> {
        self.views.get_mut(key)
    }

    /// Delete a view by name
    ///
    /// Returns the deleted view if it existed.
    pub fn delete(&mut self, key: &str) -> Option<SceneView> {
        self.view_order.retain(|k| k != key);
        if self.current.as_deref() == Some(key) {
            self.current = None;
        }
        self.views.remove(key)
    }

    /// Check if a view exists
    pub fn contains(&self, key: &str) -> bool {
        self.views.contains_key(key)
    }

    /// Get the list of view names in order
    pub fn list(&self) -> &[String] {
        &self.view_order
    }

    /// Get an iterator over all views in order
    pub fn iter(&self) -> impl Iterator<Item = (&str, &SceneView)> {
        self.view_order
            .iter()
            .filter_map(|k| self.views.get(k).map(|v| (k.as_str(), v)))
    }

    /// Get the currently active view name
    pub fn current(&self) -> Option<&str> {
        self.current.as_deref()
    }

    /// Get the next view in the list (wrapping around)
    pub fn next(&self) -> Option<&str> {
        if self.view_order.is_empty() {
            return None;
        }

        let current_idx = self
            .current
            .as_ref()
            .and_then(|c| self.view_order.iter().position(|k| k == c))
            .unwrap_or(self.view_order.len() - 1);

        let next_idx = (current_idx + 1) % self.view_order.len();
        Some(&self.view_order[next_idx])
    }

    /// Get the previous view in the list (wrapping around)
    pub fn prev(&self) -> Option<&str> {
        if self.view_order.is_empty() {
            return None;
        }

        let current_idx = self
            .current
            .as_ref()
            .and_then(|c| self.view_order.iter().position(|k| k == c))
            .unwrap_or(0);

        let prev_idx = if current_idx == 0 {
            self.view_order.len() - 1
        } else {
            current_idx - 1
        };

        Some(&self.view_order[prev_idx])
    }

    /// Recall the next view
    pub fn recall_next(&mut self, camera: &mut Camera, animate: f32) -> SceneResult<()> {
        let key = self
            .next()
            .ok_or_else(|| SceneError::ViewNotFound("no views".to_string()))?
            .to_string();
        self.recall(&key, camera, animate)
    }

    /// Recall the previous view
    pub fn recall_prev(&mut self, camera: &mut Camera, animate: f32) -> SceneResult<()> {
        let key = self
            .prev()
            .ok_or_else(|| SceneError::ViewNotFound("no views".to_string()))?
            .to_string();
        self.recall(&key, camera, animate)
    }

    /// Clear all views
    pub fn clear(&mut self) {
        self.views.clear();
        self.view_order.clear();
        self.current = None;
    }

    /// Rename a view
    pub fn rename(&mut self, old_key: &str, new_key: &str) -> SceneResult<()> {
        if !self.views.contains_key(old_key) {
            return Err(SceneError::ViewNotFound(old_key.to_string()));
        }

        if self.views.contains_key(new_key) {
            return Err(SceneError::ViewExists(new_key.to_string()));
        }

        if let Some(view) = self.views.remove(old_key) {
            self.views.insert(new_key.to_string(), view);

            // Update order
            for k in &mut self.view_order {
                if k == old_key {
                    *k = new_key.to_string();
                }
            }

            // Update current
            if self.current.as_deref() == Some(old_key) {
                self.current = Some(new_key.to_string());
            }
        }

        Ok(())
    }

    /// Reorder views
    pub fn set_order(&mut self, order: &[&str]) {
        let mut new_order = Vec::new();
        let mut used = std::collections::HashSet::new();

        // Add views in the specified order
        for key in order {
            if self.views.contains_key(*key) && !used.contains(*key) {
                new_order.push((*key).to_string());
                used.insert(*key);
            }
        }

        // Add any remaining views at the end
        for key in &self.view_order {
            if !used.contains(key.as_str()) {
                new_order.push(key.clone());
            }
        }

        self.view_order = new_order;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_view_store_recall() {
        let mut manager = ViewManager::new();
        let camera = Camera::new();

        manager.store("view1", &camera);
        assert_eq!(manager.len(), 1);
        assert!(manager.get("view1").is_some());
    }

    #[test]
    fn test_view_navigation() {
        let mut manager = ViewManager::new();
        let camera = Camera::new();

        manager.store("v1", &camera);
        manager.store("v2", &camera);
        manager.store("v3", &camera);

        assert_eq!(manager.next(), Some("v1"));

        manager.current = Some("v1".to_string());
        assert_eq!(manager.next(), Some("v2"));
        assert_eq!(manager.prev(), Some("v3")); // wraps around
    }

    #[test]
    fn test_view_delete() {
        let mut manager = ViewManager::new();
        let camera = Camera::new();

        manager.store("v1", &camera);
        manager.store("v2", &camera);

        manager.delete("v1");
        assert_eq!(manager.len(), 1);
        assert!(manager.get("v1").is_none());
        assert!(manager.get("v2").is_some());
    }

    #[test]
    fn test_view_rename() {
        let mut manager = ViewManager::new();
        let camera = Camera::new();

        manager.store("old", &camera);
        manager.rename("old", "new").unwrap();

        assert!(manager.get("old").is_none());
        assert!(manager.get("new").is_some());
    }

    #[test]
    fn test_view_clear() {
        let mut manager = ViewManager::new();
        let camera = Camera::new();

        manager.store("v1", &camera);
        manager.store("v2", &camera);

        manager.clear();
        assert!(manager.is_empty());
        assert_eq!(manager.list().len(), 0);
    }
}
