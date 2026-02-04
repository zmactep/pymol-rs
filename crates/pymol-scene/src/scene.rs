//! Scene snapshot management
//!
//! This module provides named scene snapshots for storing and recalling
//! visualization states, similar to PyMOL's `scene` command.

use ahash::AHashMap;
use bitflags::bitflags;
use pymol_color::ColorIndex;
use pymol_mol::RepMask;
use pymol_settings::UniqueId;

use crate::camera::{Camera, SceneView};
use crate::error::{SceneError, SceneResult};
use crate::object::{DirtyFlags, ObjectRegistry, ObjectType};

bitflags! {
    /// Flags indicating what to store in a scene
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub struct SceneStoreMask: u32 {
        /// Store camera view
        const VIEW = 0x01;
        /// Store object enabled/disabled state
        const ACTIVE = 0x02;
        /// Store colors
        const COLOR = 0x04;
        /// Store representation visibility
        const REP = 0x08;
        /// Store frame/state index
        const FRAME = 0x10;
        /// Store thumbnail image
        const THUMBNAIL = 0x20;
        /// Store everything
        const ALL = 0x3F;
    }
}

impl Default for SceneStoreMask {
    fn default() -> Self {
        SceneStoreMask::ALL
    }
}

/// Per-atom stored properties in a scene
#[derive(Debug, Clone)]
pub struct SceneAtomData {
    /// Atom color
    pub color: ColorIndex,
    /// Visible representations
    pub visible_reps: RepMask,
}

/// Per-atom stored properties within an object
#[derive(Debug, Clone)]
pub struct ScenePerAtomData {
    /// Atom color (i32 matching Atom.color)
    pub color: i32,
    /// Visible representations
    pub visible_reps: RepMask,
}

/// Per-object stored properties in a scene
#[derive(Debug, Clone)]
pub struct SceneObjectData {
    /// Object enabled state
    pub enabled: bool,
    /// Object color
    pub color: ColorIndex,
    /// Visible representations
    pub visible_reps: RepMask,
    /// Current state index
    pub current_state: usize,
    /// Per-atom data (indexed parallel to molecule atoms)
    /// Only populated for molecule objects when COLOR or REP flags are set
    pub per_atom_data: Vec<ScenePerAtomData>,
}

/// A stored scene snapshot
///
/// A scene captures the visualization state at a particular moment,
/// including camera position, object visibility, colors, and representations.
#[derive(Debug, Clone)]
pub struct Scene {
    /// Scene name/key
    pub name: String,
    /// What was stored in this scene
    pub storemask: SceneStoreMask,
    /// Camera view state
    pub view: SceneView,
    /// Global frame/state index
    pub frame: usize,
    /// Message to display when scene is recalled
    pub message: String,
    /// Optional thumbnail image (PNG data)
    pub thumbnail: Option<Vec<u8>>,
    /// Per-atom data (keyed by unique_id)
    pub atom_data: AHashMap<UniqueId, SceneAtomData>,
    /// Per-object data (keyed by object name)
    pub object_data: AHashMap<String, SceneObjectData>,
}

impl Scene {
    /// Create a new empty scene
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            storemask: SceneStoreMask::default(),
            view: SceneView::default(),
            frame: 0,
            message: String::new(),
            thumbnail: None,
            atom_data: AHashMap::new(),
            object_data: AHashMap::new(),
        }
    }

    /// Create a scene by capturing current state
    pub fn capture(
        name: impl Into<String>,
        storemask: SceneStoreMask,
        camera: &Camera,
        registry: &ObjectRegistry,
    ) -> Self {
        let mut scene = Self::new(name);
        scene.storemask = storemask;

        // Capture view
        if storemask.contains(SceneStoreMask::VIEW) {
            scene.view = camera.current_view();
        }

        // Capture object states
        if storemask.intersects(SceneStoreMask::ACTIVE | SceneStoreMask::COLOR | SceneStoreMask::REP | SceneStoreMask::FRAME) {
            for obj in registry.iter() {
                let state = obj.state();
                
                // Capture per-atom data for molecule objects
                let per_atom_data = if storemask.intersects(SceneStoreMask::COLOR | SceneStoreMask::REP)
                    && obj.object_type() == ObjectType::Molecule
                {
                    if let Some(mol_obj) = registry.get_molecule(obj.name()) {
                        mol_obj.molecule().atoms().map(|atom| {
                            ScenePerAtomData {
                                color: atom.repr.colors.base,
                                visible_reps: atom.repr.visible_reps,
                            }
                        }).collect()
                    } else {
                        Vec::new()
                    }
                } else {
                    Vec::new()
                };
                
                let obj_data = SceneObjectData {
                    enabled: state.enabled,
                    color: state.color,
                    visible_reps: state.visible_reps,
                    current_state: obj.current_state(),
                    per_atom_data,
                };
                scene.object_data.insert(obj.name().to_string(), obj_data);
            }
        }

        scene
    }

    /// Apply this scene to the camera and registry
    pub fn apply(
        &self,
        camera: &mut Camera,
        registry: &mut ObjectRegistry,
        animate: bool,
        duration: f32,
    ) {
        // Apply view
        if self.storemask.contains(SceneStoreMask::VIEW) {
            if animate && duration > 0.0 {
                camera.animate_to(self.view.clone(), duration);
            } else {
                camera.set_view(self.view.clone());
            }
        }

        // Collect names first to avoid borrow conflicts
        let names: Vec<_> = self.object_data.keys().cloned().collect();

        // Apply object states
        for name in &names {
            let obj_data = &self.object_data[name];
            
            // First, apply per-atom data if this is a molecule
            if !obj_data.per_atom_data.is_empty() {
                if let Some(mol_obj) = registry.get_molecule_mut(name) {
                    let atoms = mol_obj.molecule_mut().atoms_slice_mut();
                    let per_atom = &obj_data.per_atom_data;
                    
                    // Only apply if atom counts match (safety check)
                    if atoms.len() == per_atom.len() {
                        for (atom, stored) in atoms.iter_mut().zip(per_atom.iter()) {
                            if self.storemask.contains(SceneStoreMask::COLOR) {
                                atom.repr.colors.base = stored.color;
                            }
                            if self.storemask.contains(SceneStoreMask::REP) {
                                atom.repr.visible_reps = stored.visible_reps;
                            }
                        }
                    }
                    
                    // Mark the molecule as dirty so it rebuilds representations
                    mol_obj.invalidate(DirtyFlags::COLOR | DirtyFlags::REPS);
                }
            }
            
            // Then apply object-level state
            if let Some(obj) = registry.get_mut(name) {
                let state = obj.state_mut();

                if self.storemask.contains(SceneStoreMask::ACTIVE) {
                    state.enabled = obj_data.enabled;
                }

                if self.storemask.contains(SceneStoreMask::COLOR) {
                    state.color = obj_data.color;
                }

                if self.storemask.contains(SceneStoreMask::REP) {
                    state.visible_reps = obj_data.visible_reps;
                }

                if self.storemask.contains(SceneStoreMask::FRAME) {
                    obj.set_current_state(obj_data.current_state);
                }
            }
        }
    }
}

/// Scene manager for storing and recalling named scenes
#[derive(Debug, Default)]
pub struct SceneManager {
    /// Stored scenes by name
    scenes: AHashMap<String, Scene>,
    /// Ordered list of scene names
    scene_order: Vec<String>,
    /// Currently active scene (if any)
    current: Option<String>,
}

impl SceneManager {
    /// Create a new empty scene manager
    pub fn new() -> Self {
        Self::default()
    }

    /// Get the number of stored scenes
    pub fn len(&self) -> usize {
        self.scenes.len()
    }

    /// Check if there are no stored scenes
    pub fn is_empty(&self) -> bool {
        self.scenes.is_empty()
    }

    /// Store a new scene or update an existing one
    pub fn store(
        &mut self,
        key: &str,
        storemask: SceneStoreMask,
        camera: &Camera,
        registry: &ObjectRegistry,
    ) {
        let scene = Scene::capture(key, storemask, camera, registry);

        // Update order if new scene
        if !self.scenes.contains_key(key) {
            self.scene_order.push(key.to_string());
        }

        self.scenes.insert(key.to_string(), scene);
    }

    /// Recall a scene by name
    pub fn recall(
        &mut self,
        key: &str,
        camera: &mut Camera,
        registry: &mut ObjectRegistry,
        animate: bool,
        duration: f32,
    ) -> SceneResult<()> {
        let scene = self
            .scenes
            .get(key)
            .ok_or_else(|| SceneError::SceneNotFound(key.to_string()))?;

        scene.apply(camera, registry, animate, duration);
        self.current = Some(key.to_string());

        Ok(())
    }

    /// Get a scene by name
    pub fn get(&self, key: &str) -> Option<&Scene> {
        self.scenes.get(key)
    }

    /// Get a mutable scene by name
    pub fn get_mut(&mut self, key: &str) -> Option<&mut Scene> {
        self.scenes.get_mut(key)
    }

    /// Delete a scene by name
    pub fn delete(&mut self, key: &str) -> Option<Scene> {
        self.scene_order.retain(|k| k != key);
        if self.current.as_deref() == Some(key) {
            self.current = None;
        }
        self.scenes.remove(key)
    }

    /// Get the list of scene names in order
    pub fn list(&self) -> &[String] {
        &self.scene_order
    }

    /// Get an iterator over all scenes
    pub fn iter(&self) -> impl Iterator<Item = (&str, &Scene)> {
        self.scene_order
            .iter()
            .filter_map(|k| self.scenes.get(k).map(|s| (k.as_str(), s)))
    }

    /// Get the currently active scene name
    pub fn current(&self) -> Option<&str> {
        self.current.as_deref()
    }

    /// Get the next scene in the list (wrapping around)
    pub fn next(&self) -> Option<&str> {
        if self.scene_order.is_empty() {
            return None;
        }

        let current_idx = self
            .current
            .as_ref()
            .and_then(|c| self.scene_order.iter().position(|k| k == c))
            .unwrap_or(self.scene_order.len() - 1);

        let next_idx = (current_idx + 1) % self.scene_order.len();
        Some(&self.scene_order[next_idx])
    }

    /// Get the previous scene in the list (wrapping around)
    pub fn prev(&self) -> Option<&str> {
        if self.scene_order.is_empty() {
            return None;
        }

        let current_idx = self
            .current
            .as_ref()
            .and_then(|c| self.scene_order.iter().position(|k| k == c))
            .unwrap_or(0);

        let prev_idx = if current_idx == 0 {
            self.scene_order.len() - 1
        } else {
            current_idx - 1
        };

        Some(&self.scene_order[prev_idx])
    }

    /// Recall the next scene
    pub fn recall_next(
        &mut self,
        camera: &mut Camera,
        registry: &mut ObjectRegistry,
        animate: bool,
        duration: f32,
    ) -> SceneResult<()> {
        let key = self
            .next()
            .ok_or_else(|| SceneError::SceneNotFound("no scenes".to_string()))?
            .to_string();
        self.recall(&key, camera, registry, animate, duration)
    }

    /// Recall the previous scene
    pub fn recall_prev(
        &mut self,
        camera: &mut Camera,
        registry: &mut ObjectRegistry,
        animate: bool,
        duration: f32,
    ) -> SceneResult<()> {
        let key = self
            .prev()
            .ok_or_else(|| SceneError::SceneNotFound("no scenes".to_string()))?
            .to_string();
        self.recall(&key, camera, registry, animate, duration)
    }

    /// Clear all scenes
    pub fn clear(&mut self) {
        self.scenes.clear();
        self.scene_order.clear();
        self.current = None;
    }

    /// Rename a scene
    pub fn rename(&mut self, old_key: &str, new_key: &str) -> SceneResult<()> {
        if !self.scenes.contains_key(old_key) {
            return Err(SceneError::SceneNotFound(old_key.to_string()));
        }

        if let Some(mut scene) = self.scenes.remove(old_key) {
            scene.name = new_key.to_string();
            self.scenes.insert(new_key.to_string(), scene);

            // Update order
            for k in &mut self.scene_order {
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

    /// Reorder scenes
    pub fn set_order(&mut self, order: &[&str]) {
        let mut new_order = Vec::new();
        let mut used = std::collections::HashSet::new();

        // Add scenes in the specified order
        for key in order {
            if self.scenes.contains_key(*key) && !used.contains(*key) {
                new_order.push((*key).to_string());
                used.insert(*key);
            }
        }

        // Add any remaining scenes at the end
        for key in &self.scene_order {
            if !used.contains(key.as_str()) {
                new_order.push(key.clone());
            }
        }

        self.scene_order = new_order;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::camera::Camera;
    use crate::object::ObjectRegistry;

    #[test]
    fn test_scene_store_recall() {
        let mut manager = SceneManager::new();
        let camera = Camera::new();
        let registry = ObjectRegistry::new();

        manager.store("scene1", SceneStoreMask::ALL, &camera, &registry);
        assert_eq!(manager.len(), 1);
        assert!(manager.get("scene1").is_some());
    }

    #[test]
    fn test_scene_navigation() {
        let mut manager = SceneManager::new();
        let camera = Camera::new();
        let registry = ObjectRegistry::new();

        manager.store("s1", SceneStoreMask::ALL, &camera, &registry);
        manager.store("s2", SceneStoreMask::ALL, &camera, &registry);
        manager.store("s3", SceneStoreMask::ALL, &camera, &registry);

        assert_eq!(manager.next(), Some("s1"));

        manager.current = Some("s1".to_string());
        assert_eq!(manager.next(), Some("s2"));
        assert_eq!(manager.prev(), Some("s3")); // wraps around
    }

    #[test]
    fn test_scene_delete() {
        let mut manager = SceneManager::new();
        let camera = Camera::new();
        let registry = ObjectRegistry::new();

        manager.store("s1", SceneStoreMask::ALL, &camera, &registry);
        manager.store("s2", SceneStoreMask::ALL, &camera, &registry);

        manager.delete("s1");
        assert_eq!(manager.len(), 1);
        assert!(manager.get("s1").is_none());
        assert!(manager.get("s2").is_some());
    }

    #[test]
    fn test_scene_rename() {
        let mut manager = SceneManager::new();
        let camera = Camera::new();
        let registry = ObjectRegistry::new();

        manager.store("old", SceneStoreMask::ALL, &camera, &registry);
        manager.rename("old", "new").unwrap();

        assert!(manager.get("old").is_none());
        assert!(manager.get("new").is_some());
    }

    #[test]
    fn test_storemask() {
        let mask = SceneStoreMask::VIEW | SceneStoreMask::ACTIVE;
        assert!(mask.contains(SceneStoreMask::VIEW));
        assert!(mask.contains(SceneStoreMask::ACTIVE));
        assert!(!mask.contains(SceneStoreMask::COLOR));
    }
}
