//! Object management for scene objects
//!
//! This module provides:
//! - [`Object`] trait for all scene objects
//! - [`ObjectType`] enum for object type identification
//! - [`ObjectState`] for per-object visual state
//! - [`ObjectRegistry`] for managing named objects

mod cgo;
mod group;
mod label;
mod map;
mod molecule;
mod surface;

pub use cgo::CgoObject;
pub use group::GroupObject;
pub use label::{Label, LabelAnchor, LabelObject};
pub use map::MapObject;
pub use molecule::{DirtyFlags, MoleculeObject};
pub use surface::SurfaceObject;

use ahash::AHashMap;
use lin_alg::f32::{Mat4, Vec3};
use pymol_color::ColorIndex;
use pymol_mol::RepMask;
use pymol_settings::GlobalSettings;

use crate::error::{SceneError, SceneResult};

/// Object type enumeration
///
/// Identifies the type of a scene object for type-safe downcasting and
/// type-specific operations.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ObjectType {
    /// Molecular object (atoms, bonds, coordinates)
    Molecule,
    /// Electron density map
    Map,
    /// Triangle mesh
    Mesh,
    /// Molecular surface
    Surface,
    /// Compiled Graphics Object (custom primitives)
    Cgo,
    /// Group of objects
    Group,
    /// Measurement (distances, angles, etc.)
    Measurement,
    /// Volume rendering
    Volume,
    /// Alignment object
    Alignment,
}

impl std::fmt::Display for ObjectType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ObjectType::Molecule => write!(f, "molecule"),
            ObjectType::Map => write!(f, "map"),
            ObjectType::Mesh => write!(f, "mesh"),
            ObjectType::Surface => write!(f, "surface"),
            ObjectType::Cgo => write!(f, "cgo"),
            ObjectType::Group => write!(f, "group"),
            ObjectType::Measurement => write!(f, "measurement"),
            ObjectType::Volume => write!(f, "volume"),
            ObjectType::Alignment => write!(f, "alignment"),
        }
    }
}

/// Per-object visual state
///
/// Contains the common visual properties shared by all object types.
#[derive(Debug, Clone)]
pub struct ObjectState {
    /// Whether the object is enabled/visible
    pub enabled: bool,
    /// Default object color
    pub color: ColorIndex,
    /// Visible representations (for applicable objects)
    pub visible_reps: RepMask,
    /// 4x4 transformation matrix (TTT in PyMOL)
    ///
    /// This matrix is applied to all coordinates when rendering.
    pub transform: Mat4,
}

impl Default for ObjectState {
    fn default() -> Self {
        Self {
            enabled: true,
            color: ColorIndex::default(),
            visible_reps: RepMask::default(),
            transform: Mat4::new_identity(),
        }
    }
}

impl ObjectState {
    /// Create a new object state with default values
    pub fn new() -> Self {
        Self::default()
    }

    /// Create a disabled object state
    pub fn disabled() -> Self {
        Self {
            enabled: false,
            ..Default::default()
        }
    }

    /// Set the transformation matrix
    pub fn set_transform(&mut self, transform: Mat4) {
        self.transform = transform;
    }

    /// Reset the transformation to identity
    pub fn reset_transform(&mut self) {
        self.transform = Mat4::new_identity();
    }

    /// Check if a representation is visible
    pub fn rep_visible(&self, rep: u32) -> bool {
        self.visible_reps.is_visible(rep)
    }

    /// Show a representation
    pub fn show_rep(&mut self, rep: u32) {
        self.visible_reps.set_visible(rep);
    }

    /// Hide a representation
    pub fn hide_rep(&mut self, rep: u32) {
        self.visible_reps.set_hidden(rep);
    }
}

/// Trait for all scene objects
///
/// This trait defines the common interface for all objects that can be
/// added to a scene, including molecules, maps, surfaces, etc.
pub trait Object: Send + Sync {
    /// Get the object name
    fn name(&self) -> &str;

    /// Get the object type
    fn object_type(&self) -> ObjectType;

    /// Get the visual state
    fn state(&self) -> &ObjectState;

    /// Get mutable access to the visual state
    fn state_mut(&mut self) -> &mut ObjectState;

    /// Get the bounding box (min, max) for the current coordinate state
    ///
    /// Returns `None` if the object has no geometry.
    fn extent(&self) -> Option<(Vec3, Vec3)>;

    /// Get the number of coordinate states
    ///
    /// For single-state objects, returns 1.
    fn n_states(&self) -> usize {
        1
    }

    /// Get the current state index
    fn current_state(&self) -> usize {
        0
    }

    /// Set the current state index
    fn set_current_state(&mut self, _state: usize) -> bool {
        false
    }

    /// Get object-level settings override (if any)
    fn settings(&self) -> Option<&GlobalSettings> {
        None
    }

    /// Get mutable object-level settings
    fn settings_mut(&mut self) -> Option<&mut GlobalSettings> {
        None
    }

    /// Check if this object is enabled
    fn is_enabled(&self) -> bool {
        self.state().enabled
    }

    /// Enable the object
    fn enable(&mut self) {
        self.state_mut().enabled = true;
    }

    /// Disable the object
    fn disable(&mut self) {
        self.state_mut().enabled = false;
    }

    /// Get the center of the object
    fn center(&self) -> Option<Vec3> {
        let (min, max) = self.extent()?;
        Some(Vec3::new(
            (min.x + max.x) * 0.5,
            (min.y + max.y) * 0.5,
            (min.z + max.z) * 0.5,
        ))
    }
}

/// Object registry for managing named objects
///
/// The registry stores all scene objects by name and maintains render order.
pub struct ObjectRegistry {
    /// Objects stored by name
    objects: AHashMap<String, Box<dyn Object>>,
    /// Ordered list of object names for rendering
    render_order: Vec<String>,
    /// Counter for generating unique names
    next_id: u32,
}

impl Default for ObjectRegistry {
    fn default() -> Self {
        Self::new()
    }
}

impl ObjectRegistry {
    /// Create a new empty registry
    pub fn new() -> Self {
        Self {
            objects: AHashMap::new(),
            render_order: Vec::new(),
            next_id: 1,
        }
    }

    /// Get the number of objects
    pub fn len(&self) -> usize {
        self.objects.len()
    }

    /// Check if the registry is empty
    pub fn is_empty(&self) -> bool {
        self.objects.is_empty()
    }

    /// Add an object to the registry
    ///
    /// Returns the name used for the object. If an object with the same name
    /// already exists, it will be replaced.
    pub fn add<O: Object + 'static>(&mut self, obj: O) -> &str {
        let name = obj.name().to_string();
        self.render_order.retain(|n| n != &name);
        self.render_order.push(name.clone());
        self.objects.insert(name, Box::new(obj));
        self.render_order.last().map(|s| s.as_str()).unwrap()
    }

    /// Add an object with an auto-generated name
    pub fn add_with_name<O: Object + 'static>(&mut self, mut obj: O, base_name: &str) -> String
    where
        O: ObjectWithName,
    {
        let name = self.generate_unique_name(base_name);
        obj.set_name(name.clone());
        self.add(obj);
        name
    }

    /// Generate a unique name based on a base name
    fn generate_unique_name(&mut self, base: &str) -> String {
        if !self.objects.contains_key(base) {
            return base.to_string();
        }

        loop {
            let name = format!("{}_{}", base, self.next_id);
            self.next_id += 1;
            if !self.objects.contains_key(&name) {
                return name;
            }
        }
    }

    /// Remove an object by name
    pub fn remove(&mut self, name: &str) -> Option<Box<dyn Object>> {
        self.render_order.retain(|n| n != name);
        self.objects.remove(name)
    }

    /// Get an object by name
    pub fn get(&self, name: &str) -> Option<&dyn Object> {
        self.objects.get(name).map(|b| b.as_ref())
    }

    /// Get mutable access to an object by name
    pub fn get_mut(&mut self, name: &str) -> Option<&mut dyn Object> {
        match self.objects.get_mut(name) {
            Some(boxed) => Some(boxed.as_mut()),
            None => None,
        }
    }

    /// Get a molecule object by name
    pub fn get_molecule(&self, name: &str) -> Option<&MoleculeObject> {
        let obj = self.objects.get(name)?;
        if obj.object_type() == ObjectType::Molecule {
            // Safety: we just checked the type
            let ptr = obj.as_ref() as *const dyn Object;
            Some(unsafe { &*(ptr as *const MoleculeObject) })
        } else {
            None
        }
    }

    /// Get mutable access to a molecule object by name
    pub fn get_molecule_mut(&mut self, name: &str) -> Option<&mut MoleculeObject> {
        let obj = self.objects.get_mut(name)?;
        if obj.object_type() == ObjectType::Molecule {
            // Safety: we just checked the type
            let ptr = obj.as_mut() as *mut dyn Object;
            Some(unsafe { &mut *(ptr as *mut MoleculeObject) })
        } else {
            None
        }
    }

    /// Get a group object by name
    pub fn get_group(&self, name: &str) -> Option<&GroupObject> {
        let obj = self.objects.get(name)?;
        if obj.object_type() == ObjectType::Group {
            // Safety: we just checked the type
            let ptr = obj.as_ref() as *const dyn Object;
            Some(unsafe { &*(ptr as *const GroupObject) })
        } else {
            None
        }
    }

    /// Get mutable access to a group object by name
    pub fn get_group_mut(&mut self, name: &str) -> Option<&mut GroupObject> {
        let obj = self.objects.get_mut(name)?;
        if obj.object_type() == ObjectType::Group {
            // Safety: we just checked the type
            let ptr = obj.as_mut() as *mut dyn Object;
            Some(unsafe { &mut *(ptr as *mut GroupObject) })
        } else {
            None
        }
    }

    /// Get a surface object by name
    pub fn get_surface(&self, name: &str) -> Option<&SurfaceObject> {
        let obj = self.objects.get(name)?;
        if obj.object_type() == ObjectType::Surface {
            // Safety: we just checked the type
            let ptr = obj.as_ref() as *const dyn Object;
            Some(unsafe { &*(ptr as *const SurfaceObject) })
        } else {
            None
        }
    }

    /// Get mutable access to a surface object by name
    pub fn get_surface_mut(&mut self, name: &str) -> Option<&mut SurfaceObject> {
        let obj = self.objects.get_mut(name)?;
        if obj.object_type() == ObjectType::Surface {
            // Safety: we just checked the type
            let ptr = obj.as_mut() as *mut dyn Object;
            Some(unsafe { &mut *(ptr as *mut SurfaceObject) })
        } else {
            None
        }
    }

    /// Get a CGO object by name
    pub fn get_cgo(&self, name: &str) -> Option<&CgoObject> {
        let obj = self.objects.get(name)?;
        if obj.object_type() == ObjectType::Cgo {
            // Safety: we just checked the type
            let ptr = obj.as_ref() as *const dyn Object;
            Some(unsafe { &*(ptr as *const CgoObject) })
        } else {
            None
        }
    }

    /// Get mutable access to a CGO object by name
    pub fn get_cgo_mut(&mut self, name: &str) -> Option<&mut CgoObject> {
        let obj = self.objects.get_mut(name)?;
        if obj.object_type() == ObjectType::Cgo {
            // Safety: we just checked the type
            let ptr = obj.as_mut() as *mut dyn Object;
            Some(unsafe { &mut *(ptr as *mut CgoObject) })
        } else {
            None
        }
    }

    /// Get a map object by name
    pub fn get_map(&self, name: &str) -> Option<&MapObject> {
        let obj = self.objects.get(name)?;
        if obj.object_type() == ObjectType::Map {
            // Safety: we just checked the type
            let ptr = obj.as_ref() as *const dyn Object;
            Some(unsafe { &*(ptr as *const MapObject) })
        } else {
            None
        }
    }

    /// Get mutable access to a map object by name
    pub fn get_map_mut(&mut self, name: &str) -> Option<&mut MapObject> {
        let obj = self.objects.get_mut(name)?;
        if obj.object_type() == ObjectType::Map {
            // Safety: we just checked the type
            let ptr = obj.as_mut() as *mut dyn Object;
            Some(unsafe { &mut *(ptr as *mut MapObject) })
        } else {
            None
        }
    }

    /// Get the combined bounding box of all children in a group
    ///
    /// Recursively computes the extent including nested groups.
    pub fn group_extent(&self, group_name: &str) -> Option<(Vec3, Vec3)> {
        let group = self.get_group(group_name)?;
        let children = group.children().to_vec();

        let mut min = Vec3::new(f32::MAX, f32::MAX, f32::MAX);
        let mut max = Vec3::new(f32::MIN, f32::MIN, f32::MIN);
        let mut has_extent = false;

        for child_name in &children {
            // Check if child is a group (recursive)
            if let Some(_) = self.get_group(child_name) {
                if let Some((child_min, child_max)) = self.group_extent(child_name) {
                    min.x = min.x.min(child_min.x);
                    min.y = min.y.min(child_min.y);
                    min.z = min.z.min(child_min.z);
                    max.x = max.x.max(child_max.x);
                    max.y = max.y.max(child_max.y);
                    max.z = max.z.max(child_max.z);
                    has_extent = true;
                }
            } else if let Some(obj) = self.get(child_name) {
                if let Some((child_min, child_max)) = obj.extent() {
                    min.x = min.x.min(child_min.x);
                    min.y = min.y.min(child_min.y);
                    min.z = min.z.min(child_min.z);
                    max.x = max.x.max(child_max.x);
                    max.y = max.y.max(child_max.y);
                    max.z = max.z.max(child_max.z);
                    has_extent = true;
                }
            }
        }

        if has_extent {
            Some((min, max))
        } else {
            None
        }
    }

    /// Enable or disable a group and all its children recursively
    pub fn set_group_enabled(&mut self, group_name: &str, enabled: bool) -> SceneResult<()> {
        // First, collect the children names
        let children = {
            let group = self
                .get_group(group_name)
                .ok_or_else(|| SceneError::ObjectNotFound(group_name.to_string()))?;
            group.children().to_vec()
        };

        // Enable/disable the group itself
        if let Some(obj) = self.objects.get_mut(group_name) {
            obj.state_mut().enabled = enabled;
        }

        // Recursively enable/disable children
        for child_name in children {
            // Check if child is a group
            if self.get_group(&child_name).is_some() {
                // Recursive call for nested groups
                let _ = self.set_group_enabled(&child_name, enabled);
            } else if let Some(obj) = self.objects.get_mut(&child_name) {
                obj.state_mut().enabled = enabled;
            }
        }

        Ok(())
    }

    /// Check if an object exists
    pub fn contains(&self, name: &str) -> bool {
        self.objects.contains_key(name)
    }

    /// Get all object names
    pub fn names(&self) -> impl Iterator<Item = &str> {
        self.render_order.iter().map(|s| s.as_str())
    }

    /// Get all objects in render order
    pub fn iter(&self) -> impl Iterator<Item = &dyn Object> {
        self.render_order
            .iter()
            .filter_map(|name| self.objects.get(name).map(|b| b.as_ref()))
    }

    /// Get all objects mutably as boxed references
    ///
    /// Note: Returns iterator over Box references, not trait object references,
    /// to avoid lifetime issues with trait objects.
    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut Box<dyn Object>> + '_ {
        self.objects.values_mut()
    }

    /// Enable an object
    pub fn enable(&mut self, name: &str, enabled: bool) -> SceneResult<()> {
        let obj = self
            .objects
            .get_mut(name)
            .ok_or_else(|| SceneError::ObjectNotFound(name.to_string()))?;
        obj.state_mut().enabled = enabled;
        Ok(())
    }

    /// Disable all objects
    pub fn disable_all(&mut self) {
        for obj in self.objects.values_mut() {
            obj.state_mut().enabled = false;
        }
    }

    /// Enable all objects
    pub fn enable_all(&mut self) {
        for obj in self.objects.values_mut() {
            obj.state_mut().enabled = true;
        }
    }

    /// Get enabled objects in render order
    pub fn enabled_objects(&self) -> impl Iterator<Item = &dyn Object> {
        self.render_order.iter().filter_map(|name| {
            self.objects.get(name).and_then(|obj| {
                if obj.state().enabled {
                    Some(obj.as_ref())
                } else {
                    None
                }
            })
        })
    }

    /// Get object names matching a pattern
    ///
    /// Currently supports simple wildcard matching with `*`.
    pub fn matching<'a>(&'a self, pattern: &'a str) -> Vec<&'a str> {
        if pattern == "*" || pattern == "all" {
            return self.render_order.iter().map(|s| s.as_str()).collect();
        }

        // Simple prefix/suffix matching
        if let Some(prefix) = pattern.strip_suffix('*') {
            return self
                .render_order
                .iter()
                .filter(|name| name.starts_with(prefix))
                .map(|s| s.as_str())
                .collect();
        }

        if let Some(suffix) = pattern.strip_prefix('*') {
            return self
                .render_order
                .iter()
                .filter(|name| name.ends_with(suffix))
                .map(|s| s.as_str())
                .collect();
        }

        // Exact match
        if self.objects.contains_key(pattern) {
            // Return a reference to our stored name, not the input
            self.render_order
                .iter()
                .find(|n| n.as_str() == pattern)
                .map(|s| vec![s.as_str()])
                .unwrap_or_default()
        } else {
            vec![]
        }
    }

    /// Set the render order
    pub fn set_render_order(&mut self, names: &[&str]) {
        // Keep only names that exist and preserve any not in the list at the end
        let mut new_order = Vec::new();
        let mut used = std::collections::HashSet::new();

        for name in names {
            if self.objects.contains_key(*name) && !used.contains(*name) {
                new_order.push((*name).to_string());
                used.insert(*name);
            }
        }

        // Add any remaining objects at the end
        for name in &self.render_order {
            if !used.contains(name.as_str()) {
                new_order.push(name.clone());
            }
        }

        self.render_order = new_order;
    }

    /// Move an object to the front of the render order
    pub fn bring_to_front(&mut self, name: &str) {
        if let Some(pos) = self.render_order.iter().position(|n| n == name) {
            let name = self.render_order.remove(pos);
            self.render_order.push(name);
        }
    }

    /// Move an object to the back of the render order
    pub fn send_to_back(&mut self, name: &str) {
        if let Some(pos) = self.render_order.iter().position(|n| n == name) {
            let name = self.render_order.remove(pos);
            self.render_order.insert(0, name);
        }
    }

    /// Compute the bounding box of all enabled objects
    pub fn extent(&self) -> Option<(Vec3, Vec3)> {
        let mut min = Vec3::new(f32::MAX, f32::MAX, f32::MAX);
        let mut max = Vec3::new(f32::MIN, f32::MIN, f32::MIN);
        let mut has_extent = false;

        for obj in self.enabled_objects() {
            if let Some((obj_min, obj_max)) = obj.extent() {
                min.x = min.x.min(obj_min.x);
                min.y = min.y.min(obj_min.y);
                min.z = min.z.min(obj_min.z);
                max.x = max.x.max(obj_max.x);
                max.y = max.y.max(obj_max.y);
                max.z = max.z.max(obj_max.z);
                has_extent = true;
            }
        }

        if has_extent {
            Some((min, max))
        } else {
            None
        }
    }

    /// Get the center of all enabled objects
    pub fn center(&self) -> Option<Vec3> {
        let (min, max) = self.extent()?;
        Some(Vec3::new(
            (min.x + max.x) * 0.5,
            (min.y + max.y) * 0.5,
            (min.z + max.z) * 0.5,
        ))
    }

    /// Rename an object
    pub fn rename(&mut self, old_name: &str, new_name: &str) -> SceneResult<()> {
        if !self.objects.contains_key(old_name) {
            return Err(SceneError::ObjectNotFound(old_name.to_string()));
        }
        if self.objects.contains_key(new_name) {
            return Err(SceneError::ObjectExists(new_name.to_string()));
        }

        if let Some(obj) = self.objects.remove(old_name) {
            self.objects.insert(new_name.to_string(), obj);

            // Update render order
            for name in &mut self.render_order {
                if name == old_name {
                    *name = new_name.to_string();
                }
            }
        }

        Ok(())
    }

    /// Clear all objects
    pub fn clear(&mut self) {
        self.objects.clear();
        self.render_order.clear();
    }
}

/// Trait for objects that can have their name set
pub trait ObjectWithName {
    /// Set the object name
    fn set_name(&mut self, name: String);
}

#[cfg(test)]
mod tests {
    use super::*;

    // Mock object for testing
    struct MockObject {
        name: String,
        state: ObjectState,
    }

    impl MockObject {
        fn new(name: &str) -> Self {
            Self {
                name: name.to_string(),
                state: ObjectState::default(),
            }
        }
    }

    impl Object for MockObject {
        fn name(&self) -> &str {
            &self.name
        }

        fn object_type(&self) -> ObjectType {
            ObjectType::Molecule
        }

        fn state(&self) -> &ObjectState {
            &self.state
        }

        fn state_mut(&mut self) -> &mut ObjectState {
            &mut self.state
        }

        fn extent(&self) -> Option<(Vec3, Vec3)> {
            None
        }
    }

    #[test]
    fn test_registry_add_remove() {
        let mut registry = ObjectRegistry::new();

        registry.add(MockObject::new("obj1"));
        registry.add(MockObject::new("obj2"));

        assert_eq!(registry.len(), 2);
        assert!(registry.contains("obj1"));
        assert!(registry.contains("obj2"));

        registry.remove("obj1");
        assert_eq!(registry.len(), 1);
        assert!(!registry.contains("obj1"));
    }

    #[test]
    fn test_registry_render_order() {
        let mut registry = ObjectRegistry::new();

        registry.add(MockObject::new("a"));
        registry.add(MockObject::new("b"));
        registry.add(MockObject::new("c"));

        let names: Vec<_> = registry.names().collect();
        assert_eq!(names, vec!["a", "b", "c"]);

        registry.bring_to_front("a");
        let names: Vec<_> = registry.names().collect();
        assert_eq!(names, vec!["b", "c", "a"]);

        registry.send_to_back("a");
        let names: Vec<_> = registry.names().collect();
        assert_eq!(names, vec!["a", "b", "c"]);
    }

    #[test]
    fn test_registry_matching() {
        let mut registry = ObjectRegistry::new();

        registry.add(MockObject::new("protein"));
        registry.add(MockObject::new("protein_chain_A"));
        registry.add(MockObject::new("ligand"));

        let matches = registry.matching("protein*");
        assert_eq!(matches.len(), 2);

        let matches = registry.matching("*chain_A");
        assert_eq!(matches.len(), 1);

        let matches = registry.matching("*");
        assert_eq!(matches.len(), 3);
    }

    #[test]
    fn test_enable_disable() {
        let mut registry = ObjectRegistry::new();
        registry.add(MockObject::new("obj1"));
        registry.add(MockObject::new("obj2"));

        registry.disable_all();
        assert!(registry.enabled_objects().next().is_none());

        registry.enable("obj1", true).unwrap();
        assert_eq!(registry.enabled_objects().count(), 1);
    }
}
