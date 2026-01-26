//! High-level picking API for atom/object selection
//!
//! Provides methods for picking atoms and objects via mouse click,
//! using either GPU-based color coding or CPU-based ray casting.

use lin_alg::f32::Vec3;
use pymol_mol::AtomIndex;
use pymol_render::picking::{PickResult, PickingConfig, PickingId, Ray};

use crate::object::{ObjectRegistry, ObjectType};

/// A pick hit with object and atom information
#[derive(Debug, Clone)]
pub struct PickHit {
    /// Name of the object that was hit
    pub object_name: String,
    /// Type of the object
    pub object_type: ObjectType,
    /// Atom index if an atom was hit
    pub atom_index: Option<AtomIndex>,
    /// World-space position of the hit
    pub position: Vec3,
    /// Distance from camera
    pub distance: f32,
}

impl PickHit {
    /// Check if this hit an atom
    pub fn is_atom(&self) -> bool {
        self.atom_index.is_some()
    }
}

/// Picker for interactive selection
///
/// Provides methods for picking objects and atoms in the scene.
pub struct Picker {
    /// Configuration
    config: PickingConfig,
    /// Last pick result (for multi-pick operations)
    last_result: Option<PickHit>,
}

impl Default for Picker {
    fn default() -> Self {
        Self::new()
    }
}

impl Picker {
    /// Create a new picker with default configuration
    pub fn new() -> Self {
        Self {
            config: PickingConfig::default(),
            last_result: None,
        }
    }

    /// Create a picker with custom configuration
    pub fn with_config(config: PickingConfig) -> Self {
        Self {
            config,
            last_result: None,
        }
    }

    /// Get the current configuration
    pub fn config(&self) -> &PickingConfig {
        &self.config
    }

    /// Set the configuration
    pub fn set_config(&mut self, config: PickingConfig) {
        self.config = config;
    }

    /// Get the last pick result
    pub fn last_result(&self) -> Option<&PickHit> {
        self.last_result.as_ref()
    }

    /// Clear the last pick result
    pub fn clear_last(&mut self) {
        self.last_result = None;
    }

    /// Pick using CPU ray-casting
    ///
    /// This is a fallback method that doesn't require GPU readback.
    /// It's slower but works without specialized shaders.
    ///
    /// # Arguments
    /// * `ray` - The picking ray in world space
    /// * `registry` - The object registry to pick from
    ///
    /// # Returns
    /// The closest hit, if any
    pub fn pick_ray(&mut self, ray: &Ray, registry: &ObjectRegistry) -> Option<PickHit> {
        let mut closest: Option<PickHit> = None;
        let mut closest_distance = f32::MAX;

        // Iterate through all enabled objects
        for obj in registry.enabled_objects() {
            let name = obj.name().to_string();
            let obj_type = obj.object_type();

            // For molecules, test against atom spheres
            if let Some(mol_obj) = registry.get_molecule(&name) {
                let molecule = mol_obj.molecule();
                let state_idx = mol_obj.display_state();

                if let Some(coord_set) = molecule.get_coord_set(state_idx) {
                    for (atom_idx, coord) in coord_set.iter_with_atoms() {
                        if let Some(atom) = molecule.get_atom(atom_idx) {
                            let center = [coord.x, coord.y, coord.z];
                            let radius = atom.effective_vdw();

                            if let Some(t) = ray.intersect_sphere(center, radius) {
                                if t < closest_distance {
                                    let hit_pos = ray.at(t);
                                    closest = Some(PickHit {
                                        object_name: name.clone(),
                                        object_type: obj_type,
                                        atom_index: Some(atom_idx),
                                        position: Vec3::new(hit_pos[0], hit_pos[1], hit_pos[2]),
                                        distance: t,
                                    });
                                    closest_distance = t;
                                }
                            }
                        }
                    }
                }
            }
            // For other object types, test against bounding box
            else if let Some((min, max)) = obj.extent() {
                // Simple AABB test (could be improved with actual geometry)
                if let Some(t) = Self::ray_aabb_intersection(ray, min, max) {
                    if t < closest_distance {
                        let hit_pos = ray.at(t);
                        closest = Some(PickHit {
                            object_name: name,
                            object_type: obj_type,
                            atom_index: None,
                            position: Vec3::new(hit_pos[0], hit_pos[1], hit_pos[2]),
                            distance: t,
                        });
                        closest_distance = t;
                    }
                }
            }
        }

        self.last_result = closest.clone();
        closest
    }

    /// Convert GPU picking result to a PickHit
    ///
    /// # Arguments
    /// * `result` - The GPU picking result
    /// * `registry` - The object registry
    /// * `object_names` - Mapping from object IDs to names
    pub fn from_gpu_result(
        &mut self,
        result: &PickResult,
        registry: &ObjectRegistry,
        object_names: &[String],
    ) -> Option<PickHit> {
        if !result.is_hit() {
            self.last_result = None;
            return None;
        }

        let object_id = result.object_id as usize;
        if object_id >= object_names.len() {
            self.last_result = None;
            return None;
        }

        let object_name = &object_names[object_id];
        let obj = registry.get(object_name)?;

        let hit = PickHit {
            object_name: object_name.clone(),
            object_type: obj.object_type(),
            atom_index: if result.is_atom_hit() {
                Some(AtomIndex::from(result.atom_index as usize))
            } else {
                None
            },
            position: result
                .position
                .map(|p| Vec3::new(p[0], p[1], p[2]))
                .unwrap_or_else(|| Vec3::new(0.0, 0.0, 0.0)),
            distance: result.depth,
        };

        self.last_result = Some(hit.clone());
        Some(hit)
    }

    /// Ray-AABB intersection test
    fn ray_aabb_intersection(ray: &Ray, min: Vec3, max: Vec3) -> Option<f32> {
        let mut tmin = f32::NEG_INFINITY;
        let mut tmax = f32::INFINITY;

        let min = [min.x, min.y, min.z];
        let max = [max.x, max.y, max.z];

        for i in 0..3 {
            if ray.direction[i].abs() > 1e-6 {
                let t1 = (min[i] - ray.origin[i]) / ray.direction[i];
                let t2 = (max[i] - ray.origin[i]) / ray.direction[i];

                let t_near = t1.min(t2);
                let t_far = t1.max(t2);

                tmin = tmin.max(t_near);
                tmax = tmax.min(t_far);

                if tmin > tmax {
                    return None;
                }
            } else {
                // Ray is parallel to slab
                if ray.origin[i] < min[i] || ray.origin[i] > max[i] {
                    return None;
                }
            }
        }

        if tmax < 0.0 {
            None
        } else if tmin < 0.0 {
            Some(tmax)
        } else {
            Some(tmin)
        }
    }
}

/// Build a mapping from object IDs to names for GPU picking
#[allow(dead_code)]
pub fn build_object_id_map(registry: &ObjectRegistry) -> Vec<String> {
    registry.names().map(|s| s.to_string()).collect()
}

/// Get the object ID for a named object
#[allow(dead_code)]
pub fn get_object_id(object_names: &[String], name: &str) -> Option<u32> {
    object_names
        .iter()
        .position(|n| n == name)
        .map(|i| i as u32)
}

/// Create a picking ID for an object
#[allow(dead_code)]
pub fn object_picking_id(object_names: &[String], name: &str) -> Option<PickingId> {
    get_object_id(object_names, name).map(PickingId::object)
}

/// Create a picking ID for an atom
#[allow(dead_code)]
pub fn atom_picking_id(
    object_names: &[String],
    object_name: &str,
    atom_index: AtomIndex,
) -> Option<PickingId> {
    get_object_id(object_names, object_name)
        .map(|obj_id| PickingId::new(obj_id, usize::from(atom_index) as u32))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_picker_creation() {
        let picker = Picker::new();
        assert!(picker.last_result().is_none());
    }

    #[test]
    fn test_build_object_id_map() {
        // This would need a mock registry to test properly
        let registry = ObjectRegistry::new();
        let map = build_object_id_map(&registry);
        assert!(map.is_empty());
    }

    #[test]
    fn test_ray_aabb_intersection() {
        let ray = Ray::new([0.0, 0.0, 0.0], [0.0, 0.0, -1.0]);
        let min = Vec3::new(-1.0, -1.0, -5.0);
        let max = Vec3::new(1.0, 1.0, -3.0);

        let t = Picker::ray_aabb_intersection(&ray, min, max);
        assert!(t.is_some());
        assert!((t.unwrap() - 3.0).abs() < 0.001);
    }

    #[test]
    fn test_ray_aabb_miss() {
        let ray = Ray::new([0.0, 0.0, 0.0], [0.0, 0.0, -1.0]);
        let min = Vec3::new(10.0, 10.0, -5.0);
        let max = Vec3::new(11.0, 11.0, -3.0);

        let t = Picker::ray_aabb_intersection(&ray, min, max);
        assert!(t.is_none());
    }

    #[test]
    fn test_pick_hit() {
        let hit = PickHit {
            object_name: "protein".to_string(),
            object_type: ObjectType::Molecule,
            atom_index: Some(AtomIndex::from(42usize)),
            position: Vec3::new(1.0, 2.0, 3.0),
            distance: 5.0,
        };

        assert!(hit.is_atom());
        assert_eq!(hit.object_name, "protein");
    }
}
