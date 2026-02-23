//! High-level picking API for atom/object selection
//!
//! Provides methods for picking atoms and objects via mouse click,
//! using either GPU-based color coding or CPU-based ray casting.

use lin_alg::f32::Vec3;
use pymol_mol::{AtomIndex, ObjectMolecule};
use pymol_render::picking::{PickResult, PickingConfig, Ray};
use pymol_select::SelectionResult;

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

/// Expand a pick hit to a selection based on `mouse_selection_mode`.
///
/// Modes: 0=atoms, 1=residues, 2=chains, 3=segments, 4=objects,
///        5=molecules, 6=C-alphas.
pub fn expand_pick_to_selection(
    hit: &PickHit,
    mode: i32,
    molecule: &ObjectMolecule,
) -> SelectionResult {
    let atom_count = molecule.atom_count();

    // Objects / molecules: select everything
    if mode == 4 || mode == 5 {
        return SelectionResult::all(atom_count);
    }

    // Need a hit atom for the remaining modes
    let hit_idx = match hit.atom_index {
        Some(idx) => idx,
        None => return SelectionResult::all(atom_count),
    };

    // Atoms: single atom
    if mode == 0 {
        return SelectionResult::from_indices(atom_count, std::iter::once(hit_idx));
    }

    // Get the reference atom for grouping
    let ref_atom = match molecule.get_atom(hit_idx) {
        Some(a) => a,
        None => return SelectionResult::none(atom_count),
    };

    // C-alphas: select the CA atom of the hit residue
    if mode == 6 {
        let indices = molecule.atoms_indexed().filter_map(|(idx, atom)| {
            if &*atom.name == "CA"
                && atom.residue.key.chain == ref_atom.residue.key.chain
                && atom.residue.key.resv == ref_atom.residue.key.resv
                && atom.residue.key.inscode == ref_atom.residue.key.inscode
            {
                Some(idx)
            } else {
                None
            }
        });
        return SelectionResult::from_indices(atom_count, indices);
    }

    let indices = molecule.atoms_indexed().filter_map(|(idx, atom)| {
        let matches = match mode {
            1 => {
                // Residues: same chain + resv + inscode
                atom.residue.key.chain == ref_atom.residue.key.chain
                    && atom.residue.key.resv == ref_atom.residue.key.resv
                    && atom.residue.key.inscode == ref_atom.residue.key.inscode
            }
            2 => {
                // Chains: same chain
                atom.residue.key.chain == ref_atom.residue.key.chain
            }
            3 => {
                // Segments: same segi
                atom.residue.segi == ref_atom.residue.segi
            }
            _ => false,
        };
        if matches { Some(idx) } else { None }
    });

    SelectionResult::from_indices(atom_count, indices)
}

/// Generate a PyMOL selection expression string for a pick hit.
///
/// The expression depends on `mouse_selection_mode`:
/// - 0 (atoms): `model OBJ and index N`
/// - 1 (residues): `model OBJ and chain C and resi R[inscode]`
/// - 2 (chains): `model OBJ and chain C`
/// - 3 (segments): `model OBJ and segi S`
/// - 4 (objects): `model OBJ`
/// - 5 (molecules): `bymolecule (model OBJ and index N)` (connected component)
/// - 6 (C-alphas): `model OBJ and chain C and resi R[inscode] and name CA`
pub fn pick_expression_for_hit(
    hit: &PickHit,
    mode: i32,
    molecule: &ObjectMolecule,
) -> Option<String> {
    let obj = &hit.object_name;

    // Object mode: just the model name
    if mode == 4 {
        return Some(format!("model {obj}"));
    }

    let atom_idx = hit.atom_index?;
    let atom = molecule.get_atom(atom_idx)?;
    let idx = atom_idx.as_usize();

    match mode {
        0 => Some(format!("model {obj} and index {idx}")),
        5 => Some(format!("bymolecule (model {obj} and index {idx})")),
        1 => {
            let resi = format_resi(atom.residue.key.resv, atom.residue.key.inscode);
            Some(format!("model {obj} and chain {} and resi {resi}", atom.residue.key.chain))
        }
        2 => Some(format!("model {obj} and chain {}", atom.residue.key.chain)),
        3 => Some(format!("model {obj} and segi {}", atom.residue.segi)),
        6 => {
            let resi = format_resi(atom.residue.key.resv, atom.residue.key.inscode);
            Some(format!(
                "model {obj} and chain {} and resi {resi} and name CA",
                atom.residue.key.chain
            ))
        }
        _ => None,
    }
}

/// Format a residue number with optional insertion code.
fn format_resi(resv: i32, inscode: char) -> String {
    if inscode != ' ' && inscode != '\0' {
        format!("{resv}{inscode}")
    } else {
        resv.to_string()
    }
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

    /// Build a small molecule with 2 residues across 2 chains for testing.
    ///
    /// Layout: chain A / ALA 1 (atoms 0,1), chain A / GLY 2 (atom 2),
    ///         chain B / ALA 1 (atom 3) with segi "S2"
    fn test_molecule() -> ObjectMolecule {
        use pymol_mol::{AtomBuilder, Element, MoleculeBuilder};

        MoleculeBuilder::new("test")
            .add_atom(
                AtomBuilder::new().name("N").element(Element::Nitrogen)
                    .chain("A").resn("ALA").resv(1).build(),
                Vec3::new(0.0, 0.0, 0.0),
            )
            .add_atom(
                AtomBuilder::new().name("CA").element(Element::Carbon)
                    .chain("A").resn("ALA").resv(1).build(),
                Vec3::new(1.0, 0.0, 0.0),
            )
            .add_atom(
                AtomBuilder::new().name("N").element(Element::Nitrogen)
                    .chain("A").resn("GLY").resv(2).build(),
                Vec3::new(2.0, 0.0, 0.0),
            )
            .add_atom(
                AtomBuilder::new().name("CA").element(Element::Carbon)
                    .chain("B").resn("ALA").resv(1).segi("S2").build(),
                Vec3::new(3.0, 0.0, 0.0),
            )
            .build()
    }

    fn make_hit(atom_index: usize) -> PickHit {
        PickHit {
            object_name: "test".to_string(),
            object_type: ObjectType::Molecule,
            atom_index: Some(AtomIndex::from(atom_index)),
            position: Vec3::new(0.0, 0.0, 0.0),
            distance: 1.0,
        }
    }

    #[test]
    fn test_expand_atoms_mode() {
        let mol = test_molecule();
        let hit = make_hit(0);
        let sel = expand_pick_to_selection(&hit, 0, &mol);
        assert_eq!(sel.count(), 1);
        assert!(sel.contains_index(0));
    }

    #[test]
    fn test_expand_residues_mode() {
        let mol = test_molecule();
        let hit = make_hit(0); // atom 0 is in ALA 1 chain A
        let sel = expand_pick_to_selection(&hit, 1, &mol);
        // Should select atoms 0 and 1 (both ALA 1 chain A)
        assert_eq!(sel.count(), 2);
        assert!(sel.contains_index(0));
        assert!(sel.contains_index(1));
        assert!(!sel.contains_index(2));
    }

    #[test]
    fn test_expand_chains_mode() {
        let mol = test_molecule();
        let hit = make_hit(0); // chain A
        let sel = expand_pick_to_selection(&hit, 2, &mol);
        // Should select atoms 0, 1, 2 (all chain A)
        assert_eq!(sel.count(), 3);
        assert!(!sel.contains_index(3));
    }

    #[test]
    fn test_expand_segments_mode() {
        let mol = test_molecule();
        // Hit atom 3 which has segi "S2"
        let hit = make_hit(3);
        let sel = expand_pick_to_selection(&hit, 3, &mol);
        // Only atom 3 has segi "S2"
        assert_eq!(sel.count(), 1);
        assert!(sel.contains_index(3));
    }

    #[test]
    fn test_expand_objects_mode() {
        let mol = test_molecule();
        let hit = make_hit(0);
        let sel = expand_pick_to_selection(&hit, 4, &mol);
        assert_eq!(sel.count(), 4); // all atoms
    }

    #[test]
    fn test_expand_molecules_mode() {
        let mol = test_molecule();
        let hit = make_hit(0);
        let sel = expand_pick_to_selection(&hit, 5, &mol);
        assert_eq!(sel.count(), 4); // all atoms
    }

    #[test]
    fn test_expand_c_alphas_mode() {
        let mol = test_molecule();
        // Hit atom 0 (N in ALA 1, chain A) → select only CA of that residue (atom 1)
        let hit = make_hit(0);
        let sel = expand_pick_to_selection(&hit, 6, &mol);
        assert_eq!(sel.count(), 1);
        assert!(sel.contains_index(1));
        assert!(!sel.contains_index(3)); // CA in chain B — different residue

        // Hit atom 3 (CA in ALA 1, chain B) → select itself
        let hit2 = make_hit(3);
        let sel2 = expand_pick_to_selection(&hit2, 6, &mol);
        assert_eq!(sel2.count(), 1);
        assert!(sel2.contains_index(3));
    }

    // ========================================================================
    // pick_expression_for_hit tests
    // ========================================================================

    #[test]
    fn test_pick_expr_atom_mode() {
        let mol = test_molecule();
        let hit = make_hit(2); // atom 2 = N in GLY 2, chain A
        let expr = pick_expression_for_hit(&hit, 0, &mol).unwrap();
        assert_eq!(expr, "model test and index 2");
    }

    #[test]
    fn test_pick_expr_residue_mode() {
        let mol = test_molecule();
        let hit = make_hit(0); // atom 0 = N in ALA 1, chain A
        let expr = pick_expression_for_hit(&hit, 1, &mol).unwrap();
        assert_eq!(expr, "model test and chain A and resi 1");
    }

    #[test]
    fn test_pick_expr_chain_mode() {
        let mol = test_molecule();
        let hit = make_hit(3); // atom 3 = CA in ALA 1, chain B
        let expr = pick_expression_for_hit(&hit, 2, &mol).unwrap();
        assert_eq!(expr, "model test and chain B");
    }

    #[test]
    fn test_pick_expr_segment_mode() {
        let mol = test_molecule();
        let hit = make_hit(3); // atom 3 has segi "S2"
        let expr = pick_expression_for_hit(&hit, 3, &mol).unwrap();
        assert_eq!(expr, "model test and segi S2");
    }

    #[test]
    fn test_pick_expr_object_mode() {
        let mol = test_molecule();
        let hit = make_hit(0);
        let expr = pick_expression_for_hit(&hit, 4, &mol).unwrap();
        assert_eq!(expr, "model test");
    }

    #[test]
    fn test_pick_expr_molecule_mode() {
        let mol = test_molecule();
        let hit = make_hit(1); // atom 1
        let expr = pick_expression_for_hit(&hit, 5, &mol).unwrap();
        assert_eq!(expr, "bymolecule (model test and index 1)");
    }

    #[test]
    fn test_pick_expr_calpha_mode() {
        let mol = test_molecule();
        let hit = make_hit(0); // atom 0 = N in ALA 1, chain A
        let expr = pick_expression_for_hit(&hit, 6, &mol).unwrap();
        assert_eq!(expr, "model test and chain A and resi 1 and name CA");
    }

    #[test]
    fn test_pick_expr_no_atom_hit() {
        let mol = test_molecule();
        // Hit without atom_index (non-molecule object)
        let hit = PickHit {
            object_name: "test".to_string(),
            object_type: ObjectType::Molecule,
            atom_index: None,
            position: Vec3::new(0.0, 0.0, 0.0),
            distance: 1.0,
        };
        // Atom-level modes return None when no atom is hit
        assert!(pick_expression_for_hit(&hit, 0, &mol).is_none());
        // Object mode still works without an atom
        assert_eq!(pick_expression_for_hit(&hit, 4, &mol).unwrap(), "model test");
    }
}
