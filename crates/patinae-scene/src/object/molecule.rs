//! Molecule object wrapper
//!
//! Wraps `ObjectMolecule` with render state and cached representations.

use lin_alg::f32::Vec3;
use patinae_mol::{AtomIndex, CoordSet, DirtyFlags, ObjectMolecule, RepMask};
use patinae_settings::ObjectOverrides;

use super::{Object, ObjectState, ObjectType};

/// A molecular object with render state.
///
/// Wraps `ObjectMolecule` and manages visual state (enabled, colour,
/// visible reps), per-object setting overrides, dirty flags. The actual
/// GPU representations live on `patinae_render::RenderState`; this
/// struct is pure host-side data.
pub struct MoleculeObject {
    /// The underlying molecular data
    molecule: ObjectMolecule,
    /// Visual state
    state: ObjectState,
    /// Current coordinate state index being displayed
    display_state: usize,
    /// Dirty flags indicating what needs rebuilding
    dirty: DirtyFlags,
    /// Per-object settings overrides
    overrides: Option<ObjectOverrides>,
    /// Surface quality (-4 to 4, default 0)
    surface_quality: i32,
}

impl MoleculeObject {
    /// Create a new molecule object
    pub fn new(mut molecule: ObjectMolecule) -> Self {
        let mut state = ObjectState::default();
        let mut object_reps = RepMask::NONE;
        for atom in molecule.atoms_mut() {
            let flags = atom.state.flags;
            atom.repr.visible_reps = if flags.is_solvent() {
                RepMask::DOTS
            } else if atom.element.is_metal() || flags.is_inorganic() {
                RepMask::SPHERES
            } else if flags.is_organic() {
                RepMask::STICKS
            } else {
                RepMask::CARTOON
            };
            object_reps = object_reps.union(atom.repr.visible_reps);
        }
        state.visible_reps = object_reps;
        state.draw_reps = object_reps;

        Self {
            molecule,
            state,
            display_state: 0,
            dirty: DirtyFlags::ALL,
            surface_quality: 0,
            overrides: None,
        }
    }

    /// Restore a molecule object from a snapshot without resetting per-atom representations.
    ///
    /// Unlike `new()`, this preserves all atom-level state (visible_reps, colors, etc.)
    /// exactly as stored in the molecule.
    pub fn from_raw(mut molecule: ObjectMolecule) -> Self {
        molecule.rebuild_atom_bonds();
        Self {
            molecule,
            state: ObjectState::default(),
            display_state: 0,
            dirty: DirtyFlags::ALL,
            surface_quality: 0,
            overrides: None,
        }
    }

    /// Create a named molecule object preserving per-atom representation state.
    ///
    /// Like `from_raw()`, this does NOT reset per-atom `visible_reps`. The object-level
    /// `visible_reps` is computed as the union of all per-atom reps so the renderer
    /// enables all needed representation pipelines.
    pub fn from_raw_with_name(mut molecule: ObjectMolecule, name: &str) -> Self {
        molecule.name = name.to_string();
        let mut obj_reps = RepMask::default();
        for atom in molecule.atoms() {
            obj_reps = obj_reps.union(atom.repr.visible_reps);
        }
        let state = ObjectState {
            visible_reps: obj_reps,
            draw_reps: obj_reps,
            ..Default::default()
        };
        Self {
            molecule,
            state,
            display_state: 0,
            dirty: DirtyFlags::ALL,
            surface_quality: 0,
            overrides: None,
        }
    }

    /// Create a molecule object with a specific name
    pub fn with_name(molecule: ObjectMolecule, name: &str) -> Self {
        let mut obj = Self::new(molecule);
        obj.molecule.name = name.to_string();
        obj
    }

    /// Get a reference to the underlying molecule data
    pub fn molecule(&self) -> &ObjectMolecule {
        &self.molecule
    }

    /// Get a mutable reference to the underlying molecule data
    ///
    /// Note: This will mark the object as dirty.
    pub fn molecule_mut(&mut self) -> &mut ObjectMolecule {
        self.dirty = DirtyFlags::ALL;
        &mut self.molecule
    }

    /// Get mutable molecule data while marking only the supplied dirty bits.
    pub fn molecule_mut_with_dirty(&mut self, dirty: DirtyFlags) -> &mut ObjectMolecule {
        self.dirty |= dirty;
        &mut self.molecule
    }

    /// Get the display state index
    pub fn display_state(&self) -> usize {
        self.display_state
    }

    /// Get the coordinate set for the state currently displayed by this object.
    pub fn display_coord_set(&self) -> Option<&CoordSet> {
        self.molecule.get_coord_set(self.display_state)
    }

    /// Get one atom coordinate from the state currently displayed by this object.
    pub fn display_coord(&self, atom: AtomIndex) -> Option<Vec3> {
        self.display_coord_set()?.get_atom_coord(atom)
    }

    /// Set the display state index
    pub fn set_display_state(&mut self, state: usize) -> bool {
        if state < self.molecule.state_count() {
            if self.display_state != state {
                self.display_state = state;
                self.dirty |= DirtyFlags::COORDS;
            }
            true
        } else {
            false
        }
    }

    /// Invalidate cached representations
    pub fn invalidate(&mut self, flags: DirtyFlags) {
        self.dirty |= flags;
    }

    /// Clear all dirty flags. Called by the host after it has consumed the
    /// flags (e.g. patinae-render finished its sync walk).
    pub fn clear_dirty(&mut self) {
        self.dirty = DirtyFlags::empty();
    }

    /// Invalidate all representations so they get rebuilt
    ///
    /// This is a convenience method for marking all representations as dirty,
    /// typically used when settings like transparency change.
    pub fn invalidate_representations(&mut self) {
        self.dirty |= DirtyFlags::REPS;
    }

    /// Check if the molecule needs rebuilding
    pub fn is_dirty(&self) -> bool {
        !self.dirty.is_empty()
    }

    /// Get the dirty flags
    pub fn dirty_flags(&self) -> DirtyFlags {
        self.dirty
    }

    /// Get the visible representations mask
    pub fn visible_reps(&self) -> RepMask {
        self.state.visible_reps
    }

    /// Get the object-level representation draw mask.
    pub fn draw_reps(&self) -> RepMask {
        self.state.draw_reps
    }

    /// Get draw-mask reps that can be safely restored without atom rewrites.
    pub fn draw_mask_restorable_reps(&self) -> RepMask {
        self.state.draw_mask_restorable_reps
    }

    /// Check whether a draw-mask-only restore is known to be valid.
    pub fn can_restore_draw_mask(&self, rep: RepMask) -> bool {
        self.state.can_restore_draw_mask(rep)
    }

    /// Set the visible representations mask (e.g., copying from another object)
    pub fn set_visible_reps(&mut self, reps: RepMask) {
        self.state.visible_reps = reps;
        self.state.draw_reps = reps;
        self.state.draw_mask_restorable_reps = RepMask::NONE;
        self.dirty |= DirtyFlags::REPS;
    }

    /// Set the object-level representation draw mask.
    pub fn set_draw_reps(&mut self, reps: RepMask) {
        self.state.draw_reps = reps.intersection(self.state.visible_reps);
        self.dirty |= DirtyFlags::DRAW_MASK;
    }

    /// Show a representation
    ///
    /// This sets visibility both at the object level and for all atoms.
    pub fn show(&mut self, rep: RepMask) {
        if rep.can_toggle_with_draw_mask()
            && self.state.visible_reps.is_visible(rep)
            && !self.state.draw_reps.is_visible(rep)
            && self.state.can_restore_draw_mask(rep)
        {
            self.state.show_draw_rep(rep);
            self.dirty |= DirtyFlags::DRAW_MASK;
            return;
        }

        self.state.visible_reps.set_visible(rep);
        self.state.draw_reps.set_visible(rep);
        self.state.clear_draw_mask_restorable(rep);
        for atom in self.molecule.atoms_mut() {
            atom.repr.visible_reps.set_visible(rep);
        }
        self.dirty |= if rep.can_toggle_with_draw_mask() {
            DirtyFlags::VISIBILITY
        } else {
            DirtyFlags::REPS
        };
    }

    /// Hide a representation
    ///
    /// This hides the representation both at the object level and for all atoms.
    pub fn hide(&mut self, rep: RepMask) {
        if rep.can_toggle_with_draw_mask() {
            self.state.hide_draw_rep(rep);
            self.dirty |= DirtyFlags::DRAW_MASK;
            return;
        }

        self.state.visible_reps.set_hidden(rep);
        self.state.draw_reps.set_hidden(rep);
        self.state.clear_draw_mask_restorable(rep);
        for atom in self.molecule.atoms_mut() {
            atom.repr.visible_reps.set_hidden(rep);
        }
        self.dirty |= DirtyFlags::REPS;
    }

    /// Toggle a representation
    ///
    /// This toggles visibility both at the object level and for all atoms.
    pub fn toggle(&mut self, rep: RepMask) {
        if self.state.draw_reps.is_visible(rep) {
            self.hide(rep);
        } else {
            self.show(rep);
        }
    }

    /// Hide all representations
    pub fn hide_all(&mut self) {
        self.state.visible_reps = RepMask::NONE;
        self.state.draw_reps = RepMask::NONE;
        self.state.draw_mask_restorable_reps = RepMask::NONE;
        self.dirty |= DirtyFlags::REPS;
    }

    /// Set surface quality (-4 to 4, default 0)
    ///
    /// Lower values are faster but coarser, higher values are smoother but slower.
    /// - -4: Very coarse (fastest)
    /// - 0: Default quality
    /// - 4: Very fine (slowest)
    pub fn set_surface_quality(&mut self, quality: i32) {
        let quality = quality.clamp(-4, 4);
        if self.surface_quality != quality {
            self.surface_quality = quality;
            self.dirty.insert(DirtyFlags::REPS);
        }
    }

    /// Get current surface quality
    pub fn surface_quality(&self) -> i32 {
        self.surface_quality
    }

    /// Get per-object overrides (create if needed)
    pub fn get_or_create_overrides(&mut self) -> &mut ObjectOverrides {
        if self.overrides.is_none() {
            self.overrides = Some(ObjectOverrides::default());
        }
        self.overrides.as_mut().unwrap()
    }

    /// Check if cartoon representation is visible
    pub fn is_cartoon_visible(&self) -> bool {
        self.state.draw_reps.is_visible(RepMask::CARTOON)
    }

    /// Check if surface representation is visible
    pub fn is_surface_visible(&self) -> bool {
        self.state.draw_reps.is_visible(RepMask::SURFACE)
    }

    /// Check if mesh representation is visible
    pub fn is_mesh_visible(&self) -> bool {
        self.state.draw_reps.is_visible(RepMask::MESH)
    }

    /// Collect label data for screen-space rendering.
    ///
    /// Returns (world_position, label_text) for each atom that has LABELS visible
    /// and a non-empty label string.
    pub fn collect_labels(&self) -> Vec<(Vec3, &str)> {
        if !self.state.draw_reps.is_visible(RepMask::LABELS) {
            return Vec::new();
        }

        let coord_set = match self.molecule.get_coord_set(self.display_state) {
            Some(cs) => cs,
            None => return Vec::new(),
        };

        let mut labels = Vec::new();
        for (atom_idx, coord) in coord_set.iter_with_atoms() {
            let atom = match self.molecule.get_atom(atom_idx) {
                Some(a) => a,
                None => continue,
            };
            if atom.repr.visible_reps.is_visible(RepMask::LABELS) && !atom.repr.label.is_empty() {
                labels.push((coord, atom.repr.label.as_str()));
            }
        }
        labels
    }
}

impl Object for MoleculeObject {
    fn name(&self) -> &str {
        &self.molecule.name
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
        self.molecule.bounding_box(self.display_state)
    }

    fn n_states(&self) -> usize {
        self.molecule.state_count().max(1)
    }

    fn current_state(&self) -> usize {
        self.display_state
    }

    fn set_current_state(&mut self, state: usize) -> bool {
        self.set_display_state(state)
    }

    fn overrides(&self) -> Option<&ObjectOverrides> {
        self.overrides.as_ref()
    }

    fn overrides_mut(&mut self) -> Option<&mut ObjectOverrides> {
        self.overrides.as_mut()
    }

    fn get_or_create_overrides(&mut self) -> &mut ObjectOverrides {
        self.get_or_create_overrides()
    }

    fn set_name(&mut self, name: String) {
        self.molecule.name = name;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use patinae_mol::{Atom, AtomResidue, BondOrder, Element};
    use std::sync::Arc;

    fn create_test_molecule() -> ObjectMolecule {
        let mut mol = ObjectMolecule::new("test");
        let c1 = mol.add_atom(Atom::new("C1", Element::Carbon));
        let c2 = mol.add_atom(Atom::new("C2", Element::Carbon));
        mol.add_bond(c1, c2, BondOrder::Single).unwrap();

        let cs =
            patinae_mol::CoordSet::from_vec3(&[Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.5, 0.0, 0.0)]);
        mol.add_coord_set(cs);

        mol
    }

    fn collect_surface_visibility_mask(
        mol: &ObjectMolecule,
        coord_set: &patinae_mol::CoordSet,
    ) -> Vec<u32> {
        let mut mask = vec![0; coord_set.len().div_ceil(32)];
        for (atom_idx, _) in coord_set.iter_with_atoms() {
            let i = atom_idx.as_usize();
            if mol
                .get_atom(atom_idx)
                .is_some_and(|atom| atom.repr.visible_reps.is_visible(RepMask::SURFACE))
            {
                mask[i / 32] |= 1 << (i % 32);
            }
        }
        mask
    }

    #[test]
    fn test_molecule_object_creation() {
        let mol = create_test_molecule();
        let obj = MoleculeObject::new(mol);

        assert_eq!(obj.name(), "test");
        assert_eq!(obj.object_type(), ObjectType::Molecule);
        assert!(obj.state().enabled);
    }

    #[test]
    fn default_object_rep_mask_matches_atom_rep_union() {
        let mol = create_test_molecule();
        let obj = MoleculeObject::new(mol);

        assert_eq!(obj.visible_reps(), RepMask::CARTOON);
        assert_eq!(obj.draw_reps(), RepMask::CARTOON);
        assert!(obj
            .molecule()
            .atoms()
            .all(|atom| atom.repr.visible_reps == RepMask::CARTOON));
    }

    #[test]
    fn test_dirty_flags() {
        let mol = create_test_molecule();
        let mut obj = MoleculeObject::new(mol);

        assert!(obj.is_dirty());

        // Accessing molecule_mut should set dirty
        let _ = obj.molecule_mut();
        assert!(obj.dirty_flags().contains(DirtyFlags::ALL));
    }

    #[test]
    fn test_show_hide() {
        let mol = create_test_molecule();
        let mut obj = MoleculeObject::new(mol);

        obj.hide(RepMask::LINES);
        assert!(!obj.state().visible_reps.is_visible(RepMask::LINES));

        obj.show(RepMask::SPHERES);
        assert!(obj.state().visible_reps.is_visible(RepMask::SPHERES));
    }

    #[test]
    fn whole_object_cartoon_hide_uses_draw_mask_without_mutating_atoms() {
        let mol = create_test_molecule();
        let mut obj = MoleculeObject::new(mol);
        obj.clear_dirty();

        obj.hide(RepMask::CARTOON);

        assert!(obj.visible_reps().is_visible(RepMask::CARTOON));
        assert!(!obj.draw_reps().is_visible(RepMask::CARTOON));
        assert!(obj.draw_mask_restorable_reps().is_visible(RepMask::CARTOON));
        assert!(obj
            .molecule()
            .atoms()
            .all(|atom| atom.repr.visible_reps.is_visible(RepMask::CARTOON)));
        assert_eq!(obj.dirty_flags(), DirtyFlags::DRAW_MASK);

        obj.clear_dirty();
        obj.show(RepMask::CARTOON);

        assert!(obj.visible_reps().is_visible(RepMask::CARTOON));
        assert!(obj.draw_reps().is_visible(RepMask::CARTOON));
        assert!(obj.draw_mask_restorable_reps().is_visible(RepMask::CARTOON));
        assert_eq!(obj.dirty_flags(), DirtyFlags::DRAW_MASK);
    }

    #[test]
    fn first_whole_object_cartoon_show_materializes_atom_bits_once() {
        let mol = create_test_molecule();
        let mut obj = MoleculeObject::new(mol);
        obj.state_mut().visible_reps.set_hidden(RepMask::CARTOON);
        obj.state_mut().draw_reps.set_hidden(RepMask::CARTOON);
        for atom in obj.molecule_mut_with_dirty(DirtyFlags::empty()).atoms_mut() {
            atom.repr.visible_reps.set_hidden(RepMask::CARTOON);
        }
        obj.clear_dirty();

        obj.show(RepMask::CARTOON);

        assert!(obj.visible_reps().is_visible(RepMask::CARTOON));
        assert!(obj.draw_reps().is_visible(RepMask::CARTOON));
        assert!(obj
            .molecule()
            .atoms()
            .all(|atom| atom.repr.visible_reps.is_visible(RepMask::CARTOON)));
        assert_eq!(obj.dirty_flags(), DirtyFlags::VISIBILITY);
    }

    #[test]
    fn test_metal_inside_organic_ligand_defaults_to_solid_sphere() {
        let mut mol = ObjectMolecule::new("heme");
        let residue = Arc::new(AtomResidue::from_parts("A", "HEM", 142, ' ', ""));

        let mut carbon = Atom::new("C1", Element::Carbon);
        carbon.residue = residue.clone();
        carbon.state.hetatm = true;
        mol.add_atom(carbon);

        let mut iron = Atom::new("FE", Element::Iron);
        iron.residue = residue;
        iron.state.hetatm = true;
        mol.add_atom(iron);

        mol.add_coord_set(patinae_mol::CoordSet::from_vec3(&[
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.5, 0.0, 0.0),
        ]));
        mol.classify_atoms();

        let obj = MoleculeObject::new(mol);
        let atoms = obj.molecule().atoms_slice();

        assert!(atoms[0].repr.visible_reps.is_visible(RepMask::STICKS));
        assert!(atoms[1].repr.visible_reps.is_visible(RepMask::SPHERES));
        assert_eq!(obj.visible_reps(), RepMask::STICKS.union(RepMask::SPHERES));
        assert_eq!(obj.draw_reps(), obj.visible_reps());
        assert_eq!(atoms[1].repr.sphere_transparency, None);
    }

    #[test]
    fn test_collect_surface_visibility_mask_packs_bits_lsb_first() {
        let mut mol = ObjectMolecule::new("vis");
        for i in 0..40 {
            mol.add_atom(Atom::new(format!("C{}", i), Element::Carbon));
        }
        let coords: Vec<Vec3> = (0..40).map(|i| Vec3::new(i as f32, 0.0, 0.0)).collect();
        mol.add_coord_set(patinae_mol::CoordSet::from_vec3(&coords));

        // Show surface only for atoms 0, 5, 32, 39.
        for &i in &[0u32, 5, 32, 39] {
            mol.get_atom_mut(patinae_mol::AtomIndex(i))
                .unwrap()
                .repr
                .visible_reps
                .set_visible(RepMask::SURFACE);
        }

        let coord_set = mol.get_coord_set(0).unwrap();
        let mask = collect_surface_visibility_mask(&mol, coord_set);
        // 40 atoms → ceil(40/32) = 2 words.
        assert_eq!(mask.len(), 2);
        assert_eq!(mask[0], (1 << 0) | (1 << 5), "word 0: atoms 0 and 5");
        assert_eq!(mask[1], (1 << 0) | (1 << 7), "word 1: atoms 32 and 39");
    }

    #[test]
    fn test_dirty_flags_include_transparency_and_visibility() {
        // Sanity check: ALL covers the new bits so existing call sites that
        // mark ALL keep working.
        assert!(DirtyFlags::ALL.contains(DirtyFlags::TRANSPARENCY));
        assert!(DirtyFlags::ALL.contains(DirtyFlags::VISIBILITY));
        assert!(DirtyFlags::ALL.contains(DirtyFlags::DRAW_MASK));
        assert!(DirtyFlags::TRANSPARENCY.bits() != DirtyFlags::COLOR.bits());
        assert!(DirtyFlags::VISIBILITY.bits() != DirtyFlags::REPS.bits());
        assert!(DirtyFlags::DRAW_MASK.bits() != DirtyFlags::VISIBILITY.bits());
    }

    #[test]
    fn test_display_state() {
        let mut mol = ObjectMolecule::new("test");
        mol.add_atom(Atom::new("C", Element::Carbon));

        // Add multiple coordinate sets
        mol.add_coord_set(patinae_mol::CoordSet::from_vec3(&[Vec3::new(
            0.0, 0.0, 0.0,
        )]));
        mol.add_coord_set(patinae_mol::CoordSet::from_vec3(&[Vec3::new(
            1.0, 1.0, 1.0,
        )]));

        let mut obj = MoleculeObject::new(mol);

        assert_eq!(obj.n_states(), 2);
        assert_eq!(obj.display_state(), 0);

        assert!(obj.set_display_state(1));
        assert_eq!(obj.display_state(), 1);

        // Invalid state
        assert!(!obj.set_display_state(5));
    }
}
