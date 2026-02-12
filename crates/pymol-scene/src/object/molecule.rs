//! Molecule object wrapper
//!
//! Wraps `ObjectMolecule` with render state and cached representations.

use bitflags::bitflags;
use lin_alg::f32::Vec3;
use pymol_mol::{ObjectMolecule, RepMask};
use pymol_render::{
    CartoonRep, ColorResolver, DotRep, LineRep, RenderContext, Representation, RibbonRep,
    SelectionIndicatorRep, SphereRep, StickRep, SurfaceRep, WireSurfaceRep,
};
use pymol_select::SelectionResult;
use pymol_settings::{GlobalSettings, SettingResolver};

use super::{Object, ObjectState, ObjectType, ObjectWithName};

bitflags! {
    /// Flags indicating which aspects of a molecule need rebuilding
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub struct DirtyFlags: u32 {
        /// Coordinates have changed
        const COORDS = 0x01;
        /// Colors have changed
        const COLOR = 0x02;
        /// Representation settings have changed
        const REPS = 0x04;
        /// Atom selection/visibility has changed
        const SELECTION = 0x08;
        /// Everything needs rebuilding
        const ALL = 0x0F;
    }
}

/// Cached representations for a molecule
///
/// These are lazily built and cached until invalidated.
struct RepresentationCache {
    spheres: Option<SphereRep>,
    sticks: Option<StickRep>,
    lines: Option<LineRep>,
    dots: Option<DotRep>,
    cartoon: Option<CartoonRep>,
    ribbon: Option<RibbonRep>,
    surface: Option<SurfaceRep>,
    mesh: Option<WireSurfaceRep>,
    /// Selection indicator representation (rendered last, on top of everything)
    selection_indicator: Option<SelectionIndicatorRep>,
}

impl Default for RepresentationCache {
    fn default() -> Self {
        Self {
            spheres: None,
            sticks: None,
            lines: None,
            dots: None,
            cartoon: None,
            ribbon: None,
            surface: None,
            mesh: None,
            selection_indicator: None,
        }
    }
}

impl RepresentationCache {
    /// Clear all cached representations
    fn clear(&mut self) {
        self.spheres = None;
        self.sticks = None;
        self.lines = None;
        self.dots = None;
        self.cartoon = None;
        self.ribbon = None;
        self.surface = None;
        self.mesh = None;
        self.selection_indicator = None;
    }
}

/// A molecular object with render state
///
/// Wraps `ObjectMolecule` and manages:
/// - Visual state (enabled, color, visible reps)
/// - Cached GPU representations
/// - Dirty flags for invalidation
pub struct MoleculeObject {
    /// The underlying molecular data
    molecule: ObjectMolecule,
    /// Visual state
    state: ObjectState,
    /// Current coordinate state index being displayed
    display_state: usize,
    /// Cached representations
    representations: RepresentationCache,
    /// Dirty flags indicating what needs rebuilding
    dirty: DirtyFlags,
    /// Per-object settings override
    settings: Option<GlobalSettings>,
    /// Surface quality (-4 to 4, default 0)
    surface_quality: i32,
}

impl MoleculeObject {
    /// Create a new molecule object
    pub fn new(molecule: ObjectMolecule) -> Self {
        let mut state = ObjectState::default();
        // Default to lines representation
        state.visible_reps.set_visible(RepMask::LINES);

        Self {
            molecule,
            state,
            display_state: 0,
            representations: RepresentationCache::default(),
            dirty: DirtyFlags::ALL,
            surface_quality: 0,
            settings: None,
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

    /// Get the display state index
    pub fn display_state(&self) -> usize {
        self.display_state
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

    /// Show a representation
    ///
    /// This sets visibility both at the object level and for all atoms.
    pub fn show(&mut self, rep: RepMask) {
        self.state.visible_reps.set_visible(rep);
        // Also update all atoms to have this rep visible
        for atom in self.molecule.atoms_mut() {
            atom.repr.visible_reps.set_visible(rep);
        }
        self.dirty |= DirtyFlags::REPS;
    }

    /// Hide a representation
    ///
    /// This hides the representation both at the object level and for all atoms.
    pub fn hide(&mut self, rep: RepMask) {
        self.state.visible_reps.set_hidden(rep);
        // Also update all atoms to hide this rep
        for atom in self.molecule.atoms_mut() {
            atom.repr.visible_reps.set_hidden(rep);
        }
        self.dirty |= DirtyFlags::REPS;
    }

    /// Toggle a representation
    ///
    /// This toggles visibility both at the object level and for all atoms.
    pub fn toggle(&mut self, rep: RepMask) {
        self.state.visible_reps.toggle(rep);
        // Also update all atoms' visibility to match
        if self.state.visible_reps.is_visible(rep) {
            for atom in self.molecule.atoms_mut() {
                atom.repr.visible_reps.set_visible(rep);
            }
        } else {
            for atom in self.molecule.atoms_mut() {
                atom.repr.visible_reps.set_hidden(rep);
            }
        }
        self.dirty |= DirtyFlags::REPS;
    }

    /// Hide all representations
    pub fn hide_all(&mut self) {
        self.state.visible_reps = RepMask::NONE;
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
            if let Some(surface) = &mut self.representations.surface {
                surface.set_quality(quality);
            }
            self.dirty.insert(DirtyFlags::REPS);
        }
    }

    /// Get current surface quality
    pub fn surface_quality(&self) -> i32 {
        self.surface_quality
    }

    /// Set the selection indicator for this molecule
    ///
    /// This renders pink/magenta indicators at the selected atom positions.
    ///
    /// # Arguments
    /// * `selection` - The selection result indicating which atoms to show
    /// * `context` - The render context for uploading GPU data
    pub fn set_selection_indicator(
        &mut self,
        selection: &SelectionResult,
        context: &RenderContext,
    ) {
        self.set_selection_indicator_with_size(selection, context, None);
    }

    /// Set the selection indicator for this molecule with a custom size
    ///
    /// This renders pink/magenta indicators at the selected atom positions.
    ///
    /// # Arguments
    /// * `selection` - The selection result indicating which atoms to show
    /// * `context` - The render context for uploading GPU data
    /// * `size` - Optional custom size for the indicator dots
    pub fn set_selection_indicator_with_size(
        &mut self,
        selection: &SelectionResult,
        context: &RenderContext,
        size: Option<f32>,
    ) {
        // Get the current coordinate set
        let coord_set = match self.molecule.get_coord_set(self.display_state) {
            Some(cs) => cs,
            None => {
                log::debug!("No coord set available for selection indicator");
                return;
            }
        };

        // Create or get the selection indicator representation
        let indicator = self
            .representations
            .selection_indicator
            .get_or_insert_with(SelectionIndicatorRep::new);

        // Apply custom size if provided - scale up for visibility in the shader
        // The dot shader scales by 0.01, so we need larger values
        let effective_size = size.unwrap_or(pymol_render::DEFAULT_INDICATOR_SIZE) * 10.0;
        indicator.set_size(effective_size);

        // Build and upload the indicator
        indicator.build_for_selection(&self.molecule, coord_set, selection);

        log::debug!(
            "Built selection indicator with {} instances, size={}",
            indicator.instance_count(),
            indicator.size()
        );

        indicator.upload(context.device(), context.queue());
    }

    /// Clear the selection indicator for this molecule
    pub fn clear_selection_indicator(&mut self) {
        if let Some(ref mut indicator) = self.representations.selection_indicator {
            indicator.clear();
        }
    }

    /// Check if this molecule has a selection indicator
    pub fn has_selection_indicator(&self) -> bool {
        self.representations
            .selection_indicator
            .as_ref()
            .map(|i| !i.is_empty())
            .unwrap_or(false)
    }

    /// Get per-object settings (create if needed)
    pub fn get_or_create_settings(&mut self) -> &mut GlobalSettings {
        if self.settings.is_none() {
            self.settings = Some(GlobalSettings::new());
        }
        self.settings.as_mut().unwrap()
    }

    /// Build/rebuild representations for rendering
    ///
    /// This should be called before rendering if `is_dirty()` returns true.
    pub fn prepare_render(
        &mut self,
        context: &RenderContext,
        color_resolver: &ColorResolver,
        global_settings: &GlobalSettings,
    ) {
        if self.dirty.is_empty() {
            return;
        }

        // Get the current coordinate set
        let coord_set = match self.molecule.get_coord_set(self.display_state) {
            Some(cs) => cs,
            None => return,
        };

        // Create setting resolver with object settings
        let settings = if let Some(ref obj_settings) = self.settings {
            SettingResolver::with_object(global_settings, obj_settings)
        } else {
            SettingResolver::global(global_settings)
        };

        let vis = &self.state.visible_reps;

        // Build spheres if visible
        if vis.is_visible(RepMask::SPHERES) {
            let spheres = self.representations.spheres.get_or_insert_with(SphereRep::new);
            if self.dirty.intersects(DirtyFlags::COORDS | DirtyFlags::COLOR | DirtyFlags::REPS) {
                spheres.build(&self.molecule, coord_set, color_resolver, &settings);
                spheres.upload(context.device(), context.queue());
            }
        }

        // Build sticks if visible
        if vis.is_visible(RepMask::STICKS) {
            let sticks = self.representations.sticks.get_or_insert_with(StickRep::new);
            if self.dirty.intersects(DirtyFlags::COORDS | DirtyFlags::COLOR | DirtyFlags::REPS) {
                sticks.build(&self.molecule, coord_set, color_resolver, &settings);
                sticks.upload(context.device(), context.queue());
            }
        }

        // Build lines if visible
        if vis.is_visible(RepMask::LINES) {
            let lines = self.representations.lines.get_or_insert_with(LineRep::new);
            if self.dirty.intersects(DirtyFlags::COORDS | DirtyFlags::COLOR | DirtyFlags::REPS) {
                lines.build(&self.molecule, coord_set, color_resolver, &settings);
                lines.upload(context.device(), context.queue());
            }
        }

        // Build dots if visible
        if vis.is_visible(RepMask::DOTS) {
            let dots = self.representations.dots.get_or_insert_with(DotRep::new);
            if self.dirty.intersects(DirtyFlags::COORDS | DirtyFlags::COLOR | DirtyFlags::REPS) {
                dots.build(&self.molecule, coord_set, color_resolver, &settings);
                dots.upload(context.device(), context.queue());
            }
        }

        // Build cartoon if visible
        if vis.is_visible(RepMask::CARTOON) {
            let cartoon = self.representations.cartoon.get_or_insert_with(CartoonRep::new);
            if self.dirty.intersects(DirtyFlags::COORDS | DirtyFlags::COLOR | DirtyFlags::REPS) {
                cartoon.build(&self.molecule, coord_set, color_resolver, &settings);
                cartoon.upload(context.device(), context.queue());
            }
        }

        // Build ribbon if visible
        if vis.is_visible(RepMask::RIBBON) {
            let ribbon = self.representations.ribbon.get_or_insert_with(RibbonRep::new);
            if self.dirty.intersects(DirtyFlags::COORDS | DirtyFlags::COLOR | DirtyFlags::REPS) {
                ribbon.build(&self.molecule, coord_set, color_resolver, &settings);
                ribbon.upload(context.device(), context.queue());
            }
        }

        // Build surface if visible
        if vis.is_visible(RepMask::SURFACE) {
            let quality = self.surface_quality;
            let surface = self.representations.surface.get_or_insert_with(|| {
                let mut s = SurfaceRep::new();
                s.set_quality(quality);
                s
            });
            if self.dirty.intersects(DirtyFlags::COORDS | DirtyFlags::COLOR | DirtyFlags::REPS) {
                surface.set_quality(quality); // Ensure quality is up to date
                surface.build(&self.molecule, coord_set, color_resolver, &settings);
                surface.upload(context.device(), context.queue());
            }
        }

        // Build mesh (wireframe surface) if visible
        if vis.is_visible(RepMask::MESH) {
            let quality = self.surface_quality;
            let mesh = self.representations.mesh.get_or_insert_with(|| {
                let mut m = WireSurfaceRep::new();
                m.set_quality(quality);
                m
            });
            if self.dirty.intersects(DirtyFlags::COORDS | DirtyFlags::COLOR | DirtyFlags::REPS) {
                mesh.set_quality(quality);
                mesh.build(&self.molecule, coord_set, color_resolver, &settings);
                mesh.upload(context.device(), context.queue());
            }
        }

        // Clear dirty flags
        self.dirty = DirtyFlags::empty();
    }

    /// Render all visible representations
    pub fn render<'a>(&'a self, render_pass: &mut wgpu::RenderPass<'a>, context: &'a RenderContext) {
        if !self.state.enabled {
            return;
        }

        let vis = &self.state.visible_reps;

        // Render in appropriate order (opaque first, then transparent)

        // Lines
        if vis.is_visible(RepMask::LINES) {
            if let Some(ref lines) = self.representations.lines {
                if !lines.is_empty() {
                    let pipeline = context.line_pipeline(pymol_render::pipeline::BlendMode::Opaque);
                    render_pass.set_pipeline(&pipeline);
                    render_pass.set_bind_group(0, context.uniform_bind_group(), &[]);
                    lines.render(render_pass);
                }
            }
        }

        // Sticks (cylinders + sphere caps)
        if vis.is_visible(RepMask::STICKS) {
            if let Some(ref sticks) = self.representations.sticks {
                if !sticks.is_empty() {
                    // Render cylinders
                    let cylinder_pipeline = context.cylinder_pipeline(pymol_render::pipeline::BlendMode::Opaque);
                    render_pass.set_pipeline(&cylinder_pipeline);
                    render_pass.set_bind_group(0, context.uniform_bind_group(), &[]);
                    render_pass.set_vertex_buffer(0, context.billboard_vertex_buffer().slice(..));
                    render_pass.set_index_buffer(
                        context.quad_index_buffer().slice(..),
                        wgpu::IndexFormat::Uint16,
                    );
                    sticks.render(render_pass);

                    // Render sphere caps at atom positions
                    if sticks.sphere_count() > 0 {
                        if let Some(sphere_buffer) = sticks.sphere_buffer() {
                            let sphere_pipeline = context.sphere_pipeline(pymol_render::pipeline::BlendMode::Opaque);
                            render_pass.set_pipeline(&sphere_pipeline);
                            render_pass.set_bind_group(0, context.uniform_bind_group(), &[]);
                            render_pass.set_vertex_buffer(0, context.billboard_vertex_buffer().slice(..));
                            render_pass.set_index_buffer(
                                context.quad_index_buffer().slice(..),
                                wgpu::IndexFormat::Uint16,
                            );
                            render_pass.set_vertex_buffer(1, sphere_buffer.slice(..));
                            render_pass.draw_indexed(0..6, 0, 0..sticks.sphere_count());
                        }
                    }
                }
            }
        }

        // Spheres
        if vis.is_visible(RepMask::SPHERES) {
            if let Some(ref spheres) = self.representations.spheres {
                if !spheres.is_empty() {
                    let pipeline = context.sphere_pipeline(pymol_render::pipeline::BlendMode::Opaque);
                    render_pass.set_pipeline(&pipeline);
                    render_pass.set_bind_group(0, context.uniform_bind_group(), &[]);
                    render_pass.set_vertex_buffer(0, context.billboard_vertex_buffer().slice(..));
                    render_pass.set_index_buffer(
                        context.quad_index_buffer().slice(..),
                        wgpu::IndexFormat::Uint16,
                    );
                    spheres.render(render_pass);
                }
            }
        }

        // Dots
        if vis.is_visible(RepMask::DOTS) {
            if let Some(ref dots) = self.representations.dots {
                if !dots.is_empty() {
                    let pipeline = context.dot_pipeline(pymol_render::pipeline::BlendMode::Opaque);
                    render_pass.set_pipeline(&pipeline);
                    render_pass.set_bind_group(0, context.uniform_bind_group(), &[]);
                    render_pass.set_vertex_buffer(0, context.billboard_vertex_buffer().slice(..));
                    render_pass.set_index_buffer(
                        context.quad_index_buffer().slice(..),
                        wgpu::IndexFormat::Uint16,
                    );
                    dots.render(render_pass);
                }
            }
        }

        // Cartoon (mesh-based)
        if vis.is_visible(RepMask::CARTOON) {
            if let Some(ref cartoon) = self.representations.cartoon {
                if !cartoon.is_empty() {
                    let pipeline = context.mesh_pipeline(pymol_render::pipeline::BlendMode::Opaque);
                    render_pass.set_pipeline(&pipeline);
                    render_pass.set_bind_group(0, context.uniform_bind_group(), &[]);
                    cartoon.render(render_pass);
                }
            }
        }

        // Ribbon (mesh-based, uniform tube representation)
        if vis.is_visible(RepMask::RIBBON) {
            if let Some(ref ribbon) = self.representations.ribbon {
                if !ribbon.is_empty() {
                    let pipeline = context.mesh_pipeline(pymol_render::pipeline::BlendMode::Opaque);
                    render_pass.set_pipeline(&pipeline);
                    render_pass.set_bind_group(0, context.uniform_bind_group(), &[]);
                    ribbon.render(render_pass);
                }
            }
        }

        // Surface
        if vis.is_visible(RepMask::SURFACE) {
            if let Some(ref surface) = self.representations.surface {
                if !surface.is_empty() {
                    // Use transparent blend mode if the surface has any transparency
                    let blend_mode = if surface.has_transparency() {
                        pymol_render::pipeline::BlendMode::Transparent
                    } else {
                        pymol_render::pipeline::BlendMode::Opaque
                    };
                    let pipeline = context.mesh_pipeline(blend_mode);
                    render_pass.set_pipeline(&pipeline);
                    render_pass.set_bind_group(0, context.uniform_bind_group(), &[]);
                    surface.render(render_pass);
                }
            }
        }

        // Mesh (wireframe surface)
        if vis.is_visible(RepMask::MESH) {
            if let Some(ref mesh) = self.representations.mesh {
                if !mesh.is_empty() {
                    let pipeline = context.line_pipeline(pymol_render::pipeline::BlendMode::Opaque);
                    render_pass.set_pipeline(&pipeline);
                    render_pass.set_bind_group(0, context.uniform_bind_group(), &[]);
                    mesh.render(render_pass);
                }
            }
        }

        // Selection indicator (rendered last, on top of everything)
        if let Some(ref indicator) = self.representations.selection_indicator {
            if !indicator.is_empty() {
                log::debug!(
                    "Rendering selection indicator with {} instances",
                    indicator.instance_count()
                );
                let pipeline = context.dot_pipeline(pymol_render::pipeline::BlendMode::Opaque);
                render_pass.set_pipeline(&pipeline);
                render_pass.set_bind_group(0, context.uniform_bind_group(), &[]);
                render_pass.set_vertex_buffer(0, context.billboard_vertex_buffer().slice(..));
                render_pass.set_index_buffer(
                    context.quad_index_buffer().slice(..),
                    wgpu::IndexFormat::Uint16,
                );
                indicator.render(render_pass);
            }
        }
    }

    /// Clear all cached representations (forces rebuild)
    pub fn clear_cache(&mut self) {
        self.representations.clear();
        self.dirty = DirtyFlags::ALL;
    }

    /// Get cartoon mesh data for raytracing
    ///
    /// Returns (vertices, indices) if cartoon representation is visible and built.
    /// Each vertex contains position, normal, and color.
    pub fn get_cartoon_mesh(&self) -> Option<(&[pymol_render::MeshVertex], &[u32])> {
        if !self.state.visible_reps.is_visible(RepMask::CARTOON) {
            return None;
        }
        self.representations.cartoon.as_ref().map(|c| (c.vertices(), c.indices()))
    }

    /// Get surface mesh data for raytracing
    ///
    /// Returns (vertices, indices) if surface representation is visible and built.
    pub fn get_surface_mesh(&self) -> Option<(&[pymol_render::MeshVertex], &[u32])> {
        if !self.state.visible_reps.is_visible(RepMask::SURFACE) {
            return None;
        }
        self.representations.surface.as_ref().map(|s| (s.vertices(), s.indices()))
    }

    /// Get mesh wireframe edge data for raytracing
    ///
    /// Returns line vertex pairs (each consecutive pair is an edge) if mesh
    /// representation is visible and built. These are rendered as thin cylinders
    /// (sausage primitives) by the raytracer.
    pub fn get_mesh_edges(&self) -> Option<&[pymol_render::LineVertex]> {
        if !self.state.visible_reps.is_visible(RepMask::MESH) {
            return None;
        }
        self.representations
            .mesh
            .as_ref()
            .filter(|m| !m.is_empty())
            .map(|m| m.vertices())
    }

    /// Get ribbon mesh data for raytracing
    ///
    /// Returns (vertices, indices) if ribbon representation is visible and built.
    pub fn get_ribbon_mesh(&self) -> Option<(&[pymol_render::MeshVertex], &[u32])> {
        if !self.state.visible_reps.is_visible(RepMask::RIBBON) {
            return None;
        }
        self.representations.ribbon.as_ref().map(|r| (r.vertices(), r.indices()))
    }

    /// Check if cartoon representation is visible
    pub fn is_cartoon_visible(&self) -> bool {
        self.state.visible_reps.is_visible(RepMask::CARTOON)
    }

    /// Check if surface representation is visible
    pub fn is_surface_visible(&self) -> bool {
        self.state.visible_reps.is_visible(RepMask::SURFACE)
    }

    /// Check if mesh representation is visible
    pub fn is_mesh_visible(&self) -> bool {
        self.state.visible_reps.is_visible(RepMask::MESH)
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

    fn settings(&self) -> Option<&GlobalSettings> {
        self.settings.as_ref()
    }

    fn settings_mut(&mut self) -> Option<&mut GlobalSettings> {
        self.settings.as_mut()
    }
}

impl ObjectWithName for MoleculeObject {
    fn set_name(&mut self, name: String) {
        self.molecule.name = name;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use pymol_mol::{Atom, BondOrder, Element};

    fn create_test_molecule() -> ObjectMolecule {
        let mut mol = ObjectMolecule::new("test");
        let c1 = mol.add_atom(Atom::new("C1", Element::Carbon));
        let c2 = mol.add_atom(Atom::new("C2", Element::Carbon));
        mol.add_bond(c1, c2, BondOrder::Single).unwrap();

        let cs = pymol_mol::CoordSet::from_vec3(&[
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.5, 0.0, 0.0),
        ]);
        mol.add_coord_set(cs);

        mol
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
    fn test_display_state() {
        let mut mol = ObjectMolecule::new("test");
        mol.add_atom(Atom::new("C", Element::Carbon));

        // Add multiple coordinate sets
        mol.add_coord_set(pymol_mol::CoordSet::from_vec3(&[Vec3::new(0.0, 0.0, 0.0)]));
        mol.add_coord_set(pymol_mol::CoordSet::from_vec3(&[Vec3::new(1.0, 1.0, 1.0)]));

        let mut obj = MoleculeObject::new(mol);

        assert_eq!(obj.n_states(), 2);
        assert_eq!(obj.display_state(), 0);

        assert!(obj.set_display_state(1));
        assert_eq!(obj.display_state(), 1);

        // Invalid state
        assert!(!obj.set_display_state(5));
    }
}
