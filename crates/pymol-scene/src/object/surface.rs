//! Standalone surface object
//!
//! Surfaces are stored as independent objects, not tied to molecules.
//! This allows surfaces to persist after molecules are modified or deleted.

use lin_alg::f32::Vec3;
use pymol_render::{MeshRep, MeshVertex, RenderContext, Representation, SurfaceType};

use super::{Object, ObjectState, ObjectType};

/// A standalone surface object
///
/// Unlike the surface representation on MoleculeObject, SurfaceObject stores
/// the triangulated surface as an independent object. This allows:
/// - Persistence after source molecule changes
/// - Independent coloring and transparency
/// - Export/import of surface data
/// - Scene storage of surface state
pub struct SurfaceObject {
    /// Surface name
    name: String,
    /// Visual state (enabled, color, etc.)
    state: ObjectState,
    /// Mesh vertices (position, normal, color)
    vertices: Vec<MeshVertex>,
    /// Triangle indices
    indices: Vec<u32>,
    /// Reference to the source molecule name (if any)
    source_molecule: Option<String>,
    /// Surface type used to generate this surface
    surface_type: SurfaceType,
    /// Cached MeshRep for GPU rendering
    cached_rep: Option<MeshRep>,
    /// Whether the GPU buffers need updating
    dirty: bool,
}

impl SurfaceObject {
    /// Create a new empty surface object
    pub fn new(name: &str) -> Self {
        Self {
            name: name.to_string(),
            state: ObjectState::default(),
            vertices: Vec::new(),
            indices: Vec::new(),
            source_molecule: None,
            surface_type: SurfaceType::default(),
            cached_rep: None,
            dirty: true,
        }
    }

    /// Create a surface object from pre-computed mesh data
    pub fn from_mesh(
        name: &str,
        vertices: Vec<MeshVertex>,
        indices: Vec<u32>,
    ) -> Self {
        Self {
            name: name.to_string(),
            state: ObjectState::default(),
            vertices,
            indices,
            source_molecule: None,
            surface_type: SurfaceType::default(),
            cached_rep: None,
            dirty: true,
        }
    }

    /// Create a surface object from mesh data with source reference
    pub fn from_molecule_surface(
        name: &str,
        vertices: Vec<MeshVertex>,
        indices: Vec<u32>,
        source_molecule: &str,
        surface_type: SurfaceType,
    ) -> Self {
        Self {
            name: name.to_string(),
            state: ObjectState::default(),
            vertices,
            indices,
            source_molecule: Some(source_molecule.to_string()),
            surface_type,
            cached_rep: None,
            dirty: true,
        }
    }

    /// Get the vertices
    pub fn vertices(&self) -> &[MeshVertex] {
        &self.vertices
    }

    /// Get mutable access to vertices (marks dirty)
    pub fn vertices_mut(&mut self) -> &mut Vec<MeshVertex> {
        self.dirty = true;
        &mut self.vertices
    }

    /// Get the indices
    pub fn indices(&self) -> &[u32] {
        &self.indices
    }

    /// Get mutable access to indices (marks dirty)
    pub fn indices_mut(&mut self) -> &mut Vec<u32> {
        self.dirty = true;
        &mut self.indices
    }

    /// Set mesh data
    pub fn set_mesh(&mut self, vertices: Vec<MeshVertex>, indices: Vec<u32>) {
        self.vertices = vertices;
        self.indices = indices;
        self.dirty = true;
    }

    /// Clear the surface mesh
    pub fn clear(&mut self) {
        self.vertices.clear();
        self.indices.clear();
        self.dirty = true;
    }

    /// Get the source molecule name (if any)
    pub fn source_molecule(&self) -> Option<&str> {
        self.source_molecule.as_deref()
    }

    /// Set the source molecule reference
    pub fn set_source_molecule(&mut self, name: Option<String>) {
        self.source_molecule = name;
    }

    /// Get the surface type
    pub fn surface_type(&self) -> SurfaceType {
        self.surface_type
    }

    /// Set the surface type
    pub fn set_surface_type(&mut self, surface_type: SurfaceType) {
        self.surface_type = surface_type;
    }

    /// Get the number of triangles
    pub fn triangle_count(&self) -> usize {
        self.indices.len() / 3
    }

    /// Get the number of vertices
    pub fn vertex_count(&self) -> usize {
        self.vertices.len()
    }

    /// Check if the surface is empty
    pub fn is_empty(&self) -> bool {
        self.vertices.is_empty() || self.indices.is_empty()
    }

    /// Check if GPU buffers need updating
    pub fn is_dirty(&self) -> bool {
        self.dirty
    }

    /// Mark as dirty (needs GPU upload)
    pub fn invalidate(&mut self) {
        self.dirty = true;
    }

    /// Upload surface data to GPU
    pub fn prepare_render(&mut self, context: &RenderContext) {
        if !self.dirty {
            return;
        }

        let rep = self.cached_rep.get_or_insert_with(MeshRep::new);

        // Clone vertices and indices for set_mesh
        let vertices = self.vertices.clone();
        let indices = self.indices.clone();
        rep.set_mesh(vertices, indices);
        rep.upload(context.device(), context.queue());

        self.dirty = false;
    }

    /// Render the surface
    pub fn render<'a>(&'a self, render_pass: &mut wgpu::RenderPass<'a>, context: &'a RenderContext) {
        if !self.state.enabled || self.is_empty() {
            return;
        }

        if let Some(ref rep) = self.cached_rep {
            if !rep.is_empty() {
                let pipeline = context.mesh_pipeline(pymol_render::pipeline::BlendMode::Opaque);
                render_pass.set_pipeline(&pipeline);
                render_pass.set_bind_group(0, context.uniform_bind_group(), &[]);
                rep.render(render_pass);
            }
        }
    }

    /// Compute the bounding box of the surface
    fn compute_extent(&self) -> Option<(Vec3, Vec3)> {
        if self.vertices.is_empty() {
            return None;
        }

        let mut min = Vec3::new(f32::MAX, f32::MAX, f32::MAX);
        let mut max = Vec3::new(f32::MIN, f32::MIN, f32::MIN);

        for vertex in &self.vertices {
            min.x = min.x.min(vertex.position[0]);
            min.y = min.y.min(vertex.position[1]);
            min.z = min.z.min(vertex.position[2]);
            max.x = max.x.max(vertex.position[0]);
            max.y = max.y.max(vertex.position[1]);
            max.z = max.z.max(vertex.position[2]);
        }

        Some((min, max))
    }
}

impl Object for SurfaceObject {
    fn name(&self) -> &str {
        &self.name
    }

    fn object_type(&self) -> ObjectType {
        ObjectType::Surface
    }

    fn state(&self) -> &ObjectState {
        &self.state
    }

    fn state_mut(&mut self) -> &mut ObjectState {
        &mut self.state
    }

    fn extent(&self) -> Option<(Vec3, Vec3)> {
        self.compute_extent()
    }

    fn n_states(&self) -> usize {
        1
    }

    fn current_state(&self) -> usize {
        0
    }

    fn set_current_state(&mut self, _state: usize) -> bool {
        false
    }

    fn set_name(&mut self, name: String) {
        self.name = name;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_test_triangle() -> (Vec<MeshVertex>, Vec<u32>) {
        let vertices = vec![
            MeshVertex {
                position: [0.0, 0.0, 0.0],
                normal: [0.0, 0.0, 1.0],
                color: [1.0, 0.0, 0.0, 1.0],
            },
            MeshVertex {
                position: [1.0, 0.0, 0.0],
                normal: [0.0, 0.0, 1.0],
                color: [0.0, 1.0, 0.0, 1.0],
            },
            MeshVertex {
                position: [0.0, 1.0, 0.0],
                normal: [0.0, 0.0, 1.0],
                color: [0.0, 0.0, 1.0, 1.0],
            },
        ];
        let indices = vec![0, 1, 2];
        (vertices, indices)
    }

    #[test]
    fn test_surface_creation() {
        let surface = SurfaceObject::new("test_surface");
        assert_eq!(surface.name(), "test_surface");
        assert!(surface.is_empty());
        assert!(surface.is_enabled());
        assert_eq!(surface.object_type(), ObjectType::Surface);
    }

    #[test]
    fn test_surface_from_mesh() {
        let (vertices, indices) = create_test_triangle();
        let surface = SurfaceObject::from_mesh("test", vertices, indices);

        assert_eq!(surface.vertex_count(), 3);
        assert_eq!(surface.triangle_count(), 1);
        assert!(!surface.is_empty());
    }

    #[test]
    fn test_surface_from_molecule_surface() {
        let (vertices, indices) = create_test_triangle();
        let surface = SurfaceObject::from_molecule_surface(
            "test",
            vertices,
            indices,
            "my_protein",
            SurfaceType::SolventAccessible,
        );

        assert_eq!(surface.source_molecule(), Some("my_protein"));
        assert_eq!(surface.surface_type(), SurfaceType::SolventAccessible);
    }

    #[test]
    fn test_surface_extent() {
        let (vertices, indices) = create_test_triangle();
        let surface = SurfaceObject::from_mesh("test", vertices, indices);

        let (min, max) = surface.extent().expect("Should have extent");
        assert_eq!(min.x, 0.0);
        assert_eq!(min.y, 0.0);
        assert_eq!(max.x, 1.0);
        assert_eq!(max.y, 1.0);
    }

    #[test]
    fn test_surface_clear() {
        let (vertices, indices) = create_test_triangle();
        let mut surface = SurfaceObject::from_mesh("test", vertices, indices);

        assert!(!surface.is_empty());
        surface.clear();
        assert!(surface.is_empty());
    }

    #[test]
    fn test_surface_dirty_flag() {
        let mut surface = SurfaceObject::new("test");
        assert!(surface.is_dirty());

        // Setting mesh marks dirty
        let (vertices, indices) = create_test_triangle();
        surface.set_mesh(vertices, indices);
        assert!(surface.is_dirty());
    }

    #[test]
    fn test_empty_surface_extent() {
        let surface = SurfaceObject::new("test");
        assert!(surface.extent().is_none());
    }
}
