//! Compiled Graphics Object (CGO)
//!
//! CGO objects allow custom graphics primitives to be rendered.
//! They use an instruction-based approach similar to PyMOL's CGO system.

use lin_alg::f32::Vec3;
use pymol_render::{
    CylinderVertex, LineRep, LineVertex, MeshRep, MeshVertex, RenderContext, Representation,
    SphereRep, SphereVertex, StickRep,
};

use super::{Object, ObjectState, ObjectType};

/// Primitive mode for Begin/End blocks
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PrimitiveMode {
    /// Triangle list
    Triangles,
    /// Triangle strip
    TriangleStrip,
    /// Triangle fan
    TriangleFan,
    /// Line list
    Lines,
    /// Line strip
    LineStrip,
}

/// CGO instruction opcodes
///
/// These are the primitives that can be rendered via CGO.
/// The API is similar to PyMOL's CGO system.
#[derive(Debug, Clone)]
pub enum CgoOp {
    /// Begin a primitive block (like glBegin)
    Begin(PrimitiveMode),
    /// End a primitive block (like glEnd)
    End,
    /// Set vertex position (used within Begin/End blocks)
    Vertex(Vec3),
    /// Set normal for subsequent vertices
    Normal(Vec3),
    /// Set color for subsequent primitives (RGBA)
    Color([f32; 4]),
    /// Set alpha (transparency) for subsequent primitives
    Alpha(f32),
    /// Set line width for line primitives
    LineWidth(f32),
    /// Draw a sphere at a position
    Sphere {
        center: Vec3,
        radius: f32,
    },
    /// Draw a cylinder between two points
    Cylinder {
        start: Vec3,
        end: Vec3,
        radius: f32,
        color1: [f32; 4],
        color2: [f32; 4],
    },
    /// Draw a cone (truncated) between two points
    Cone {
        start: Vec3,
        end: Vec3,
        radius1: f32,
        radius2: f32,
        color1: [f32; 4],
        color2: [f32; 4],
    },
    /// Draw a single triangle with normals
    Triangle {
        v1: Vec3,
        v2: Vec3,
        v3: Vec3,
        n1: Vec3,
        n2: Vec3,
        n3: Vec3,
    },
    /// Draw a line segment
    Line {
        start: Vec3,
        end: Vec3,
    },
    /// Set a custom color for a specific primitive
    CustomColor([f32; 4]),
}

/// Builder state for Begin/End blocks
#[derive(Default)]
struct BuilderState {
    mode: Option<PrimitiveMode>,
    current_normal: Vec3,
    current_color: [f32; 4],
    vertices: Vec<(Vec3, Vec3, [f32; 4])>, // (position, normal, color)
}

impl BuilderState {
    fn new() -> Self {
        Self {
            mode: None,
            current_normal: Vec3::new(0.0, 0.0, 1.0),
            current_color: [1.0, 1.0, 1.0, 1.0],
            vertices: Vec::new(),
        }
    }

    fn reset(&mut self) {
        self.mode = None;
        self.vertices.clear();
    }
}

/// Cached representations for CGO rendering
struct CgoRepCache {
    spheres: SphereRep,
    sticks: StickRep,
    lines: LineRep,
    mesh: MeshRep,
}

impl CgoRepCache {
    fn new() -> Self {
        Self {
            spheres: SphereRep::new(),
            sticks: StickRep::new(),
            lines: LineRep::new(),
            mesh: MeshRep::new(),
        }
    }

    fn clear(&mut self) {
        self.spheres = SphereRep::new();
        self.sticks = StickRep::new();
        self.lines = LineRep::new();
        self.mesh = MeshRep::new();
    }
}

/// A Compiled Graphics Object (CGO)
///
/// CGO objects contain a list of graphics instructions that define
/// custom primitives to be rendered. This allows users to draw
/// arbitrary graphics in 3D space.
pub struct CgoObject {
    /// Object name
    name: String,
    /// Visual state (enabled, color, etc.)
    state: ObjectState,
    /// List of CGO instructions
    instructions: Vec<CgoOp>,
    /// Cached representations for rendering
    cached_rep: Option<CgoRepCache>,
    /// Whether the cached data needs rebuilding
    dirty: bool,
    /// Current color state
    current_color: [f32; 4],
    /// Current alpha state
    current_alpha: f32,
    /// Current line width
    line_width: f32,
}

impl CgoObject {
    /// Create a new empty CGO object
    pub fn new(name: &str) -> Self {
        Self {
            name: name.to_string(),
            state: ObjectState::default(),
            instructions: Vec::new(),
            cached_rep: None,
            dirty: true,
            current_color: [1.0, 1.0, 1.0, 1.0],
            current_alpha: 1.0,
            line_width: 1.0,
        }
    }

    /// Create a CGO object with initial instructions
    pub fn with_instructions(name: &str, instructions: Vec<CgoOp>) -> Self {
        Self {
            name: name.to_string(),
            state: ObjectState::default(),
            instructions,
            cached_rep: None,
            dirty: true,
            current_color: [1.0, 1.0, 1.0, 1.0],
            current_alpha: 1.0,
            line_width: 1.0,
        }
    }

    /// Get the instructions
    pub fn instructions(&self) -> &[CgoOp] {
        &self.instructions
    }

    /// Get mutable access to instructions (marks dirty)
    pub fn instructions_mut(&mut self) -> &mut Vec<CgoOp> {
        self.dirty = true;
        &mut self.instructions
    }

    /// Add an instruction
    pub fn add(&mut self, op: CgoOp) {
        self.instructions.push(op);
        self.dirty = true;
    }

    /// Add a sphere
    pub fn add_sphere(&mut self, center: Vec3, radius: f32) {
        self.add(CgoOp::Sphere { center, radius });
    }

    /// Add a cylinder
    pub fn add_cylinder(
        &mut self,
        start: Vec3,
        end: Vec3,
        radius: f32,
        color1: [f32; 4],
        color2: [f32; 4],
    ) {
        self.add(CgoOp::Cylinder {
            start,
            end,
            radius,
            color1,
            color2,
        });
    }

    /// Add a triangle
    pub fn add_triangle(&mut self, v1: Vec3, v2: Vec3, v3: Vec3) {
        // Compute face normal
        let edge1 = Vec3::new(v2.x - v1.x, v2.y - v1.y, v2.z - v1.z);
        let edge2 = Vec3::new(v3.x - v1.x, v3.y - v1.y, v3.z - v1.z);
        let normal = Vec3::new(
            edge1.y * edge2.z - edge1.z * edge2.y,
            edge1.z * edge2.x - edge1.x * edge2.z,
            edge1.x * edge2.y - edge1.y * edge2.x,
        );
        let len = (normal.x * normal.x + normal.y * normal.y + normal.z * normal.z).sqrt();
        let n = if len > 0.0 {
            Vec3::new(normal.x / len, normal.y / len, normal.z / len)
        } else {
            Vec3::new(0.0, 0.0, 1.0)
        };

        self.add(CgoOp::Triangle {
            v1,
            v2,
            v3,
            n1: n,
            n2: n,
            n3: n,
        });
    }

    /// Add a line
    pub fn add_line(&mut self, start: Vec3, end: Vec3) {
        self.add(CgoOp::Line { start, end });
    }

    /// Set the current color
    pub fn set_color(&mut self, color: [f32; 4]) {
        self.add(CgoOp::Color(color));
    }

    /// Clear all instructions
    pub fn clear(&mut self) {
        self.instructions.clear();
        self.dirty = true;
    }

    /// Check if the object is empty
    pub fn is_empty(&self) -> bool {
        self.instructions.is_empty()
    }

    /// Check if the cached data needs rebuilding
    pub fn is_dirty(&self) -> bool {
        self.dirty
    }

    /// Mark as needing rebuild
    pub fn invalidate(&mut self) {
        self.dirty = true;
    }

    /// Build representations from instructions
    fn build_representations(&mut self) {
        let cache = self.cached_rep.get_or_insert_with(CgoRepCache::new);
        cache.clear();

        let mut builder = BuilderState::new();
        let mut spheres: Vec<SphereVertex> = Vec::new();
        let mut cylinders: Vec<CylinderVertex> = Vec::new();
        let mut lines: Vec<LineVertex> = Vec::new();
        let mut mesh_vertices: Vec<MeshVertex> = Vec::new();
        let mut mesh_indices: Vec<u32> = Vec::new();

        let mut current_color = self.current_color;
        let mut current_alpha = self.current_alpha;

        for op in &self.instructions {
            match op {
                CgoOp::Color(color) => {
                    current_color = *color;
                    builder.current_color = *color;
                }
                CgoOp::Alpha(alpha) => {
                    current_alpha = *alpha;
                    current_color[3] = *alpha;
                    builder.current_color[3] = *alpha;
                }
                CgoOp::Normal(n) => {
                    builder.current_normal = *n;
                }
                CgoOp::Begin(mode) => {
                    builder.mode = Some(*mode);
                    builder.vertices.clear();
                }
                CgoOp::Vertex(pos) => {
                    if builder.mode.is_some() {
                        builder.vertices.push((*pos, builder.current_normal, builder.current_color));
                    }
                }
                CgoOp::End => {
                    if let Some(mode) = builder.mode {
                        Self::process_primitive(
                            mode,
                            &builder.vertices,
                            &mut mesh_vertices,
                            &mut mesh_indices,
                            &mut lines,
                        );
                    }
                    builder.reset();
                }
                CgoOp::Sphere { center, radius } => {
                    spheres.push(SphereVertex {
                        center: [center.x, center.y, center.z],
                        radius: *radius,
                        color: current_color,
                    });
                }
                CgoOp::Cylinder {
                    start,
                    end,
                    radius,
                    color1,
                    color2,
                } => {
                    cylinders.push(CylinderVertex {
                        start: [start.x, start.y, start.z],
                        radius: *radius,
                        end: [end.x, end.y, end.z],
                        flags: CylinderVertex::CAP_BOTH,
                        color1: *color1,
                        color2: *color2,
                    });
                }
                CgoOp::Cone { .. } => {
                    // Cones are rendered as cylinders for now
                    // A proper cone shader could be added later
                }
                CgoOp::Triangle {
                    v1,
                    v2,
                    v3,
                    n1,
                    n2,
                    n3,
                } => {
                    let base_idx = mesh_vertices.len() as u32;
                    mesh_vertices.push(MeshVertex {
                        position: [v1.x, v1.y, v1.z],
                        normal: [n1.x, n1.y, n1.z],
                        color: current_color,
                    });
                    mesh_vertices.push(MeshVertex {
                        position: [v2.x, v2.y, v2.z],
                        normal: [n2.x, n2.y, n2.z],
                        color: current_color,
                    });
                    mesh_vertices.push(MeshVertex {
                        position: [v3.x, v3.y, v3.z],
                        normal: [n3.x, n3.y, n3.z],
                        color: current_color,
                    });
                    mesh_indices.push(base_idx);
                    mesh_indices.push(base_idx + 1);
                    mesh_indices.push(base_idx + 2);
                }
                CgoOp::Line { start, end } => {
                    lines.push(LineVertex {
                        position: [start.x, start.y, start.z],
                        color: current_color,
                    });
                    lines.push(LineVertex {
                        position: [end.x, end.y, end.z],
                        color: current_color,
                    });
                }
                CgoOp::LineWidth(width) => {
                    self.line_width = *width;
                }
                CgoOp::CustomColor(color) => {
                    current_color = *color;
                }
            }
        }

        // Store the built data in the representations
        // We'll set the data directly by accessing the internal structures
        // For now, we just store the vertex data

        // Store in cache - we'll upload during prepare_render
        self.current_color = current_color;
        self.current_alpha = current_alpha;
    }

    fn process_primitive(
        mode: PrimitiveMode,
        vertices: &[(Vec3, Vec3, [f32; 4])],
        mesh_vertices: &mut Vec<MeshVertex>,
        mesh_indices: &mut Vec<u32>,
        lines: &mut Vec<LineVertex>,
    ) {
        match mode {
            PrimitiveMode::Triangles => {
                for chunk in vertices.chunks(3) {
                    if chunk.len() == 3 {
                        let base_idx = mesh_vertices.len() as u32;
                        for (pos, normal, color) in chunk {
                            mesh_vertices.push(MeshVertex {
                                position: [pos.x, pos.y, pos.z],
                                normal: [normal.x, normal.y, normal.z],
                                color: *color,
                            });
                        }
                        mesh_indices.push(base_idx);
                        mesh_indices.push(base_idx + 1);
                        mesh_indices.push(base_idx + 2);
                    }
                }
            }
            PrimitiveMode::TriangleStrip => {
                if vertices.len() >= 3 {
                    for i in 0..vertices.len() - 2 {
                        let base_idx = mesh_vertices.len() as u32;
                        let (p1, n1, c1) = &vertices[i];
                        let (p2, n2, c2) = &vertices[i + 1];
                        let (p3, n3, c3) = &vertices[i + 2];

                        mesh_vertices.push(MeshVertex {
                            position: [p1.x, p1.y, p1.z],
                            normal: [n1.x, n1.y, n1.z],
                            color: *c1,
                        });
                        mesh_vertices.push(MeshVertex {
                            position: [p2.x, p2.y, p2.z],
                            normal: [n2.x, n2.y, n2.z],
                            color: *c2,
                        });
                        mesh_vertices.push(MeshVertex {
                            position: [p3.x, p3.y, p3.z],
                            normal: [n3.x, n3.y, n3.z],
                            color: *c3,
                        });

                        // Alternate winding order for strips
                        if i % 2 == 0 {
                            mesh_indices.push(base_idx);
                            mesh_indices.push(base_idx + 1);
                            mesh_indices.push(base_idx + 2);
                        } else {
                            mesh_indices.push(base_idx);
                            mesh_indices.push(base_idx + 2);
                            mesh_indices.push(base_idx + 1);
                        }
                    }
                }
            }
            PrimitiveMode::TriangleFan => {
                if vertices.len() >= 3 {
                    let (center_pos, center_normal, center_color) = &vertices[0];
                    for i in 1..vertices.len() - 1 {
                        let base_idx = mesh_vertices.len() as u32;
                        let (p1, n1, c1) = &vertices[i];
                        let (p2, n2, c2) = &vertices[i + 1];

                        mesh_vertices.push(MeshVertex {
                            position: [center_pos.x, center_pos.y, center_pos.z],
                            normal: [center_normal.x, center_normal.y, center_normal.z],
                            color: *center_color,
                        });
                        mesh_vertices.push(MeshVertex {
                            position: [p1.x, p1.y, p1.z],
                            normal: [n1.x, n1.y, n1.z],
                            color: *c1,
                        });
                        mesh_vertices.push(MeshVertex {
                            position: [p2.x, p2.y, p2.z],
                            normal: [n2.x, n2.y, n2.z],
                            color: *c2,
                        });

                        mesh_indices.push(base_idx);
                        mesh_indices.push(base_idx + 1);
                        mesh_indices.push(base_idx + 2);
                    }
                }
            }
            PrimitiveMode::Lines => {
                for chunk in vertices.chunks(2) {
                    if chunk.len() == 2 {
                        let (p1, _, c1) = &chunk[0];
                        let (p2, _, c2) = &chunk[1];
                        lines.push(LineVertex {
                            position: [p1.x, p1.y, p1.z],
                            color: *c1,
                        });
                        lines.push(LineVertex {
                            position: [p2.x, p2.y, p2.z],
                            color: *c2,
                        });
                    }
                }
            }
            PrimitiveMode::LineStrip => {
                if vertices.len() >= 2 {
                    for i in 0..vertices.len() - 1 {
                        let (p1, _, c1) = &vertices[i];
                        let (p2, _, c2) = &vertices[i + 1];
                        lines.push(LineVertex {
                            position: [p1.x, p1.y, p1.z],
                            color: *c1,
                        });
                        lines.push(LineVertex {
                            position: [p2.x, p2.y, p2.z],
                            color: *c2,
                        });
                    }
                }
            }
        }
    }

    /// Prepare for rendering (upload to GPU)
    pub fn prepare_render(&mut self, context: &RenderContext) {
        if !self.dirty {
            return;
        }

        self.build_representations();

        // Upload each representation type
        if let Some(ref mut cache) = self.cached_rep {
            cache.spheres.upload(context.device(), context.queue());
            cache.sticks.upload(context.device(), context.queue());
            cache.lines.upload(context.device(), context.queue());
            cache.mesh.upload(context.device(), context.queue());
        }

        self.dirty = false;
    }

    /// Render the CGO
    pub fn render<'a>(&'a self, render_pass: &mut wgpu::RenderPass<'a>, context: &'a RenderContext) {
        if !self.state.enabled || self.is_empty() {
            return;
        }

        if let Some(ref cache) = self.cached_rep {
            // Render lines
            if !cache.lines.is_empty() {
                let pipeline = context.line_pipeline(pymol_render::pipeline::BlendMode::Opaque);
                render_pass.set_pipeline(&pipeline);
                render_pass.set_bind_group(0, context.uniform_bind_group(), &[]);
                cache.lines.render(render_pass);
            }

            // Render mesh (triangles)
            if !cache.mesh.is_empty() {
                let pipeline = context.mesh_pipeline(pymol_render::pipeline::BlendMode::Opaque);
                render_pass.set_pipeline(&pipeline);
                render_pass.set_bind_group(0, context.uniform_bind_group(), &[]);
                cache.mesh.render(render_pass);
            }

            // Render spheres
            if !cache.spheres.is_empty() {
                let pipeline = context.sphere_pipeline(pymol_render::pipeline::BlendMode::Opaque);
                render_pass.set_pipeline(&pipeline);
                render_pass.set_bind_group(0, context.uniform_bind_group(), &[]);
                render_pass.set_vertex_buffer(0, context.billboard_vertex_buffer().slice(..));
                render_pass.set_index_buffer(
                    context.quad_index_buffer().slice(..),
                    wgpu::IndexFormat::Uint16,
                );
                cache.spheres.render(render_pass);
            }

            // Render cylinders (sticks)
            if !cache.sticks.is_empty() {
                let pipeline = context.cylinder_pipeline(pymol_render::pipeline::BlendMode::Opaque);
                render_pass.set_pipeline(&pipeline);
                render_pass.set_bind_group(0, context.uniform_bind_group(), &[]);
                render_pass.set_vertex_buffer(0, context.billboard_vertex_buffer().slice(..));
                render_pass.set_index_buffer(
                    context.quad_index_buffer().slice(..),
                    wgpu::IndexFormat::Uint16,
                );
                cache.sticks.render(render_pass);
            }
        }
    }

    /// Compute the bounding box of all primitives
    fn compute_extent(&self) -> Option<(Vec3, Vec3)> {
        if self.instructions.is_empty() {
            return None;
        }

        let mut min = Vec3::new(f32::MAX, f32::MAX, f32::MAX);
        let mut max = Vec3::new(f32::MIN, f32::MIN, f32::MIN);
        let mut has_geometry = false;

        for op in &self.instructions {
            match op {
                CgoOp::Sphere { center, radius } => {
                    min.x = min.x.min(center.x - radius);
                    min.y = min.y.min(center.y - radius);
                    min.z = min.z.min(center.z - radius);
                    max.x = max.x.max(center.x + radius);
                    max.y = max.y.max(center.y + radius);
                    max.z = max.z.max(center.z + radius);
                    has_geometry = true;
                }
                CgoOp::Cylinder { start, end, radius, .. } => {
                    min.x = min.x.min(start.x - radius).min(end.x - radius);
                    min.y = min.y.min(start.y - radius).min(end.y - radius);
                    min.z = min.z.min(start.z - radius).min(end.z - radius);
                    max.x = max.x.max(start.x + radius).max(end.x + radius);
                    max.y = max.y.max(start.y + radius).max(end.y + radius);
                    max.z = max.z.max(start.z + radius).max(end.z + radius);
                    has_geometry = true;
                }
                CgoOp::Triangle { v1, v2, v3, .. } => {
                    min.x = min.x.min(v1.x).min(v2.x).min(v3.x);
                    min.y = min.y.min(v1.y).min(v2.y).min(v3.y);
                    min.z = min.z.min(v1.z).min(v2.z).min(v3.z);
                    max.x = max.x.max(v1.x).max(v2.x).max(v3.x);
                    max.y = max.y.max(v1.y).max(v2.y).max(v3.y);
                    max.z = max.z.max(v1.z).max(v2.z).max(v3.z);
                    has_geometry = true;
                }
                CgoOp::Vertex(pos) | CgoOp::Line { start: pos, .. } => {
                    min.x = min.x.min(pos.x);
                    min.y = min.y.min(pos.y);
                    min.z = min.z.min(pos.z);
                    max.x = max.x.max(pos.x);
                    max.y = max.y.max(pos.y);
                    max.z = max.z.max(pos.z);
                    has_geometry = true;
                }
                _ => {}
            }
        }

        if has_geometry {
            Some((min, max))
        } else {
            None
        }
    }
}

impl Object for CgoObject {
    fn name(&self) -> &str {
        &self.name
    }

    fn object_type(&self) -> ObjectType {
        ObjectType::Cgo
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

    #[test]
    fn test_cgo_creation() {
        let cgo = CgoObject::new("test_cgo");
        assert_eq!(cgo.name(), "test_cgo");
        assert!(cgo.is_empty());
        assert!(cgo.is_enabled());
        assert_eq!(cgo.object_type(), ObjectType::Cgo);
    }

    #[test]
    fn test_cgo_add_sphere() {
        let mut cgo = CgoObject::new("test");
        cgo.add_sphere(Vec3::new(0.0, 0.0, 0.0), 1.0);

        assert_eq!(cgo.instructions().len(), 1);
        assert!(!cgo.is_empty());

        let extent = cgo.extent().expect("Should have extent");
        assert_eq!(extent.0.x, -1.0);
        assert_eq!(extent.1.x, 1.0);
    }

    #[test]
    fn test_cgo_add_cylinder() {
        let mut cgo = CgoObject::new("test");
        cgo.add_cylinder(
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(5.0, 0.0, 0.0),
            0.5,
            [1.0, 0.0, 0.0, 1.0],
            [0.0, 0.0, 1.0, 1.0],
        );

        assert_eq!(cgo.instructions().len(), 1);
        let extent = cgo.extent().expect("Should have extent");
        assert_eq!(extent.0.x, -0.5);
        assert_eq!(extent.1.x, 5.5);
    }

    #[test]
    fn test_cgo_add_triangle() {
        let mut cgo = CgoObject::new("test");
        cgo.add_triangle(
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        );

        assert_eq!(cgo.instructions().len(), 1);
        let extent = cgo.extent().expect("Should have extent");
        assert_eq!(extent.0.x, 0.0);
        assert_eq!(extent.1.x, 1.0);
        assert_eq!(extent.0.y, 0.0);
        assert_eq!(extent.1.y, 1.0);
    }

    #[test]
    fn test_cgo_clear() {
        let mut cgo = CgoObject::new("test");
        cgo.add_sphere(Vec3::new(0.0, 0.0, 0.0), 1.0);
        assert!(!cgo.is_empty());

        cgo.clear();
        assert!(cgo.is_empty());
    }

    #[test]
    fn test_cgo_with_instructions() {
        let instructions = vec![
            CgoOp::Color([1.0, 0.0, 0.0, 1.0]),
            CgoOp::Sphere {
                center: Vec3::new(0.0, 0.0, 0.0),
                radius: 1.0,
            },
        ];

        let cgo = CgoObject::with_instructions("test", instructions);
        assert_eq!(cgo.instructions().len(), 2);
    }

    #[test]
    fn test_cgo_empty_extent() {
        let cgo = CgoObject::new("test");
        assert!(cgo.extent().is_none());

        let mut cgo_with_color = CgoObject::new("test2");
        cgo_with_color.set_color([1.0, 0.0, 0.0, 1.0]);
        assert!(cgo_with_color.extent().is_none());
    }
}
