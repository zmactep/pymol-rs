//! Renderer-owned GPU artifact descriptions.
//!
//! These types describe buffers that already exist inside `RenderState`.
//! They intentionally do not contain ray-tracing concepts; plugin hosts can
//! turn borrowed buffer refs into opaque ABI handles for any trusted GPU
//! plugin that wants to run compute work next to the renderer.

use crate::picking::RepKind;

/// Layout version for renderer GPU artifacts exposed to plugins.
pub const RENDER_ARTIFACT_LAYOUT_VERSION: u32 = 2;

/// Semantic role of a renderer-owned GPU buffer.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum RenderArtifactBufferRole {
    FrameUniforms,
    SceneAtoms,
    SceneCoords,
    SceneBonds,
    SceneColorLut,
    SceneMaskLut,
    SceneMarkerLut,
    SceneCsrOffsets,
    SceneCsrIndices,
    SceneObjectTable,
    SphereInstances,
    StickInstances,
    LineInstances,
    StdVertices,
    InstanceCount,
    IndirectDraw,
}

/// Primitive topology represented by a render artifact.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum RenderArtifactPrimitiveTopology {
    SphereInstances,
    CylinderInstances,
    LineInstances,
    TriangleList,
    LineList,
}

/// Borrowed renderer GPU buffer plus neutral layout metadata.
#[derive(Debug, Clone, Copy)]
pub struct RenderArtifactBufferRef<'a> {
    pub role: RenderArtifactBufferRole,
    pub buffer: &'a wgpu::Buffer,
    pub size: u64,
    pub stride: u64,
    pub element_count: u64,
}

impl<'a> RenderArtifactBufferRef<'a> {
    pub(crate) fn new(
        role: RenderArtifactBufferRole,
        buffer: &'a wgpu::Buffer,
        stride: u64,
        element_count: u64,
    ) -> Self {
        Self {
            role,
            buffer,
            size: buffer.size(),
            stride,
            element_count,
        }
    }
}

/// One displayed representation backed by renderer GPU artifacts.
#[derive(Debug, Clone, Copy)]
pub struct RenderArtifactRep {
    pub object_id: u32,
    pub rep_kind: RepKind,
    pub topology: RenderArtifactPrimitiveTopology,
    pub geometry_buffer_index: usize,
    pub count_buffer_index: Option<usize>,
    pub indirect_buffer_index: Option<usize>,
    pub element_count: u64,
    pub max_element_count: u64,
    pub atom_offset: u32,
    pub atom_count: u32,
    pub material_rgba: [f32; 4],
    pub transparency: f32,
}

/// Command-scoped snapshot of renderer GPU artifacts.
#[derive(Debug)]
pub struct RenderArtifactSnapshot<'a> {
    pub layout_version: u32,
    pub scene_generation: u64,
    pub scene_bounds_min: [f32; 3],
    pub scene_bounds_max: [f32; 3],
    pub cull_pass_initialized: bool,
    pub buffers: Vec<RenderArtifactBufferRef<'a>>,
    pub reps: Vec<RenderArtifactRep>,
}
