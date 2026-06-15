//! Portable GPU runtime descriptors for plugin command callbacks.
//!
//! The host owns all `wgpu` objects. Plugins receive only command-scoped
//! opaque handles and serializable descriptors.

use serde::{Deserialize, Serialize};

/// Command-scoped opaque GPU resource handle.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct GpuHandle {
    pub id: u64,
    pub kind: GpuHandleKind,
    pub generation: u64,
}

/// Kind tag for an opaque GPU handle.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum GpuHandleKind {
    Buffer,
    Texture,
    ShaderModule,
    BindGroupLayout,
    PipelineLayout,
    BindGroup,
    ComputePipeline,
}

/// Portable subset of device limits useful to GPU plugins.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct GpuDeviceLimits {
    pub max_buffer_size: u64,
    pub max_storage_buffer_binding_size: u64,
    pub max_compute_workgroups_per_dimension: u32,
    pub max_compute_invocations_per_workgroup: u32,
    pub max_compute_workgroup_size_x: u32,
    pub max_compute_workgroup_size_y: u32,
    pub max_compute_workgroup_size_z: u32,
}

/// Portable buffer usage flags understood by the plugin GPU runtime.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct GpuBufferUsage {
    pub bits: u32,
}

impl GpuBufferUsage {
    pub const MAP_READ: Self = Self { bits: 1 << 0 };
    pub const COPY_SRC: Self = Self { bits: 1 << 1 };
    pub const COPY_DST: Self = Self { bits: 1 << 2 };
    pub const INDEX: Self = Self { bits: 1 << 3 };
    pub const VERTEX: Self = Self { bits: 1 << 4 };
    pub const UNIFORM: Self = Self { bits: 1 << 5 };
    pub const STORAGE: Self = Self { bits: 1 << 6 };
    pub const INDIRECT: Self = Self { bits: 1 << 7 };
    pub const QUERY_RESOLVE: Self = Self { bits: 1 << 8 };

    pub const fn union(self, other: Self) -> Self {
        Self {
            bits: self.bits | other.bits,
        }
    }

    pub const fn contains(self, other: Self) -> bool {
        (self.bits & other.bits) == other.bits
    }
}

/// Portable shader stage flags understood by the plugin GPU runtime.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct GpuShaderStages {
    pub bits: u32,
}

impl GpuShaderStages {
    pub const VERTEX: Self = Self { bits: 1 << 0 };
    pub const FRAGMENT: Self = Self { bits: 1 << 1 };
    pub const COMPUTE: Self = Self { bits: 1 << 2 };

    pub const fn union(self, other: Self) -> Self {
        Self {
            bits: self.bits | other.bits,
        }
    }
}

/// Descriptor for a plugin-created GPU buffer.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct GpuBufferDescriptor {
    pub label: Option<String>,
    pub size: u64,
    pub usage: GpuBufferUsage,
}

/// Buffer binding type for a bind-group layout entry.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum GpuBufferBindingType {
    Uniform,
    StorageReadOnly,
    StorageReadWrite,
}

/// Binding type for a bind-group layout entry.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum GpuBindingType {
    Buffer {
        ty: GpuBufferBindingType,
        has_dynamic_offset: bool,
        min_binding_size: Option<u64>,
    },
}

/// One bind-group layout entry.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct GpuBindGroupLayoutEntry {
    pub binding: u32,
    pub visibility: GpuShaderStages,
    pub ty: GpuBindingType,
}

/// Descriptor for a plugin-created bind-group layout.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct GpuBindGroupLayoutDescriptor {
    pub label: Option<String>,
    pub entries: Vec<GpuBindGroupLayoutEntry>,
}

/// Descriptor for a plugin-created pipeline layout.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct GpuPipelineLayoutDescriptor {
    pub label: Option<String>,
    pub bind_group_layouts: Vec<GpuHandle>,
}

/// Descriptor for a plugin-created WGSL shader module.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct GpuShaderModuleDescriptor {
    pub label: Option<String>,
    pub wgsl: String,
}

/// Descriptor for a plugin-created compute pipeline.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct GpuComputePipelineDescriptor {
    pub label: Option<String>,
    pub layout: GpuHandle,
    pub module: GpuHandle,
    pub entry_point: String,
}

/// Cache outcome for a persistent GPU resource request.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum GpuCacheStatus {
    /// The host reused a previously-created resource.
    Hit,
    /// The host created and stored a new resource.
    Miss,
}

/// Command-scoped lease for a persistent cached GPU resource.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct GpuCachedHandle {
    pub handle: GpuHandle,
    pub status: GpuCacheStatus,
}

/// Persistent GPU cache counters for the active plugin.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct GpuCacheStats {
    pub hits: u64,
    pub misses: u64,
    pub entries: u64,
}

/// Buffer resource bound into a bind group.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct GpuBufferBinding {
    pub buffer: GpuHandle,
    pub offset: u64,
    pub size: Option<u64>,
}

/// Resource bound into a bind group.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum GpuBindingResource {
    Buffer(GpuBufferBinding),
}

/// One bind-group entry.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct GpuBindGroupEntry {
    pub binding: u32,
    pub resource: GpuBindingResource,
}

/// Descriptor for a plugin-created bind group.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct GpuBindGroupDescriptor {
    pub label: Option<String>,
    pub layout: GpuHandle,
    pub entries: Vec<GpuBindGroupEntry>,
}

/// One command recorded into a host-owned GPU batch.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum GpuBatchCommand {
    /// Write bytes into a buffer before subsequent batch commands run.
    WriteBuffer {
        buffer: GpuHandle,
        offset: u64,
        data: Vec<u8>,
    },
    /// Dispatch one compute pipeline.
    DispatchCompute {
        pipeline: GpuHandle,
        bind_groups: Vec<GpuHandle>,
        workgroups: [u32; 3],
    },
    /// Copy bytes between buffers inside the batch command encoder.
    CopyBufferToBuffer {
        source: GpuHandle,
        source_offset: u64,
        destination: GpuHandle,
        destination_offset: u64,
        size: u64,
    },
    /// Copy bytes to a temporary readback buffer after previous commands.
    ReadBuffer {
        buffer: GpuHandle,
        offset: u64,
        size: u64,
    },
}

/// Ordered GPU command batch submitted by a plugin.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct GpuSubmitBatch {
    pub label: Option<String>,
    pub commands: Vec<GpuBatchCommand>,
}

/// Result of a submitted GPU command batch.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct GpuBatchResult {
    pub readbacks: Vec<Vec<u8>>,
}

/// Semantic role of a renderer-owned GPU buffer.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
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

/// Displayed representation kind for renderer artifact snapshots.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum RenderArtifactRepKind {
    Sphere,
    Stick,
    Line,
    Cartoon,
    Ribbon,
    Surface,
    Mesh,
    Dot,
    Ellipsoid,
    Other,
}

/// Primitive topology represented by renderer artifact buffers.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum RenderArtifactPrimitiveTopology {
    SphereInstances,
    CylinderInstances,
    LineInstances,
    TriangleList,
    LineList,
}

/// Renderer-owned buffer exposed as a command-scoped GPU handle.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct RenderArtifactBufferDescriptor {
    pub handle: GpuHandle,
    pub role: RenderArtifactBufferRole,
    pub size: u64,
    pub stride: u64,
    pub element_count: u64,
}

/// One displayed representation in a render artifact snapshot.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RenderArtifactRepDescriptor {
    pub object_id: u32,
    pub rep_kind: RenderArtifactRepKind,
    pub topology: RenderArtifactPrimitiveTopology,
    pub geometry: GpuHandle,
    pub count: Option<GpuHandle>,
    pub indirect: Option<GpuHandle>,
    pub element_count: u64,
    pub max_element_count: u64,
    pub atom_offset: u32,
    pub atom_count: u32,
    pub material_rgba: [f32; 4],
    pub transparency: f32,
}

/// Command-scoped snapshot of renderer GPU artifacts.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RenderArtifactSnapshotDescriptor {
    pub snapshot_id: u64,
    pub layout_version: u32,
    pub scene_generation: u64,
    pub scene_bounds_min: [f32; 3],
    pub scene_bounds_max: [f32; 3],
    pub cull_pass_initialized: bool,
    pub device_limits: GpuDeviceLimits,
    pub buffers: Vec<RenderArtifactBufferDescriptor>,
    pub reps: Vec<RenderArtifactRepDescriptor>,
}
