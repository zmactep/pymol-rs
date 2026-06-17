//! Portable GPU runtime descriptors for plugin command callbacks.
//!
//! The host owns all `wgpu` objects. Plugins receive only command-scoped
//! opaque handles and serializable descriptors.
//!
//! This module is the public, portable contract between plugin crates and the
//! host runtime. Every descriptor is data-only so it can cross the plugin wire
//! boundary without exposing Rust pointers, `wgpu::Device`, `wgpu::Queue`, or
//! renderer internals.
//!
//! # Ownership And Lifetimes
//!
//! A `GpuHandle` is a lease, not ownership of a GPU object. Handles created by
//! ordinary GPU calls are valid only for the command callback that received
//! them. Handles returned by cached-create calls are still command-scoped
//! leases; only the host-side cached object persists across callbacks.
//!
//! Buffers, textures, texture views, samplers, bind groups, readbacks, and
//! renderer artifact handles are intentionally command-scoped. Cached handles
//! are limited to expensive immutable state such as shader modules, layouts,
//! compute pipelines, and render pipelines.
//!
//! # Batch Ordering
//!
//! `GpuSubmitBatch` commands execute in vector order inside one host-owned
//! command encoder. A nested render pass is one ordered batch command, and the
//! commands inside that pass execute in their own vector order before the host
//! records the next outer batch command.
//!
//! Readbacks are returned in the order their `ReadBuffer` commands appear in
//! the batch after all prior writes, dispatches, draws, and copies have
//! completed.
//!
//! # Validation And Errors
//!
//! The host validates handle kind, generation, plugin ownership, declared usage
//! flags, descriptor compatibility, texture copy ranges, render attachment
//! formats, and portable v1 feature limits before calling `wgpu`. Validation
//! failures are reported as portable error strings through the plugin API.
//!
//! Debug labels are for diagnostics only. Persistent cache fingerprints ignore
//! labels and include only semantic creation fields plus the host cache layout
//! version and device generation.
//!
//! # Renderer Artifacts
//!
//! `RenderArtifactSnapshotDescriptor` describes renderer-owned buffers without
//! exposing renderer types. Artifact handles are command-scoped, snapshots may
//! become stale after renderer sync or device recreation, and plugins must use
//! `layout_version` plus role-specific stride/count metadata before consuming a
//! buffer.

use serde::{Deserialize, Serialize};

/// Command-scoped opaque GPU resource handle.
///
/// The host resolves the numeric `id` only inside the current command runtime.
/// `kind` and `generation` let the host reject stale handles and handles used
/// as the wrong resource type before touching `wgpu`.
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
    TextureView,
    Sampler,
    ShaderModule,
    BindGroupLayout,
    PipelineLayout,
    BindGroup,
    ComputePipeline,
    RenderPipeline,
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
    pub buffer_binding_array: bool,
    pub storage_resource_binding_array: bool,
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

/// Portable texture usage flags understood by the plugin GPU runtime.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct GpuTextureUsage {
    pub bits: u32,
}

impl GpuTextureUsage {
    pub const COPY_SRC: Self = Self { bits: 1 << 0 };
    pub const COPY_DST: Self = Self { bits: 1 << 1 };
    pub const TEXTURE_BINDING: Self = Self { bits: 1 << 2 };
    pub const STORAGE_BINDING: Self = Self { bits: 1 << 3 };
    pub const RENDER_ATTACHMENT: Self = Self { bits: 1 << 4 };

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

/// Extent of a texture copy or allocation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct GpuExtent3d {
    pub width: u32,
    pub height: u32,
    pub depth_or_array_layers: u32,
}

/// Origin of a texture copy region.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct GpuOrigin3d {
    pub x: u32,
    pub y: u32,
    pub z: u32,
}

/// Texture allocation dimension.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum GpuTextureDimension {
    D1,
    D2,
    D3,
}

/// Texture view dimension used in bindings.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum GpuTextureViewDimension {
    D1,
    D2,
    D2Array,
    Cube,
    CubeArray,
    D3,
}

/// Texture aspect addressed by a copy or view.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum GpuTextureAspect {
    All,
    StencilOnly,
    DepthOnly,
}

/// Portable subset of texture formats exposed to plugins.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum GpuTextureFormat {
    R8Unorm,
    R8Uint,
    R8Sint,
    R16Uint,
    R16Sint,
    R16Float,
    R32Uint,
    R32Sint,
    R32Float,
    Rg8Unorm,
    Rg8Uint,
    Rg8Sint,
    Rg16Float,
    Rg32Float,
    Rgba8Unorm,
    Rgba8UnormSrgb,
    Rgba8Uint,
    Rgba8Sint,
    Bgra8Unorm,
    Bgra8UnormSrgb,
    Rgba16Float,
    Rgba32Float,
    Depth32Float,
}

/// Descriptor for a plugin-created GPU texture.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct GpuTextureDescriptor {
    pub label: Option<String>,
    pub size: GpuExtent3d,
    pub mip_level_count: u32,
    pub sample_count: u32,
    pub dimension: GpuTextureDimension,
    pub format: GpuTextureFormat,
    pub usage: GpuTextureUsage,
}

/// Descriptor for a plugin-created texture view.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct GpuTextureViewDescriptor {
    pub label: Option<String>,
    pub texture: GpuHandle,
    pub format: Option<GpuTextureFormat>,
    pub dimension: Option<GpuTextureViewDimension>,
    pub usage: Option<GpuTextureUsage>,
    pub aspect: GpuTextureAspect,
    pub base_mip_level: u32,
    pub mip_level_count: Option<u32>,
    pub base_array_layer: u32,
    pub array_layer_count: Option<u32>,
}

/// Addressing mode for a sampler axis.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum GpuAddressMode {
    ClampToEdge,
    Repeat,
    MirrorRepeat,
    ClampToBorder,
}

/// Filtering mode for sampler lookups.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum GpuFilterMode {
    Nearest,
    Linear,
}

/// Comparison function for sampler depth comparisons.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum GpuCompareFunction {
    Never,
    Less,
    Equal,
    LessEqual,
    Greater,
    NotEqual,
    GreaterEqual,
    Always,
}

/// Descriptor for a plugin-created sampler.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct GpuSamplerDescriptor {
    pub label: Option<String>,
    pub address_mode_u: GpuAddressMode,
    pub address_mode_v: GpuAddressMode,
    pub address_mode_w: GpuAddressMode,
    pub mag_filter: GpuFilterMode,
    pub min_filter: GpuFilterMode,
    pub mipmap_filter: GpuFilterMode,
    pub lod_min_clamp: f32,
    pub lod_max_clamp: f32,
    pub compare: Option<GpuCompareFunction>,
    pub anisotropy_clamp: u16,
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
    Sampler(GpuSamplerBindingType),
    Texture {
        sample_type: GpuTextureSampleType,
        view_dimension: GpuTextureViewDimension,
        multisampled: bool,
    },
    StorageTexture {
        access: GpuStorageTextureAccess,
        format: GpuTextureFormat,
        view_dimension: GpuTextureViewDimension,
    },
}

/// Sampler binding type for a bind-group layout entry.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum GpuSamplerBindingType {
    Filtering,
    NonFiltering,
    Comparison,
}

/// Texture sample type for a bind-group layout entry.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum GpuTextureSampleType {
    Float { filterable: bool },
    Depth,
    Sint,
    Uint,
}

/// Storage texture access for a bind-group layout entry.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum GpuStorageTextureAccess {
    WriteOnly,
    ReadOnly,
    ReadWrite,
    Atomic,
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

/// Vertex buffer stepping mode for a render pipeline.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum GpuVertexStepMode {
    Vertex,
    Instance,
}

/// Portable subset of vertex formats for plugin render pipelines.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum GpuVertexFormat {
    Uint8,
    Uint8x2,
    Uint8x4,
    Sint8,
    Sint8x2,
    Sint8x4,
    Unorm8,
    Unorm8x2,
    Unorm8x4,
    Snorm8,
    Snorm8x2,
    Snorm8x4,
    Uint16,
    Uint16x2,
    Uint16x4,
    Sint16,
    Sint16x2,
    Sint16x4,
    Unorm16,
    Unorm16x2,
    Unorm16x4,
    Snorm16,
    Snorm16x2,
    Snorm16x4,
    Float16,
    Float16x2,
    Float16x4,
    Float32,
    Float32x2,
    Float32x3,
    Float32x4,
    Uint32,
    Uint32x2,
    Uint32x3,
    Uint32x4,
    Sint32,
    Sint32x2,
    Sint32x3,
    Sint32x4,
    Unorm10_10_10_2,
    Unorm8x4Bgra,
}

/// One vertex attribute in a render pipeline vertex buffer layout.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct GpuVertexAttribute {
    pub format: GpuVertexFormat,
    pub offset: u64,
    pub shader_location: u32,
}

/// One vertex buffer layout for a render pipeline.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct GpuVertexBufferLayout {
    pub array_stride: u64,
    pub step_mode: GpuVertexStepMode,
    pub attributes: Vec<GpuVertexAttribute>,
}

/// Vertex shader state for a render pipeline.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct GpuVertexState {
    pub module: GpuHandle,
    pub entry_point: String,
    pub buffers: Vec<GpuVertexBufferLayout>,
}

/// Primitive topology for a render pipeline.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum GpuPrimitiveTopology {
    PointList,
    LineList,
    LineStrip,
    TriangleList,
    TriangleStrip,
}

/// Index format for indexed render commands.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum GpuIndexFormat {
    Uint16,
    Uint32,
}

/// Front-face winding for render pipeline culling.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum GpuFrontFace {
    Ccw,
    Cw,
}

/// Cullable face for render pipeline culling.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum GpuFace {
    Front,
    Back,
}

/// Primitive assembly state for a render pipeline.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct GpuPrimitiveState {
    pub topology: GpuPrimitiveTopology,
    pub strip_index_format: Option<GpuIndexFormat>,
    pub front_face: GpuFrontFace,
    pub cull_mode: Option<GpuFace>,
}

/// Color write mask for a render pipeline color target.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct GpuColorWriteMask {
    pub bits: u32,
}

impl GpuColorWriteMask {
    pub const RED: Self = Self { bits: 1 << 0 };
    pub const GREEN: Self = Self { bits: 1 << 1 };
    pub const BLUE: Self = Self { bits: 1 << 2 };
    pub const ALPHA: Self = Self { bits: 1 << 3 };
    pub const COLOR: Self = Self {
        bits: Self::RED.bits | Self::GREEN.bits | Self::BLUE.bits,
    };
    pub const ALL: Self = Self {
        bits: Self::COLOR.bits | Self::ALPHA.bits,
    };

    pub const fn union(self, other: Self) -> Self {
        Self {
            bits: self.bits | other.bits,
        }
    }

    pub const fn contains(self, other: Self) -> bool {
        (self.bits & other.bits) == other.bits
    }
}

/// Preset blend state for a render pipeline color target.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum GpuBlendState {
    Replace,
    AlphaBlending,
    PremultipliedAlphaBlending,
}

/// Color target state for a fragment shader output.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct GpuColorTargetState {
    pub format: GpuTextureFormat,
    pub blend: Option<GpuBlendState>,
    pub write_mask: GpuColorWriteMask,
}

/// Fragment shader state for a render pipeline.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct GpuFragmentState {
    pub module: GpuHandle,
    pub entry_point: String,
    pub targets: Vec<Option<GpuColorTargetState>>,
}

/// Depth state for a render pipeline.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct GpuDepthStencilState {
    pub format: GpuTextureFormat,
    pub depth_write_enabled: bool,
    pub depth_compare: GpuCompareFunction,
}

/// Multisample state for a render pipeline.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct GpuMultisampleState {
    pub count: u32,
    pub mask: u64,
    pub alpha_to_coverage_enabled: bool,
}

/// Descriptor for a plugin-created render pipeline.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct GpuRenderPipelineDescriptor {
    pub label: Option<String>,
    pub layout: GpuHandle,
    pub vertex: GpuVertexState,
    pub primitive: GpuPrimitiveState,
    pub depth_stencil: Option<GpuDepthStencilState>,
    pub multisample: GpuMultisampleState,
    pub fragment: Option<GpuFragmentState>,
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
///
/// A cache hit or miss changes only host-side creation cost. The returned
/// `handle` follows the same command-scoped lifetime rules as any other
/// `GpuHandle`, while the cached object remains host-owned for later commands
/// from the same plugin and device generation.
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
    TextureView(GpuHandle),
    Sampler(GpuHandle),
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

/// Buffer layout for a texture copy.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct GpuTexelCopyBufferLayout {
    pub offset: u64,
    pub bytes_per_row: Option<u32>,
    pub rows_per_image: Option<u32>,
}

/// Texture subresource addressed by a texture copy.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct GpuTexelCopyTextureInfo {
    pub texture: GpuHandle,
    pub mip_level: u32,
    pub origin: GpuOrigin3d,
    pub aspect: GpuTextureAspect,
}

/// Clear color used by a render pass color attachment.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct GpuColor {
    pub r: f64,
    pub g: f64,
    pub b: f64,
    pub a: f64,
}

/// Color attachment load operation.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum GpuColorLoadOp {
    Load,
    Clear(GpuColor),
}

/// Depth attachment load operation.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum GpuDepthLoadOp {
    Load,
    Clear(f32),
}

/// Attachment store operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum GpuStoreOp {
    Store,
    Discard,
}

/// Color attachment operations for a render pass.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct GpuColorOperations {
    pub load: GpuColorLoadOp,
    pub store: GpuStoreOp,
}

/// Depth attachment operations for a render pass.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct GpuDepthOperations {
    pub load: GpuDepthLoadOp,
    pub store: GpuStoreOp,
}

/// Color attachment for an offscreen render pass.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct GpuRenderPassColorAttachment {
    pub view: GpuHandle,
    pub ops: GpuColorOperations,
}

/// Depth attachment for an offscreen render pass.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct GpuRenderPassDepthStencilAttachment {
    pub view: GpuHandle,
    pub depth_ops: Option<GpuDepthOperations>,
}

/// Descriptor for an offscreen plugin render pass.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct GpuRenderPassDescriptor {
    pub label: Option<String>,
    pub color_attachments: Vec<Option<GpuRenderPassColorAttachment>>,
    pub depth_stencil_attachment: Option<GpuRenderPassDepthStencilAttachment>,
}

/// One command recorded inside a plugin render pass.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum GpuRenderPassCommand {
    SetPipeline {
        pipeline: GpuHandle,
    },
    SetBindGroup {
        index: u32,
        bind_group: GpuHandle,
    },
    SetVertexBuffer {
        slot: u32,
        buffer: GpuHandle,
        offset: u64,
        size: Option<u64>,
    },
    SetIndexBuffer {
        buffer: GpuHandle,
        format: GpuIndexFormat,
        offset: u64,
        size: Option<u64>,
    },
    SetViewport {
        x: f32,
        y: f32,
        width: f32,
        height: f32,
        min_depth: f32,
        max_depth: f32,
    },
    SetScissorRect {
        x: u32,
        y: u32,
        width: u32,
        height: u32,
    },
    Draw {
        vertex_start: u32,
        vertex_count: u32,
        instance_start: u32,
        instance_count: u32,
    },
    DrawIndexed {
        index_start: u32,
        index_count: u32,
        base_vertex: i32,
        instance_start: u32,
        instance_count: u32,
    },
    DrawIndirect {
        buffer: GpuHandle,
        offset: u64,
    },
    DrawIndexedIndirect {
        buffer: GpuHandle,
        offset: u64,
    },
}

/// One command recorded into a host-owned GPU batch.
///
/// The outer command vector is an ordering contract. A render pass command is
/// recorded as one nested pass, which keeps compute, render, copy, and readback
/// work ordered without exposing an open-ended begin/end pass stream.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
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
    /// Dispatch one compute pipeline from a GPU-written dispatch buffer.
    DispatchComputeIndirect {
        pipeline: GpuHandle,
        bind_groups: Vec<GpuHandle>,
        indirect_buffer: GpuHandle,
        indirect_offset: u64,
    },
    /// Copy bytes between buffers inside the batch command encoder.
    CopyBufferToBuffer {
        source: GpuHandle,
        source_offset: u64,
        destination: GpuHandle,
        destination_offset: u64,
        size: u64,
    },
    /// Copy texels from a buffer into a texture.
    CopyBufferToTexture {
        source: GpuHandle,
        source_layout: GpuTexelCopyBufferLayout,
        destination: GpuTexelCopyTextureInfo,
        size: GpuExtent3d,
    },
    /// Copy texels from a texture into a buffer.
    CopyTextureToBuffer {
        source: GpuTexelCopyTextureInfo,
        destination: GpuHandle,
        destination_layout: GpuTexelCopyBufferLayout,
        size: GpuExtent3d,
    },
    /// Copy texels between textures.
    CopyTextureToTexture {
        source: GpuTexelCopyTextureInfo,
        destination: GpuTexelCopyTextureInfo,
        size: GpuExtent3d,
    },
    /// Copy bytes to a temporary readback buffer after previous commands.
    ReadBuffer {
        buffer: GpuHandle,
        offset: u64,
        size: u64,
    },
    /// Record one offscreen render pass at this batch position.
    RenderPass {
        descriptor: GpuRenderPassDescriptor,
        commands: Vec<GpuRenderPassCommand>,
    },
    /// Set the profiling scope for subsequent batch commands.
    SetProfileScope { name: String },
}

/// Ordered GPU command batch submitted by a plugin.
///
/// The host records all commands into a fresh command encoder, submits once,
/// waits for completion when requested or when readbacks are requested, and
/// returns readback payloads in command order. Batch submission never transfers
/// ownership of a plugin handle to the plugin.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct GpuSubmitBatch {
    pub label: Option<String>,
    pub wait_for_completion: bool,
    pub commands: Vec<GpuBatchCommand>,
}

/// Result of a submitted GPU command batch.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct GpuBatchResult {
    pub readbacks: Vec<Vec<u8>>,
}

/// Semantic role of a renderer-owned GPU buffer.
///
/// Roles identify layout, not ownership. Plugins must combine the role with
/// `stride`, `element_count`, and the snapshot `layout_version` before reading
/// from a renderer artifact buffer.
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
///
/// `element_count` is the direct draw/dispatch count when no count or indirect
/// buffer is present. `max_element_count` is the backing capacity. If `count`
/// or `indirect` is present, that buffer is authoritative for GPU-side culling
/// or draw count and `element_count` may be zero.
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
///
/// A snapshot is a point-in-time list of renderer buffers already synchronized
/// for drawing. Artifact handles expire with the command runtime; a later
/// renderer sync, culling pass, or device change can make old metadata stale.
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
