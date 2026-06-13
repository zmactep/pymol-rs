//! GPU marching cubes.
//!
//! One compute dispatch per grid; each thread classifies one cube of 8
//! corners and emits triangles via `atomicAdd` on a shared vertex counter.
//! Output is a flat `StdVertex` storage buffer the size of the next stage's
//! draw call (no index buffer — every triangle has three unique vertices,
//! matching the CPU oracle).
//!
//! ## Bind group layout (group 0)
//!
//! | Binding | Resource                                         | Access  |
//! |--------:|--------------------------------------------------|---------|
//! |       0 | `McParams` uniform                               | uniform |
//! |       1 | `field_3d` rgba16float storage texture           | read    |
//! |       2 | `owner_3d` r32uint storage texture               | read    |
//! |       3 | `mc_tables` storage buffer (packed Bourke tables)| read    |
//! |       4 | `vertex_count` `atomic<u32>` storage buffer      | rw      |
//! |       5 | `vertices` `StdVertex[]` storage buffer          | rw      |
//!
//! `vertex_count` MUST be reset to 0 by the host before each dispatch
//! ([`SurfaceMcCompute::reset_counter`]). After the dispatch the host reads
//! the counter back to know how many vertices were emitted.

use bytemuck::{Pod, Zeroable};
use wgpu::util::DeviceExt;

use crate::representations::surface::mc_tables::{pack_for_gpu, GpuMcTablesLayout};
use crate::shader_source;

pub const WORKGROUP: u32 = 4;

/// Mirrors `McParams` in `surface_mc.wgsl`. 80 B.
#[repr(C)]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
pub struct McParams {
    pub bbox_min: [f32; 3],
    pub voxel_size: f32,
    pub dims: [u32; 3],
    pub iso: f32,
    pub max_vertices: u32,
    /// 0 ⇒ inside means `value >= iso` (SAS density).
    /// 1 ⇒ inside means `value <= iso` (SES SDF).
    pub invert_inside: u32,
    /// 0 ⇒ emit triangles (3 vertices per `n_tri` in TriangleList layout).
    /// 1 ⇒ emit line segments (6 vertices per `n_tri` in LineList layout)
    ///     using the per-triangle edge permutation `[0,1, 1,2, 2,0]`.
    /// Used by mesh wireframe rendering — same density grid + same MC
    /// classification, different output topology.
    pub emit_lines: u32,
    pub _pad1: u32,
    /// World-space cube-center bounds that this dispatch is allowed to emit.
    /// Untiled surfaces use the full grid domain.
    pub emit_core_min: [f32; 3],
    pub _pad2: f32,
    pub emit_core_max: [f32; 3],
    pub _pad3: f32,
}

impl McParams {
    pub const SIZE: u64 = std::mem::size_of::<Self>() as u64;
}

pub struct SurfaceMcCompute {
    pub pipeline: wgpu::ComputePipeline,
    pub bind_group_layout: wgpu::BindGroupLayout,
    /// Packed Bourke tables. Created once, reused across dispatches.
    pub tables_buf: wgpu::Buffer,
}

impl SurfaceMcCompute {
    pub fn new(device: &wgpu::Device) -> Self {
        let bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("patinae.surface_mc.bgl"),
            entries: &[
                // 0: McParams uniform
                wgpu::BindGroupLayoutEntry {
                    binding: 0,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Uniform,
                        has_dynamic_offset: false,
                        min_binding_size: wgpu::BufferSize::new(McParams::SIZE),
                    },
                    count: None,
                },
                // 1: field_3d (rgba16float, read)
                wgpu::BindGroupLayoutEntry {
                    binding: 1,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::StorageTexture {
                        access: wgpu::StorageTextureAccess::ReadOnly,
                        format: wgpu::TextureFormat::Rgba16Float,
                        view_dimension: wgpu::TextureViewDimension::D3,
                    },
                    count: None,
                },
                // 2: owner_3d (r32uint, read)
                wgpu::BindGroupLayoutEntry {
                    binding: 2,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::StorageTexture {
                        access: wgpu::StorageTextureAccess::ReadOnly,
                        format: wgpu::TextureFormat::R32Uint,
                        view_dimension: wgpu::TextureViewDimension::D3,
                    },
                    count: None,
                },
                // 3: mc_tables (storage<u32>, read)
                wgpu::BindGroupLayoutEntry {
                    binding: 3,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Storage { read_only: true },
                        has_dynamic_offset: false,
                        min_binding_size: wgpu::BufferSize::new(4),
                    },
                    count: None,
                },
                // 4: vertex_count (atomic<u32>)
                wgpu::BindGroupLayoutEntry {
                    binding: 4,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Storage { read_only: false },
                        has_dynamic_offset: false,
                        min_binding_size: wgpu::BufferSize::new(4),
                    },
                    count: None,
                },
                // 5: vertices (StdVertex[])
                wgpu::BindGroupLayoutEntry {
                    binding: 5,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Storage { read_only: false },
                        has_dynamic_offset: false,
                        min_binding_size: wgpu::BufferSize::new(24),
                    },
                    count: None,
                },
            ],
        });

        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("patinae.surface_mc.pl"),
            bind_group_layouts: &[&bind_group_layout],
            immediate_size: 0,
        });

        let module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("patinae.surface_mc.wgsl"),
            source: wgpu::ShaderSource::Wgsl(
                shader_source::expand(shader_source::SURFACE_MC_WGSL).into(),
            ),
        });

        let pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("patinae.surface_mc.pipeline"),
            layout: Some(&pipeline_layout),
            module: &module,
            entry_point: Some("cs_main"),
            compilation_options: wgpu::PipelineCompilationOptions::default(),
            cache: None,
        });

        let tables = pack_for_gpu();
        debug_assert_eq!(tables.len(), GpuMcTablesLayout::TOTAL_WORDS);
        let tables_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("patinae.surface_mc.tables"),
            contents: bytemuck::cast_slice(&tables),
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
        });

        Self {
            pipeline,
            bind_group_layout,
            tables_buf,
        }
    }

    /// Build all per-dispatch buffers + the bind group. The caller owns
    /// `vertex_buf` and `count_buf` so they can be re-bound elsewhere
    /// (e.g. as a vertex buffer in the translucent pass).
    pub fn build_inputs(
        &self,
        device: &wgpu::Device,
        params: McParams,
        field_view: &wgpu::TextureView,
        owner_view: &wgpu::TextureView,
        vertex_buf: &wgpu::Buffer,
        count_buf: &wgpu::Buffer,
    ) -> SurfaceMcInputs {
        let params_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("patinae.surface_mc.params"),
            contents: bytemuck::bytes_of(&params),
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        });
        let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("patinae.surface_mc.bg"),
            layout: &self.bind_group_layout,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: params_buf.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: wgpu::BindingResource::TextureView(field_view),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: wgpu::BindingResource::TextureView(owner_view),
                },
                wgpu::BindGroupEntry {
                    binding: 3,
                    resource: self.tables_buf.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 4,
                    resource: count_buf.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 5,
                    resource: vertex_buf.as_entire_binding(),
                },
            ],
        });
        SurfaceMcInputs {
            params_buf,
            bind_group,
            dims: params.dims,
        }
    }

    /// Zero the vertex counter. Must be called before each dispatch.
    pub fn reset_counter(&self, queue: &wgpu::Queue, count_buf: &wgpu::Buffer) {
        queue.write_buffer(count_buf, 0, bytemuck::bytes_of(&0u32));
    }

    pub fn dispatch(&self, encoder: &mut wgpu::CommandEncoder, inputs: &SurfaceMcInputs) {
        let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
            label: Some("patinae.surface_mc.pass"),
            timestamp_writes: None,
        });
        pass.set_pipeline(&self.pipeline);
        pass.set_bind_group(0, &inputs.bind_group, &[]);
        // The shader processes (dims-1)³ cubes; round up by WORKGROUP. Cubes
        // beyond `dims-1` early-return, so over-dispatching is safe.
        pass.dispatch_workgroups(
            inputs.dims[0].div_ceil(WORKGROUP),
            inputs.dims[1].div_ceil(WORKGROUP),
            inputs.dims[2].div_ceil(WORKGROUP),
        );
    }
}

pub struct SurfaceMcInputs {
    pub params_buf: wgpu::Buffer,
    pub bind_group: wgpu::BindGroup,
    pub dims: [u32; 3],
}

/// Convenience: allocate a vertex buffer big enough for `max_vertices`
/// `StdVertex` records (24 B each). Reuse across frames as long as the cap
/// is unchanged.
pub fn create_vertex_buffer(device: &wgpu::Device, max_vertices: u32) -> wgpu::Buffer {
    device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("patinae.surface_mc.vertices"),
        size: (max_vertices as u64) * 24,
        usage: wgpu::BufferUsages::STORAGE
            | wgpu::BufferUsages::VERTEX
            | wgpu::BufferUsages::COPY_SRC,
        mapped_at_creation: false,
    })
}

/// Convenience: allocate the atomic counter buffer (4 B). Single u32.
pub fn create_counter_buffer(device: &wgpu::Device) -> wgpu::Buffer {
    device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("patinae.surface_mc.counter"),
        size: 4,
        usage: wgpu::BufferUsages::STORAGE
            | wgpu::BufferUsages::COPY_SRC
            | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mc_params_layout() {
        // 5×16 = 80 B (matches the WGSL struct).
        assert_eq!(std::mem::size_of::<McParams>(), 80);
        assert_eq!(std::mem::align_of::<McParams>(), 4);
    }
}
