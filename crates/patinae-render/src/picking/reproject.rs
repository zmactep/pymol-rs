//! Picking-texture reprojection compute pass.
//!
//! Warps the previous frame's `picking` texture (Rg32Uint @ half viewport)
//! into the current frame's view using the previous frame's depth +
//! the two `view_proj` matrices. Two-pass kernel
//! (`shaders/compute/picking_reproject.wgsl`):
//!
//! 1. `cs_min_depth` — each source pixel `atomicMin`'s its reprojected
//!    depth into `best_depth[target_idx]`, claiming the nearest source
//!    pixel for that target.
//! 2. `cs_write`     — each source pixel re-reads `best_depth`; if it
//!    matches, it `textureStore`s its packed picking value to the
//!    current frame's picking texture.
//!
//! Pixels whose reprojected NDC falls outside the current frustum are
//! dropped silently — they would have been "no atom" anyway. A scene_dirty
//! re-record corrects any accumulated drift.

use std::num::NonZeroU64;

use bytemuck::{Pod, Zeroable};
use wgpu::util::DeviceExt;

use crate::context::RenderContext;
use crate::frame::PICKING_FORMAT;
use crate::memory::{buffer_usage, GpuMemoryUsage};
use crate::shader_source::{self, PICKING_REPROJECT_WGSL};

/// Mirror of the WGSL `ReprojectParams` uniform.
#[repr(C)]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
pub struct ReprojectParams {
    /// Inverse of the view-proj matrix that rendered the `picking_prev`
    /// texture (the last full re-record).
    pub view_proj_prev_inv: [[f32; 4]; 4],
    /// Current frame's view-proj — where we reproject into.
    pub view_proj_cur: [[f32; 4]; 4],
}

impl ReprojectParams {
    pub const SIZE: u64 = std::mem::size_of::<Self>() as u64;
}

/// Compute pipelines + bind-group layout + best-depth storage for the
/// reprojection pass. Owned by `RenderState`.
pub struct PickingReproject {
    pub bind_group_layout: wgpu::BindGroupLayout,
    pub pipeline_min_depth: wgpu::ComputePipeline,
    pub pipeline_write: wgpu::ComputePipeline,
    pub params_buffer: wgpu::Buffer,
    /// `array<atomic<u32>>` sized `picking_w * picking_h`. Reallocated on
    /// resize via [`Self::ensure_best_depth_capacity`].
    pub best_depth_buffer: wgpu::Buffer,
    best_depth_capacity_words: u32,
}

impl PickingReproject {
    pub fn new(ctx: &RenderContext, initial_pick_w: u32, initial_pick_h: u32) -> Self {
        let device = &ctx.device;

        let module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("patinae.picking.reproject.shader"),
            source: wgpu::ShaderSource::Wgsl(shader_source::expand(PICKING_REPROJECT_WGSL).into()),
        });

        let bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("patinae.picking.reproject.layout"),
            entries: &[
                // 0: ReprojectParams uniform
                wgpu::BindGroupLayoutEntry {
                    binding: 0,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Uniform,
                        has_dynamic_offset: false,
                        min_binding_size: NonZeroU64::new(ReprojectParams::SIZE),
                    },
                    count: None,
                },
                // 1: picking_in (prev frame)
                wgpu::BindGroupLayoutEntry {
                    binding: 1,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Texture {
                        sample_type: wgpu::TextureSampleType::Uint,
                        view_dimension: wgpu::TextureViewDimension::D2,
                        multisampled: false,
                    },
                    count: None,
                },
                // 2: depth_in (prev frame's picking_depth, read as Float)
                wgpu::BindGroupLayoutEntry {
                    binding: 2,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Texture {
                        sample_type: wgpu::TextureSampleType::Float { filterable: false },
                        view_dimension: wgpu::TextureViewDimension::D2,
                        multisampled: false,
                    },
                    count: None,
                },
                // 3: best_depth storage<read_write>
                wgpu::BindGroupLayoutEntry {
                    binding: 3,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Storage { read_only: false },
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
                // 4: picking_out (current frame, write-only storage texture)
                wgpu::BindGroupLayoutEntry {
                    binding: 4,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::StorageTexture {
                        access: wgpu::StorageTextureAccess::WriteOnly,
                        format: PICKING_FORMAT,
                        view_dimension: wgpu::TextureViewDimension::D2,
                    },
                    count: None,
                },
            ],
        });

        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("patinae.picking.reproject.pipeline_layout"),
            bind_group_layouts: &[&bind_group_layout],
            immediate_size: 0,
        });

        let pipeline_min_depth = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("patinae.picking.reproject.cs_min_depth"),
            layout: Some(&pipeline_layout),
            module: &module,
            entry_point: Some("cs_min_depth"),
            compilation_options: Default::default(),
            cache: None,
        });

        let pipeline_write = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("patinae.picking.reproject.cs_write"),
            layout: Some(&pipeline_layout),
            module: &module,
            entry_point: Some("cs_write"),
            compilation_options: Default::default(),
            cache: None,
        });

        let params_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("patinae.picking.reproject.params"),
            contents: bytemuck::bytes_of(&ReprojectParams {
                view_proj_prev_inv: identity4(),
                view_proj_cur: identity4(),
            }),
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        });

        let words = initial_pick_w.max(1) * initial_pick_h.max(1);
        let best_depth_buffer = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("patinae.picking.reproject.best_depth"),
            size: (words as u64) * 4,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        Self {
            bind_group_layout,
            pipeline_min_depth,
            pipeline_write,
            params_buffer,
            best_depth_buffer,
            best_depth_capacity_words: words,
        }
    }

    /// Ensure `best_depth` covers the current picking dimensions. Reallocates
    /// only on growth — callers can hand sticky-larger sizes (resize up,
    /// then down) without thrash.
    pub fn ensure_best_depth_capacity(&mut self, device: &wgpu::Device, pick_w: u32, pick_h: u32) {
        let words = pick_w.max(1) * pick_h.max(1);
        if words <= self.best_depth_capacity_words {
            return;
        }
        self.best_depth_buffer = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("patinae.picking.reproject.best_depth"),
            size: (words as u64) * 4,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        self.best_depth_capacity_words = words;
    }

    /// Write `ReprojectParams` into the GPU uniform.
    pub fn upload_params(&self, queue: &wgpu::Queue, params: ReprojectParams) {
        queue.write_buffer(&self.params_buffer, 0, bytemuck::bytes_of(&params));
    }

    /// Estimated GPU bytes allocated by the reprojection buffers.
    pub fn memory_usage(&self) -> GpuMemoryUsage {
        let mut usage = GpuMemoryUsage::default();
        usage.add(buffer_usage(&self.params_buffer));
        usage.add(buffer_usage(&self.best_depth_buffer));
        usage
    }

    /// Build a bind group from the current frame's resources. Cheap to
    /// rebuild every dispatch — done on the call site rather than cached
    /// to avoid invalidation on resize.
    pub fn make_bind_group(
        &self,
        device: &wgpu::Device,
        picking_prev: &wgpu::TextureView,
        depth_prev: &wgpu::TextureView,
        picking_out: &wgpu::TextureView,
    ) -> wgpu::BindGroup {
        device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("patinae.picking.reproject.bind_group"),
            layout: &self.bind_group_layout,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: self.params_buffer.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: wgpu::BindingResource::TextureView(picking_prev),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: wgpu::BindingResource::TextureView(depth_prev),
                },
                wgpu::BindGroupEntry {
                    binding: 3,
                    resource: self.best_depth_buffer.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 4,
                    resource: wgpu::BindingResource::TextureView(picking_out),
                },
            ],
        })
    }

    /// Record both compute passes into `encoder`. Caller is responsible
    /// for clearing `best_depth_buffer` to `0xFFFFFFFF` (INF) before this
    /// runs — see [`Self::clear_best_depth`].
    pub fn record(
        &self,
        encoder: &mut wgpu::CommandEncoder,
        bind_group: &wgpu::BindGroup,
        pick_w: u32,
        pick_h: u32,
    ) {
        let wg_x = pick_w.div_ceil(8);
        let wg_y = pick_h.div_ceil(8);
        {
            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("patinae.picking.reproject.min_depth"),
                timestamp_writes: None,
            });
            pass.set_pipeline(&self.pipeline_min_depth);
            pass.set_bind_group(0, bind_group, &[]);
            pass.dispatch_workgroups(wg_x, wg_y, 1);
        }
        {
            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("patinae.picking.reproject.write"),
                timestamp_writes: None,
            });
            pass.set_pipeline(&self.pipeline_write);
            pass.set_bind_group(0, bind_group, &[]);
            pass.dispatch_workgroups(wg_x, wg_y, 1);
        }
    }

    /// Reset `best_depth_buffer` to the `INF_DEPTH` sentinel (`0xFFFFFFFF`).
    /// Done via a small staging write — for half-res 256² that's 256 KB,
    /// negligible. Larger viewports might want a clear-via-compute kernel.
    pub fn clear_best_depth(&self, queue: &wgpu::Queue, pick_w: u32, pick_h: u32) {
        let words = (pick_w * pick_h) as usize;
        let bytes: Vec<u8> = vec![0xFF; words * 4];
        queue.write_buffer(&self.best_depth_buffer, 0, &bytes);
    }
}

fn identity4() -> [[f32; 4]; 4] {
    [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ]
}
