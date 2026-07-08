//! Viewport sphere-count compute pass for camera-aware sphere LOD.
//!
//! The pass reads the same `SceneStore` atom/coord data as the sphere build
//! shader and atomically counts sphere atoms inside the current camera
//! frustum. `SphereRep` uses the asynchronous 4-byte result to decide whether
//! a large scene should keep its sampled instance buffer or rebuild all
//! spheres for a zoomed-in fragment.

use std::num::NonZeroU64;

use crate::compute::cull::CullParams;
use crate::scene_store::SceneStoreLayout;
use crate::shader_source;

pub const WORKGROUP: u32 = 64;

pub struct SphereLodCountPipeline {
    pub pipeline: wgpu::ComputePipeline,
    pub count_layout: wgpu::BindGroupLayout,
}

impl SphereLodCountPipeline {
    pub fn new(device: &wgpu::Device, scene_layout: &SceneStoreLayout) -> Self {
        let count_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("patinae.sphere_lod_count.layout"),
            entries: &[
                wgpu::BindGroupLayoutEntry {
                    binding: 0,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Uniform,
                        has_dynamic_offset: false,
                        min_binding_size: NonZeroU64::new(CullParams::SIZE),
                    },
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 1,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Storage { read_only: false },
                        has_dynamic_offset: false,
                        min_binding_size: NonZeroU64::new(4),
                    },
                    count: None,
                },
            ],
        });
        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("patinae.sphere_lod_count.pipeline_layout"),
            bind_group_layouts: &[Some(&scene_layout.bind_group_layout), Some(&count_layout)],
            immediate_size: 0,
        });
        let module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("patinae.sphere_lod_count.shader"),
            source: wgpu::ShaderSource::Wgsl(
                shader_source::expand(shader_source::SPHERE_LOD_COUNT_WGSL).into(),
            ),
        });
        let pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("patinae.sphere_lod_count.pipeline"),
            layout: Some(&pipeline_layout),
            module: &module,
            entry_point: Some("cs_main"),
            compilation_options: wgpu::PipelineCompilationOptions::default(),
            cache: None,
        });
        Self {
            pipeline,
            count_layout,
        }
    }

    pub fn make_bind_group(
        &self,
        device: &wgpu::Device,
        params_buf: &wgpu::Buffer,
        count_buf: &wgpu::Buffer,
    ) -> wgpu::BindGroup {
        device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("patinae.sphere_lod_count.bg"),
            layout: &self.count_layout,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: params_buf.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: count_buf.as_entire_binding(),
                },
            ],
        })
    }

    pub fn dispatch(
        &self,
        encoder: &mut wgpu::CommandEncoder,
        scene_bg: &wgpu::BindGroup,
        obj_dynamic_offset: u32,
        count_bg: &wgpu::BindGroup,
        atom_count: u32,
    ) {
        if atom_count == 0 {
            return;
        }
        let groups = atom_count.div_ceil(WORKGROUP);
        let (wg_x, wg_y) = super::split_1d_dispatch(groups);
        let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
            label: Some("patinae.sphere_lod_count.dispatch"),
            timestamp_writes: None,
        });
        pass.set_pipeline(&self.pipeline);
        pass.set_bind_group(0, scene_bg, &[obj_dynamic_offset]);
        pass.set_bind_group(1, count_bg, &[]);
        pass.dispatch_workgroups(wg_x, wg_y, 1);
    }
}
