//! GPU vdW signed distance field (first half of SES).
//!
//! Mirror of [`super::surface_density`] for SES: same acceleration buffers,
//! but the output is r32float (signed distance) instead of rgba16float
//! (density + gradient).
//!
//! ## Bind groups
//!
//! - **Group 0** — `SceneStore` (with dynamic offset). Atoms + coords
//!   come from here.
//! - **Group 1** — per-rep build state:
//!   | Binding | Resource                                    |
//!   |--------:|---------------------------------------------|
//!   |       0 | `SdfParams` uniform                         |
//!   |       1 | `sdf_3d` r32float storage texture           |
//!   |       2 | `owner_3d` r32uint storage texture          |
//!   |       3 | cell offsets storage                        |
//!   |       4 | cell-local atom index storage               |
//!   |       5 | `SurfaceAccelParams` uniform                |

use bytemuck::{Pod, Zeroable};
use wgpu::util::DeviceExt;

use super::surface_density::{SurfaceAccelParams, WORKGROUP};
use crate::scene_store::SceneStoreLayout;
use crate::shader_source;

/// Mirrors `SdfParams` in `surface_vdw_sdf.wgsl`. 32 B.
#[repr(C)]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
pub struct SdfParams {
    pub bbox_min: [f32; 3],
    pub voxel_size: f32,
    pub dims: [u32; 3],
    pub _pad_a: u32,
}

impl SdfParams {
    pub const SIZE: u64 = std::mem::size_of::<Self>() as u64;
}

pub struct SurfaceVdwSdfCompute {
    pub pipeline: wgpu::ComputePipeline,
    pub build_layout: wgpu::BindGroupLayout,
}

impl SurfaceVdwSdfCompute {
    pub fn new(device: &wgpu::Device, scene_layout: &SceneStoreLayout) -> Self {
        let build_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("patinae.surface_vdw_sdf.build_layout"),
            entries: &[
                wgpu::BindGroupLayoutEntry {
                    binding: 0,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Uniform,
                        has_dynamic_offset: false,
                        min_binding_size: wgpu::BufferSize::new(SdfParams::SIZE),
                    },
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 1,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::StorageTexture {
                        access: wgpu::StorageTextureAccess::WriteOnly,
                        format: wgpu::TextureFormat::R32Float,
                        view_dimension: wgpu::TextureViewDimension::D3,
                    },
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 2,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::StorageTexture {
                        access: wgpu::StorageTextureAccess::WriteOnly,
                        format: wgpu::TextureFormat::R32Uint,
                        view_dimension: wgpu::TextureViewDimension::D3,
                    },
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 3,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Storage { read_only: true },
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 4,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Storage { read_only: true },
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 5,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Uniform,
                        has_dynamic_offset: false,
                        min_binding_size: wgpu::BufferSize::new(SurfaceAccelParams::SIZE),
                    },
                    count: None,
                },
            ],
        });

        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("patinae.surface_vdw_sdf.pl"),
            bind_group_layouts: &[Some(&scene_layout.bind_group_layout), Some(&build_layout)],
            immediate_size: 0,
        });

        let module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("patinae.surface_vdw_sdf.wgsl"),
            source: wgpu::ShaderSource::Wgsl(shader_source::SURFACE_VDW_SDF_WGSL.into()),
        });

        let pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("patinae.surface_vdw_sdf.pipeline"),
            layout: Some(&pipeline_layout),
            module: &module,
            entry_point: Some("cs_main"),
            compilation_options: wgpu::PipelineCompilationOptions::default(),
            cache: None,
        });

        Self {
            pipeline,
            build_layout,
        }
    }

    pub fn build_inputs(
        &self,
        device: &wgpu::Device,
        args: SurfaceVdwSdfBuildArgs<'_>,
    ) -> SurfaceVdwSdfInputs {
        let params_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("patinae.surface_vdw_sdf.params"),
            contents: bytemuck::bytes_of(&args.params),
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        });
        let accel_params_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("patinae.surface_vdw_sdf.accel_params"),
            contents: bytemuck::bytes_of(&args.accel_params),
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        });
        let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("patinae.surface_vdw_sdf.bg"),
            layout: &self.build_layout,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: params_buf.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: wgpu::BindingResource::TextureView(args.sdf_view),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: wgpu::BindingResource::TextureView(args.owner_view),
                },
                wgpu::BindGroupEntry {
                    binding: 3,
                    resource: args.cell_offsets.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 4,
                    resource: args.cell_indices.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 5,
                    resource: accel_params_buf.as_entire_binding(),
                },
            ],
        });
        SurfaceVdwSdfInputs {
            params_buf,
            accel_params_buf,
            bind_group,
            dims: args.params.dims,
        }
    }

    pub fn dispatch(
        &self,
        encoder: &mut wgpu::CommandEncoder,
        scene_bg: &wgpu::BindGroup,
        obj_dynamic_offset: u32,
        inputs: &SurfaceVdwSdfInputs,
    ) {
        let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
            label: Some("patinae.surface_vdw_sdf.pass"),
            timestamp_writes: None,
        });
        pass.set_pipeline(&self.pipeline);
        pass.set_bind_group(0, scene_bg, &[obj_dynamic_offset]);
        pass.set_bind_group(1, &inputs.bind_group, &[]);
        pass.dispatch_workgroups(
            inputs.dims[0].div_ceil(WORKGROUP),
            inputs.dims[1].div_ceil(WORKGROUP),
            inputs.dims[2].div_ceil(WORKGROUP),
        );
    }
}

/// Borrow-bundle of inputs for [`SurfaceVdwSdfCompute::build_inputs`].
/// See [`crate::compute::surface_density::SurfaceDensityBuildArgs`] for
/// the rationale (avoids a 7-arg signature).
pub struct SurfaceVdwSdfBuildArgs<'a> {
    pub params: SdfParams,
    pub sdf_view: &'a wgpu::TextureView,
    pub owner_view: &'a wgpu::TextureView,
    pub cell_offsets: &'a wgpu::Buffer,
    pub cell_indices: &'a wgpu::Buffer,
    pub accel_params: SurfaceAccelParams,
}

pub struct SurfaceVdwSdfInputs {
    pub params_buf: wgpu::Buffer,
    pub accel_params_buf: wgpu::Buffer,
    pub bind_group: wgpu::BindGroup,
    pub dims: [u32; 3],
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sdf_params_layout() {
        assert_eq!(std::mem::size_of::<SdfParams>(), 32);
    }
}
