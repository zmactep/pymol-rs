//! GPU SAS density splat.
//!
//! Owns the compute pipeline and bind-group layout for
//! `shaders/compute/surface_density.wgsl`. The shader is one thread per
//! voxel and gathers atoms from a surface-local spatial binning buffer.
//!
//! ## Bind groups
//!
//! - **Group 0** — `SceneStore` (with dynamic offset for per-object
//!   `obj` uniform). Atoms + coords come from here; the kernel reads
//!   `[obj.atom_offset .. obj.atom_offset + obj.atom_count)`.
//! - **Group 1** — per-rep build state:
//!   | Binding | Resource                                    |
//!   |--------:|---------------------------------------------|
//!   |       0 | `DensityParams` uniform                     |
//!   |       1 | `density_3d` rgba16float storage texture    |
//!   |       2 | `owner_3d` r32uint storage texture          |
//!   |       3 | cell offsets storage                        |
//!   |       4 | cell-local atom index storage               |
//!   |       5 | `SurfaceAccelParams` uniform                |

use bytemuck::{Pod, Zeroable};
use wgpu::util::DeviceExt;

use crate::scene_store::SceneStoreLayout;
use crate::shader_source;

/// Workgroup edge — total threads per group is `WORKGROUP * WORKGROUP * WORKGROUP`.
pub const WORKGROUP: u32 = 4;

/// Mirrors `SurfaceAccelParams` in the surface producer WGSL files. 32 B.
///
/// `cell_offsets` is indexed by `x + y * dims.x + z * dims.x * dims.y`
/// and has `cell_count + 1` entries. `cell_indices[offsets[c]..offsets[c+1])`
/// stores object-local atom ids for that cell.
#[repr(C)]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
pub struct SurfaceAccelParams {
    pub origin: [f32; 3],
    pub cell_size: f32,
    pub dims: [u32; 3],
    pub _pad: u32,
}

impl SurfaceAccelParams {
    pub const SIZE: u64 = std::mem::size_of::<Self>() as u64;
}

/// Mirrors `DensityParams` in `surface_density.wgsl`. 48 B.
#[repr(C)]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
pub struct DensityParams {
    pub bbox_min: [f32; 3],
    pub voxel_size: f32,
    pub dims: [u32; 3],
    pub _pad_a: u32,
    pub probe_radius: f32,
    pub alpha: f32,
    pub _pad0: f32,
    pub _pad1: f32,
}

impl DensityParams {
    pub const SIZE: u64 = std::mem::size_of::<Self>() as u64;
}

/// Compute pipeline + cached bind-group layout.
pub struct SurfaceDensityCompute {
    pub pipeline: wgpu::ComputePipeline,
    pub build_layout: wgpu::BindGroupLayout,
}

impl SurfaceDensityCompute {
    pub fn new(device: &wgpu::Device, scene_layout: &SceneStoreLayout) -> Self {
        let build_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("patinae.surface_density.build_layout"),
            entries: &[
                wgpu::BindGroupLayoutEntry {
                    binding: 0,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Uniform,
                        has_dynamic_offset: false,
                        min_binding_size: wgpu::BufferSize::new(DensityParams::SIZE),
                    },
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 1,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::StorageTexture {
                        access: wgpu::StorageTextureAccess::WriteOnly,
                        format: wgpu::TextureFormat::Rgba16Float,
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
            label: Some("patinae.surface_density.pl"),
            bind_group_layouts: &[&scene_layout.bind_group_layout, &build_layout],
            immediate_size: 0,
        });

        let module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("patinae.surface_density.wgsl"),
            source: wgpu::ShaderSource::Wgsl(shader_source::SURFACE_DENSITY_WGSL.into()),
        });

        let pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("patinae.surface_density.pipeline"),
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

    /// Build the per-rep build bind group + uniform buffer. SceneStore
    /// bind group is supplied at dispatch time with the per-object dynamic
    /// offset, so it isn't bundled here.
    pub fn build_inputs(
        &self,
        device: &wgpu::Device,
        args: SurfaceDensityBuildArgs<'_>,
    ) -> SurfaceDensityInputs {
        let params_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("patinae.surface_density.params"),
            contents: bytemuck::bytes_of(&args.params),
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        });
        let accel_params_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("patinae.surface_density.accel_params"),
            contents: bytemuck::bytes_of(&args.accel_params),
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        });
        let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("patinae.surface_density.bg"),
            layout: &self.build_layout,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: params_buf.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: wgpu::BindingResource::TextureView(args.density_view),
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
        SurfaceDensityInputs {
            params_buf,
            accel_params_buf,
            bind_group,
            dims: args.params.dims,
        }
    }

    /// Record a single dispatch covering the full grid.
    pub fn dispatch(
        &self,
        encoder: &mut wgpu::CommandEncoder,
        scene_bg: &wgpu::BindGroup,
        obj_dynamic_offset: u32,
        inputs: &SurfaceDensityInputs,
    ) {
        let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
            label: Some("patinae.surface_density.pass"),
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

/// Borrow-bundle of inputs for [`SurfaceDensityCompute::build_inputs`].
/// Avoids a 7-arg signature; each field is a `Copy` value or a borrow
/// into resources the caller already holds.
pub struct SurfaceDensityBuildArgs<'a> {
    pub params: DensityParams,
    pub density_view: &'a wgpu::TextureView,
    pub owner_view: &'a wgpu::TextureView,
    pub cell_offsets: &'a wgpu::Buffer,
    pub cell_indices: &'a wgpu::Buffer,
    pub accel_params: SurfaceAccelParams,
}

/// Bag of GPU resources for one density dispatch.
pub struct SurfaceDensityInputs {
    pub params_buf: wgpu::Buffer,
    pub accel_params_buf: wgpu::Buffer,
    pub bind_group: wgpu::BindGroup,
    pub dims: [u32; 3],
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn density_params_layout() {
        // 3×16 = 48 bytes. Mirrors the WGSL struct layout — break this and the
        // shader will silently misalign uniforms.
        assert_eq!(std::mem::size_of::<DensityParams>(), 48);
        assert_eq!(std::mem::align_of::<DensityParams>(), 4);
    }

    #[test]
    fn accel_params_layout() {
        assert_eq!(std::mem::size_of::<SurfaceAccelParams>(), 32);
        assert_eq!(std::mem::align_of::<SurfaceAccelParams>(), 4);
    }
}
