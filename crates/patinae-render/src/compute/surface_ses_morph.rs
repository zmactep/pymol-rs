//! Connolly SES probe-erosion morph.
//!
//! Reads the analytical vdW SDF written by [`super::surface_vdw_sdf`] and
//! writes the SES scalar field consumed by [`super::surface_mc`] (with
//! `invert_inside = 1`, `iso = 0`).
//!
//! See `shaders/compute/surface_ses_morph.wgsl` for the math; the short story
//! is `g_ses(x) = max_{y in ball(x, probe)} f_vdw(y) - probe`. Negative inside
//! the molecular volume, zero on the SES surface, positive in solvent.
//!
//! ## Bind group layout (group 0)
//!
//! | Binding | Resource                                     | Stage usage    |
//! |--------:|----------------------------------------------|----------------|
//! |       0 | `MorphParams` uniform                        | compute        |
//! |       1 | `sdf_in` r32float storage texture            | compute, read  |
//! |       2 | `field_out` rgba16float storage texture      | compute, write |

use bytemuck::{Pod, Zeroable};
use wgpu::util::DeviceExt;

use super::surface_density::WORKGROUP;
use crate::shader_source;

/// Mirrors `MorphParams` in `surface_ses_morph.wgsl`. 32 B
/// (`vec3<u32>` + f32 + f32 + u32 + 2×u32 padding).
#[repr(C)]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
pub struct MorphParams {
    pub dims: [u32; 3],
    pub voxel_size: f32,
    pub probe_radius: f32,
    pub stencil_radius_voxels: u32,
    pub _pad0: u32,
    pub _pad1: u32,
}

impl MorphParams {
    pub const SIZE: u64 = std::mem::size_of::<Self>() as u64;
}

/// Compute the integer ball-stencil radius (in voxels) that fully contains a
/// sphere of radius `probe_radius` centred on a voxel. `ceil(probe / voxel)`
/// — the WGSL pass uses `<=` against `r²` so the boundary cubes are sampled.
pub fn stencil_radius_voxels(probe_radius: f32, voxel_size: f32) -> u32 {
    if voxel_size <= 0.0 {
        return 0;
    }
    (probe_radius / voxel_size).ceil().max(0.0) as u32
}

pub struct SurfaceSesMorphCompute {
    pub pipeline: wgpu::ComputePipeline,
    pub bind_group_layout: wgpu::BindGroupLayout,
}

impl SurfaceSesMorphCompute {
    pub fn new(device: &wgpu::Device) -> Self {
        let bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("patinae.surface_ses_morph.bgl"),
            entries: &[
                wgpu::BindGroupLayoutEntry {
                    binding: 0,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Uniform,
                        has_dynamic_offset: false,
                        min_binding_size: wgpu::BufferSize::new(MorphParams::SIZE),
                    },
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 1,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    // r32float as a regular sampled texture — `textureLoad`
                    // only, no sampler. `Float { filterable: false }` is
                    // baseline-supported (Float { filterable: true } would
                    // need the `Float32Filterable` feature).
                    ty: wgpu::BindingType::Texture {
                        sample_type: wgpu::TextureSampleType::Float { filterable: false },
                        view_dimension: wgpu::TextureViewDimension::D3,
                        multisampled: false,
                    },
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 2,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::StorageTexture {
                        access: wgpu::StorageTextureAccess::WriteOnly,
                        format: wgpu::TextureFormat::Rgba16Float,
                        view_dimension: wgpu::TextureViewDimension::D3,
                    },
                    count: None,
                },
            ],
        });

        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("patinae.surface_ses_morph.pl"),
            bind_group_layouts: &[&bind_group_layout],
            immediate_size: 0,
        });

        let module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("patinae.surface_ses_morph.wgsl"),
            source: wgpu::ShaderSource::Wgsl(shader_source::SURFACE_SES_MORPH_WGSL.into()),
        });

        let pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("patinae.surface_ses_morph.pipeline"),
            layout: Some(&pipeline_layout),
            module: &module,
            entry_point: Some("cs_main"),
            compilation_options: wgpu::PipelineCompilationOptions::default(),
            cache: None,
        });

        Self {
            pipeline,
            bind_group_layout,
        }
    }

    pub fn build_inputs(
        &self,
        device: &wgpu::Device,
        params: MorphParams,
        sdf_view: &wgpu::TextureView,
        field_view: &wgpu::TextureView,
    ) -> SurfaceSesMorphInputs {
        let params_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("patinae.surface_ses_morph.params"),
            contents: bytemuck::bytes_of(&params),
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        });
        let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("patinae.surface_ses_morph.bg"),
            layout: &self.bind_group_layout,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: params_buf.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: wgpu::BindingResource::TextureView(sdf_view),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: wgpu::BindingResource::TextureView(field_view),
                },
            ],
        });
        SurfaceSesMorphInputs {
            params_buf,
            bind_group,
            dims: params.dims,
        }
    }

    pub fn dispatch(&self, encoder: &mut wgpu::CommandEncoder, inputs: &SurfaceSesMorphInputs) {
        let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
            label: Some("patinae.surface_ses_morph.pass"),
            timestamp_writes: None,
        });
        pass.set_pipeline(&self.pipeline);
        pass.set_bind_group(0, &inputs.bind_group, &[]);
        pass.dispatch_workgroups(
            inputs.dims[0].div_ceil(WORKGROUP),
            inputs.dims[1].div_ceil(WORKGROUP),
            inputs.dims[2].div_ceil(WORKGROUP),
        );
    }
}

pub struct SurfaceSesMorphInputs {
    pub params_buf: wgpu::Buffer,
    pub bind_group: wgpu::BindGroup,
    pub dims: [u32; 3],
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn morph_params_layout() {
        // vec3<u32> (12) + f32 (4) + f32 (4) + u32 (4) + 8 padding = 32 B.
        assert_eq!(std::mem::size_of::<MorphParams>(), 32);
    }

    #[test]
    fn stencil_radius_rounds_up() {
        assert_eq!(stencil_radius_voxels(1.4, 0.5), 3); // ceil(2.8) = 3
        assert_eq!(stencil_radius_voxels(1.4, 1.0), 2); // ceil(1.4) = 2
        assert_eq!(stencil_radius_voxels(0.5, 1.0), 1); // ceil(0.5) = 1
        assert_eq!(stencil_radius_voxels(0.0, 0.5), 0);
        assert_eq!(stencil_radius_voxels(1.0, 0.0), 0);
    }
}
