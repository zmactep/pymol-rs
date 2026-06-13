//! Compute pipeline that builds the stick instance buffer on the GPU.
//!
//! One thread per object-local bond. Reads [`crate::scene_store::SceneStore`] (group 0) at
//! the per-object dynamic offset; reads per-rep params + writes the
//! instance buffer + indirect args (group 1). Multi-bonds emit 2 or 3
//! instances per thread via repeated `atomicAdd`s.

use bytemuck::{Pod, Zeroable};

use crate::compute::{
    build_3binding_build_layout, make_compute_build_pipeline, ComputeBuildPipeline,
};
use crate::scene_store::SceneStoreLayout;
use crate::shader_source;

pub const WORKGROUP: u32 = 64;

/// Mirror of `StickBuildParams` in `build_stick.wgsl`. 16 B.
#[repr(C, align(16))]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
pub struct StickBuildParams {
    pub radius: f32,
    pub valence_scale: f32,
    pub valence_enabled: u32,
    pub sample_shift: u32,
}

impl StickBuildParams {
    pub const SIZE: u64 = std::mem::size_of::<Self>() as u64;
}

const _: () = assert!(std::mem::size_of::<StickBuildParams>() == 16);

pub use crate::compute::indirect_seed;

pub struct StickBuildPipeline {
    pub pipeline: wgpu::ComputePipeline,
    pub build_layout: wgpu::BindGroupLayout,
}

impl StickBuildPipeline {
    pub fn new(device: &wgpu::Device, scene_layout: &SceneStoreLayout) -> Self {
        // StickInstance = 48 B (vec4 + vec4 + vec2 + vec2 padding).
        let build_layout = build_3binding_build_layout(
            device,
            "patinae.stick_build.layout",
            StickBuildParams::SIZE,
            48,
        );
        let ComputeBuildPipeline {
            pipeline,
            build_layout,
        } = make_compute_build_pipeline(
            device,
            scene_layout,
            "patinae.stick_build",
            shader_source::BUILD_STICK_WGSL,
            build_layout,
        );
        Self {
            pipeline,
            build_layout,
        }
    }

    pub fn make_bind_group(
        &self,
        device: &wgpu::Device,
        params_buf: &wgpu::Buffer,
        raw_instance_buf: &wgpu::Buffer,
        raw_count_buf: &wgpu::Buffer,
    ) -> wgpu::BindGroup {
        device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("patinae.stick_build.bg"),
            layout: &self.build_layout,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: params_buf.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: raw_instance_buf.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: raw_count_buf.as_entire_binding(),
                },
            ],
        })
    }

    pub fn dispatch(
        &self,
        encoder: &mut wgpu::CommandEncoder,
        scene_bg: &wgpu::BindGroup,
        obj_dynamic_offset: u32,
        build_bg: &wgpu::BindGroup,
        bond_count: u32,
    ) {
        if bond_count == 0 {
            return;
        }
        let groups = bond_count.div_ceil(WORKGROUP);
        let (wg_x, wg_y) = super::split_1d_dispatch(groups);
        let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
            label: Some("patinae.stick_build.dispatch"),
            timestamp_writes: None,
        });
        pass.set_pipeline(&self.pipeline);
        pass.set_bind_group(0, scene_bg, &[obj_dynamic_offset]);
        pass.set_bind_group(1, build_bg, &[]);
        pass.dispatch_workgroups(wg_x, wg_y, 1);
    }
}
