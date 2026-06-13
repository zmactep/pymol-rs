//! Compute pipeline that builds the dot atom-instance buffer on the GPU. One
//! thread per atom; the vertex shader expands each visible atom into procedural
//! surface samples.

use bytemuck::{Pod, Zeroable};

use crate::compute::{
    build_3binding_build_layout, make_compute_build_pipeline, ComputeBuildPipeline,
};
use crate::scene_store::SceneStoreLayout;
use crate::shader_source;

pub const WORKGROUP: u32 = 64;

#[repr(C, align(16))]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
pub struct DotBuildParams {
    pub instance_capacity: u32,
    pub _pad0: u32,
    pub _pad1: u32,
    pub _pad2: u32,
}

impl DotBuildParams {
    pub const SIZE: u64 = std::mem::size_of::<Self>() as u64;
}

pub use crate::compute::indirect_seed;

pub struct DotBuildPipeline {
    pub pipeline: wgpu::ComputePipeline,
    pub build_layout: wgpu::BindGroupLayout,
}

impl DotBuildPipeline {
    pub fn new(device: &wgpu::Device, scene_layout: &SceneStoreLayout) -> Self {
        // DotAtomInstance = 32 B (center/group + vdW radius/pad).
        let build_layout = build_3binding_build_layout(
            device,
            "patinae.dot_build.layout",
            DotBuildParams::SIZE,
            32,
        );
        let ComputeBuildPipeline {
            pipeline,
            build_layout,
        } = make_compute_build_pipeline(
            device,
            scene_layout,
            "patinae.dot_build",
            shader_source::BUILD_DOT_WGSL,
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
            label: Some("patinae.dot_build.bg"),
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
        atom_count: u32,
    ) {
        if atom_count == 0 {
            return;
        }
        let groups = atom_count.div_ceil(WORKGROUP);
        let (wg_x, wg_y) = super::split_1d_dispatch(groups);
        let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
            label: Some("patinae.dot_build.dispatch"),
            timestamp_writes: None,
        });
        pass.set_pipeline(&self.pipeline);
        pass.set_bind_group(0, scene_bg, &[obj_dynamic_offset]);
        pass.set_bind_group(1, build_bg, &[]);
        pass.dispatch_workgroups(wg_x, wg_y, 1);
    }
}
