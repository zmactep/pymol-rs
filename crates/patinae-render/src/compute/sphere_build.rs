//! Compute pipeline that builds the sphere instance buffer on the GPU.
//!
//! Replaces the CPU-side atom iteration in [`crate::representations::sphere::SphereRep`]. Reads
//! [`crate::scene_store::SceneStore`] (group 0) at the per-object dynamic offset; reads
//! per-rep params + writes the **raw** instance buffer + the
//! `raw_count` atomic counter (group 1). A separate cull pass
//! ([`crate::compute::cull`]) reads the raw buffer, applies
//! frustum culling, and writes the surviving subset into
//! a compacted buffer + the indirect-draw args. Build no longer touches
//! the indirect args — that ownership belongs to cull.
//!
//! Lifecycle:
//! - Pipeline + BGLs are built once (owned by `RenderContext` or
//!   `RenderState`).
//! - Per-rep state (`params_buf`, `raw_count_buf`, `bind_group`) is
//!   created by `SphereRep` against the layout exposed here.
//! - Each rebuild wave: host writes `0u32` into `raw_count_buf`, then
//!   dispatches `(atom_count + WORKGROUP - 1) / WORKGROUP` workgroups.
//!   The kernel `atomicAdd`s emitted instance slots into
//!   `raw_count[0]`.

use bytemuck::{Pod, Zeroable};

use crate::compute::{
    build_3binding_build_layout, make_compute_build_pipeline, ComputeBuildPipeline,
};
use crate::scene_store::SceneStoreLayout;
use crate::shader_source;

pub const WORKGROUP: u32 = 64;

/// Mirror of `SphereBuildParams` in `build_sphere.wgsl`. 16 B.
#[repr(C, align(16))]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
pub struct SphereBuildParams {
    pub sphere_scale: f32,
    pub sample_shift: u32,
    pub _pad0: u32,
    pub _pad1: u32,
}

impl SphereBuildParams {
    pub const SIZE: u64 = std::mem::size_of::<Self>() as u64;
}

pub use crate::compute::indirect_seed;

pub struct SphereBuildPipeline {
    pub pipeline: wgpu::ComputePipeline,
    pub build_layout: wgpu::BindGroupLayout,
}

impl SphereBuildPipeline {
    pub fn new(device: &wgpu::Device, scene_layout: &SceneStoreLayout) -> Self {
        // SphereInstance = 32 B (WGSL struct alignment).
        let build_layout = build_3binding_build_layout(
            device,
            "patinae.sphere_build.layout",
            SphereBuildParams::SIZE,
            32,
        );
        let ComputeBuildPipeline {
            pipeline,
            build_layout,
        } = make_compute_build_pipeline(
            device,
            scene_layout,
            "patinae.sphere_build",
            shader_source::BUILD_SPHERE_WGSL,
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
            label: Some("patinae.sphere_build.bg"),
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
            label: Some("patinae.sphere_build.dispatch"),
            timestamp_writes: None,
        });
        pass.set_pipeline(&self.pipeline);
        pass.set_bind_group(0, scene_bg, &[obj_dynamic_offset]);
        pass.set_bind_group(1, build_bg, &[]);
        pass.dispatch_workgroups(wg_x, wg_y, 1);
    }
}
