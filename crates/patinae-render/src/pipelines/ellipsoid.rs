//! Ellipsoid color pipeline — atom-impostor billboard + ray-quadric FS,
//! WBOIT out.
//!
//! Bind groups:
//! - 0: `FrameUniforms` (`RenderContext`)
//! - 1: shadow placeholder (`RenderContext`)
//! - 2: [`SceneStoreLayout`] — scene-wide LUTs + `obj` dynamic-uniform
//! - 3: [`EllipsoidParamsLayout`] — per-rep alpha multiplier

use std::num::NonZeroU64;

use bytemuck::{Pod, Zeroable};

use crate::context::RenderContext;
use crate::pipelines::build_draw_pair;
use crate::representations::ellipsoid::EllipsoidInstance;
use crate::scene_store::SceneStoreLayout;
use crate::shader_source::ELLIPSOID_WGSL;

/// Mirror of `EllipsoidParams` in `shaders/representations/ellipsoid.wgsl`.
#[repr(C, align(16))]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
pub struct EllipsoidParams {
    pub alpha_mul: f32,
    pub _pad0: f32,
    pub _pad1: f32,
    pub _pad2: f32,
}

impl EllipsoidParams {
    pub const SIZE: u64 = std::mem::size_of::<Self>() as u64;
}

impl Default for EllipsoidParams {
    fn default() -> Self {
        Self {
            alpha_mul: 1.0,
            _pad0: 0.0,
            _pad1: 0.0,
            _pad2: 0.0,
        }
    }
}

pub struct EllipsoidParamsLayout {
    pub bind_group_layout: wgpu::BindGroupLayout,
}

impl EllipsoidParamsLayout {
    pub fn new(device: &wgpu::Device) -> Self {
        let bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("patinae.ellipsoid.params_layout"),
            entries: &[wgpu::BindGroupLayoutEntry {
                binding: 0,
                visibility: wgpu::ShaderStages::VERTEX | wgpu::ShaderStages::FRAGMENT,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Uniform,
                    has_dynamic_offset: false,
                    min_binding_size: NonZeroU64::new(EllipsoidParams::SIZE),
                },
                count: None,
            }],
        });
        Self { bind_group_layout }
    }
}

pub struct EllipsoidPipeline {
    pub pipeline: wgpu::RenderPipeline,
    pub pipeline_opaque: wgpu::RenderPipeline,
}

impl EllipsoidPipeline {
    pub fn new(
        ctx: &RenderContext,
        scene_layout: &SceneStoreLayout,
        ellipsoid_layout: &EllipsoidParamsLayout,
    ) -> Self {
        let pair = build_draw_pair(
            ctx,
            "patinae.ellipsoid",
            ELLIPSOID_WGSL,
            &[
                &ctx.frame.bind_group_layout,
                &ctx.lighting.bind_group_layout,
                &scene_layout.bind_group_layout,
                &ellipsoid_layout.bind_group_layout,
            ],
            EllipsoidInstance::vertex_layout(),
            wgpu::PrimitiveTopology::TriangleList,
        );
        Self {
            pipeline: pair.translucent,
            pipeline_opaque: pair.opaque,
        }
    }
}
