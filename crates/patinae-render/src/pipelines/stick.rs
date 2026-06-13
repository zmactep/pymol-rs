//! Stick (capsule) color pipeline — bond-impostor billboard + ray-capsule FS.
//!
//! Bind groups:
//! - 0: `FrameUniforms` (`RenderContext`)
//! - 1: shadow placeholder (`RenderContext`)
//! - 2: [`SceneStoreLayout`] — scene-wide LUTs + `obj` dynamic-uniform
//! - 3: [`StickParamsLayout`] — per-rep alpha multiplier
//!
//! The WBOIT path keeps `depth_write = false` and writes through
//! `write_translucent`.

use std::num::NonZeroU64;

use bytemuck::{Pod, Zeroable};

use crate::context::RenderContext;
use crate::pipelines::build_draw_pair;
use crate::representations::stick::StickInstance;
use crate::scene_store::SceneStoreLayout;
use crate::shader_source::STICK_WGSL;

/// Mirror of `StickParams` in `shaders/representations/stick.wgsl`. 16 B.
#[repr(C, align(16))]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
pub struct StickParams {
    pub alpha_mul: f32,
    pub _pad0: f32,
    pub _pad1: f32,
    pub _pad2: f32,
}

impl StickParams {
    pub const SIZE: u64 = std::mem::size_of::<Self>() as u64;
}

impl Default for StickParams {
    fn default() -> Self {
        Self {
            alpha_mul: 1.0,
            _pad0: 0.0,
            _pad1: 0.0,
            _pad2: 0.0,
        }
    }
}

pub struct StickParamsLayout {
    pub bind_group_layout: wgpu::BindGroupLayout,
}

impl StickParamsLayout {
    pub fn new(device: &wgpu::Device) -> Self {
        let bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("patinae.stick.params_layout"),
            entries: &[wgpu::BindGroupLayoutEntry {
                binding: 0,
                visibility: wgpu::ShaderStages::VERTEX | wgpu::ShaderStages::FRAGMENT,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Uniform,
                    has_dynamic_offset: false,
                    min_binding_size: NonZeroU64::new(StickParams::SIZE),
                },
                count: None,
            }],
        });
        Self { bind_group_layout }
    }
}

pub struct StickPipeline {
    pub pipeline: wgpu::RenderPipeline,
    pub pipeline_opaque: wgpu::RenderPipeline,
}

impl StickPipeline {
    pub fn new(
        ctx: &RenderContext,
        scene_layout: &SceneStoreLayout,
        stick_layout: &StickParamsLayout,
    ) -> Self {
        let pair = build_draw_pair(
            ctx,
            "patinae.stick",
            STICK_WGSL,
            &[
                &ctx.frame.bind_group_layout,
                &ctx.lighting.bind_group_layout,
                &scene_layout.bind_group_layout,
                &stick_layout.bind_group_layout,
            ],
            StickInstance::vertex_layout(),
            wgpu::PrimitiveTopology::TriangleList,
        );
        Self {
            pipeline: pair.translucent,
            pipeline_opaque: pair.opaque,
        }
    }
}
