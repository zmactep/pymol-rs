//! Cartoon / ribbon color pipeline — `StdVertex` triangle list, scene-store
//! group 2, per-rep alpha multiplier on group 3.
//!
//! Bind groups:
//! - 0: `FrameUniforms` (`RenderContext`)
//! - 1: shadow placeholder (`RenderContext`)
//! - 2: [`SceneStoreLayout`] — scene-wide LUTs + `obj` dynamic-uniform
//! - 3: [`CartoonParamsLayout`] — per-rep alpha multiplier
//!
//! Reuses `StdVertex` from the mesh module. `cartoon.wgsl` reads atom
//! state from the scene-wide `SceneStore` via `INCLUDE_SCENE`.

use std::num::NonZeroU64;

use bytemuck::{Pod, Zeroable};

use crate::context::RenderContext;
use crate::pipelines::build_draw_pair;
use crate::representations::mesh::StdVertex;
use crate::scene_store::SceneStoreLayout;
use crate::shader_source::CARTOON_WGSL;

/// Mirror of `CartoonParams` in `shaders/representations/cartoon.wgsl`.
#[repr(C, align(16))]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
pub struct CartoonParams {
    pub alpha_mul: f32,
    pub color_slot: u32,
    pub _pad1: u32,
    pub _pad2: u32,
}

impl CartoonParams {
    pub const SIZE: u64 = std::mem::size_of::<Self>() as u64;
}

impl Default for CartoonParams {
    fn default() -> Self {
        Self {
            alpha_mul: 1.0,
            color_slot: 0,
            _pad1: 0,
            _pad2: 0,
        }
    }
}

pub struct CartoonParamsLayout {
    pub bind_group_layout: wgpu::BindGroupLayout,
}

impl CartoonParamsLayout {
    pub fn new(device: &wgpu::Device) -> Self {
        let bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("patinae.cartoon.params_layout"),
            entries: &[wgpu::BindGroupLayoutEntry {
                binding: 0,
                visibility: wgpu::ShaderStages::VERTEX | wgpu::ShaderStages::FRAGMENT,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Uniform,
                    has_dynamic_offset: false,
                    min_binding_size: NonZeroU64::new(CartoonParams::SIZE),
                },
                count: None,
            }],
        });
        Self { bind_group_layout }
    }
}

pub struct CartoonPipeline {
    pub pipeline: wgpu::RenderPipeline,
    pub pipeline_opaque: wgpu::RenderPipeline,
}

impl CartoonPipeline {
    pub fn new(
        ctx: &RenderContext,
        scene_layout: &SceneStoreLayout,
        cartoon_layout: &CartoonParamsLayout,
    ) -> Self {
        let pair = build_draw_pair(
            ctx,
            "patinae.cartoon",
            CARTOON_WGSL,
            &[
                &ctx.frame.bind_group_layout,
                &ctx.lighting.bind_group_layout,
                &scene_layout.bind_group_layout,
                &cartoon_layout.bind_group_layout,
            ],
            StdVertex::vertex_layout(),
            wgpu::PrimitiveTopology::TriangleList,
        );
        Self {
            pipeline: pair.translucent,
            pipeline_opaque: pair.opaque,
        }
    }
}
