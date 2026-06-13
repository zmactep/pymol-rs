//! Surface color pipeline — `StdVertex` triangle list emitted by the MC
//! compute pipeline. Scene-store group 2, per-rep alpha multiplier on
//! group 3.
//!
//! Bind groups:
//! - 0: `FrameUniforms` (`RenderContext`)
//! - 1: shadow placeholder (`RenderContext`)
//! - 2: [`SceneStoreLayout`] — scene-wide LUTs + `obj` dynamic-uniform
//! - 3: [`SurfaceParamsLayout`] — per-rep alpha multiplier
//!
//! Reuses `StdVertex` from the mesh module. `surface.wgsl` reads atom
//! state from the scene-wide `SceneStore` via `INCLUDE_SCENE`.

use std::num::NonZeroU64;

use bytemuck::{Pod, Zeroable};

use crate::context::RenderContext;
use crate::pipelines::build_draw_pair;
use crate::representations::mesh::StdVertex;
use crate::scene_store::SceneStoreLayout;
use crate::shader_source::SURFACE_WGSL;

/// Mirror of `SurfaceParams` in `shaders/representations/surface.wgsl`.
#[repr(C, align(16))]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
pub struct SurfaceParams {
    pub alpha_mul: f32,
    pub color_smoothing: f32,
    pub color_smoothing_threshold: f32,
    pub _pad2: f32,
}

impl SurfaceParams {
    pub const SIZE: u64 = std::mem::size_of::<Self>() as u64;
}

impl Default for SurfaceParams {
    fn default() -> Self {
        Self {
            alpha_mul: 1.0,
            color_smoothing: 1.0,
            color_smoothing_threshold: 0.05,
            _pad2: 0.0,
        }
    }
}

pub struct SurfaceParamsLayout {
    pub bind_group_layout: wgpu::BindGroupLayout,
}

impl SurfaceParamsLayout {
    pub fn new(device: &wgpu::Device) -> Self {
        let bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("patinae.surface.params_layout"),
            entries: &[wgpu::BindGroupLayoutEntry {
                binding: 0,
                visibility: wgpu::ShaderStages::VERTEX | wgpu::ShaderStages::FRAGMENT,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Uniform,
                    has_dynamic_offset: false,
                    min_binding_size: NonZeroU64::new(SurfaceParams::SIZE),
                },
                count: None,
            }],
        });
        Self { bind_group_layout }
    }
}

pub struct SurfacePipeline {
    pub pipeline: wgpu::RenderPipeline,
    pub pipeline_opaque: wgpu::RenderPipeline,
}

impl SurfacePipeline {
    pub fn new(
        ctx: &RenderContext,
        scene_layout: &SceneStoreLayout,
        surface_layout: &SurfaceParamsLayout,
    ) -> Self {
        let pair = build_draw_pair(
            ctx,
            "patinae.surface",
            SURFACE_WGSL,
            &[
                &ctx.frame.bind_group_layout,
                &ctx.lighting.bind_group_layout,
                &scene_layout.bind_group_layout,
                &surface_layout.bind_group_layout,
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
