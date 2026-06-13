//! Generic map contour pipelines.
//!
//! Map contours are not atom-owned molecular representations, so these
//! pipelines bind only frame, lighting, and a per-map material/model block.

use std::num::NonZeroU64;

use bytemuck::{Pod, Zeroable};

use crate::context::RenderContext;
use crate::map_contour::MapVertex;
use crate::pipelines::build_draw_pair;
use crate::shader_source::MAP_WGSL;

/// Uniform payload for one rendered map contour.
#[repr(C, align(16))]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
pub struct MapParams {
    pub color: [f32; 4],
    pub model: [[f32; 4]; 4],
}

impl MapParams {
    pub const SIZE: u64 = std::mem::size_of::<Self>() as u64;
}

/// Bind-group layout for map contour uniforms.
pub struct MapParamsLayout {
    pub bind_group_layout: wgpu::BindGroupLayout,
}

impl MapParamsLayout {
    pub fn new(device: &wgpu::Device) -> Self {
        let bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("patinae.map.params_layout"),
            entries: &[wgpu::BindGroupLayoutEntry {
                binding: 0,
                visibility: wgpu::ShaderStages::VERTEX | wgpu::ShaderStages::FRAGMENT,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Uniform,
                    has_dynamic_offset: false,
                    min_binding_size: NonZeroU64::new(MapParams::SIZE),
                },
                count: None,
            }],
        });
        Self { bind_group_layout }
    }
}

/// Render pipelines for map line and triangle contours.
pub struct MapPipeline {
    pub line: wgpu::RenderPipeline,
    pub line_opaque: wgpu::RenderPipeline,
    pub triangles: wgpu::RenderPipeline,
    pub triangles_opaque: wgpu::RenderPipeline,
}

impl MapPipeline {
    pub fn new(ctx: &RenderContext, map_layout: &MapParamsLayout) -> Self {
        let line_pair = build_draw_pair(
            ctx,
            "patinae.map.line",
            MAP_WGSL,
            &[
                &ctx.frame.bind_group_layout,
                &ctx.lighting.bind_group_layout,
                &map_layout.bind_group_layout,
            ],
            MapVertex::vertex_layout(),
            wgpu::PrimitiveTopology::LineList,
        );
        let triangle_pair = build_draw_pair(
            ctx,
            "patinae.map.triangle",
            MAP_WGSL,
            &[
                &ctx.frame.bind_group_layout,
                &ctx.lighting.bind_group_layout,
                &map_layout.bind_group_layout,
            ],
            MapVertex::vertex_layout(),
            wgpu::PrimitiveTopology::TriangleList,
        );
        Self {
            line: line_pair.translucent,
            line_opaque: line_pair.opaque,
            triangles: triangle_pair.translucent,
            triangles_opaque: triangle_pair.opaque,
        }
    }
}
