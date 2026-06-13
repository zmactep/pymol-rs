//! Dot color pipeline — one atom instance, procedural surface-sample billboards.

use bytemuck::{Pod, Zeroable};
use wgpu::util::DeviceExt;

use crate::context::RenderContext;
use crate::pipelines::{build_draw_pair, build_fast_overlay_pipeline};
use crate::representations::dot::DotAtomInstance;
use crate::scene_store::SceneStoreLayout;
use crate::shader_source::DOT_WGSL;

pub const DOT_SAMPLE_COUNTS: [u32; 6] = [8, 16, 32, 64, 128, 256];
pub const DOT_DIRECTION_COUNT: u32 = 8 + 16 + 32 + 64 + 128 + 256;

#[repr(C, align(16))]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
pub struct DotDrawParams {
    pub samples_per_atom: u32,
    pub dir_offset: u32,
    pub radius_px: f32,
    pub _pad0: u32,
}

impl DotDrawParams {
    pub const SIZE: u64 = std::mem::size_of::<Self>() as u64;
}

pub struct DotParamsLayout {
    pub bind_group_layout: wgpu::BindGroupLayout,
    directions_buffer: wgpu::Buffer,
}

impl DotParamsLayout {
    pub fn new(device: &wgpu::Device) -> Self {
        let bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("patinae.dot.params.layout"),
            entries: &[
                wgpu::BindGroupLayoutEntry {
                    binding: 0,
                    visibility: wgpu::ShaderStages::VERTEX,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Uniform,
                        has_dynamic_offset: false,
                        min_binding_size: std::num::NonZeroU64::new(DotDrawParams::SIZE),
                    },
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 1,
                    visibility: wgpu::ShaderStages::VERTEX,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Uniform,
                        has_dynamic_offset: false,
                        min_binding_size: std::num::NonZeroU64::new(
                            DOT_DIRECTION_COUNT as u64 * 16,
                        ),
                    },
                    count: None,
                },
            ],
        });
        let directions = dot_direction_table();
        let directions_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("patinae.dot.directions"),
            contents: bytemuck::cast_slice(&directions),
            usage: wgpu::BufferUsages::UNIFORM,
        });
        Self {
            bind_group_layout,
            directions_buffer,
        }
    }

    pub fn make_bind_group(
        &self,
        device: &wgpu::Device,
        params_buffer: &wgpu::Buffer,
    ) -> wgpu::BindGroup {
        device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("patinae.dot.params.bg"),
            layout: &self.bind_group_layout,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: params_buffer.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: self.directions_buffer.as_entire_binding(),
                },
            ],
        })
    }
}

pub fn dot_direction_offset(samples: u32) -> u32 {
    match samples {
        8 => 0,
        16 => 8,
        32 => 24,
        64 => 56,
        128 => 120,
        256 => 248,
        _ => 0,
    }
}

fn dot_direction_table() -> Vec<[f32; 4]> {
    let mut dirs = Vec::with_capacity(DOT_DIRECTION_COUNT as usize);
    for &count in &DOT_SAMPLE_COUNTS {
        for slot in 0..count {
            dirs.push(golden_spiral_direction(slot, count));
        }
    }
    dirs
}

fn golden_spiral_direction(slot: u32, count: u32) -> [f32; 4] {
    const GOLDEN_ANGLE: f32 = 2.399_963_1;
    let n = count.max(1) as f32;
    let denom = (n - 1.0).max(1.0);
    let i = slot as f32;
    let y = 1.0 - (i / denom) * 2.0;
    let r = (1.0 - y * y).max(0.0).sqrt();
    let theta = GOLDEN_ANGLE * i;
    [theta.cos() * r, y, theta.sin() * r, 0.0]
}

pub struct DotPipeline {
    pub pipeline: wgpu::RenderPipeline,
    pub pipeline_opaque: wgpu::RenderPipeline,
    pub pipeline_fast_overlay: wgpu::RenderPipeline,
}

impl DotPipeline {
    pub fn new(
        ctx: &RenderContext,
        scene_layout: &SceneStoreLayout,
        dot_layout: &DotParamsLayout,
    ) -> Self {
        let layouts = [
            &ctx.frame.bind_group_layout,
            &ctx.lighting.bind_group_layout,
            &scene_layout.bind_group_layout,
            &dot_layout.bind_group_layout,
        ];
        let pair = build_draw_pair(
            ctx,
            "patinae.dot",
            DOT_WGSL,
            &layouts,
            DotAtomInstance::vertex_layout(),
            wgpu::PrimitiveTopology::TriangleList,
        );
        let pipeline_fast_overlay = build_fast_overlay_pipeline(
            ctx,
            "patinae.dot.fast_overlay",
            DOT_WGSL,
            &layouts,
            DotAtomInstance::vertex_layout(),
            wgpu::PrimitiveTopology::TriangleList,
        );
        Self {
            pipeline: pair.translucent,
            pipeline_opaque: pair.opaque,
            pipeline_fast_overlay,
        }
    }
}
