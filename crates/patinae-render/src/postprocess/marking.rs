//! Selection and hover overlay.
//!
//! The overlay is a small postprocess chain:
//! 1. resolve the visible id texture + marker LUT into a binary mask,
//! 2. composite a bright magenta selection tint over the scene colour.

use bytemuck::{Pod, Zeroable};
use wgpu::util::DeviceExt;

use crate::context::RenderContext;
use crate::frame::MARKING_MASK_FORMAT;
use crate::shader_source;

pub const MAX_MARKING_OBJECTS: usize = 4096;

#[repr(C, align(16))]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
pub struct MarkingParams {
    /// `(1 / width, 1 / height, rim_px, _)`.
    pub inv_size_radii: [f32; 4],
}

pub struct MarkingBindGroups {
    pub mask: wgpu::BindGroup,
    pub composite: wgpu::BindGroup,
}

pub struct MarkingPass {
    pub mask_pipeline: wgpu::RenderPipeline,
    pub composite_pipeline: wgpu::RenderPipeline,
    pub mask_layout: wgpu::BindGroupLayout,
    pub composite_layout: wgpu::BindGroupLayout,
    pub params_buffer: wgpu::Buffer,
    pub object_offsets_buffer: wgpu::Buffer,
}

struct FullscreenPass<'a> {
    target: &'a wgpu::TextureView,
    pipeline: &'a wgpu::RenderPipeline,
    bg: &'a wgpu::BindGroup,
    load: wgpu::LoadOp<wgpu::Color>,
    timestamp_writes: Option<wgpu::RenderPassTimestampWrites<'a>>,
    label: &'static str,
}

impl MarkingPass {
    pub fn new(ctx: &RenderContext) -> Self {
        let device = &ctx.device;

        let params_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("patinae.marking.params"),
            contents: bytemuck::bytes_of(&MarkingParams {
                inv_size_radii: [1.0, 1.0, 2.0, 10.0],
            }),
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        });
        let object_offsets_buffer = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("patinae.marking.object_offsets"),
            size: (MAX_MARKING_OBJECTS * std::mem::size_of::<u32>()) as u64,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        let mask_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("patinae.marking.mask.bgl"),
            entries: &[
                wgpu::BindGroupLayoutEntry {
                    binding: 0,
                    visibility: wgpu::ShaderStages::FRAGMENT,
                    ty: wgpu::BindingType::Texture {
                        sample_type: wgpu::TextureSampleType::Uint,
                        view_dimension: wgpu::TextureViewDimension::D2,
                        multisampled: false,
                    },
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 1,
                    visibility: wgpu::ShaderStages::FRAGMENT,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Storage { read_only: true },
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 2,
                    visibility: wgpu::ShaderStages::FRAGMENT,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Storage { read_only: true },
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
            ],
        });

        let composite_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("patinae.marking.composite.bgl"),
            entries: &[
                wgpu::BindGroupLayoutEntry {
                    binding: 0,
                    visibility: wgpu::ShaderStages::FRAGMENT,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Uniform,
                        has_dynamic_offset: false,
                        min_binding_size: wgpu::BufferSize::new(
                            std::mem::size_of::<MarkingParams>() as u64,
                        ),
                    },
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 1,
                    visibility: wgpu::ShaderStages::FRAGMENT,
                    ty: wgpu::BindingType::Texture {
                        sample_type: wgpu::TextureSampleType::Float { filterable: true },
                        view_dimension: wgpu::TextureViewDimension::D2,
                        multisampled: false,
                    },
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 2,
                    visibility: wgpu::ShaderStages::FRAGMENT,
                    ty: wgpu::BindingType::Texture {
                        sample_type: wgpu::TextureSampleType::Float { filterable: false },
                        view_dimension: wgpu::TextureViewDimension::D2,
                        multisampled: false,
                    },
                    count: None,
                },
            ],
        });

        let mask_pipeline = make_pipeline(
            device,
            "patinae.marking.mask.pipeline",
            &mask_layout,
            MARKING_MASK_FORMAT,
            shader_source::MARKING_MASK_WGSL,
        );
        let composite_pipeline = make_pipeline(
            device,
            "patinae.marking.composite.pipeline",
            &composite_layout,
            ctx.color_format,
            shader_source::MARKING_COMPOSITE_WGSL,
        );

        Self {
            mask_pipeline,
            composite_pipeline,
            mask_layout,
            composite_layout,
            params_buffer,
            object_offsets_buffer,
        }
    }

    pub fn upload_params(&self, queue: &wgpu::Queue, target: (u32, u32), selection_width: f32) {
        let width = target.0.max(1) as f32;
        let height = target.1.max(1) as f32;
        let style_width = selection_width.clamp(0.5, 20.0);
        let rim_px = (style_width * 2.0).clamp(1.5, 4.0);
        queue.write_buffer(
            &self.params_buffer,
            0,
            bytemuck::bytes_of(&MarkingParams {
                inv_size_radii: [1.0 / width, 1.0 / height, rim_px, 0.0],
            }),
        );
    }

    pub fn upload_object_offsets<I>(&self, queue: &wgpu::Queue, offsets: I)
    where
        I: IntoIterator<Item = (u32, u32)>,
    {
        let mut table = [0u32; MAX_MARKING_OBJECTS];
        for (object_id, atom_offset) in offsets {
            let idx = object_id as usize;
            if idx < table.len() {
                table[idx] = atom_offset;
            }
        }
        queue.write_buffer(&self.object_offsets_buffer, 0, bytemuck::cast_slice(&table));
    }

    pub fn make_bind_groups(
        &self,
        device: &wgpu::Device,
        overlay_id: &wgpu::TextureView,
        marker_lut: &wgpu::Buffer,
        scene_color: &wgpu::TextureView,
        mask: &wgpu::TextureView,
    ) -> MarkingBindGroups {
        let mask_bg = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("patinae.marking.mask.bg"),
            layout: &self.mask_layout,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: wgpu::BindingResource::TextureView(overlay_id),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: self.object_offsets_buffer.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: marker_lut.as_entire_binding(),
                },
            ],
        });
        let composite = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("patinae.marking.composite.bg"),
            layout: &self.composite_layout,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: self.params_buffer.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: wgpu::BindingResource::TextureView(scene_color),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: wgpu::BindingResource::TextureView(mask),
                },
            ],
        });
        MarkingBindGroups {
            mask: mask_bg,
            composite,
        }
    }

    pub fn record_mask<'a>(
        &'a self,
        encoder: &'a mut wgpu::CommandEncoder,
        mask: &'a wgpu::TextureView,
        bg: &'a wgpu::BindGroup,
        timestamp_writes: Option<wgpu::RenderPassTimestampWrites<'a>>,
    ) {
        self.record_fullscreen(
            encoder,
            FullscreenPass {
                target: mask,
                pipeline: &self.mask_pipeline,
                bg,
                load: wgpu::LoadOp::Clear(wgpu::Color::TRANSPARENT),
                timestamp_writes,
                label: "patinae.marking.mask_pass",
            },
        );
    }

    pub fn record_composite<'a>(
        &'a self,
        encoder: &'a mut wgpu::CommandEncoder,
        target: &'a wgpu::TextureView,
        bg: &'a wgpu::BindGroup,
    ) {
        self.record_fullscreen(
            encoder,
            FullscreenPass {
                target,
                pipeline: &self.composite_pipeline,
                bg,
                load: wgpu::LoadOp::Clear(wgpu::Color::TRANSPARENT),
                timestamp_writes: None,
                label: "patinae.marking.composite_pass",
            },
        );
    }

    fn record_fullscreen<'a>(
        &'a self,
        encoder: &'a mut wgpu::CommandEncoder,
        pass_desc: FullscreenPass<'a>,
    ) {
        let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
            label: Some(pass_desc.label),
            color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                view: pass_desc.target,
                depth_slice: None,
                resolve_target: None,
                ops: wgpu::Operations {
                    load: pass_desc.load,
                    store: wgpu::StoreOp::Store,
                },
            })],
            depth_stencil_attachment: None,
            timestamp_writes: pass_desc.timestamp_writes,
            occlusion_query_set: None,
            multiview_mask: None,
        });
        pass.set_pipeline(pass_desc.pipeline);
        pass.set_bind_group(0, pass_desc.bg, &[]);
        pass.draw(0..3, 0..1);
    }
}

fn make_pipeline(
    device: &wgpu::Device,
    label: &'static str,
    layout: &wgpu::BindGroupLayout,
    target_format: wgpu::TextureFormat,
    source: &str,
) -> wgpu::RenderPipeline {
    let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
        label: Some(label),
        bind_group_layouts: &[Some(layout)],
        immediate_size: 0,
    });
    let module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
        label: Some(label),
        source: wgpu::ShaderSource::Wgsl(source.into()),
    });
    device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
        label: Some(label),
        layout: Some(&pipeline_layout),
        vertex: wgpu::VertexState {
            module: &module,
            entry_point: Some("vs_main"),
            compilation_options: Default::default(),
            buffers: &[],
        },
        fragment: Some(wgpu::FragmentState {
            module: &module,
            entry_point: Some("fs_main"),
            compilation_options: Default::default(),
            targets: &[Some(wgpu::ColorTargetState {
                format: target_format,
                blend: None,
                write_mask: wgpu::ColorWrites::ALL,
            })],
        }),
        primitive: wgpu::PrimitiveState {
            topology: wgpu::PrimitiveTopology::TriangleList,
            ..Default::default()
        },
        depth_stencil: None,
        multisample: wgpu::MultisampleState::default(),
        multiview_mask: None,
        cache: None,
    })
}
