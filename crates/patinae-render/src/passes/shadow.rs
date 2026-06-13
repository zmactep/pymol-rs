//! Directional shadow pass resources.

use wgpu::util::DeviceExt;

use crate::passes::lighting::{create_depth_texture, DEFAULT_SHADOW_MAP_SIZE};
use crate::uniforms::FrameUniforms;

pub struct DirectionalShadowPass {
    pub depth_texture: wgpu::Texture,
    pub depth_view: wgpu::TextureView,
    pub frame_bind_group: wgpu::BindGroup,
    frame_buffer: wgpu::Buffer,
    pub size: u32,
}

impl DirectionalShadowPass {
    pub fn new(device: &wgpu::Device, frame_layout: &wgpu::BindGroupLayout) -> Self {
        let frame_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("patinae.shadow.frame_uniforms"),
            contents: bytemuck::bytes_of(&FrameUniforms::default()),
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        });
        let frame_bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("patinae.shadow.frame_bind_group"),
            layout: frame_layout,
            entries: &[wgpu::BindGroupEntry {
                binding: 0,
                resource: frame_buffer.as_entire_binding(),
            }],
        });
        let (depth_texture, depth_view) = create_depth_texture(device, DEFAULT_SHADOW_MAP_SIZE);
        Self {
            depth_texture,
            depth_view,
            frame_bind_group,
            frame_buffer,
            size: DEFAULT_SHADOW_MAP_SIZE,
        }
    }

    pub fn ensure_size(&mut self, device: &wgpu::Device, requested: u32) {
        let size = requested.clamp(64, 4096).next_power_of_two();
        if size == self.size {
            return;
        }
        let (depth_texture, depth_view) = create_depth_texture(device, size);
        self.depth_texture = depth_texture;
        self.depth_view = depth_view;
        self.size = size;
    }

    pub fn upload_frame(&self, queue: &wgpu::Queue, frame: &FrameUniforms) {
        queue.write_buffer(&self.frame_buffer, 0, bytemuck::bytes_of(frame));
    }
}
