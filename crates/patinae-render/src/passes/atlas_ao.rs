//! Skripkin multi-direction atlas AO pass resources and cache state.

use wgpu::util::DeviceExt;

use crate::memory::{buffer_usage, estimate_texture_2d_bytes, GpuMemoryUsage};
use crate::passes::lighting::{create_depth_texture, MAX_ATLAS_DIRECTIONS};
use crate::uniforms::FrameUniforms;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AtlasAoKey {
    pub bounds_center: [u32; 3],
    pub bounds_radius: u32,
    pub shadow_signature: u64,
    pub directions: u32,
    pub map_size: u32,
    pub bias: u32,
}

#[derive(Debug, Clone, Copy, Default)]
pub struct AtlasAoCache {
    key: Option<AtlasAoKey>,
    pub rebuilt: u32,
    pub reused: u32,
}

impl AtlasAoCache {
    pub fn should_rebuild(&mut self, key: AtlasAoKey) -> bool {
        if self.key == Some(key) {
            self.reused = self.reused.saturating_add(1);
            false
        } else {
            self.key = Some(key);
            self.rebuilt = self.rebuilt.saturating_add(1);
            true
        }
    }

    pub fn invalidate(&mut self) {
        self.key = None;
    }
}

pub struct AtlasAoPass {
    pub depth_texture: wgpu::Texture,
    pub depth_view: wgpu::TextureView,
    frame_bind_group_layout: wgpu::BindGroupLayout,
    atlas_frame_buffers: Vec<wgpu::Buffer>,
    atlas_frame_bind_groups: Vec<wgpu::BindGroup>,
    pub size: u32,
    pub tile_size: u32,
    pub atlas_grid: u32,
    pub cache: AtlasAoCache,
}

impl AtlasAoPass {
    pub fn new(device: &wgpu::Device, frame_layout: &wgpu::BindGroupLayout) -> Self {
        let (depth_texture, depth_view) = create_depth_texture(device, 32);
        Self {
            depth_texture,
            depth_view,
            frame_bind_group_layout: frame_layout.clone(),
            atlas_frame_buffers: Vec::new(),
            atlas_frame_bind_groups: Vec::new(),
            size: 32,
            tile_size: 32,
            atlas_grid: 1,
            cache: AtlasAoCache::default(),
        }
    }

    pub fn ensure_atlas(
        &mut self,
        device: &wgpu::Device,
        direction_count: u32,
        requested_tile_size: u32,
    ) -> (u32, u32) {
        let count = direction_count.max(1).min(MAX_ATLAS_DIRECTIONS as u32);
        let grid = atlas_grid(count);
        let max_dim = device.limits().max_texture_dimension_2d.max(1);
        let max_tile = (max_dim / grid).max(1);
        let requested = requested_tile_size.clamp(32, 4096).next_power_of_two();
        let tile_size = requested.min(max_tile).max(1);
        let atlas_size = grid * tile_size;

        if atlas_size != self.size || tile_size != self.tile_size || grid != self.atlas_grid {
            let (depth_texture, depth_view) = create_depth_texture(device, atlas_size);
            self.depth_texture = depth_texture;
            self.depth_view = depth_view;
            self.size = atlas_size;
            self.tile_size = tile_size;
            self.atlas_grid = grid;
            self.cache.invalidate();
        }

        while self.atlas_frame_buffers.len() < count as usize {
            let buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("patinae.atlas_ao.frame_uniforms"),
                contents: bytemuck::bytes_of(&FrameUniforms::default()),
                usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            });
            let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some("patinae.atlas_ao.frame_bind_group"),
                layout: &self.frame_bind_group_layout,
                entries: &[wgpu::BindGroupEntry {
                    binding: 0,
                    resource: buffer.as_entire_binding(),
                }],
            });
            self.atlas_frame_buffers.push(buffer);
            self.atlas_frame_bind_groups.push(bind_group);
        }

        (grid, tile_size)
    }

    pub fn upload_atlas_frame(&self, queue: &wgpu::Queue, index: usize, frame: &FrameUniforms) {
        if let Some(buffer) = self.atlas_frame_buffers.get(index) {
            queue.write_buffer(buffer, 0, bytemuck::bytes_of(frame));
        }
    }

    pub fn atlas_frame_bind_group(&self, index: usize) -> Option<&wgpu::BindGroup> {
        self.atlas_frame_bind_groups.get(index)
    }

    pub(crate) fn memory_usage(&self) -> GpuMemoryUsage {
        let mut usage = GpuMemoryUsage::default();
        usage.add(GpuMemoryUsage::allocation(estimate_texture_2d_bytes(
            self.size,
            self.size,
            crate::passes::lighting::OCCLUSION_DEPTH_FORMAT,
        )));
        for buffer in &self.atlas_frame_buffers {
            usage.add(buffer_usage(buffer));
        }
        usage
    }
}

fn atlas_grid(direction_count: u32) -> u32 {
    let mut grid = 1u32;
    while grid.saturating_mul(grid) < direction_count {
        grid += 1;
    }
    grid
}
