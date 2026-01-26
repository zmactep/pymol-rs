//! GPU buffer management utilities

use wgpu::util::DeviceExt;

/// Create a vertex buffer from a slice of vertex data
#[allow(dead_code)]
pub fn create_vertex_buffer<T: bytemuck::Pod>(
    device: &wgpu::Device,
    label: &str,
    data: &[T],
) -> wgpu::Buffer {
    device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some(label),
        contents: bytemuck::cast_slice(data),
        usage: wgpu::BufferUsages::VERTEX | wgpu::BufferUsages::COPY_DST,
    })
}

/// Create an index buffer from a slice of indices
#[allow(dead_code)]
pub fn create_index_buffer(device: &wgpu::Device, label: &str, data: &[u32]) -> wgpu::Buffer {
    device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some(label),
        contents: bytemuck::cast_slice(data),
        usage: wgpu::BufferUsages::INDEX | wgpu::BufferUsages::COPY_DST,
    })
}

/// Create a uniform buffer from data
pub fn create_uniform_buffer<T: bytemuck::Pod>(
    device: &wgpu::Device,
    label: &str,
    data: &T,
) -> wgpu::Buffer {
    device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some(label),
        contents: bytemuck::bytes_of(data),
        usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
    })
}

/// Update an existing buffer with new data
#[allow(dead_code)]
pub fn update_buffer<T: bytemuck::Pod>(queue: &wgpu::Queue, buffer: &wgpu::Buffer, data: &[T]) {
    queue.write_buffer(buffer, 0, bytemuck::cast_slice(data));
}

/// Update a uniform buffer with new data
pub fn update_uniform_buffer<T: bytemuck::Pod>(
    queue: &wgpu::Queue,
    buffer: &wgpu::Buffer,
    data: &T,
) {
    queue.write_buffer(buffer, 0, bytemuck::bytes_of(data));
}

/// A managed buffer that tracks its capacity and can grow
pub struct GrowableBuffer {
    buffer: Option<wgpu::Buffer>,
    capacity: usize,
    usage: wgpu::BufferUsages,
    label: String,
}

impl GrowableBuffer {
    /// Create a new growable buffer
    pub fn new(label: impl Into<String>, usage: wgpu::BufferUsages) -> Self {
        Self {
            buffer: None,
            capacity: 0,
            usage,
            label: label.into(),
        }
    }

    /// Ensure the buffer has at least the given capacity (in bytes)
    pub fn ensure_capacity(&mut self, device: &wgpu::Device, required_bytes: usize) {
        if self.capacity < required_bytes {
            // Grow by at least 2x or to required size
            let new_capacity = (self.capacity * 2).max(required_bytes).max(256);
            self.buffer = Some(device.create_buffer(&wgpu::BufferDescriptor {
                label: Some(&self.label),
                size: new_capacity as u64,
                usage: self.usage | wgpu::BufferUsages::COPY_DST,
                mapped_at_creation: false,
            }));
            self.capacity = new_capacity;
        }
    }

    /// Write data to the buffer, growing if necessary
    pub fn write<T: bytemuck::Pod>(
        &mut self,
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        data: &[T],
    ) {
        let bytes = bytemuck::cast_slice(data);
        self.ensure_capacity(device, bytes.len());
        if let Some(buffer) = &self.buffer {
            queue.write_buffer(buffer, 0, bytes);
        }
    }

    /// Get a reference to the underlying buffer
    pub fn buffer(&self) -> Option<&wgpu::Buffer> {
        self.buffer.as_ref()
    }
}
