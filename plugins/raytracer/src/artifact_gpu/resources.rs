use patinae_plugin::prelude::ViewerLike;
use patinae_scene::{
    GpuBindGroupEntry, GpuBindGroupLayoutEntry, GpuBindingResource, GpuBindingType,
    GpuBufferBinding, GpuBufferBindingType, GpuBufferDescriptor, GpuBufferUsage, GpuDeviceLimits,
    GpuHandle, GpuShaderStages,
};

use super::layout::EMPTY_STORAGE_BYTES;

pub(super) fn storage_bytes_for<T>(count: u32, label: &str) -> Result<u64, String> {
    let element_size = std::mem::size_of::<T>() as u64;
    u64::from(count.max(1))
        .checked_mul(element_size)
        .ok_or_else(|| format!("{label} buffer size overflow"))
}

pub(super) fn storage_bytes_for_device<T>(
    count: u32,
    limits: &GpuDeviceLimits,
    label: &str,
) -> Result<u64, String> {
    let bytes = storage_bytes_for::<T>(count, label)?;
    checked_storage_buffer_size(bytes, limits, label)
}

pub(super) fn checked_storage_bytes(
    elements: u64,
    stride: u64,
    label: &str,
) -> Result<u64, String> {
    elements
        .checked_mul(stride)
        .ok_or_else(|| format!("{label} binding size overflow"))
}

pub(super) fn checked_storage_buffer_size(
    size: u64,
    limits: &GpuDeviceLimits,
    label: &str,
) -> Result<u64, String> {
    let effective_size = size.max(EMPTY_STORAGE_BYTES);
    let limit = limits
        .max_buffer_size
        .min(limits.max_storage_buffer_binding_size);
    if effective_size > limit {
        return Err(format!(
            "{label} buffer size {effective_size} exceeds GPU storage buffer limit {limit}"
        ));
    }
    Ok(size)
}

pub(super) fn direct_draw_args_buffer(
    viewer: &mut dyn ViewerLike,
    label: &str,
    args: [u32; 4],
) -> Result<GpuHandle, String> {
    create_buffer(
        viewer,
        label,
        std::mem::size_of::<[u32; 4]>() as u64,
        storage_usage(),
        Some(bytemuck::cast_slice(&args).to_vec()),
    )
}

pub(super) fn create_buffer(
    viewer: &mut dyn ViewerLike,
    label: &str,
    size: u64,
    usage: GpuBufferUsage,
    initial_data: Option<Vec<u8>>,
) -> Result<GpuHandle, String> {
    viewer.gpu_create_buffer(
        GpuBufferDescriptor {
            label: Some(label.to_string()),
            size: size.max(EMPTY_STORAGE_BYTES),
            usage,
        },
        initial_data,
    )
}

pub(super) fn storage_usage() -> GpuBufferUsage {
    usage(&[
        GpuBufferUsage::STORAGE,
        GpuBufferUsage::COPY_SRC,
        GpuBufferUsage::COPY_DST,
    ])
}

pub(super) fn uniform_usage() -> GpuBufferUsage {
    usage(&[GpuBufferUsage::UNIFORM, GpuBufferUsage::COPY_DST])
}

fn usage(parts: &[GpuBufferUsage]) -> GpuBufferUsage {
    parts
        .iter()
        .copied()
        .fold(GpuBufferUsage { bits: 0 }, GpuBufferUsage::union)
}

pub(super) fn storage_layout(binding: u32, ty: GpuBufferBindingType) -> GpuBindGroupLayoutEntry {
    GpuBindGroupLayoutEntry {
        binding,
        visibility: GpuShaderStages::COMPUTE,
        ty: GpuBindingType::Buffer {
            ty,
            has_dynamic_offset: false,
            min_binding_size: None,
        },
    }
}

pub(super) fn buffer_entry(
    binding: u32,
    buffer: GpuHandle,
    offset: u64,
    size: Option<u64>,
) -> GpuBindGroupEntry {
    GpuBindGroupEntry {
        binding,
        resource: GpuBindingResource::Buffer(GpuBufferBinding {
            buffer,
            offset,
            size,
        }),
    }
}
