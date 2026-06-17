use patinae_plugin::prelude::ViewerLike;
use patinae_scene::{
    GpuBatchCommand, GpuBindGroupDescriptor, GpuBindGroupLayoutDescriptor, GpuBufferBindingType,
    GpuComputePipelineDescriptor, GpuDeviceLimits, GpuHandle, GpuPipelineLayoutDescriptor,
    GpuShaderModuleDescriptor,
};

use super::dispatch::validate_workgroups;
use super::layout::{ArtifactPrimitiveMetadata, FinalizeStreamingMetadataParams, PrimitiveCounts};
use super::resources::{
    buffer_entry, checked_storage_buffer_size, create_buffer, storage_indirect_usage,
    storage_layout, storage_usage, uniform_usage,
};
use crate::shaders;

pub(super) const DISPATCH_INDIRECT_ARGS_BYTES: u64 = std::mem::size_of::<[u32; 3]>() as u64;

#[derive(Clone, Copy)]
pub(super) struct FinalizeStreamingMetadataInput {
    pub(super) metadata: GpuHandle,
    pub(super) visible_counter: GpuHandle,
    pub(super) visible_counter_bytes: u64,
    pub(super) leaf_dispatch_args: GpuHandle,
    pub(super) params: FinalizeStreamingMetadataParams,
    pub(super) max_dispatch_dimension: u32,
}

pub(super) fn create_primitive_metadata(
    viewer: &mut dyn ViewerLike,
    label: &str,
    counts: PrimitiveCounts,
    triangle_capacity: u32,
    device_limits: &GpuDeviceLimits,
) -> Result<GpuHandle, String> {
    let mut metadata = ArtifactPrimitiveMetadata::from_counts(counts)?;
    metadata.triangle_capacity = triangle_capacity;
    let bytes = checked_storage_buffer_size(
        std::mem::size_of::<ArtifactPrimitiveMetadata>() as u64,
        device_limits,
        label,
    )?;
    create_buffer(
        viewer,
        label,
        bytes,
        storage_usage(),
        Some(bytemuck::bytes_of(&metadata).to_vec()),
    )
}

pub(super) fn create_dispatch_indirect_args(
    viewer: &mut dyn ViewerLike,
    label: &str,
) -> Result<GpuHandle, String> {
    create_buffer(
        viewer,
        label,
        DISPATCH_INDIRECT_ARGS_BYTES,
        storage_indirect_usage(),
        Some(vec![0; DISPATCH_INDIRECT_ARGS_BYTES as usize]),
    )
}

pub(super) fn finalize_streaming_metadata(
    viewer: &mut dyn ViewerLike,
    input: FinalizeStreamingMetadataInput,
    batch_commands: &mut Vec<GpuBatchCommand>,
) -> Result<(), String> {
    validate_workgroups(
        [1, 1, 1],
        input.max_dispatch_dimension,
        "streaming metadata dispatch",
    )?;
    let shader = viewer
        .gpu_create_cached_shader_module(GpuShaderModuleDescriptor {
            label: Some("ray.artifact.finalize_streaming_metadata.shader".to_string()),
            wgsl: shaders::artifact_finalize_streaming_metadata(),
        })?
        .handle;
    let layout = viewer
        .gpu_create_cached_bind_group_layout(GpuBindGroupLayoutDescriptor {
            label: Some("ray.artifact.finalize_streaming_metadata.layout".to_string()),
            entries: vec![
                storage_layout(0, GpuBufferBindingType::StorageReadWrite),
                storage_layout(1, GpuBufferBindingType::StorageReadWrite),
                storage_layout(2, GpuBufferBindingType::Uniform),
                storage_layout(3, GpuBufferBindingType::StorageReadWrite),
            ],
        })?
        .handle;
    let pipeline_layout = viewer
        .gpu_create_cached_pipeline_layout(GpuPipelineLayoutDescriptor {
            label: Some("ray.artifact.finalize_streaming_metadata.pipeline_layout".to_string()),
            bind_group_layouts: vec![layout],
        })?
        .handle;
    let pipeline = viewer
        .gpu_create_cached_compute_pipeline(GpuComputePipelineDescriptor {
            label: Some("ray.artifact.finalize_streaming_metadata.pipeline".to_string()),
            layout: pipeline_layout,
            module: shader,
            entry_point: "finalize_streaming_metadata".to_string(),
        })?
        .handle;
    let params_buffer = create_buffer(
        viewer,
        "ray.artifact.finalize_streaming_metadata.params",
        std::mem::size_of::<FinalizeStreamingMetadataParams>() as u64,
        uniform_usage(),
        Some(bytemuck::bytes_of(&input.params).to_vec()),
    )?;
    let bind_group = viewer.gpu_create_bind_group(GpuBindGroupDescriptor {
        label: Some("ray.artifact.finalize_streaming_metadata.bind_group".to_string()),
        layout,
        entries: vec![
            buffer_entry(
                0,
                input.metadata,
                0,
                Some(std::mem::size_of::<ArtifactPrimitiveMetadata>() as u64),
            ),
            buffer_entry(
                1,
                input.visible_counter,
                0,
                Some(input.visible_counter_bytes),
            ),
            buffer_entry(
                2,
                params_buffer,
                0,
                Some(std::mem::size_of::<FinalizeStreamingMetadataParams>() as u64),
            ),
            buffer_entry(
                3,
                input.leaf_dispatch_args,
                0,
                Some(DISPATCH_INDIRECT_ARGS_BYTES),
            ),
        ],
    })?;
    super::set_profile_scope(batch_commands, super::PROFILE_SCOPE_METADATA_FINALIZE);
    batch_commands.push(GpuBatchCommand::DispatchCompute {
        pipeline,
        bind_groups: vec![bind_group],
        workgroups: [1, 1, 1],
    });
    Ok(())
}
