use patinae_plugin::prelude::ViewerLike;
use patinae_scene::{
    GpuBatchCommand, GpuBindGroupDescriptor, GpuBindGroupLayoutDescriptor, GpuBufferBindingType,
    GpuComputePipelineDescriptor, GpuDeviceLimits, GpuHandle, GpuPipelineLayoutDescriptor,
    GpuShaderModuleDescriptor, GpuSubmitBatch,
};

use super::dispatch::dispatch_grid_for;
use super::layout::ArtifactVisibleTriangleParams;
use super::reps::{ArtifactPlan, TriangleArtifactRep};
use super::resources::{
    buffer_entry, checked_storage_buffer_size, checked_storage_bytes, create_buffer,
    storage_layout, storage_usage, uniform_usage,
};
use crate::gpu::RaytraceParams;
use crate::shaders;

pub(super) struct SurfaceVisibility {
    pub(super) counter_buffer: GpuHandle,
    pub(super) counts: Vec<u32>,
}

pub(super) fn count_visible_surface_triangles(
    viewer: &mut dyn ViewerLike,
    plan: &ArtifactPlan<'_>,
    params: &RaytraceParams,
    device_limits: &GpuDeviceLimits,
    max_dispatch_dimension: u32,
) -> Result<Option<SurfaceVisibility>, String> {
    let surface_rep_count = plan
        .triangle_reps
        .iter()
        .filter(|rep| rep.uses_surface_visibility_culling())
        .count();
    if surface_rep_count == 0 {
        return Ok(None);
    }

    let counter_bytes = checked_storage_buffer_size(
        checked_storage_bytes(
            surface_rep_count as u64,
            std::mem::size_of::<u32>() as u64,
            "surface visibility counters",
        )?,
        device_limits,
        "surface visibility counters",
    )?;
    let counter_buffer = create_buffer(
        viewer,
        "ray.artifact.surface_visible_counts",
        counter_bytes,
        storage_usage(),
        Some(vec![0; counter_bytes as usize]),
    )?;

    let shader = viewer
        .gpu_create_cached_shader_module(GpuShaderModuleDescriptor {
            label: Some("ray.artifact.count_visible_surface_triangles.shader".to_string()),
            wgsl: shaders::artifact_count_visible_triangles(),
        })?
        .handle;
    let layout = viewer
        .gpu_create_cached_bind_group_layout(GpuBindGroupLayoutDescriptor {
            label: Some("ray.artifact.count_visible_surface_triangles.layout".to_string()),
            entries: vec![
                storage_layout(0, GpuBufferBindingType::StorageReadOnly),
                storage_layout(1, GpuBufferBindingType::StorageReadWrite),
                storage_layout(2, GpuBufferBindingType::Uniform),
            ],
        })?
        .handle;
    let pipeline_layout = viewer
        .gpu_create_cached_pipeline_layout(GpuPipelineLayoutDescriptor {
            label: Some("ray.artifact.count_visible_surface_triangles.pipeline_layout".to_string()),
            bind_group_layouts: vec![layout],
        })?
        .handle;
    let pipeline = viewer
        .gpu_create_cached_compute_pipeline(GpuComputePipelineDescriptor {
            label: Some("ray.artifact.count_visible_surface_triangles.pipeline".to_string()),
            layout: pipeline_layout,
            module: shader,
            entry_point: "count_visible_triangles".to_string(),
        })?
        .handle;
    let params_buffer = create_buffer(
        viewer,
        "ray.artifact.count_visible_surface_triangles.params",
        std::mem::size_of::<ArtifactVisibleTriangleParams>() as u64,
        uniform_usage(),
        None,
    )?;

    let mut batch_commands = Vec::new();
    let mut counter_index = 0_u32;
    for rep in plan
        .triangle_reps
        .iter()
        .filter(|rep| rep.uses_surface_visibility_culling())
    {
        let grid = dispatch_grid_for(rep.source_triangle_count(), max_dispatch_dimension)?;
        let cull_params =
            visible_triangle_params(params, rep, 0, counter_index, grid.invocation_width);
        batch_commands.push(GpuBatchCommand::WriteBuffer {
            buffer: params_buffer,
            offset: 0,
            data: bytemuck::bytes_of(&cull_params).to_vec(),
        });
        let bind_group = viewer.gpu_create_bind_group(GpuBindGroupDescriptor {
            label: Some("ray.artifact.count_visible_surface_triangles.bind_group".to_string()),
            layout,
            entries: vec![
                buffer_entry(0, rep.rep.geometry, 0, Some(rep.geometry_binding_size)),
                buffer_entry(1, counter_buffer, 0, Some(counter_bytes)),
                buffer_entry(
                    2,
                    params_buffer,
                    0,
                    Some(std::mem::size_of::<ArtifactVisibleTriangleParams>() as u64),
                ),
            ],
        })?;
        batch_commands.push(GpuBatchCommand::DispatchCompute {
            pipeline,
            bind_groups: vec![bind_group],
            workgroups: grid.workgroups,
        });
        counter_index = counter_index
            .checked_add(1)
            .ok_or_else(|| "surface visibility counter index overflow".to_string())?;
    }
    batch_commands.push(GpuBatchCommand::ReadBuffer {
        buffer: counter_buffer,
        offset: 0,
        size: counter_bytes,
    });

    let result = viewer.gpu_submit_batch(GpuSubmitBatch {
        label: Some("ray.artifact.surface_visibility.batch".to_string()),
        commands: batch_commands,
    })?;
    if result.readbacks.len() != 1 {
        return Err(format!(
            "surface visibility batch returned {} readbacks, expected 1",
            result.readbacks.len()
        ));
    }
    let counts = decode_visible_triangle_counts(&result.readbacks[0], surface_rep_count)?;
    Ok(Some(SurfaceVisibility {
        counter_buffer,
        counts,
    }))
}

pub(super) fn visible_triangle_params(
    params: &RaytraceParams,
    rep: &TriangleArtifactRep<'_>,
    output_triangle_count: u32,
    counter_index: u32,
    dispatch_width: u32,
) -> ArtifactVisibleTriangleParams {
    ArtifactVisibleTriangleParams {
        view_matrix: params.view_matrix,
        proj_matrix: params.proj_matrix,
        source_vertex_count: rep.vertex_count,
        triangle_offset: rep.triangle_offset,
        output_triangle_count,
        atom_offset: rep.rep.atom_offset,
        rep_slot: rep.rep_slot,
        transparency: rep.rep.transparency,
        dispatch_width,
        counter_index,
    }
}

pub(super) fn decode_visible_triangle_counts(
    bytes: &[u8],
    count: usize,
) -> Result<Vec<u32>, String> {
    let expected_bytes = count
        .checked_mul(std::mem::size_of::<u32>())
        .ok_or_else(|| "surface visibility readback size overflow".to_string())?;
    if bytes.len() != expected_bytes {
        return Err(format!(
            "surface visibility readback returned {} bytes, expected {expected_bytes}",
            bytes.len()
        ));
    }

    let mut counts = Vec::with_capacity(count);
    for chunk in bytes.chunks_exact(std::mem::size_of::<u32>()) {
        let bytes: [u8; 4] = chunk
            .try_into()
            .map_err(|_| "surface visibility counter chunk size mismatch".to_string())?;
        counts.push(u32::from_le_bytes(bytes));
    }
    Ok(counts)
}
