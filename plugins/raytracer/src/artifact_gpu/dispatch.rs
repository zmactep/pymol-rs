use patinae_plugin::prelude::ViewerLike;
use patinae_scene::{
    GpuBatchCommand, GpuBindGroupDescriptor, GpuBindGroupLayoutDescriptor, GpuBufferBindingType,
    GpuComputePipelineDescriptor, GpuHandle, GpuPipelineLayoutDescriptor,
    GpuShaderModuleDescriptor, GpuSubmitBatch,
};

use super::layout::{
    DownsampleInput, DownsampleParams, RaytraceDispatchInput, MAX_OUTPUT_READBACK_BYTES,
    WORKGROUP_SIZE,
};
use super::resources::{buffer_entry, create_buffer, storage_layout, storage_usage, uniform_usage};
use super::rt_profile_enabled;
use crate::gpu::uniforms::RaytraceUniforms;
use crate::shaders;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) struct DispatchGrid {
    pub(super) workgroups: [u32; 3],
    pub(super) invocation_width: u32,
}

pub(super) fn dispatch_raytrace(
    viewer: &mut dyn ViewerLike,
    input: RaytraceDispatchInput<'_>,
    batch_commands: &mut Vec<GpuBatchCommand>,
) -> Result<Vec<u8>, String> {
    let params = input.params;
    let supersample = params.antialias.max(1);
    let render_width = params
        .width
        .checked_mul(supersample)
        .ok_or_else(|| "ray output width overflow".to_string())?;
    let render_height = params
        .height
        .checked_mul(supersample)
        .ok_or_else(|| "ray output height overflow".to_string())?;
    let render_bytes = rgba_bytes(render_width, render_height, "ray render output")?;
    let final_bytes = rgba_bytes(params.width, params.height, "ray final output")?;
    if final_bytes > MAX_OUTPUT_READBACK_BYTES {
        return Err(format!(
            "ray output readback {} bytes exceeds the 64 MiB command ABI limit",
            final_bytes
        ));
    }

    let output_buffer = create_buffer(
        viewer,
        "ray.artifact.output_rgba",
        render_bytes,
        storage_usage(),
        None,
    )?;
    let uniform_buffer = create_buffer(
        viewer,
        "ray.artifact.raytrace.uniforms",
        std::mem::size_of::<RaytraceUniforms>() as u64,
        uniform_usage(),
        Some(
            bytemuck::bytes_of(&RaytraceUniforms::from_counts(
                params,
                input.counts.spheres,
                input.counts.cylinders,
                input.counts.capsules,
                input.counts.triangles,
                input.bvh_node_count,
            ))
            .to_vec(),
        ),
    )?;

    let shader = viewer
        .gpu_create_cached_shader_module(GpuShaderModuleDescriptor {
            label: Some("ray.artifact.raytrace.shader".to_string()),
            wgsl: shaders::artifact_raytrace_output_buffer(),
        })?
        .handle;
    let layout = viewer
        .gpu_create_cached_bind_group_layout(GpuBindGroupLayoutDescriptor {
            label: Some("ray.artifact.raytrace.layout".to_string()),
            entries: vec![
                storage_layout(0, GpuBufferBindingType::Uniform),
                storage_layout(1, GpuBufferBindingType::StorageReadOnly),
                storage_layout(2, GpuBufferBindingType::StorageReadOnly),
                storage_layout(3, GpuBufferBindingType::StorageReadOnly),
                storage_layout(4, GpuBufferBindingType::StorageReadOnly),
                storage_layout(5, GpuBufferBindingType::StorageReadOnly),
                storage_layout(6, GpuBufferBindingType::StorageReadOnly),
                storage_layout(7, GpuBufferBindingType::StorageReadWrite),
            ],
        })?
        .handle;
    let pipeline_layout = viewer
        .gpu_create_cached_pipeline_layout(GpuPipelineLayoutDescriptor {
            label: Some("ray.artifact.raytrace.pipeline_layout".to_string()),
            bind_group_layouts: vec![layout],
        })?
        .handle;
    let pipeline = viewer
        .gpu_create_cached_compute_pipeline(GpuComputePipelineDescriptor {
            label: Some("ray.artifact.raytrace.pipeline".to_string()),
            layout: pipeline_layout,
            module: shader,
            entry_point: "main".to_string(),
        })?
        .handle;
    let bind_group = viewer.gpu_create_bind_group(GpuBindGroupDescriptor {
        label: Some("ray.artifact.raytrace.bind_group".to_string()),
        layout,
        entries: vec![
            buffer_entry(
                0,
                uniform_buffer,
                0,
                Some(std::mem::size_of::<RaytraceUniforms>() as u64),
            ),
            buffer_entry(1, input.primitives.spheres, 0, None),
            buffer_entry(2, input.primitives.cylinders, 0, None),
            buffer_entry(3, input.primitives.capsules, 0, None),
            buffer_entry(4, input.primitives.triangles, 0, None),
            buffer_entry(5, input.bvh.nodes, 0, None),
            buffer_entry(6, input.bvh.indices, 0, None),
            buffer_entry(7, output_buffer, 0, None),
        ],
    })?;

    let ray_workgroups = [render_width.div_ceil(8), render_height.div_ceil(8), 1];
    validate_workgroups(ray_workgroups, input.max_dispatch_dimension, "ray dispatch")?;
    batch_commands.push(GpuBatchCommand::DispatchCompute {
        pipeline,
        bind_groups: vec![bind_group],
        workgroups: ray_workgroups,
    });

    let (final_buffer, downsample_mode) = if supersample > 1 {
        (
            append_downsample(
                viewer,
                DownsampleInput {
                    source_buffer: output_buffer,
                    src_width: render_width,
                    dst_width: params.width,
                    dst_height: params.height,
                    factor: supersample,
                    max_dispatch_dimension: input.max_dispatch_dimension,
                },
                batch_commands,
            )?,
            format!("gpu_{}x", supersample),
        )
    } else {
        (output_buffer, "none".to_string())
    };

    batch_commands.push(GpuBatchCommand::ReadBuffer {
        buffer: final_buffer,
        offset: 0,
        size: final_bytes,
    });

    let batch_start = std::time::Instant::now();
    let batch = GpuSubmitBatch {
        label: Some("ray.artifact.batch".to_string()),
        commands: std::mem::take(batch_commands),
    };
    let result = viewer.gpu_submit_batch(batch)?;
    if rt_profile_enabled() {
        log::info!(
            "patinae.rt_profile plugin.gpu_batch_ms={} downsample_mode={} final_readback_bytes={}",
            batch_start.elapsed().as_millis(),
            downsample_mode,
            final_bytes
        );
    }

    let mut readbacks = result.readbacks;
    if readbacks.len() != 1 {
        return Err(format!(
            "ray GPU batch returned {} readbacks, expected 1",
            readbacks.len()
        ));
    }
    Ok(readbacks.remove(0))
}

fn append_downsample(
    viewer: &mut dyn ViewerLike,
    input: DownsampleInput,
    batch_commands: &mut Vec<GpuBatchCommand>,
) -> Result<GpuHandle, String> {
    let final_bytes = rgba_bytes(input.dst_width, input.dst_height, "ray downsample output")?;
    let output_buffer = create_buffer(
        viewer,
        "ray.artifact.downsample_rgba",
        final_bytes,
        storage_usage(),
        None,
    )?;
    let params = DownsampleParams {
        src_width: input.src_width,
        dst_width: input.dst_width,
        dst_height: input.dst_height,
        factor: input.factor,
    };
    let params_buffer = create_buffer(
        viewer,
        "ray.artifact.downsample.params",
        std::mem::size_of::<DownsampleParams>() as u64,
        uniform_usage(),
        Some(bytemuck::bytes_of(&params).to_vec()),
    )?;
    let shader = viewer
        .gpu_create_cached_shader_module(GpuShaderModuleDescriptor {
            label: Some("ray.artifact.downsample.shader".to_string()),
            wgsl: shaders::ARTIFACT_DOWNSAMPLE.to_string(),
        })?
        .handle;
    let layout = viewer
        .gpu_create_cached_bind_group_layout(GpuBindGroupLayoutDescriptor {
            label: Some("ray.artifact.downsample.layout".to_string()),
            entries: vec![
                storage_layout(0, GpuBufferBindingType::StorageReadOnly),
                storage_layout(1, GpuBufferBindingType::StorageReadWrite),
                storage_layout(2, GpuBufferBindingType::Uniform),
            ],
        })?
        .handle;
    let pipeline_layout = viewer
        .gpu_create_cached_pipeline_layout(GpuPipelineLayoutDescriptor {
            label: Some("ray.artifact.downsample.pipeline_layout".to_string()),
            bind_group_layouts: vec![layout],
        })?
        .handle;
    let pipeline = viewer
        .gpu_create_cached_compute_pipeline(GpuComputePipelineDescriptor {
            label: Some("ray.artifact.downsample.pipeline".to_string()),
            layout: pipeline_layout,
            module: shader,
            entry_point: "main".to_string(),
        })?
        .handle;
    let bind_group = viewer.gpu_create_bind_group(GpuBindGroupDescriptor {
        label: Some("ray.artifact.downsample.bind_group".to_string()),
        layout,
        entries: vec![
            buffer_entry(0, input.source_buffer, 0, None),
            buffer_entry(1, output_buffer, 0, None),
            buffer_entry(
                2,
                params_buffer,
                0,
                Some(std::mem::size_of::<DownsampleParams>() as u64),
            ),
        ],
    })?;
    let workgroups = [input.dst_width.div_ceil(8), input.dst_height.div_ceil(8), 1];
    validate_workgroups(
        workgroups,
        input.max_dispatch_dimension,
        "ray downsample dispatch",
    )?;
    batch_commands.push(GpuBatchCommand::DispatchCompute {
        pipeline,
        bind_groups: vec![bind_group],
        workgroups,
    });
    Ok(output_buffer)
}

fn rgba_bytes(width: u32, height: u32, label: &str) -> Result<u64, String> {
    u64::from(width)
        .checked_mul(u64::from(height))
        .and_then(|pixels| pixels.checked_mul(4))
        .ok_or_else(|| format!("{label} buffer size overflow"))
}

fn workgroups_for(items: u32) -> u32 {
    items.max(1).div_ceil(WORKGROUP_SIZE)
}

pub(super) fn dispatch_grid_for(
    items: u32,
    max_workgroups_per_dimension: u32,
) -> Result<DispatchGrid, String> {
    if max_workgroups_per_dimension == 0 {
        return Err("GPU device reports zero max compute workgroups per dimension".to_string());
    }
    let total_workgroups = workgroups_for(items);
    let workgroups_x = total_workgroups.min(max_workgroups_per_dimension);
    let workgroups_y = total_workgroups.div_ceil(workgroups_x);
    if workgroups_y > max_workgroups_per_dimension {
        return Err(format!(
            "GPU dispatch requires {} workgroups, exceeding {}x{} tiling capacity",
            total_workgroups, max_workgroups_per_dimension, max_workgroups_per_dimension
        ));
    }
    let invocation_width = workgroups_x
        .checked_mul(WORKGROUP_SIZE)
        .ok_or_else(|| "GPU dispatch invocation width overflow".to_string())?;
    Ok(DispatchGrid {
        workgroups: [workgroups_x, workgroups_y, 1],
        invocation_width,
    })
}

pub(super) fn validate_workgroups(
    workgroups: [u32; 3],
    max_workgroups_per_dimension: u32,
    label: &str,
) -> Result<(), String> {
    for (axis, count) in ["x", "y", "z"].into_iter().zip(workgroups) {
        if count > max_workgroups_per_dimension {
            return Err(format!(
                "{label} workgroups.{axis}={} exceeds GPU limit {}",
                count, max_workgroups_per_dimension
            ));
        }
    }
    Ok(())
}
