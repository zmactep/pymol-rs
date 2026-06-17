use patinae_plugin::prelude::ViewerLike;
use patinae_scene::{
    GpuBatchCommand, GpuBindGroupDescriptor, GpuBindGroupLayoutDescriptor, GpuBufferBindingType,
    GpuComputePipelineDescriptor, GpuHandle, GpuPipelineLayoutDescriptor,
    GpuShaderModuleDescriptor, GpuSubmitBatch,
};

use super::layout::{
    ArtifactPrimitiveMetadata, DownsampleInput, DownsampleParams, RaytraceDispatchInput,
    MAX_OUTPUT_READBACK_BYTES, WORKGROUP_SIZE,
};
use super::resources::{buffer_entry, create_buffer, storage_layout, storage_usage, uniform_usage};
use super::rt_profile_counts_enabled;
use crate::gpu::uniforms::RaytraceUniforms;
use crate::shaders;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) struct DispatchGrid {
    pub(super) workgroups: [u32; 3],
    pub(super) invocation_width: u32,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum RaytraceDispatchTarget {
    CpuReadback,
    ViewportGpu,
}

pub(super) enum RaytraceDispatchOutput {
    CpuImage {
        data: Vec<u8>,
        profile_lines: Vec<String>,
    },
    ViewportGpu {
        profile_lines: Vec<String>,
    },
}

pub(super) struct RaytraceReadbacks {
    pub(super) metadata: Option<ArtifactPrimitiveMetadata>,
    pub(super) pixels: Option<Vec<u8>>,
}

pub(super) fn dispatch_raytrace(
    viewer: &mut dyn ViewerLike,
    input: RaytraceDispatchInput<'_>,
    batch_commands: &mut Vec<GpuBatchCommand>,
    target: RaytraceDispatchTarget,
) -> Result<RaytraceDispatchOutput, String> {
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
    if target == RaytraceDispatchTarget::CpuReadback && final_bytes > MAX_OUTPUT_READBACK_BYTES {
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
                storage_layout(8, GpuBufferBindingType::StorageReadOnly),
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
            buffer_entry(8, input.metadata, 0, None),
        ],
    })?;

    let ray_workgroups = [render_width.div_ceil(8), render_height.div_ceil(8), 1];
    validate_workgroups(ray_workgroups, input.max_dispatch_dimension, "ray dispatch")?;
    super::set_profile_scope(batch_commands, super::PROFILE_SCOPE_RAYTRACE);
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

    let metadata_readback_enabled = rt_profile_counts_enabled();
    if metadata_readback_enabled {
        batch_commands.push(GpuBatchCommand::ReadBuffer {
            buffer: input.metadata,
            offset: 0,
            size: std::mem::size_of::<ArtifactPrimitiveMetadata>() as u64,
        });
    }
    if target == RaytraceDispatchTarget::CpuReadback {
        batch_commands.push(GpuBatchCommand::ReadBuffer {
            buffer: final_buffer,
            offset: 0,
            size: final_bytes,
        });
    }

    let batch_start = std::time::Instant::now();
    let batch_command_count = batch_commands.len();
    let batch = GpuSubmitBatch {
        label: Some("ray.artifact.batch".to_string()),
        wait_for_completion: target == RaytraceDispatchTarget::ViewportGpu,
        commands: std::mem::take(batch_commands),
    };
    let result = viewer.gpu_submit_batch(batch)?;
    let mut profile_lines = Vec::new();
    super::push_profile_line(
        &mut profile_lines,
        format!(
            "patinae.rt_profile plugin.gpu_batch_ms={} batch_commands={} ray_workgroups={}x{} downsample_mode={} final_readback_bytes={} metadata_readback_bytes={}",
            batch_start.elapsed().as_millis(),
            batch_command_count,
            ray_workgroups[0],
            ray_workgroups[1],
            downsample_mode,
            if target == RaytraceDispatchTarget::CpuReadback { final_bytes } else { 0 },
            if metadata_readback_enabled {
                std::mem::size_of::<ArtifactPrimitiveMetadata>() as u64
            } else {
                0
            }
        ),
    );
    super::push_profile_line(
        &mut profile_lines,
        format!(
            "patinae.rt_profile plugin.ray_settings ray_shadow={} ray_transparency_shadows={} ray_max_passes={} transparency_mode={} light_count={} spec_count={} ray_trace_fog={}",
            u8::from(params.settings.ray_shadow),
            u8::from(params.settings.ray_transparency_shadows),
            params.settings.ray_max_passes,
            params.settings.transparency_mode,
            params.settings.light_count,
            params.settings.spec_count,
            u8::from(params.settings.ray_trace_fog),
        ),
    );

    let readbacks = split_raytrace_readbacks(target, metadata_readback_enabled, result.readbacks)?;
    if let Some(metadata) = readbacks.metadata {
        super::push_profile_line(
            &mut profile_lines,
            format!(
            "patinae.rt_profile plugin.gpu_metadata active_spheres={} active_cylinders={} active_capsules={} active_triangles={} visible_triangles={} active_primitive_count={} triangle_capacity={} active_triangle_capacity_pct={} visible_triangle_capacity_pct={} overflow={}",
            metadata.sphere_count,
            metadata.cylinder_count,
            metadata.capsule_count,
            metadata.triangle_count,
            metadata.visible_triangle_count,
            metadata.primitive_count,
            metadata.triangle_capacity,
            profile_percent(metadata.triangle_count, metadata.triangle_capacity),
            profile_percent(metadata.visible_triangle_count, metadata.triangle_capacity),
            metadata.overflow
            ),
        );
    }
    match target {
        RaytraceDispatchTarget::CpuReadback => {
            let Some(pixels) = readbacks.pixels else {
                return Err("ray GPU export batch returned no pixel readback".to_string());
            };
            Ok(RaytraceDispatchOutput::CpuImage {
                data: pixels,
                profile_lines,
            })
        }
        RaytraceDispatchTarget::ViewportGpu => {
            viewer.set_viewport_gpu_image_from_buffer(final_buffer, params.width, params.height)?;
            Ok(RaytraceDispatchOutput::ViewportGpu { profile_lines })
        }
    }
}

pub(super) fn split_raytrace_readbacks(
    target: RaytraceDispatchTarget,
    metadata_readback_enabled: bool,
    mut readbacks: Vec<Vec<u8>>,
) -> Result<RaytraceReadbacks, String> {
    let expected = usize::from(metadata_readback_enabled)
        + usize::from(target == RaytraceDispatchTarget::CpuReadback);
    if readbacks.len() != expected {
        return Err(format!(
            "ray GPU batch returned {} readbacks, expected {expected}",
            readbacks.len()
        ));
    }

    let metadata = if metadata_readback_enabled {
        Some(decode_primitive_metadata_readback(&readbacks.remove(0))?)
    } else {
        None
    };
    let pixels = if target == RaytraceDispatchTarget::CpuReadback {
        Some(readbacks.remove(0))
    } else {
        None
    };

    Ok(RaytraceReadbacks { metadata, pixels })
}

fn profile_percent(numerator: u32, denominator: u32) -> String {
    if denominator == 0 {
        return "0.00".to_string();
    }
    format!(
        "{:.2}",
        f64::from(numerator) * 100.0 / f64::from(denominator)
    )
}

pub(super) fn decode_primitive_metadata_readback(
    bytes: &[u8],
) -> Result<ArtifactPrimitiveMetadata, String> {
    let expected = std::mem::size_of::<ArtifactPrimitiveMetadata>();
    if bytes.len() != expected {
        return Err(format!(
            "ray GPU metadata readback returned {} bytes, expected {expected}",
            bytes.len()
        ));
    }
    Ok(bytemuck::pod_read_unaligned(bytes))
}

pub(super) fn append_downsample(
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
    super::set_profile_scope(batch_commands, super::PROFILE_SCOPE_DOWNSAMPLE);
    batch_commands.push(GpuBatchCommand::DispatchCompute {
        pipeline,
        bind_groups: vec![bind_group],
        workgroups,
    });
    Ok(output_buffer)
}

pub(super) fn rgba_bytes(width: u32, height: u32, label: &str) -> Result<u64, String> {
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
