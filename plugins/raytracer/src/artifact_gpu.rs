//! GPU artifact ray path for command-runtime handles.
//!
//! The product `ray` command must not pull displayed geometry back through the
//! ABI. This module consumes renderer-owned artifact handles, builds the ray
//! triangle and BVH buffers on the host GPU, and dispatches the ray shader via
//! the safe GPU command callbacks.
//!
//! The generic artifact contract ends at `RenderArtifactSnapshotDescriptor`.
//! Ray-specific logic starts here: supported `StdVertices` triangle-list
//! representations are converted into `GpuTriangle` records, then the plugin
//! builds a BVH and dispatches the ray shader through the portable GPU batch
//! API.
//!
//! Supported input is intentionally narrow in this path. Direct triangle-list
//! cartoon and ribbon artifacts are accepted when they provide a direct
//! `element_count` and the shared `SceneColorLut` stride is 64 bytes. A direct
//! triangle-list surface descriptor would follow the same path, but current
//! renderer-generated surface parts are indirect-counted and therefore return
//! a clear error until count/indirect support lands.
//!
//! Sphere, stick, line, mesh line-list, count-buffer, and indirect-draw
//! artifacts also return clear errors until native ray kernels for those
//! layouts are implemented.
//!
//! Work buffers and bind groups remain temporary command resources. Shader
//! modules, bind-group layouts, pipeline layouts, and compute pipelines use the
//! persistent host cache exposed by the GPU runtime.

use bytemuck::{Pod, Zeroable};
use patinae_plugin::prelude::ViewerLike;
use patinae_scene::{
    GpuBatchCommand, GpuBindGroupDescriptor, GpuBindGroupEntry, GpuBindGroupLayoutDescriptor,
    GpuBindGroupLayoutEntry, GpuBindingResource, GpuBindingType, GpuBufferBinding,
    GpuBufferBindingType, GpuBufferDescriptor, GpuBufferUsage, GpuComputePipelineDescriptor,
    GpuHandle, GpuPipelineLayoutDescriptor, GpuShaderModuleDescriptor, GpuShaderStages,
    GpuSubmitBatch, RenderArtifactBufferDescriptor, RenderArtifactBufferRole,
    RenderArtifactPrimitiveTopology, RenderArtifactRepDescriptor, RenderArtifactRepKind,
    RenderArtifactSnapshotDescriptor,
};

use crate::bvh::BvhNode;
use crate::gpu::uniforms::RaytraceUniforms;
use crate::gpu::RaytraceParams;
use crate::primitive::{GpuCylinder, GpuSphere, GpuTriangle};

const WORKGROUP_SIZE: u32 = 128;
const BVH_LEAF_SIZE: u32 = 4;
const EMPTY_STORAGE_BYTES: u64 = 16;
const COLOR_LUT_STRIDE: u64 = 64;
const MAX_OUTPUT_READBACK_BYTES: u64 = 64 * 1024 * 1024;

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct ArtifactConvertParams {
    vertex_count: u32,
    triangle_offset: u32,
    atom_offset: u32,
    rep_slot: u32,
    transparency: f32,
    dispatch_width: u32,
    _pad0: u32,
    _pad1: u32,
}

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct ArtifactBvhParams {
    triangle_count: u32,
    leaf_slots: u32,
    leaf_start: u32,
    level_start: u32,
    level_count: u32,
    dispatch_width: u32,
    _pad0: u32,
    _pad1: u32,
}

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct DownsampleParams {
    src_width: u32,
    dst_width: u32,
    dst_height: u32,
    factor: u32,
}

struct DirectTriangleRep<'a> {
    rep: &'a RenderArtifactRepDescriptor,
    triangle_offset: u32,
    triangle_count: u32,
    rep_slot: u32,
}

struct ArtifactPlan<'a> {
    color_lut: GpuHandle,
    triangle_reps: Vec<DirectTriangleRep<'a>>,
    triangle_count: u32,
}

struct BvhShape {
    leaf_slots: u32,
    leaf_start: u32,
    node_count: u32,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct DispatchGrid {
    workgroups: [u32; 3],
    invocation_width: u32,
}

/// Raytrace renderer GPU artifacts through the command-runtime GPU ABI.
///
/// This function consumes a command-scoped renderer artifact snapshot. It
/// validates supported layouts, records all conversion/BVH/raytrace work into
/// one ordered GPU batch, and returns only the final RGBA readback bytes.
pub(crate) fn raytrace_artifacts(
    viewer: &mut dyn ViewerLike,
    snapshot: &RenderArtifactSnapshotDescriptor,
    params: &RaytraceParams,
) -> Result<Vec<u8>, String> {
    if params.settings.ray_trace_mode != 0 {
        return Err(
            "native GPU artifact ray path currently supports ray_trace_mode=0 (Normal)".to_string(),
        );
    }

    let plan = plan_artifact_triangles(snapshot)?;
    if plan.triangle_count == 0 {
        return Err(
            "GPU artifact snapshot contains no raytraceable triangle primitives".to_string(),
        );
    }

    let setup_start = std::time::Instant::now();
    let mut batch_commands = Vec::new();
    let triangle_buffer = create_buffer(
        viewer,
        "ray.artifact.triangles",
        u64::from(plan.triangle_count) * std::mem::size_of::<GpuTriangle>() as u64,
        storage_usage(),
        None,
    )?;
    let max_dispatch_dimension = snapshot.device_limits.max_compute_workgroups_per_dimension;
    build_triangles(
        viewer,
        &plan,
        plan.color_lut,
        triangle_buffer,
        max_dispatch_dimension,
        &mut batch_commands,
    )?;

    let bvh_shape = bvh_shape_for(plan.triangle_count)?;
    let bvh_node_buffer = create_buffer(
        viewer,
        "ray.artifact.bvh_nodes",
        u64::from(bvh_shape.node_count) * std::mem::size_of::<BvhNode>() as u64,
        storage_usage(),
        None,
    )?;
    let bvh_index_buffer = create_buffer(
        viewer,
        "ray.artifact.bvh_indices",
        u64::from(plan.triangle_count) * std::mem::size_of::<u32>() as u64,
        storage_usage(),
        None,
    )?;
    build_bvh(
        viewer,
        triangle_buffer,
        bvh_node_buffer,
        bvh_index_buffer,
        plan.triangle_count,
        &bvh_shape,
        max_dispatch_dimension,
        &mut batch_commands,
    )?;
    if rt_profile_enabled() {
        log::info!(
            "patinae.rt_profile plugin.temp_resource_setup_ms={} triangles={} bvh_nodes={}",
            setup_start.elapsed().as_millis(),
            plan.triangle_count,
            bvh_shape.node_count
        );
    }

    let image = dispatch_raytrace(
        viewer,
        params,
        triangle_buffer,
        bvh_node_buffer,
        bvh_index_buffer,
        plan.triangle_count,
        bvh_shape.node_count,
        max_dispatch_dimension,
        &mut batch_commands,
    )?;

    Ok(image)
}

fn plan_artifact_triangles<'a>(
    snapshot: &'a RenderArtifactSnapshotDescriptor,
) -> Result<ArtifactPlan<'a>, String> {
    let color_lut = snapshot
        .buffers
        .iter()
        .find(|buffer| buffer.role == RenderArtifactBufferRole::SceneColorLut)
        .ok_or_else(|| "render artifact snapshot is missing SceneColorLut".to_string())?;
    validate_color_lut(color_lut)?;

    let mut triangle_reps = Vec::new();
    let mut triangle_count = 0_u32;
    for rep in &snapshot.reps {
        match rep.topology {
            RenderArtifactPrimitiveTopology::TriangleList => {
                let Some(rep_slot) = direct_triangle_rep_slot(rep.rep_kind) else {
                    continue;
                };
                if rep.indirect.is_some() || rep.count.is_some() || rep.element_count == 0 {
                    return Err(format!(
                        "native GPU artifact ray path requires direct element_count for {:?}",
                        rep.rep_kind
                    ));
                }
                if rep.element_count % 3 != 0 {
                    return Err(format!(
                        "{:?} artifact vertex count {} is not divisible by 3",
                        rep.rep_kind, rep.element_count
                    ));
                }
                let rep_triangles = u32::try_from(rep.element_count / 3)
                    .map_err(|_| "artifact triangle count exceeds u32".to_string())?;
                triangle_reps.push(DirectTriangleRep {
                    rep,
                    triangle_offset: triangle_count,
                    triangle_count: rep_triangles,
                    rep_slot,
                });
                triangle_count = triangle_count
                    .checked_add(rep_triangles)
                    .ok_or_else(|| "artifact triangle count overflow".to_string())?;
            }
            RenderArtifactPrimitiveTopology::SphereInstances
            | RenderArtifactPrimitiveTopology::CylinderInstances
            | RenderArtifactPrimitiveTopology::LineInstances => {
                if rep.max_element_count > 0 {
                    return Err(format!(
                        "native GPU artifact ray path does not yet support {:?} instance artifacts",
                        rep.rep_kind
                    ));
                }
            }
            RenderArtifactPrimitiveTopology::LineList => {
                if rep.max_element_count > 0 || rep.element_count > 0 {
                    return Err(format!(
                        "native GPU artifact ray path does not yet support {:?} line-list artifacts",
                        rep.rep_kind
                    ));
                }
            }
        }
    }

    Ok(ArtifactPlan {
        color_lut: color_lut.handle,
        triangle_reps,
        triangle_count,
    })
}

fn validate_color_lut(buffer: &RenderArtifactBufferDescriptor) -> Result<(), String> {
    if buffer.stride != COLOR_LUT_STRIDE {
        return Err(format!(
            "SceneColorLut stride {} does not match expected {}",
            buffer.stride, COLOR_LUT_STRIDE
        ));
    }
    Ok(())
}

fn direct_triangle_rep_slot(kind: RenderArtifactRepKind) -> Option<u32> {
    match kind {
        RenderArtifactRepKind::Cartoon => Some(4),
        RenderArtifactRepKind::Ribbon => Some(5),
        RenderArtifactRepKind::Surface => Some(6),
        _ => None,
    }
}

// Consumes the public `StdVertex` packing documented by render artifacts:
// position lanes, packed normal, group id, and flags in six u32 words.
fn build_triangles(
    viewer: &mut dyn ViewerLike,
    plan: &ArtifactPlan<'_>,
    color_lut: GpuHandle,
    triangle_buffer: GpuHandle,
    max_dispatch_dimension: u32,
    batch_commands: &mut Vec<GpuBatchCommand>,
) -> Result<(), String> {
    let shader = viewer
        .gpu_create_cached_shader_module(GpuShaderModuleDescriptor {
            label: Some("ray.artifact.build_triangles.shader".to_string()),
            wgsl: ARTIFACT_TRIANGLE_SHADER.to_string(),
        })?
        .handle;
    let layout = viewer
        .gpu_create_cached_bind_group_layout(GpuBindGroupLayoutDescriptor {
            label: Some("ray.artifact.build_triangles.layout".to_string()),
            entries: vec![
                storage_layout(0, GpuBufferBindingType::StorageReadOnly),
                storage_layout(1, GpuBufferBindingType::StorageReadOnly),
                storage_layout(2, GpuBufferBindingType::StorageReadWrite),
                storage_layout(3, GpuBufferBindingType::Uniform),
            ],
        })?
        .handle;
    let pipeline_layout = viewer
        .gpu_create_cached_pipeline_layout(GpuPipelineLayoutDescriptor {
            label: Some("ray.artifact.build_triangles.pipeline_layout".to_string()),
            bind_group_layouts: vec![layout],
        })?
        .handle;
    let pipeline = viewer
        .gpu_create_cached_compute_pipeline(GpuComputePipelineDescriptor {
            label: Some("ray.artifact.build_triangles.pipeline".to_string()),
            layout: pipeline_layout,
            module: shader,
            entry_point: "build_triangles".to_string(),
        })?
        .handle;
    let params_buffer = create_buffer(
        viewer,
        "ray.artifact.build_triangles.params",
        std::mem::size_of::<ArtifactConvertParams>() as u64,
        uniform_usage(),
        None,
    )?;

    for rep in &plan.triangle_reps {
        let grid = dispatch_grid_for(rep.triangle_count, max_dispatch_dimension)?;
        let params = ArtifactConvertParams {
            vertex_count: rep
                .rep
                .element_count
                .try_into()
                .map_err(|_| "artifact vertex count exceeds u32".to_string())?,
            triangle_offset: rep.triangle_offset,
            atom_offset: rep.rep.atom_offset,
            rep_slot: rep.rep_slot,
            transparency: rep.rep.transparency,
            dispatch_width: grid.invocation_width,
            _pad0: 0,
            _pad1: 0,
        };
        batch_commands.push(GpuBatchCommand::WriteBuffer {
            buffer: params_buffer,
            offset: 0,
            data: bytemuck::bytes_of(&params).to_vec(),
        });
        let bind_group = viewer.gpu_create_bind_group(GpuBindGroupDescriptor {
            label: Some("ray.artifact.build_triangles.bind_group".to_string()),
            layout,
            entries: vec![
                buffer_entry(0, rep.rep.geometry, 0, None),
                buffer_entry(1, color_lut, 0, None),
                buffer_entry(2, triangle_buffer, 0, None),
                buffer_entry(
                    3,
                    params_buffer,
                    0,
                    Some(std::mem::size_of::<ArtifactConvertParams>() as u64),
                ),
            ],
        })?;
        batch_commands.push(GpuBatchCommand::DispatchCompute {
            pipeline,
            bind_groups: vec![bind_group],
            workgroups: grid.workgroups,
        });
    }
    Ok(())
}

fn build_bvh(
    viewer: &mut dyn ViewerLike,
    triangle_buffer: GpuHandle,
    bvh_node_buffer: GpuHandle,
    bvh_index_buffer: GpuHandle,
    triangle_count: u32,
    shape: &BvhShape,
    max_dispatch_dimension: u32,
    batch_commands: &mut Vec<GpuBatchCommand>,
) -> Result<(), String> {
    let shader = viewer
        .gpu_create_cached_shader_module(GpuShaderModuleDescriptor {
            label: Some("ray.artifact.bvh.shader".to_string()),
            wgsl: ARTIFACT_BVH_SHADER.to_string(),
        })?
        .handle;
    let layout = viewer
        .gpu_create_cached_bind_group_layout(GpuBindGroupLayoutDescriptor {
            label: Some("ray.artifact.bvh.layout".to_string()),
            entries: vec![
                storage_layout(0, GpuBufferBindingType::StorageReadOnly),
                storage_layout(1, GpuBufferBindingType::StorageReadWrite),
                storage_layout(2, GpuBufferBindingType::StorageReadWrite),
                storage_layout(3, GpuBufferBindingType::Uniform),
            ],
        })?
        .handle;
    let pipeline_layout = viewer
        .gpu_create_cached_pipeline_layout(GpuPipelineLayoutDescriptor {
            label: Some("ray.artifact.bvh.pipeline_layout".to_string()),
            bind_group_layouts: vec![layout],
        })?
        .handle;
    let leaf_pipeline = viewer
        .gpu_create_cached_compute_pipeline(GpuComputePipelineDescriptor {
            label: Some("ray.artifact.bvh.leaves.pipeline".to_string()),
            layout: pipeline_layout,
            module: shader,
            entry_point: "build_leaves".to_string(),
        })?
        .handle;
    let internal_pipeline = viewer
        .gpu_create_cached_compute_pipeline(GpuComputePipelineDescriptor {
            label: Some("ray.artifact.bvh.internal.pipeline".to_string()),
            layout: pipeline_layout,
            module: shader,
            entry_point: "build_internal".to_string(),
        })?
        .handle;
    let params_buffer = create_buffer(
        viewer,
        "ray.artifact.bvh.params",
        std::mem::size_of::<ArtifactBvhParams>() as u64,
        uniform_usage(),
        None,
    )?;
    let bind_group = viewer.gpu_create_bind_group(GpuBindGroupDescriptor {
        label: Some("ray.artifact.bvh.bind_group".to_string()),
        layout,
        entries: vec![
            buffer_entry(0, triangle_buffer, 0, None),
            buffer_entry(1, bvh_node_buffer, 0, None),
            buffer_entry(2, bvh_index_buffer, 0, None),
            buffer_entry(
                3,
                params_buffer,
                0,
                Some(std::mem::size_of::<ArtifactBvhParams>() as u64),
            ),
        ],
    })?;

    let mut params = ArtifactBvhParams {
        triangle_count,
        leaf_slots: shape.leaf_slots,
        leaf_start: shape.leaf_start,
        level_start: shape.leaf_start,
        level_count: shape.leaf_slots,
        dispatch_width: 0,
        _pad0: 0,
        _pad1: 0,
    };
    let leaf_grid = dispatch_grid_for(shape.leaf_slots, max_dispatch_dimension)?;
    params.dispatch_width = leaf_grid.invocation_width;
    batch_commands.push(GpuBatchCommand::WriteBuffer {
        buffer: params_buffer,
        offset: 0,
        data: bytemuck::bytes_of(&params).to_vec(),
    });
    batch_commands.push(GpuBatchCommand::DispatchCompute {
        pipeline: leaf_pipeline,
        bind_groups: vec![bind_group],
        workgroups: leaf_grid.workgroups,
    });

    let mut child_level_start = shape.leaf_start;
    let mut child_level_count = shape.leaf_slots;
    while child_level_count > 1 {
        let parent_count = child_level_count / 2;
        let grid = dispatch_grid_for(parent_count, max_dispatch_dimension)?;
        params.level_start = child_level_start;
        params.level_count = parent_count;
        params.dispatch_width = grid.invocation_width;
        batch_commands.push(GpuBatchCommand::WriteBuffer {
            buffer: params_buffer,
            offset: 0,
            data: bytemuck::bytes_of(&params).to_vec(),
        });
        batch_commands.push(GpuBatchCommand::DispatchCompute {
            pipeline: internal_pipeline,
            bind_groups: vec![bind_group],
            workgroups: grid.workgroups,
        });
        child_level_start -= parent_count;
        child_level_count = parent_count;
    }

    Ok(())
}

fn dispatch_raytrace(
    viewer: &mut dyn ViewerLike,
    params: &RaytraceParams,
    triangle_buffer: GpuHandle,
    bvh_node_buffer: GpuHandle,
    bvh_index_buffer: GpuHandle,
    triangle_count: u32,
    bvh_node_count: u32,
    max_dispatch_dimension: u32,
    batch_commands: &mut Vec<GpuBatchCommand>,
) -> Result<Vec<u8>, String> {
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

    let sphere_buffer = create_buffer(
        viewer,
        "ray.artifact.empty_spheres",
        std::mem::size_of::<GpuSphere>() as u64,
        storage_usage(),
        None,
    )?;
    let cylinder_buffer = create_buffer(
        viewer,
        "ray.artifact.empty_cylinders",
        std::mem::size_of::<GpuCylinder>() as u64,
        storage_usage(),
        None,
    )?;
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
                0,
                0,
                triangle_count,
                bvh_node_count,
            ))
            .to_vec(),
        ),
    )?;

    let shader = viewer
        .gpu_create_cached_shader_module(GpuShaderModuleDescriptor {
            label: Some("ray.artifact.raytrace.shader".to_string()),
            wgsl: raytrace_output_buffer_shader(),
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
                storage_layout(6, GpuBufferBindingType::StorageReadWrite),
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
            buffer_entry(1, sphere_buffer, 0, None),
            buffer_entry(2, cylinder_buffer, 0, None),
            buffer_entry(3, triangle_buffer, 0, None),
            buffer_entry(4, bvh_node_buffer, 0, None),
            buffer_entry(5, bvh_index_buffer, 0, None),
            buffer_entry(6, output_buffer, 0, None),
        ],
    })?;

    let ray_workgroups = [render_width.div_ceil(8), render_height.div_ceil(8), 1];
    validate_workgroups(ray_workgroups, max_dispatch_dimension, "ray dispatch")?;
    batch_commands.push(GpuBatchCommand::DispatchCompute {
        pipeline,
        bind_groups: vec![bind_group],
        workgroups: ray_workgroups,
    });

    let (final_buffer, downsample_mode) = if supersample > 1 {
        (
            append_downsample(
                viewer,
                output_buffer,
                render_width,
                params.width,
                params.height,
                supersample,
                max_dispatch_dimension,
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
    source_buffer: GpuHandle,
    src_width: u32,
    dst_width: u32,
    dst_height: u32,
    factor: u32,
    max_dispatch_dimension: u32,
    batch_commands: &mut Vec<GpuBatchCommand>,
) -> Result<GpuHandle, String> {
    let final_bytes = rgba_bytes(dst_width, dst_height, "ray downsample output")?;
    let output_buffer = create_buffer(
        viewer,
        "ray.artifact.downsample_rgba",
        final_bytes,
        storage_usage(),
        None,
    )?;
    let params = DownsampleParams {
        src_width,
        dst_width,
        dst_height,
        factor,
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
            wgsl: DOWNSAMPLE_SHADER.to_string(),
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
            buffer_entry(0, source_buffer, 0, None),
            buffer_entry(1, output_buffer, 0, None),
            buffer_entry(
                2,
                params_buffer,
                0,
                Some(std::mem::size_of::<DownsampleParams>() as u64),
            ),
        ],
    })?;
    let workgroups = [dst_width.div_ceil(8), dst_height.div_ceil(8), 1];
    validate_workgroups(
        workgroups,
        max_dispatch_dimension,
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

fn bvh_shape_for(triangle_count: u32) -> Result<BvhShape, String> {
    let leaf_count = triangle_count.div_ceil(BVH_LEAF_SIZE).max(1);
    let leaf_slots = leaf_count.next_power_of_two();
    let node_count = leaf_slots
        .checked_mul(2)
        .and_then(|value| value.checked_sub(1))
        .ok_or_else(|| "BVH node count overflow".to_string())?;
    Ok(BvhShape {
        leaf_slots,
        leaf_start: leaf_slots - 1,
        node_count,
    })
}

fn workgroups_for(items: u32) -> u32 {
    items.max(1).div_ceil(WORKGROUP_SIZE)
}

fn dispatch_grid_for(
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

fn validate_workgroups(
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

fn create_buffer(
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

fn storage_usage() -> GpuBufferUsage {
    usage(&[
        GpuBufferUsage::STORAGE,
        GpuBufferUsage::COPY_SRC,
        GpuBufferUsage::COPY_DST,
    ])
}

fn uniform_usage() -> GpuBufferUsage {
    usage(&[GpuBufferUsage::UNIFORM, GpuBufferUsage::COPY_DST])
}

fn usage(parts: &[GpuBufferUsage]) -> GpuBufferUsage {
    parts
        .iter()
        .copied()
        .fold(GpuBufferUsage { bits: 0 }, GpuBufferUsage::union)
}

fn storage_layout(binding: u32, ty: GpuBufferBindingType) -> GpuBindGroupLayoutEntry {
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

fn buffer_entry(
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

fn raytrace_output_buffer_shader() -> String {
    include_str!("shaders/raytrace.wgsl")
        .replace(
            "@group(0) @binding(6) var output: texture_storage_2d<rgba8unorm, write>;",
            "@group(0) @binding(6) var<storage, read_write> output_pixels: array<u32>;",
        )
        .replace(
            "@group(0) @binding(7) var depth_output: texture_storage_2d<r32float, write>;\n",
            "",
        )
        .replace(
            "@group(0) @binding(8) var normal_output: texture_storage_2d<rgba16float, write>;\n",
            "",
        )
        .replace(
            "let node = bvh_nodes[node_idx];\n        \n        // Test AABB",
            "let node = bvh_nodes[node_idx];\n        if node.left_or_first == 0xffffffffu && node.count == 0u {\n            continue;\n        }\n        \n        // Test AABB",
        )
        .replace(
            "textureStore(output, pixel_coord, vec4<f32>(final_color, accumulated_alpha));",
            "let output_index = global_id.y * u32(uniforms.viewport.x) + global_id.x;\n    output_pixels[output_index] = pack4x8unorm(vec4<f32>(final_color, accumulated_alpha));",
        )
        .replace(
            "    // Output depth and normal for edge detection pass\n    textureStore(depth_output, pixel_coord, vec4<f32>(depth, 0.0, 0.0, 1.0));\n    textureStore(normal_output, pixel_coord, vec4<f32>(normal, 1.0));\n",
            "",
        )
}

fn rt_profile_enabled() -> bool {
    std::env::var_os("PATINAE_RT_PROFILE").is_some()
}

const DOWNSAMPLE_SHADER: &str = r#"
struct DownsampleParams {
    src_width: u32,
    dst_width: u32,
    dst_height: u32,
    factor: u32,
}

@group(0) @binding(0) var<storage, read> src_pixels: array<u32>;
@group(0) @binding(1) var<storage, read_write> dst_pixels: array<u32>;
@group(0) @binding(2) var<uniform> params: DownsampleParams;

@compute @workgroup_size(8, 8, 1)
fn main(@builtin(global_invocation_id) global_id: vec3<u32>) {
    if global_id.x >= params.dst_width || global_id.y >= params.dst_height {
        return;
    }

    var rgba = vec4<f32>(0.0);
    for (var sy = 0u; sy < params.factor; sy = sy + 1u) {
        for (var sx = 0u; sx < params.factor; sx = sx + 1u) {
            let src_x = global_id.x * params.factor + sx;
            let src_y = global_id.y * params.factor + sy;
            rgba = rgba + unpack4x8unorm(src_pixels[src_y * params.src_width + src_x]);
        }
    }

    let scale = 1.0 / f32(params.factor * params.factor);
    let dst_index = global_id.y * params.dst_width + global_id.x;
    dst_pixels[dst_index] = pack4x8unorm(rgba * scale);
}
"#;

const ARTIFACT_TRIANGLE_SHADER: &str = r#"
struct ConvertParams {
    vertex_count: u32,
    triangle_offset: u32,
    atom_offset: u32,
    rep_slot: u32,
    transparency: f32,
    dispatch_width: u32,
    _pad0: u32,
    _pad1: u32,
}

struct RepColorLutEntry {
    sphere: u32,
    stick: u32,
    line: u32,
    dot: u32,
    cartoon: u32,
    ribbon: u32,
    surface: u32,
    mesh: u32,
    ellipsoid: u32,
    _pad0: u32,
    _pad1: u32,
    _pad2: u32,
}

struct ColorLutEntry {
    base: vec4<f32>,
    reps: RepColorLutEntry,
}

struct Triangle {
    v0: vec3<f32>, _pad0: f32,
    v1: vec3<f32>, _pad1: f32,
    v2: vec3<f32>, _pad2: f32,
    n0: vec3<f32>, _pad3: f32,
    n1: vec3<f32>, _pad4: f32,
    n2: vec3<f32>, _pad5: f32,
    color: vec4<f32>,
    transparency: f32,
    _pad_a: f32, _pad_b: f32, _pad_c: f32,
}

@group(0) @binding(0) var<storage, read> std_vertex_words: array<u32>;
@group(0) @binding(1) var<storage, read> color_lut: array<ColorLutEntry>;
@group(0) @binding(2) var<storage, read_write> triangles: array<Triangle>;
@group(0) @binding(3) var<uniform> params: ConvertParams;

const REP_COLOR_INHERIT: u32 = 0xffffffffu;

fn vertex_word(vertex_index: u32, lane: u32) -> u32 {
    return std_vertex_words[vertex_index * 6u + lane];
}

fn vertex_position(vertex_index: u32) -> vec3<f32> {
    return vec3<f32>(
        bitcast<f32>(vertex_word(vertex_index, 0u)),
        bitcast<f32>(vertex_word(vertex_index, 1u)),
        bitcast<f32>(vertex_word(vertex_index, 2u)),
    );
}

fn vertex_normal_oct(vertex_index: u32) -> u32 {
    return vertex_word(vertex_index, 3u);
}

fn vertex_group_id(vertex_index: u32) -> u32 {
    return vertex_word(vertex_index, 4u);
}

fn sign_extend_16(value: u32) -> i32 {
    let low = i32(value & 0xffffu);
    if low >= 32768 {
        return low - 65536;
    }
    return low;
}

fn oct_decode(packed: u32) -> vec3<f32> {
    var x = max(f32(sign_extend_16(packed)) / 32767.0, -1.0);
    var y = max(f32(sign_extend_16(packed >> 16u)) / 32767.0, -1.0);
    let z = 1.0 - abs(x) - abs(y);
    var out: vec3<f32>;
    if z < 0.0 {
        let ox = (1.0 - abs(y)) * select(-1.0, 1.0, x >= 0.0);
        let oy = (1.0 - abs(x)) * select(-1.0, 1.0, y >= 0.0);
        out = vec3<f32>(ox, oy, z);
    } else {
        out = vec3<f32>(x, y, z);
    }
    return normalize(out);
}

fn rep_packed_color(reps: RepColorLutEntry, slot: u32) -> u32 {
    switch slot {
        case 0u: { return reps.sphere; }
        case 1u: { return reps.stick; }
        case 2u: { return reps.line; }
        case 3u: { return reps.dot; }
        case 4u: { return reps.cartoon; }
        case 5u: { return reps.ribbon; }
        case 6u: { return reps.surface; }
        case 7u: { return reps.mesh; }
        case 8u: { return reps.ellipsoid; }
        default: { return REP_COLOR_INHERIT; }
    }
}

fn unpack_rep_rgb8(packed: u32) -> vec4<f32> {
    return vec4<f32>(
        f32(packed & 0xffu) / 255.0,
        f32((packed >> 8u) & 0xffu) / 255.0,
        f32((packed >> 16u) & 0xffu) / 255.0,
        1.0,
    );
}

fn material_for_group(group_id: u32) -> vec4<f32> {
    let entry = color_lut[params.atom_offset + group_id];
    let packed = rep_packed_color(entry.reps, params.rep_slot);
    var rgba = select(unpack_rep_rgb8(packed), entry.base, packed == REP_COLOR_INHERIT);
    rgba.a = clamp(rgba.a * (1.0 - params.transparency), 0.0, 1.0);
    return rgba;
}

@compute @workgroup_size(128)
fn build_triangles(@builtin(global_invocation_id) gid: vec3<u32>) {
    let triangle_index = gid.y * params.dispatch_width + gid.x;
    let triangle_count = params.vertex_count / 3u;
    if triangle_index >= triangle_count {
        return;
    }

    let base_vertex = triangle_index * 3u;
    let v0 = base_vertex;
    let v1 = base_vertex + 1u;
    let v2 = base_vertex + 2u;
    let out_index = params.triangle_offset + triangle_index;
    let material = material_for_group(vertex_group_id(v0));

    triangles[out_index] = Triangle(
        vertex_position(v0), 0.0,
        vertex_position(v1), 0.0,
        vertex_position(v2), 0.0,
        oct_decode(vertex_normal_oct(v0)), 0.0,
        oct_decode(vertex_normal_oct(v1)), 0.0,
        oct_decode(vertex_normal_oct(v2)), 0.0,
        material,
        1.0 - material.a,
        0.0, 0.0, 0.0,
    );
}
"#;

const ARTIFACT_BVH_SHADER: &str = r#"
struct BvhParams {
    triangle_count: u32,
    leaf_slots: u32,
    leaf_start: u32,
    level_start: u32,
    level_count: u32,
    dispatch_width: u32,
    _pad0: u32,
    _pad1: u32,
}

struct Triangle {
    v0: vec3<f32>, _pad0: f32,
    v1: vec3<f32>, _pad1: f32,
    v2: vec3<f32>, _pad2: f32,
    n0: vec3<f32>, _pad3: f32,
    n1: vec3<f32>, _pad4: f32,
    n2: vec3<f32>, _pad5: f32,
    color: vec4<f32>,
    transparency: f32,
    _pad_a: f32, _pad_b: f32, _pad_c: f32,
}

struct BvhNode {
    min: vec3<f32>,
    left_or_first: u32,
    max: vec3<f32>,
    count: u32,
}

@group(0) @binding(0) var<storage, read> triangles: array<Triangle>;
@group(0) @binding(1) var<storage, read_write> bvh_nodes: array<BvhNode>;
@group(0) @binding(2) var<storage, read_write> bvh_indices: array<u32>;
@group(0) @binding(3) var<uniform> params: BvhParams;

const LEAF_SIZE: u32 = 4u;
const INF: f32 = 1.0e30;
const EMPTY_NODE: u32 = 0xffffffffu;

fn empty_node() -> BvhNode {
    return BvhNode(vec3<f32>(0.0), EMPTY_NODE, vec3<f32>(0.0), 0u);
}

fn is_empty_node(node: BvhNode) -> bool {
    return node.left_or_first == EMPTY_NODE && node.count == 0u;
}

fn grow_min(current: vec3<f32>, value: vec3<f32>) -> vec3<f32> {
    return min(current, value);
}

fn grow_max(current: vec3<f32>, value: vec3<f32>) -> vec3<f32> {
    return max(current, value);
}

@compute @workgroup_size(128)
fn build_leaves(@builtin(global_invocation_id) gid: vec3<u32>) {
    let leaf_id = gid.y * params.dispatch_width + gid.x;
    if leaf_id >= params.leaf_slots {
        return;
    }

    let node_index = params.leaf_start + leaf_id;
    let first = leaf_id * LEAF_SIZE;
    if first >= params.triangle_count {
        bvh_nodes[node_index] = empty_node();
        return;
    }

    let count = min(LEAF_SIZE, params.triangle_count - first);
    var bmin = vec3<f32>(INF);
    var bmax = vec3<f32>(-INF);
    for (var i = 0u; i < LEAF_SIZE; i = i + 1u) {
        if i < count {
            let triangle_index = first + i;
            let tri = triangles[triangle_index];
            bvh_indices[triangle_index] = (2u << 30u) | triangle_index;
            bmin = grow_min(bmin, min(tri.v0, min(tri.v1, tri.v2)));
            bmax = grow_max(bmax, max(tri.v0, max(tri.v1, tri.v2)));
        }
    }
    bvh_nodes[node_index] = BvhNode(bmin, first, bmax, count);
}

@compute @workgroup_size(128)
fn build_internal(@builtin(global_invocation_id) gid: vec3<u32>) {
    let local_id = gid.y * params.dispatch_width + gid.x;
    if local_id >= params.level_count {
        return;
    }
    let parent_start = params.level_start - params.level_count;
    let parent_index = parent_start + local_id;
    let left_index = params.level_start + local_id * 2u;
    let right_index = left_index + 1u;
    let left = bvh_nodes[left_index];
    let right = bvh_nodes[right_index];
    if is_empty_node(left) && is_empty_node(right) {
        bvh_nodes[parent_index] = empty_node();
        return;
    }
    if is_empty_node(left) {
        bvh_nodes[parent_index] = right;
        return;
    }
    if is_empty_node(right) {
        bvh_nodes[parent_index] = left;
        return;
    }
    let bmin = min(left.min, right.min);
    let bmax = max(left.max, right.max);
    bvh_nodes[parent_index] = BvhNode(bmin, left_index, bmax, 0u);
}
"#;

#[cfg(test)]
mod tests {
    use super::*;
    use patinae_scene::GpuHandleKind;

    fn handle(id: u64) -> GpuHandle {
        GpuHandle {
            id,
            kind: GpuHandleKind::Buffer,
            generation: 1,
        }
    }

    #[test]
    fn bvh_shape_pads_leaves_to_power_of_two() {
        let shape = bvh_shape_for(9).expect("shape");

        assert_eq!(shape.leaf_slots, 4);
        assert_eq!(shape.leaf_start, 3);
        assert_eq!(shape.node_count, 7);
    }

    #[test]
    fn dispatch_grid_tiles_past_single_dimension_limit() {
        let items = 90_473 * WORKGROUP_SIZE;
        let grid = dispatch_grid_for(items, 65_535).expect("dispatch grid");

        assert_eq!(grid.workgroups, [65_535, 2, 1]);
        assert_eq!(grid.invocation_width, 65_535 * WORKGROUP_SIZE);
    }

    #[test]
    fn raytrace_shader_rewrite_uses_output_buffer_and_empty_node_guard() {
        let shader = raytrace_output_buffer_shader();

        assert!(shader.contains("output_pixels"));
        assert!(shader.contains("global_id.y * u32(uniforms.viewport.x) + global_id.x"));
        assert!(shader.contains("node.left_or_first == 0xffffffffu"));
        assert!(!shader.contains("texture_storage_2d"));
        assert!(!shader.contains("textureStore("));
    }

    #[test]
    fn artifact_wgsl_modules_parse_and_validate() {
        let modules = [
            (
                "ray.artifact.triangles",
                ARTIFACT_TRIANGLE_SHADER.to_string(),
            ),
            ("ray.artifact.bvh", ARTIFACT_BVH_SHADER.to_string()),
            ("ray.artifact.raytrace", raytrace_output_buffer_shader()),
            ("ray.artifact.downsample", DOWNSAMPLE_SHADER.to_string()),
        ];

        let mut failures = Vec::new();

        for (label, source) in modules {
            let module = match naga::front::wgsl::parse_str(&source) {
                Ok(module) => module,
                Err(err) => {
                    failures.push(format!(
                        "{label}: WGSL parse failed\n{}",
                        err.emit_to_string_with_path(&source, label)
                    ));
                    continue;
                }
            };

            let mut validator = naga::valid::Validator::new(
                naga::valid::ValidationFlags::all(),
                naga::valid::Capabilities::all(),
            );
            if let Err(err) = validator.validate(&module) {
                failures.push(format!("{label}: WGSL validation failed\n{err:#?}"));
            }
        }

        assert!(
            failures.is_empty(),
            "invalid ray artifact WGSL modules:\n{}",
            failures.join("\n\n")
        );
    }

    #[test]
    fn artifact_plan_uses_scene_color_lut_and_cartoon_slot() {
        let snapshot = RenderArtifactSnapshotDescriptor {
            snapshot_id: 1,
            layout_version: 2,
            scene_generation: 7,
            scene_bounds_min: [0.0, 0.0, 0.0],
            scene_bounds_max: [1.0, 1.0, 1.0],
            cull_pass_initialized: true,
            device_limits: patinae_scene::GpuDeviceLimits {
                max_buffer_size: 1024,
                max_storage_buffer_binding_size: 1024,
                max_compute_workgroups_per_dimension: 65_535,
                max_compute_invocations_per_workgroup: 256,
                max_compute_workgroup_size_x: 256,
                max_compute_workgroup_size_y: 1,
                max_compute_workgroup_size_z: 1,
            },
            buffers: vec![
                RenderArtifactBufferDescriptor {
                    handle: handle(10),
                    role: RenderArtifactBufferRole::SceneColorLut,
                    size: 64,
                    stride: COLOR_LUT_STRIDE,
                    element_count: 1,
                },
                RenderArtifactBufferDescriptor {
                    handle: handle(11),
                    role: RenderArtifactBufferRole::StdVertices,
                    size: 72,
                    stride: 24,
                    element_count: 3,
                },
            ],
            reps: vec![RenderArtifactRepDescriptor {
                object_id: 1,
                rep_kind: RenderArtifactRepKind::Cartoon,
                topology: RenderArtifactPrimitiveTopology::TriangleList,
                geometry: handle(11),
                count: None,
                indirect: None,
                element_count: 3,
                max_element_count: 3,
                atom_offset: 5,
                atom_count: 1,
                material_rgba: [1.0, 0.0, 0.0, 1.0],
                transparency: 0.0,
            }],
        };

        let plan = plan_artifact_triangles(&snapshot).expect("plan");

        assert_eq!(plan.color_lut, handle(10));
        assert_eq!(plan.triangle_count, 1);
        assert_eq!(plan.triangle_reps[0].rep_slot, 4);
    }
}
