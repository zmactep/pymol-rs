use patinae_plugin::prelude::ViewerLike;
use patinae_scene::{
    GpuBatchCommand, GpuBindGroupDescriptor, GpuBindGroupLayoutDescriptor, GpuBufferBindingType,
    GpuComputePipelineDescriptor, GpuPipelineLayoutDescriptor, GpuShaderModuleDescriptor,
};

use super::dispatch::dispatch_grid_for;
use super::layout::{ArtifactBvhParams, BvhBuildInput, BvhShape, BVH_LEAF_SIZE};
use super::resources::{buffer_entry, create_buffer, storage_layout, uniform_usage};
use crate::shaders;

pub(super) fn build_bvh(
    viewer: &mut dyn ViewerLike,
    input: BvhBuildInput<'_>,
    batch_commands: &mut Vec<GpuBatchCommand>,
) -> Result<(), String> {
    let shader = viewer
        .gpu_create_cached_shader_module(GpuShaderModuleDescriptor {
            label: Some("ray.artifact.bvh.shader".to_string()),
            wgsl: shaders::artifact_bvh(),
        })?
        .handle;
    let layout = viewer
        .gpu_create_cached_bind_group_layout(GpuBindGroupLayoutDescriptor {
            label: Some("ray.artifact.bvh.layout".to_string()),
            entries: vec![
                storage_layout(0, GpuBufferBindingType::StorageReadOnly),
                storage_layout(1, GpuBufferBindingType::StorageReadOnly),
                storage_layout(2, GpuBufferBindingType::StorageReadOnly),
                storage_layout(3, GpuBufferBindingType::StorageReadOnly),
                storage_layout(4, GpuBufferBindingType::StorageReadWrite),
                storage_layout(5, GpuBufferBindingType::StorageReadWrite),
                storage_layout(6, GpuBufferBindingType::Uniform),
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
            buffer_entry(0, input.primitives.spheres, 0, None),
            buffer_entry(1, input.primitives.cylinders, 0, None),
            buffer_entry(2, input.primitives.capsules, 0, None),
            buffer_entry(3, input.primitives.triangles, 0, None),
            buffer_entry(4, input.bvh.nodes, 0, None),
            buffer_entry(5, input.bvh.indices, 0, None),
            buffer_entry(
                6,
                params_buffer,
                0,
                Some(std::mem::size_of::<ArtifactBvhParams>() as u64),
            ),
        ],
    })?;

    let mut params = ArtifactBvhParams {
        primitive_count: input.primitive_count,
        sphere_count: input.counts.spheres,
        cylinder_count: input.counts.cylinders,
        capsule_count: input.counts.capsules,
        triangle_count: input.counts.triangles,
        leaf_slots: input.shape.leaf_slots,
        leaf_start: input.shape.leaf_start,
        level_start: input.shape.leaf_start,
        level_count: input.shape.leaf_slots,
        dispatch_width: 0,
        _pad0: 0,
        _pad1: 0,
    };
    let leaf_grid = dispatch_grid_for(input.shape.leaf_slots, input.max_dispatch_dimension)?;
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

    let mut child_level_start = input.shape.leaf_start;
    let mut child_level_count = input.shape.leaf_slots;
    while child_level_count > 1 {
        let parent_count = child_level_count / 2;
        let grid = dispatch_grid_for(parent_count, input.max_dispatch_dimension)?;
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

pub(super) fn bvh_shape_for(triangle_count: u32) -> Result<BvhShape, String> {
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
