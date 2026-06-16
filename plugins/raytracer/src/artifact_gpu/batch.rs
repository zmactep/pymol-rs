use patinae_plugin::prelude::ViewerLike;
use patinae_scene::{
    GpuBatchCommand, GpuBindGroupDescriptor, GpuBindGroupLayoutDescriptor, GpuBufferBindingType,
    GpuComputePipelineDescriptor, GpuHandle, GpuPipelineLayoutDescriptor,
    GpuShaderModuleDescriptor,
};

use super::dispatch::dispatch_grid_for;
use super::resources::{buffer_entry, create_buffer, storage_layout, uniform_usage};
use super::{
    ArtifactCapsuleParams, ArtifactCylinderParams, ArtifactPlan, ArtifactSphereParams,
    ArtifactTriangleParams,
};
use crate::shaders;

pub(super) fn build_spheres(
    viewer: &mut dyn ViewerLike,
    plan: &ArtifactPlan<'_>,
    color_lut: GpuHandle,
    sphere_buffer: GpuHandle,
    max_dispatch_dimension: u32,
    batch_commands: &mut Vec<GpuBatchCommand>,
) -> Result<(), String> {
    if plan.sphere_reps.is_empty() {
        return Ok(());
    }
    let shader = viewer
        .gpu_create_cached_shader_module(GpuShaderModuleDescriptor {
            label: Some("ray.artifact.build_spheres.shader".to_string()),
            wgsl: shaders::artifact_spheres(),
        })?
        .handle;
    let layout = viewer
        .gpu_create_cached_bind_group_layout(GpuBindGroupLayoutDescriptor {
            label: Some("ray.artifact.build_spheres.layout".to_string()),
            entries: vec![
                storage_layout(0, GpuBufferBindingType::StorageReadOnly),
                storage_layout(1, GpuBufferBindingType::StorageReadOnly),
                storage_layout(2, GpuBufferBindingType::StorageReadOnly),
                storage_layout(3, GpuBufferBindingType::StorageReadWrite),
                storage_layout(4, GpuBufferBindingType::Uniform),
            ],
        })?
        .handle;
    let pipeline_layout = viewer
        .gpu_create_cached_pipeline_layout(GpuPipelineLayoutDescriptor {
            label: Some("ray.artifact.build_spheres.pipeline_layout".to_string()),
            bind_group_layouts: vec![layout],
        })?
        .handle;
    let pipeline = viewer
        .gpu_create_cached_compute_pipeline(GpuComputePipelineDescriptor {
            label: Some("ray.artifact.build_spheres.pipeline".to_string()),
            layout: pipeline_layout,
            module: shader,
            entry_point: "build_spheres".to_string(),
        })?
        .handle;
    let params_buffer = create_buffer(
        viewer,
        "ray.artifact.build_spheres.params",
        std::mem::size_of::<ArtifactSphereParams>() as u64,
        uniform_usage(),
        None,
    )?;

    for rep in &plan.sphere_reps {
        let grid = dispatch_grid_for(rep.instance_capacity, max_dispatch_dimension)?;
        let draw_args = rep.indirect_or_direct_args(viewer, "ray.artifact.sphere.direct_args")?;
        let params = ArtifactSphereParams {
            instance_capacity: rep.instance_capacity,
            sphere_offset: rep.sphere_offset,
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
            label: Some("ray.artifact.build_spheres.bind_group".to_string()),
            layout,
            entries: vec![
                buffer_entry(0, rep.rep.geometry, 0, Some(rep.geometry_binding_size)),
                buffer_entry(1, color_lut, 0, None),
                buffer_entry(2, draw_args, 0, None),
                buffer_entry(3, sphere_buffer, 0, None),
                buffer_entry(
                    4,
                    params_buffer,
                    0,
                    Some(std::mem::size_of::<ArtifactSphereParams>() as u64),
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

pub(super) fn build_cylinders(
    viewer: &mut dyn ViewerLike,
    plan: &ArtifactPlan<'_>,
    color_lut: GpuHandle,
    cylinder_buffer: GpuHandle,
    max_dispatch_dimension: u32,
    batch_commands: &mut Vec<GpuBatchCommand>,
) -> Result<(), String> {
    if plan.cylinder_reps.is_empty() {
        return Ok(());
    }
    let line_shader = viewer
        .gpu_create_cached_shader_module(GpuShaderModuleDescriptor {
            label: Some("ray.artifact.build_line_cylinders.shader".to_string()),
            wgsl: shaders::artifact_line_cylinders(),
        })?
        .handle;
    let layout = viewer
        .gpu_create_cached_bind_group_layout(GpuBindGroupLayoutDescriptor {
            label: Some("ray.artifact.build_cylinders.layout".to_string()),
            entries: vec![
                storage_layout(0, GpuBufferBindingType::StorageReadOnly),
                storage_layout(1, GpuBufferBindingType::StorageReadOnly),
                storage_layout(2, GpuBufferBindingType::StorageReadOnly),
                storage_layout(3, GpuBufferBindingType::StorageReadWrite),
                storage_layout(4, GpuBufferBindingType::Uniform),
            ],
        })?
        .handle;
    let pipeline_layout = viewer
        .gpu_create_cached_pipeline_layout(GpuPipelineLayoutDescriptor {
            label: Some("ray.artifact.build_cylinders.pipeline_layout".to_string()),
            bind_group_layouts: vec![layout],
        })?
        .handle;
    let line_pipeline = viewer
        .gpu_create_cached_compute_pipeline(GpuComputePipelineDescriptor {
            label: Some("ray.artifact.build_line_cylinders.pipeline".to_string()),
            layout: pipeline_layout,
            module: line_shader,
            entry_point: "build_cylinders".to_string(),
        })?
        .handle;
    let params_buffer = create_buffer(
        viewer,
        "ray.artifact.build_cylinders.params",
        std::mem::size_of::<ArtifactCylinderParams>() as u64,
        uniform_usage(),
        None,
    )?;

    for rep in &plan.cylinder_reps {
        let grid = dispatch_grid_for(rep.instance_capacity, max_dispatch_dimension)?;
        let draw_args = rep.indirect_or_direct_args(viewer, "ray.artifact.cylinder.direct_args")?;
        let params = ArtifactCylinderParams {
            instance_capacity: rep.instance_capacity,
            cylinder_offset: rep.cylinder_offset,
            atom_offset: rep.rep.atom_offset,
            rep_slot: rep.rep_slot,
            transparency: rep.rep.transparency,
            radius: rep.radius,
            dispatch_width: grid.invocation_width,
            _pad0: 0,
        };
        batch_commands.push(GpuBatchCommand::WriteBuffer {
            buffer: params_buffer,
            offset: 0,
            data: bytemuck::bytes_of(&params).to_vec(),
        });
        let bind_group = viewer.gpu_create_bind_group(GpuBindGroupDescriptor {
            label: Some("ray.artifact.build_cylinders.bind_group".to_string()),
            layout,
            entries: vec![
                buffer_entry(0, rep.rep.geometry, 0, Some(rep.geometry_binding_size)),
                buffer_entry(1, color_lut, 0, None),
                buffer_entry(2, draw_args, 0, None),
                buffer_entry(3, cylinder_buffer, 0, None),
                buffer_entry(
                    4,
                    params_buffer,
                    0,
                    Some(std::mem::size_of::<ArtifactCylinderParams>() as u64),
                ),
            ],
        })?;
        batch_commands.push(GpuBatchCommand::DispatchCompute {
            pipeline: line_pipeline,
            bind_groups: vec![bind_group],
            workgroups: grid.workgroups,
        });
    }
    Ok(())
}

pub(super) fn build_capsules(
    viewer: &mut dyn ViewerLike,
    plan: &ArtifactPlan<'_>,
    color_lut: GpuHandle,
    capsule_buffer: GpuHandle,
    max_dispatch_dimension: u32,
    batch_commands: &mut Vec<GpuBatchCommand>,
) -> Result<(), String> {
    if plan.capsule_reps.is_empty() {
        return Ok(());
    }
    let shader = viewer
        .gpu_create_cached_shader_module(GpuShaderModuleDescriptor {
            label: Some("ray.artifact.build_stick_capsules.shader".to_string()),
            wgsl: shaders::artifact_stick_capsules(),
        })?
        .handle;
    let layout = viewer
        .gpu_create_cached_bind_group_layout(GpuBindGroupLayoutDescriptor {
            label: Some("ray.artifact.build_capsules.layout".to_string()),
            entries: vec![
                storage_layout(0, GpuBufferBindingType::StorageReadOnly),
                storage_layout(1, GpuBufferBindingType::StorageReadOnly),
                storage_layout(2, GpuBufferBindingType::StorageReadOnly),
                storage_layout(3, GpuBufferBindingType::StorageReadWrite),
                storage_layout(4, GpuBufferBindingType::Uniform),
            ],
        })?
        .handle;
    let pipeline_layout = viewer
        .gpu_create_cached_pipeline_layout(GpuPipelineLayoutDescriptor {
            label: Some("ray.artifact.build_capsules.pipeline_layout".to_string()),
            bind_group_layouts: vec![layout],
        })?
        .handle;
    let pipeline = viewer
        .gpu_create_cached_compute_pipeline(GpuComputePipelineDescriptor {
            label: Some("ray.artifact.build_stick_capsules.pipeline".to_string()),
            layout: pipeline_layout,
            module: shader,
            entry_point: "build_capsules".to_string(),
        })?
        .handle;
    let params_buffer = create_buffer(
        viewer,
        "ray.artifact.build_capsules.params",
        std::mem::size_of::<ArtifactCapsuleParams>() as u64,
        uniform_usage(),
        None,
    )?;

    for rep in &plan.capsule_reps {
        let grid = dispatch_grid_for(rep.instance_capacity, max_dispatch_dimension)?;
        let draw_args = rep.indirect_or_direct_args(viewer, "ray.artifact.capsule.direct_args")?;
        let params = ArtifactCapsuleParams {
            instance_capacity: rep.instance_capacity,
            capsule_offset: rep.capsule_offset,
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
            label: Some("ray.artifact.build_capsules.bind_group".to_string()),
            layout,
            entries: vec![
                buffer_entry(0, rep.rep.geometry, 0, Some(rep.geometry_binding_size)),
                buffer_entry(1, color_lut, 0, None),
                buffer_entry(2, draw_args, 0, None),
                buffer_entry(3, capsule_buffer, 0, None),
                buffer_entry(
                    4,
                    params_buffer,
                    0,
                    Some(std::mem::size_of::<ArtifactCapsuleParams>() as u64),
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

// Consumes the public `StdVertex` packing documented by render artifacts:
// position lanes, packed normal, group id, and flags in six u32 words.
pub(super) fn build_triangles(
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
            wgsl: shaders::artifact_triangles(),
        })?
        .handle;
    let layout = viewer
        .gpu_create_cached_bind_group_layout(GpuBindGroupLayoutDescriptor {
            label: Some("ray.artifact.build_triangles.layout".to_string()),
            entries: vec![
                storage_layout(0, GpuBufferBindingType::StorageReadOnly),
                storage_layout(1, GpuBufferBindingType::StorageReadOnly),
                storage_layout(2, GpuBufferBindingType::StorageReadWrite),
                storage_layout(3, GpuBufferBindingType::StorageReadOnly),
                storage_layout(4, GpuBufferBindingType::Uniform),
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
        std::mem::size_of::<ArtifactTriangleParams>() as u64,
        uniform_usage(),
        None,
    )?;

    for rep in &plan.triangle_reps {
        let grid = dispatch_grid_for(rep.triangle_count, max_dispatch_dimension)?;
        let draw_args = rep.indirect_or_direct_args(viewer, "ray.artifact.triangle.direct_args")?;
        let params = ArtifactTriangleParams {
            vertex_capacity: rep.vertex_count,
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
                buffer_entry(0, rep.rep.geometry, 0, Some(rep.geometry_binding_size)),
                buffer_entry(1, color_lut, 0, None),
                buffer_entry(2, triangle_buffer, 0, None),
                buffer_entry(3, draw_args, 0, None),
                buffer_entry(
                    4,
                    params_buffer,
                    0,
                    Some(std::mem::size_of::<ArtifactTriangleParams>() as u64),
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
