use patinae_plugin::prelude::ViewerLike;
use patinae_scene::{
    GpuBatchCommand, GpuBindGroupDescriptor, GpuBindGroupLayoutDescriptor, GpuBufferBindingType,
    GpuComputePipelineDescriptor, GpuHandle, GpuPipelineLayoutDescriptor,
    GpuShaderModuleDescriptor,
};

use super::dispatch::dispatch_grid_for;
use super::layout::{
    ArtifactCapsuleParams, ArtifactCylinderParams, ArtifactSphereParams, ArtifactTriangleParams,
    ArtifactVisibleTriangleParams,
};
use super::reps::ArtifactPlan;
use super::resources::{buffer_entry, create_buffer, storage_layout, uniform_usage};
use crate::gpu::RaytraceParams;
use crate::shaders;

#[derive(Clone, Copy)]
struct TriangleBuildPipeline {
    pipeline: GpuHandle,
    layout: GpuHandle,
    params_buffer: GpuHandle,
}

#[derive(Clone, Copy)]
struct VisibleTriangleBuildPipeline {
    pipeline: GpuHandle,
    layout: GpuHandle,
    params_buffer: GpuHandle,
    counter_buffer: GpuHandle,
    counter_bytes: u64,
}

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
    visibility_counter_buffer: Option<GpuHandle>,
    ray_params: &RaytraceParams,
    max_dispatch_dimension: u32,
    batch_commands: &mut Vec<GpuBatchCommand>,
) -> Result<(), String> {
    if plan.triangle_reps.is_empty() {
        return Ok(());
    }
    let scene_atoms = triangle_scene_atoms(plan)?;
    let triangle_pipeline = create_triangle_build_pipeline(viewer)?;

    let visible_pipeline = if needs_visible_triangle_pipeline(plan) {
        let counter_buffer = visibility_counter_buffer.ok_or_else(|| {
            "surface triangle visibility counters are missing for visible surface reps".to_string()
        })?;
        Some(create_visible_triangle_build_pipeline(
            viewer,
            plan,
            counter_buffer,
            batch_commands,
        )?)
    } else {
        None
    };

    for rep in &plan.triangle_reps {
        if let Some(counter_index) = rep.visibility_counter_index {
            let Some(visible_pipeline) = visible_pipeline else {
                return Err(
                    "surface triangle visibility pipeline is missing for visible surface rep"
                        .to_string(),
                );
            };
            let grid = dispatch_grid_for(rep.source_triangle_count(), max_dispatch_dimension)?;
            let params = super::visibility::visible_triangle_params(
                ray_params,
                rep,
                rep.triangle_count,
                counter_index,
                grid.invocation_width,
            );
            batch_commands.push(GpuBatchCommand::WriteBuffer {
                buffer: visible_pipeline.params_buffer,
                offset: 0,
                data: bytemuck::bytes_of(&params).to_vec(),
            });
            let bind_group = viewer.gpu_create_bind_group(GpuBindGroupDescriptor {
                label: Some("ray.artifact.build_visible_triangles.bind_group".to_string()),
                layout: visible_pipeline.layout,
                entries: vec![
                    buffer_entry(0, rep.rep.geometry, 0, Some(rep.geometry_binding_size)),
                    buffer_entry(1, color_lut, 0, None),
                    buffer_entry(2, triangle_buffer, 0, None),
                    buffer_entry(3, scene_atoms, 0, None),
                    buffer_entry(
                        4,
                        visible_pipeline.counter_buffer,
                        0,
                        Some(visible_pipeline.counter_bytes),
                    ),
                    buffer_entry(
                        5,
                        visible_pipeline.params_buffer,
                        0,
                        Some(std::mem::size_of::<ArtifactVisibleTriangleParams>() as u64),
                    ),
                ],
            })?;
            batch_commands.push(GpuBatchCommand::DispatchCompute {
                pipeline: visible_pipeline.pipeline,
                bind_groups: vec![bind_group],
                workgroups: grid.workgroups,
            });
            continue;
        }

        let grid = dispatch_grid_for(rep.triangle_count, max_dispatch_dimension)?;
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
            buffer: triangle_pipeline.params_buffer,
            offset: 0,
            data: bytemuck::bytes_of(&params).to_vec(),
        });
        let bind_group = viewer.gpu_create_bind_group(GpuBindGroupDescriptor {
            label: Some("ray.artifact.build_triangles.bind_group".to_string()),
            layout: triangle_pipeline.layout,
            entries: vec![
                buffer_entry(0, rep.rep.geometry, 0, Some(rep.geometry_binding_size)),
                buffer_entry(1, color_lut, 0, None),
                buffer_entry(2, triangle_buffer, 0, None),
                buffer_entry(3, scene_atoms, 0, None),
                buffer_entry(
                    4,
                    triangle_pipeline.params_buffer,
                    0,
                    Some(std::mem::size_of::<ArtifactTriangleParams>() as u64),
                ),
            ],
        })?;
        batch_commands.push(GpuBatchCommand::DispatchCompute {
            pipeline: triangle_pipeline.pipeline,
            bind_groups: vec![bind_group],
            workgroups: grid.workgroups,
        });
    }
    Ok(())
}

fn triangle_scene_atoms(plan: &ArtifactPlan<'_>) -> Result<GpuHandle, String> {
    plan.scene_atoms.ok_or_else(|| {
        "render artifact snapshot is missing SceneAtoms for triangle artifacts".to_string()
    })
}

fn needs_visible_triangle_pipeline(plan: &ArtifactPlan<'_>) -> bool {
    plan.triangle_reps
        .iter()
        .any(|rep| rep.visibility_counter_index.is_some())
}

fn create_triangle_build_pipeline(
    viewer: &mut dyn ViewerLike,
) -> Result<TriangleBuildPipeline, String> {
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

    Ok(TriangleBuildPipeline {
        pipeline,
        layout,
        params_buffer,
    })
}

fn create_visible_triangle_build_pipeline(
    viewer: &mut dyn ViewerLike,
    plan: &ArtifactPlan<'_>,
    counter_buffer: GpuHandle,
    batch_commands: &mut Vec<GpuBatchCommand>,
) -> Result<VisibleTriangleBuildPipeline, String> {
    let counter_bytes = surface_visibility_counter_bytes(plan)?;
    batch_commands.push(GpuBatchCommand::WriteBuffer {
        buffer: counter_buffer,
        offset: 0,
        data: vec![0; counter_bytes as usize],
    });

    let shader = viewer
        .gpu_create_cached_shader_module(GpuShaderModuleDescriptor {
            label: Some("ray.artifact.build_visible_triangles.shader".to_string()),
            wgsl: shaders::artifact_visible_triangles(),
        })?
        .handle;
    let layout = viewer
        .gpu_create_cached_bind_group_layout(GpuBindGroupLayoutDescriptor {
            label: Some("ray.artifact.build_visible_triangles.layout".to_string()),
            entries: vec![
                storage_layout(0, GpuBufferBindingType::StorageReadOnly),
                storage_layout(1, GpuBufferBindingType::StorageReadOnly),
                storage_layout(2, GpuBufferBindingType::StorageReadWrite),
                storage_layout(3, GpuBufferBindingType::StorageReadOnly),
                storage_layout(4, GpuBufferBindingType::StorageReadWrite),
                storage_layout(5, GpuBufferBindingType::Uniform),
            ],
        })?
        .handle;
    let pipeline_layout = viewer
        .gpu_create_cached_pipeline_layout(GpuPipelineLayoutDescriptor {
            label: Some("ray.artifact.build_visible_triangles.pipeline_layout".to_string()),
            bind_group_layouts: vec![layout],
        })?
        .handle;
    let pipeline = viewer
        .gpu_create_cached_compute_pipeline(GpuComputePipelineDescriptor {
            label: Some("ray.artifact.build_visible_triangles.pipeline".to_string()),
            layout: pipeline_layout,
            module: shader,
            entry_point: "build_visible_triangles".to_string(),
        })?
        .handle;
    let params_buffer = create_buffer(
        viewer,
        "ray.artifact.build_visible_triangles.params",
        std::mem::size_of::<ArtifactVisibleTriangleParams>() as u64,
        uniform_usage(),
        None,
    )?;

    Ok(VisibleTriangleBuildPipeline {
        pipeline,
        layout,
        params_buffer,
        counter_buffer,
        counter_bytes,
    })
}

fn surface_visibility_counter_bytes(plan: &ArtifactPlan<'_>) -> Result<u64, String> {
    let Some(max_counter_index) = plan
        .triangle_reps
        .iter()
        .filter_map(|rep| rep.visibility_counter_index)
        .max()
    else {
        return Ok(0);
    };
    let counter_count = max_counter_index
        .checked_add(1)
        .ok_or_else(|| "surface visibility counter buffer size overflow".to_string())?;
    u64::from(counter_count)
        .checked_mul(std::mem::size_of::<u32>() as u64)
        .ok_or_else(|| "surface visibility counter buffer size overflow".to_string())
}
