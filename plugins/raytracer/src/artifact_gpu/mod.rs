//! GPU artifact ray path for command-runtime handles.
//!
//! The product `ray` command must not pull displayed geometry back through the
//! ABI. This module consumes renderer-owned artifact handles, builds the ray
//! sphere, cylinder, capsule, triangle, and BVH buffers on the host GPU, and dispatches
//! the ray shader via the safe GPU command callbacks.
//!
//! The generic artifact contract ends at `RenderArtifactSnapshotDescriptor`.
//! Ray-specific logic starts here: renderer sphere/stick/line instances and
//! supported `StdVertices` triangle-list representations are converted into
//! ray primitives, then the plugin builds a BVH and dispatches the ray shader
//! through the portable GPU batch API.
//!
//! Supported input is intentionally renderer-neutral: sphere instances become
//! `GpuSphere`, stick instances become `GpuCapsule`, line instances become
//! `GpuCylinder`, direct
//! cartoon/ribbon triangles use `element_count`, and surface triangle parts use
//! their indirect vertex count on the GPU. Mesh line-list artifacts still
//! return a clear error until a line-list ray path lands.
//!
//! Work buffers and bind groups remain temporary command resources. Shader
//! modules, bind-group layouts, pipeline layouts, and compute pipelines use the
//! persistent host cache exposed by the GPU runtime.

use patinae_plugin::prelude::ViewerLike;
use patinae_scene::RenderArtifactSnapshotDescriptor;

use crate::bvh::BvhNode;
use crate::gpu::RaytraceParams;
use crate::primitive::{GpuCapsule, GpuCylinder, GpuSphere};

mod batch;
mod bvh;
mod dispatch;
mod indirect;
mod layout;
mod plan;
mod reps;
mod resources;
mod triangles;
mod visibility;

use layout::{
    ArtifactTriangle, BvhBuffers, BvhBuildInput, PrimitiveBuffers, RaytraceDispatchInput,
};
use resources::{
    checked_storage_buffer_size, checked_storage_bytes, create_buffer, storage_bytes_for_device,
    storage_usage,
};
use triangles::{apply_visible_surface_triangle_counts, resolve_indirect_triangle_counts};

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

    let mut plan = plan::plan_artifact_primitives(snapshot)?;
    resolve_indirect_triangle_counts(viewer, &mut plan)?;
    let surface_visibility = visibility::count_visible_surface_triangles(
        viewer,
        &plan,
        params,
        &snapshot.device_limits,
        snapshot.device_limits.max_compute_workgroups_per_dimension,
    )?;
    if let Some(visibility) = &surface_visibility {
        apply_visible_surface_triangle_counts(&mut plan, &visibility.counts)?;
    }
    let primitive_count = plan.primitive_count()?;
    if primitive_count == 0 {
        return Err("GPU artifact snapshot contains no raytraceable GPU primitives".to_string());
    }

    let setup_start = std::time::Instant::now();
    let mut batch_commands = Vec::new();
    let device_limits = &snapshot.device_limits;
    let sphere_buffer = create_buffer(
        viewer,
        "ray.artifact.spheres",
        storage_bytes_for_device::<GpuSphere>(
            plan.sphere_count,
            device_limits,
            "ray artifact spheres",
        )?,
        storage_usage(),
        None,
    )?;
    let cylinder_buffer = create_buffer(
        viewer,
        "ray.artifact.cylinders",
        storage_bytes_for_device::<GpuCylinder>(
            plan.cylinder_count,
            device_limits,
            "ray artifact cylinders",
        )?,
        storage_usage(),
        None,
    )?;
    let capsule_buffer = create_buffer(
        viewer,
        "ray.artifact.capsules",
        storage_bytes_for_device::<GpuCapsule>(
            plan.capsule_count,
            device_limits,
            "ray artifact capsules",
        )?,
        storage_usage(),
        None,
    )?;
    let triangle_buffer = create_buffer(
        viewer,
        "ray.artifact.triangles",
        storage_bytes_for_device::<ArtifactTriangle>(
            plan.triangle_count,
            device_limits,
            "ray artifact triangles",
        )?,
        storage_usage(),
        None,
    )?;
    let primitive_buffers = PrimitiveBuffers {
        spheres: sphere_buffer,
        cylinders: cylinder_buffer,
        capsules: capsule_buffer,
        triangles: triangle_buffer,
    };
    let primitive_counts = plan.primitive_counts();
    let max_dispatch_dimension = snapshot.device_limits.max_compute_workgroups_per_dimension;
    batch::build_spheres(
        viewer,
        &plan,
        plan.color_lut,
        primitive_buffers.spheres,
        max_dispatch_dimension,
        &mut batch_commands,
    )?;
    batch::build_cylinders(
        viewer,
        &plan,
        plan.color_lut,
        primitive_buffers.cylinders,
        max_dispatch_dimension,
        &mut batch_commands,
    )?;
    batch::build_capsules(
        viewer,
        &plan,
        plan.color_lut,
        primitive_buffers.capsules,
        max_dispatch_dimension,
        &mut batch_commands,
    )?;
    batch::build_triangles(
        viewer,
        &plan,
        plan.color_lut,
        primitive_buffers.triangles,
        surface_visibility
            .as_ref()
            .map(|visibility| visibility.counter_buffer),
        params,
        max_dispatch_dimension,
        &mut batch_commands,
    )?;

    let bvh_shape = bvh::bvh_shape_for(primitive_count)?;
    let bvh_node_buffer = create_buffer(
        viewer,
        "ray.artifact.bvh_nodes",
        checked_storage_buffer_size(
            checked_storage_bytes(
                u64::from(bvh_shape.node_count),
                std::mem::size_of::<BvhNode>() as u64,
                "ray artifact BVH nodes",
            )?,
            device_limits,
            "ray artifact BVH nodes",
        )?,
        storage_usage(),
        None,
    )?;
    let bvh_index_buffer = create_buffer(
        viewer,
        "ray.artifact.bvh_indices",
        checked_storage_buffer_size(
            checked_storage_bytes(
                u64::from(primitive_count),
                std::mem::size_of::<u32>() as u64,
                "ray artifact BVH indices",
            )?,
            device_limits,
            "ray artifact BVH indices",
        )?,
        storage_usage(),
        None,
    )?;
    let bvh_buffers = BvhBuffers {
        nodes: bvh_node_buffer,
        indices: bvh_index_buffer,
    };
    bvh::build_bvh(
        viewer,
        BvhBuildInput {
            primitives: primitive_buffers,
            bvh: bvh_buffers,
            counts: primitive_counts,
            primitive_count,
            shape: &bvh_shape,
            max_dispatch_dimension,
        },
        &mut batch_commands,
    )?;
    if rt_profile_enabled() {
        log::info!(
            "patinae.rt_profile plugin.temp_resource_setup_ms={} spheres={} cylinders={} capsules={} triangles={} bvh_nodes={}",
            setup_start.elapsed().as_millis(),
            plan.sphere_count,
            plan.cylinder_count,
            plan.capsule_count,
            plan.triangle_count,
            bvh_shape.node_count
        );
    }

    let image = dispatch::dispatch_raytrace(
        viewer,
        RaytraceDispatchInput {
            params,
            primitives: primitive_buffers,
            bvh: bvh_buffers,
            counts: primitive_counts,
            bvh_node_count: bvh_shape.node_count,
            max_dispatch_dimension,
        },
        &mut batch_commands,
    )?;

    Ok(image)
}

fn rt_profile_enabled() -> bool {
    std::env::var_os("PATINAE_RT_PROFILE").is_some()
}

#[cfg(test)]
mod tests;
