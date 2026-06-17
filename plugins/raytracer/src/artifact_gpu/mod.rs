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
use patinae_scene::{GpuBatchCommand, RenderArtifactSnapshotDescriptor};

use crate::bvh::BvhNode;
use crate::gpu::RaytraceParams;
use crate::primitive::{GpuCapsule, GpuCylinder, GpuSphere};

mod batch;
mod bvh;
mod dispatch;
mod layout;
mod metadata;
mod plan;
mod reps;
mod resources;
mod streaming;
mod triangles;
mod visibility;

use layout::{
    ArtifactTriangle, BvhBuffers, BvhBuildInput, PrimitiveBuffers, RaytraceDispatchInput,
};
use resources::{
    checked_storage_buffer_size, checked_storage_bytes, create_buffer, storage_bytes_for_device,
    storage_usage,
};
use triangles::prepare_triangle_gpu_metadata;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum RaytraceArtifactTarget {
    CpuReadback,
    ViewportGpu,
}

pub(crate) enum RaytraceArtifactOutput {
    CpuImage {
        data: Vec<u8>,
        profile_lines: Vec<String>,
    },
    ViewportGpu {
        profile_lines: Vec<String>,
    },
}

#[derive(Default)]
struct ArtifactProfile {
    plan_ms: u128,
    triangle_metadata_ms: u128,
    temp_resource_setup_ms: u128,
    primitive_command_setup_ms: u128,
    bvh_resource_setup_ms: u128,
    bvh_command_setup_ms: u128,
    dispatch_call_ms: u128,
    total_ms: u128,
}

pub(super) const PROFILE_SCOPE_PRIMITIVE_BUILD: &str = "primitive_build";
pub(super) const PROFILE_SCOPE_SURFACE_COMPACT: &str = "surface_compact";
pub(super) const PROFILE_SCOPE_METADATA_FINALIZE: &str = "metadata_finalize";
pub(super) const PROFILE_SCOPE_BVH_BUILD: &str = "bvh_build";
pub(super) const PROFILE_SCOPE_RAYTRACE: &str = "raytrace";
pub(super) const PROFILE_SCOPE_DOWNSAMPLE: &str = "downsample";

pub(super) fn set_profile_scope(batch_commands: &mut Vec<GpuBatchCommand>, name: &'static str) {
    batch_commands.push(GpuBatchCommand::SetProfileScope {
        name: name.to_string(),
    });
}

/// Raytrace renderer GPU artifacts through the command-runtime GPU ABI.
///
/// This function consumes a command-scoped renderer artifact snapshot. It
/// validates supported layouts, records all conversion/BVH/raytrace work into
/// one ordered GPU batch, then either returns RGBA bytes for export or promotes
/// the final GPU buffer into the native viewport.
pub(crate) fn raytrace_artifacts(
    viewer: &mut dyn ViewerLike,
    snapshot: &RenderArtifactSnapshotDescriptor,
    params: &RaytraceParams,
    target: RaytraceArtifactTarget,
) -> Result<RaytraceArtifactOutput, String> {
    if params.settings.ray_trace_mode != 0 {
        return Err(
            "native GPU artifact ray path currently supports ray_trace_mode=0 (Normal)".to_string(),
        );
    }

    let total_start = std::time::Instant::now();
    let mut profile = ArtifactProfile::default();
    let stage_start = std::time::Instant::now();
    let mut plan = plan::plan_artifact_primitives(snapshot)?;
    profile.plan_ms = stage_start.elapsed().as_millis();

    let stage_start = std::time::Instant::now();
    let surface_visibility_rep_count = prepare_triangle_gpu_metadata(&mut plan)?;
    profile.triangle_metadata_ms = stage_start.elapsed().as_millis();

    let primitive_count = plan.primitive_count()?;
    if primitive_count == 0 {
        return Err("GPU artifact snapshot contains no raytraceable GPU primitives".to_string());
    }
    if streaming::choose_artifact_ray_plan(&plan, &snapshot.device_limits)
        == streaming::ArtifactRayPlan::StreamingFallback
    {
        if rt_profile_enabled() {
            log::info!(
                "patinae.rt_profile plugin.artifact_streaming_fallback triangles={} primitive_count={} surface_visibility_reps={}",
                plan.triangle_count,
                primitive_count,
                surface_visibility_rep_count
            );
        }
        return streaming::raytrace_streamed_artifacts(
            viewer,
            &plan,
            &snapshot.device_limits,
            params,
            target,
        );
    }

    let stage_start = std::time::Instant::now();
    let mut batch_commands = Vec::new();
    let mut profile_lines = Vec::new();
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
    let surface_visibility =
        visibility::create_surface_visibility_counters(viewer, &plan, &snapshot.device_limits)?;
    let surface_visibility_rep_count = surface_visibility
        .as_ref()
        .map(|visibility| visibility.rep_count)
        .unwrap_or(surface_visibility_rep_count);
    let primitive_buffers = PrimitiveBuffers {
        spheres: sphere_buffer,
        cylinders: cylinder_buffer,
        capsules: capsule_buffer,
        triangles: triangle_buffer,
    };
    profile.temp_resource_setup_ms = stage_start.elapsed().as_millis();

    let stage_start = std::time::Instant::now();
    let primitive_counts = plan.primitive_counts();
    let max_dispatch_dimension = snapshot.device_limits.max_compute_workgroups_per_dimension;
    let primitive_metadata = metadata::create_primitive_metadata(
        viewer,
        "ray.artifact.primitive_metadata",
        primitive_counts,
        plan.triangle_count,
        device_limits,
    )?;
    set_profile_scope(&mut batch_commands, PROFILE_SCOPE_PRIMITIVE_BUILD);
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
        batch::TriangleBuildInput {
            color_lut: plan.color_lut,
            triangle_buffer: primitive_buffers.triangles,
            visibility_counter_buffer: surface_visibility
                .as_ref()
                .map(|visibility| (visibility.counter_buffer, visibility.counter_bytes)),
            ray_params: params,
            max_dispatch_dimension,
        },
        &mut batch_commands,
    )?;
    profile.primitive_command_setup_ms = stage_start.elapsed().as_millis();

    let stage_start = std::time::Instant::now();
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
    profile.bvh_resource_setup_ms = stage_start.elapsed().as_millis();

    let stage_start = std::time::Instant::now();
    bvh::build_bvh(
        viewer,
        BvhBuildInput {
            primitives: primitive_buffers,
            bvh: bvh_buffers,
            metadata: primitive_metadata,
            leaf_dispatch_args: None,
            shape: &bvh_shape,
            max_dispatch_dimension,
        },
        &mut batch_commands,
    )?;
    profile.bvh_command_setup_ms = stage_start.elapsed().as_millis();
    push_profile_line(
        &mut profile_lines,
        format!(
            "patinae.rt_profile plugin.artifact_prepare plan_ms={} triangle_metadata_ms={} temp_resource_setup_ms={} primitive_command_setup_ms={} bvh_resource_setup_ms={} bvh_command_setup_ms={} spheres={} cylinders={} capsules={} triangles={} primitive_count={} bvh_nodes={} surface_visibility_reps={} pre_dispatch_commands={}",
            profile.plan_ms,
            profile.triangle_metadata_ms,
            profile.temp_resource_setup_ms,
            profile.primitive_command_setup_ms,
            profile.bvh_resource_setup_ms,
            profile.bvh_command_setup_ms,
            plan.sphere_count,
            plan.cylinder_count,
            plan.capsule_count,
            plan.triangle_count,
            primitive_count,
            bvh_shape.node_count,
            surface_visibility_rep_count,
            batch_commands.len()
        ),
    );

    let stage_start = std::time::Instant::now();
    let output = dispatch::dispatch_raytrace(
        viewer,
        RaytraceDispatchInput {
            params,
            primitives: primitive_buffers,
            bvh: bvh_buffers,
            counts: primitive_counts,
            metadata: primitive_metadata,
            bvh_node_count: bvh_shape.node_count,
            max_dispatch_dimension,
        },
        &mut batch_commands,
        match target {
            RaytraceArtifactTarget::CpuReadback => dispatch::RaytraceDispatchTarget::CpuReadback,
            RaytraceArtifactTarget::ViewportGpu => dispatch::RaytraceDispatchTarget::ViewportGpu,
        },
    )?;
    profile.dispatch_call_ms = stage_start.elapsed().as_millis();
    profile.total_ms = total_start.elapsed().as_millis();
    let output_bytes = match &output {
        dispatch::RaytraceDispatchOutput::CpuImage { data, .. } => data.len(),
        dispatch::RaytraceDispatchOutput::ViewportGpu { .. } => 0,
    };
    push_profile_line(
        &mut profile_lines,
        format!(
            "patinae.rt_profile plugin.artifact_total total_ms={} dispatch_call_ms={} output_bytes={} output_width={} output_height={} antialias={}",
            profile.total_ms,
            profile.dispatch_call_ms,
            output_bytes,
            params.width,
            params.height,
            params.antialias.max(1)
        ),
    );

    match output {
        dispatch::RaytraceDispatchOutput::CpuImage {
            data,
            profile_lines: mut dispatch_profile_lines,
        } => {
            profile_lines.append(&mut dispatch_profile_lines);
            Ok(RaytraceArtifactOutput::CpuImage {
                data,
                profile_lines,
            })
        }
        dispatch::RaytraceDispatchOutput::ViewportGpu {
            profile_lines: mut dispatch_profile_lines,
        } => {
            profile_lines.append(&mut dispatch_profile_lines);
            Ok(RaytraceArtifactOutput::ViewportGpu { profile_lines })
        }
    }
}

pub(super) fn push_profile_line(lines: &mut Vec<String>, line: String) {
    if rt_profile_enabled() {
        log::info!("{line}");
        lines.push(line);
    }
}

fn rt_profile_enabled() -> bool {
    std::env::var_os("PATINAE_RT_PROFILE").is_some()
}

fn rt_profile_counts_enabled() -> bool {
    rt_profile_enabled() && std::env::var_os("PATINAE_RT_PROFILE_COUNTS").is_some()
}

#[cfg(test)]
mod tests;
