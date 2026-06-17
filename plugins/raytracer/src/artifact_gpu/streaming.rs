use patinae_plugin::prelude::ViewerLike;
use patinae_scene::GpuDeviceLimits;

use super::batch::TriangleBuildWindow;
use super::dispatch::{self, RaytraceDispatchTarget};
use super::layout::{
    ArtifactTriangle, BvhBuffers, BvhBuildInput, FinalizeStreamingMetadataParams, PrimitiveBuffers,
    PrimitiveCounts, RaytraceDispatchInput, MAX_ENCODED_PRIMITIVE_INDEX,
};
use super::reps::ArtifactPlan;
use super::resources::{
    checked_storage_buffer_size, checked_storage_bytes, create_buffer, storage_bytes_for_device,
    storage_usage,
};
use super::{batch, bvh, metadata, visibility, RaytraceArtifactOutput, RaytraceArtifactTarget};
use crate::bvh::BvhNode;
use crate::gpu::RaytraceParams;
use crate::primitive::{GpuCapsule, GpuCylinder, GpuSphere};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum ArtifactRayPlan {
    SingleGlobal,
    StreamingFallback,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) struct StreamingChunkShape {
    pub(super) triangle_capacity: u32,
    pub(super) primitive_capacity: u32,
    pub(super) bvh_node_capacity: u32,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) struct StreamingTriangleBudget {
    pub(super) direct_triangle_count: u32,
    pub(super) skipped_direct_triangle_count: u32,
    pub(super) visible_triangle_capacity: u32,
    pub(super) include_direct_triangles: bool,
}

pub(super) fn choose_artifact_ray_plan(
    plan: &ArtifactPlan<'_>,
    limits: &GpuDeviceLimits,
) -> ArtifactRayPlan {
    if single_triangle_storage_fits(plan, limits) {
        ArtifactRayPlan::SingleGlobal
    } else {
        ArtifactRayPlan::StreamingFallback
    }
}

pub(super) fn single_triangle_storage_fits(
    plan: &ArtifactPlan<'_>,
    limits: &GpuDeviceLimits,
) -> bool {
    storage_bytes_for_device::<ArtifactTriangle>(
        plan.triangle_count,
        limits,
        "ray artifact triangles",
    )
    .is_ok()
}

pub(super) fn choose_streaming_chunk_shape(
    plan: &ArtifactPlan<'_>,
    limits: &GpuDeviceLimits,
    max_dispatch_dimension: u32,
) -> Result<StreamingChunkShape, String> {
    let has_triangle_source = plan
        .triangle_reps
        .iter()
        .any(|rep| rep.source_triangle_count() > 0);
    if !has_triangle_source {
        return Err("streaming ray fallback requires triangle artifacts".to_string());
    }

    let storage_limit = limits
        .max_buffer_size
        .min(limits.max_storage_buffer_binding_size);
    let max_by_triangle_bytes =
        (storage_limit / std::mem::size_of::<ArtifactTriangle>() as u64).max(1);
    let mut candidate = u32::try_from(max_by_triangle_bytes.min(u64::from(u32::MAX)))
        .unwrap_or(u32::MAX)
        .clamp(1, MAX_ENCODED_PRIMITIVE_INDEX);
    let analytic_count = plan
        .sphere_count
        .checked_add(plan.cylinder_count)
        .and_then(|count| count.checked_add(plan.capsule_count))
        .ok_or_else(|| "streaming ray analytic primitive count overflow".to_string())?;

    while candidate > 0 {
        let primitive_capacity = analytic_count
            .checked_add(candidate)
            .ok_or_else(|| "streaming ray primitive capacity overflow".to_string())?;
        if storage_bytes_for_device::<ArtifactTriangle>(
            candidate,
            limits,
            "streaming ray triangle chunk",
        )
        .is_ok()
        {
            let shape = bvh::bvh_shape_for(primitive_capacity)?;
            let bvh_nodes_ok = checked_storage_buffer_size(
                checked_storage_bytes(
                    u64::from(shape.node_count),
                    std::mem::size_of::<BvhNode>() as u64,
                    "streaming ray BVH nodes",
                )?,
                limits,
                "streaming ray BVH nodes",
            )
            .is_ok();
            let bvh_indices_ok = checked_storage_buffer_size(
                checked_storage_bytes(
                    u64::from(primitive_capacity),
                    std::mem::size_of::<u32>() as u64,
                    "streaming ray BVH indices",
                )?,
                limits,
                "streaming ray BVH indices",
            )
            .is_ok();
            let dispatch_ok = dispatch::dispatch_grid_for(shape.leaf_slots, max_dispatch_dimension)
                .is_ok()
                && dispatch::dispatch_grid_for(candidate, max_dispatch_dimension).is_ok();
            if bvh_nodes_ok && bvh_indices_ok && dispatch_ok {
                return Ok(StreamingChunkShape {
                    triangle_capacity: candidate,
                    primitive_capacity,
                    bvh_node_capacity: shape.node_count,
                });
            }
        }
        candidate /= 2;
    }

    Err(format!(
        "streaming ray fallback cannot allocate a one-triangle chunk below GPU storage buffer limit {storage_limit}"
    ))
}

pub(super) fn raytrace_streamed_artifacts(
    viewer: &mut dyn ViewerLike,
    plan: &ArtifactPlan<'_>,
    device_limits: &GpuDeviceLimits,
    params: &RaytraceParams,
    target: RaytraceArtifactTarget,
) -> Result<RaytraceArtifactOutput, String> {
    let max_dispatch_dimension = device_limits.max_compute_workgroups_per_dimension;
    let chunk_shape = choose_streaming_chunk_shape(plan, device_limits, max_dispatch_dimension)?;
    let mut batch_commands = Vec::new();

    let sphere_buffer = create_buffer(
        viewer,
        "ray.artifact.streaming.spheres",
        storage_bytes_for_device::<GpuSphere>(
            plan.sphere_count,
            device_limits,
            "streaming ray spheres",
        )?,
        storage_usage(),
        None,
    )?;
    let cylinder_buffer = create_buffer(
        viewer,
        "ray.artifact.streaming.cylinders",
        storage_bytes_for_device::<GpuCylinder>(
            plan.cylinder_count,
            device_limits,
            "streaming ray cylinders",
        )?,
        storage_usage(),
        None,
    )?;
    let capsule_buffer = create_buffer(
        viewer,
        "ray.artifact.streaming.capsules",
        storage_bytes_for_device::<GpuCapsule>(
            plan.capsule_count,
            device_limits,
            "streaming ray capsules",
        )?,
        storage_usage(),
        None,
    )?;
    let triangle_buffer = create_buffer(
        viewer,
        "ray.artifact.streaming.triangle_chunk",
        storage_bytes_for_device::<ArtifactTriangle>(
            chunk_shape.triangle_capacity,
            device_limits,
            "streaming ray triangle chunk",
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
    super::set_profile_scope(&mut batch_commands, super::PROFILE_SCOPE_PRIMITIVE_BUILD);
    batch::build_spheres(
        viewer,
        plan,
        plan.color_lut,
        sphere_buffer,
        max_dispatch_dimension,
        &mut batch_commands,
    )?;
    batch::build_cylinders(
        viewer,
        plan,
        plan.color_lut,
        cylinder_buffer,
        max_dispatch_dimension,
        &mut batch_commands,
    )?;
    batch::build_capsules(
        viewer,
        plan,
        plan.color_lut,
        capsule_buffer,
        max_dispatch_dimension,
        &mut batch_commands,
    )?;

    let bvh_node_buffer = create_buffer(
        viewer,
        "ray.artifact.streaming.bvh_nodes",
        checked_storage_buffer_size(
            checked_storage_bytes(
                u64::from(chunk_shape.bvh_node_capacity),
                std::mem::size_of::<BvhNode>() as u64,
                "streaming ray BVH nodes",
            )?,
            device_limits,
            "streaming ray BVH nodes",
        )?,
        storage_usage(),
        None,
    )?;
    let bvh_index_buffer = create_buffer(
        viewer,
        "ray.artifact.streaming.bvh_indices",
        checked_storage_buffer_size(
            checked_storage_bytes(
                u64::from(chunk_shape.primitive_capacity),
                std::mem::size_of::<u32>() as u64,
                "streaming ray BVH indices",
            )?,
            device_limits,
            "streaming ray BVH indices",
        )?,
        storage_usage(),
        None,
    )?;
    let bvh_buffers = BvhBuffers {
        nodes: bvh_node_buffer,
        indices: bvh_index_buffer,
    };

    let visibility = visibility::create_surface_visibility_counters(viewer, plan, device_limits)?;
    let triangle_budget = streaming_triangle_budget(plan, chunk_shape.triangle_capacity)?;
    let source_stats = streaming_source_stats(plan)?;
    if triangle_budget.skipped_direct_triangle_count > 0 {
        log::warn!(
            "native GPU artifact streaming fallback skipped {} direct triangles to reserve reusable triangle capacity for visible surfaces",
            triangle_budget.skipped_direct_triangle_count
        );
    }
    let active_counts = PrimitiveCounts {
        spheres: plan.sphere_count,
        cylinders: plan.cylinder_count,
        capsules: plan.capsule_count,
        triangles: triangle_budget.direct_triangle_count,
    };
    let primitive_metadata = metadata::create_primitive_metadata(
        viewer,
        "ray.artifact.streaming.primitive_metadata",
        active_counts,
        chunk_shape.triangle_capacity,
        device_limits,
    )?;
    let leaf_dispatch_args = if visibility.is_some() {
        Some(metadata::create_dispatch_indirect_args(
            viewer,
            "ray.artifact.streaming.bvh_leaf_dispatch_args",
        )?)
    } else {
        None
    };

    let mut direct_offset = 0_u32;
    let mut surface_window_count = 0_u32;
    let mut reset_visible_counter = true;
    for rep in &plan.triangle_reps {
        if rep.visibility_counter_index.is_none() {
            if !triangle_budget.include_direct_triangles {
                continue;
            }
            let source_total = rep.source_triangle_count();
            batch::build_triangle_window(
                viewer,
                batch::TriangleWindowBuildInput {
                    plan,
                    rep,
                    color_lut: plan.color_lut,
                    triangle_buffer,
                    visibility_counter_buffer: None,
                    ray_params: params,
                    window: TriangleBuildWindow {
                        source_triangle_start: 0,
                        source_triangle_count: source_total,
                        triangle_offset: direct_offset,
                        output_triangle_count: rep.triangle_count,
                        counter_index: None,
                        reset_counter: false,
                    },
                    max_dispatch_dimension,
                },
                &mut batch_commands,
            )?;
            direct_offset = direct_offset
                .checked_add(rep.triangle_count)
                .ok_or_else(|| "streaming ray direct triangle offset overflow".to_string())?;
            continue;
        }

        let mut source_start = 0_u32;
        let source_total = rep.source_triangle_count();
        while source_start < source_total {
            let remaining = source_total - source_start;
            let source_count = remaining.min(chunk_shape.triangle_capacity);
            batch::build_triangle_window(
                viewer,
                batch::TriangleWindowBuildInput {
                    plan,
                    rep,
                    color_lut: plan.color_lut,
                    triangle_buffer,
                    visibility_counter_buffer: visibility
                        .as_ref()
                        .map(|visibility| (visibility.counter_buffer, visibility.counter_bytes)),
                    ray_params: params,
                    window: TriangleBuildWindow {
                        source_triangle_start: source_start,
                        source_triangle_count: source_count,
                        triangle_offset: triangle_budget.direct_triangle_count,
                        output_triangle_count: triangle_budget.visible_triangle_capacity,
                        counter_index: Some(0),
                        reset_counter: reset_visible_counter,
                    },
                    max_dispatch_dimension,
                },
                &mut batch_commands,
            )?;
            reset_visible_counter = false;
            source_start = source_start
                .checked_add(source_count)
                .ok_or_else(|| "streaming ray source window overflow".to_string())?;
            surface_window_count = surface_window_count
                .checked_add(1)
                .ok_or_else(|| "streaming ray source window count overflow".to_string())?;
        }
    }

    if let Some(visibility) = &visibility {
        let leaf_dispatch_args = leaf_dispatch_args
            .ok_or_else(|| "streaming ray missing BVH leaf dispatch args".to_string())?;
        metadata::finalize_streaming_metadata(
            viewer,
            metadata::FinalizeStreamingMetadataInput {
                metadata: primitive_metadata,
                visible_counter: visibility.counter_buffer,
                visible_counter_bytes: visibility.counter_bytes,
                leaf_dispatch_args,
                params: FinalizeStreamingMetadataParams {
                    direct_triangle_count: triangle_budget.direct_triangle_count,
                    visible_triangle_capacity: triangle_budget.visible_triangle_capacity,
                    counter_index: 0,
                    bvh_leaf_size: super::layout::BVH_LEAF_SIZE,
                    max_workgroups_per_dimension: max_dispatch_dimension,
                    _pad0: 0,
                    _pad1: 0,
                    _pad2: 0,
                },
                max_dispatch_dimension,
            },
            &mut batch_commands,
        )?;
    }

    let capacity_counts = PrimitiveCounts {
        spheres: plan.sphere_count,
        cylinders: plan.cylinder_count,
        capsules: plan.capsule_count,
        triangles: chunk_shape.triangle_capacity,
    };
    let primitive_count = capacity_counts.primitive_count()?;
    let bvh_shape = bvh::bvh_shape_for(primitive_count)?;
    bvh::build_bvh(
        viewer,
        BvhBuildInput {
            primitives: primitive_buffers,
            bvh: bvh_buffers,
            metadata: primitive_metadata,
            leaf_dispatch_args,
            shape: &bvh_shape,
            max_dispatch_dimension,
        },
        &mut batch_commands,
    )?;
    let mut profile_lines = Vec::new();
    super::push_profile_line(
        &mut profile_lines,
        format!(
            "patinae.rt_profile plugin.artifact_path=streaming_fallback_active_counts bvh_leaf_dispatch={} plugin.streaming_compact_prepare_commands={} surface_reps={} surface_windows={} direct_source_triangles={} surface_source_triangles={} direct_triangles={} visible_triangle_capacity={} chunk_triangle_capacity={} primitive_capacity={} bvh_nodes={}",
            if leaf_dispatch_args.is_some() {
                "indirect_active"
            } else {
                "direct_capacity"
            },
            batch_commands.len(),
            source_stats.surface_rep_count,
            surface_window_count,
            source_stats.direct_source_triangles,
            source_stats.surface_source_triangles,
            triangle_budget.direct_triangle_count,
            triangle_budget.visible_triangle_capacity,
            chunk_shape.triangle_capacity,
            primitive_count,
            bvh_shape.node_count
        ),
    );

    match dispatch::dispatch_raytrace(
        viewer,
        RaytraceDispatchInput {
            params,
            primitives: primitive_buffers,
            bvh: bvh_buffers,
            counts: capacity_counts,
            metadata: primitive_metadata,
            bvh_node_count: bvh_shape.node_count,
            max_dispatch_dimension,
        },
        &mut batch_commands,
        match target {
            RaytraceArtifactTarget::CpuReadback => RaytraceDispatchTarget::CpuReadback,
            RaytraceArtifactTarget::ViewportGpu => RaytraceDispatchTarget::ViewportGpu,
        },
    )? {
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

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct StreamingSourceStats {
    direct_source_triangles: u32,
    surface_source_triangles: u32,
    surface_rep_count: u32,
}

fn streaming_source_stats(plan: &ArtifactPlan<'_>) -> Result<StreamingSourceStats, String> {
    let mut stats = StreamingSourceStats {
        direct_source_triangles: 0,
        surface_source_triangles: 0,
        surface_rep_count: 0,
    };
    for rep in &plan.triangle_reps {
        if rep.visibility_counter_index.is_some() {
            stats.surface_source_triangles = stats
                .surface_source_triangles
                .checked_add(rep.source_triangle_count())
                .ok_or_else(|| {
                    "streaming ray surface source triangle count overflow".to_string()
                })?;
            stats.surface_rep_count = stats
                .surface_rep_count
                .checked_add(1)
                .ok_or_else(|| "streaming ray surface rep count overflow".to_string())?;
        } else {
            stats.direct_source_triangles = stats
                .direct_source_triangles
                .checked_add(rep.source_triangle_count())
                .ok_or_else(|| "streaming ray direct source triangle count overflow".to_string())?;
        }
    }
    Ok(stats)
}

fn direct_triangle_count(plan: &ArtifactPlan<'_>) -> Result<u32, String> {
    plan.triangle_reps
        .iter()
        .filter(|rep| rep.visibility_counter_index.is_none())
        .try_fold(0_u32, |count, rep| {
            count
                .checked_add(rep.triangle_count)
                .ok_or_else(|| "streaming ray direct triangle count overflow".to_string())
        })
}

fn has_visible_surface_triangles(plan: &ArtifactPlan<'_>) -> bool {
    plan.triangle_reps
        .iter()
        .any(|rep| rep.visibility_counter_index.is_some())
}

pub(super) fn streaming_triangle_budget(
    plan: &ArtifactPlan<'_>,
    triangle_capacity: u32,
) -> Result<StreamingTriangleBudget, String> {
    let requested_direct_triangle_count = direct_triangle_count(plan)?;
    let has_visible_surfaces = has_visible_surface_triangles(plan);
    let include_direct_triangles =
        !has_visible_surfaces || requested_direct_triangle_count < triangle_capacity;
    if requested_direct_triangle_count > triangle_capacity && !has_visible_surfaces {
        return Err(format!(
            "streaming ray fallback has {} direct triangles, exceeding reusable triangle capacity {}",
            requested_direct_triangle_count, triangle_capacity
        ));
    }
    let direct_triangle_count = if include_direct_triangles {
        requested_direct_triangle_count
    } else {
        0
    };
    let visible_triangle_capacity = triangle_capacity
        .checked_sub(direct_triangle_count)
        .ok_or_else(|| "streaming ray visible triangle capacity underflow".to_string())?;
    if visible_triangle_capacity == 0 && has_visible_surfaces {
        return Err(
            "streaming ray fallback has no remaining triangle capacity for visible surfaces"
                .to_string(),
        );
    }
    Ok(StreamingTriangleBudget {
        direct_triangle_count,
        skipped_direct_triangle_count: requested_direct_triangle_count - direct_triangle_count,
        visible_triangle_capacity,
        include_direct_triangles,
    })
}
