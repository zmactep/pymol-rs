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

use bytemuck::{Pod, Zeroable};
use patinae_plugin::prelude::ViewerLike;
use patinae_scene::{GpuHandle, RenderArtifactRepDescriptor, RenderArtifactSnapshotDescriptor};

use crate::bvh::BvhNode;
use crate::gpu::RaytraceParams;
use crate::primitive::{GpuCapsule, GpuCylinder, GpuSphere, GpuTriangle};

mod batch;
mod bvh;
mod dispatch;
mod plan;
mod resources;

use resources::{create_buffer, direct_draw_args_buffer, storage_bytes_for, storage_usage};

const WORKGROUP_SIZE: u32 = 128;
const BVH_LEAF_SIZE: u32 = 4;
const EMPTY_STORAGE_BYTES: u64 = 16;
const COLOR_LUT_STRIDE: u64 = 64;
const SPHERE_INSTANCE_STRIDE: u64 = 32;
const STICK_INSTANCE_STRIDE: u64 = 48;
const LINE_INSTANCE_STRIDE: u64 = 48;
const STD_VERTEX_STRIDE: u64 = 24;
const MAX_OUTPUT_READBACK_BYTES: u64 = 64 * 1024 * 1024;
const MAX_ENCODED_PRIMITIVE_INDEX: u32 = 0x3fff_ffff;
const RAY_LINE_RADIUS: f32 = 0.035;

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct ArtifactTriangleParams {
    vertex_capacity: u32,
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
struct ArtifactSphereParams {
    instance_capacity: u32,
    sphere_offset: u32,
    atom_offset: u32,
    rep_slot: u32,
    transparency: f32,
    dispatch_width: u32,
    _pad0: u32,
    _pad1: u32,
}

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct ArtifactCylinderParams {
    instance_capacity: u32,
    cylinder_offset: u32,
    atom_offset: u32,
    rep_slot: u32,
    transparency: f32,
    radius: f32,
    dispatch_width: u32,
    _pad0: u32,
}

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct ArtifactCapsuleParams {
    instance_capacity: u32,
    capsule_offset: u32,
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
    primitive_count: u32,
    sphere_count: u32,
    cylinder_count: u32,
    capsule_count: u32,
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

struct SphereArtifactRep<'a> {
    rep: &'a RenderArtifactRepDescriptor,
    sphere_offset: u32,
    instance_capacity: u32,
    geometry_binding_size: u64,
    rep_slot: u32,
}

impl SphereArtifactRep<'_> {
    fn indirect_or_direct_args(
        &self,
        viewer: &mut dyn ViewerLike,
        label: &str,
    ) -> Result<GpuHandle, String> {
        if let Some(indirect) = self.rep.indirect {
            return Ok(indirect);
        }
        direct_draw_args_buffer(viewer, label, [0, self.instance_capacity, 0, 0])
    }
}

struct CylinderArtifactRep<'a> {
    rep: &'a RenderArtifactRepDescriptor,
    cylinder_offset: u32,
    instance_capacity: u32,
    geometry_binding_size: u64,
    rep_slot: u32,
    radius: f32,
}

impl CylinderArtifactRep<'_> {
    fn indirect_or_direct_args(
        &self,
        viewer: &mut dyn ViewerLike,
        label: &str,
    ) -> Result<GpuHandle, String> {
        if let Some(indirect) = self.rep.indirect {
            return Ok(indirect);
        }
        direct_draw_args_buffer(viewer, label, [0, self.instance_capacity, 0, 0])
    }
}

struct CapsuleArtifactRep<'a> {
    rep: &'a RenderArtifactRepDescriptor,
    capsule_offset: u32,
    instance_capacity: u32,
    geometry_binding_size: u64,
    rep_slot: u32,
}

impl CapsuleArtifactRep<'_> {
    fn indirect_or_direct_args(
        &self,
        viewer: &mut dyn ViewerLike,
        label: &str,
    ) -> Result<GpuHandle, String> {
        if let Some(indirect) = self.rep.indirect {
            return Ok(indirect);
        }
        direct_draw_args_buffer(viewer, label, [0, self.instance_capacity, 0, 0])
    }
}

struct TriangleArtifactRep<'a> {
    rep: &'a RenderArtifactRepDescriptor,
    triangle_offset: u32,
    triangle_count: u32,
    vertex_count: u32,
    geometry_binding_size: u64,
    rep_slot: u32,
}

impl TriangleArtifactRep<'_> {
    fn indirect_or_direct_args(
        &self,
        viewer: &mut dyn ViewerLike,
        label: &str,
    ) -> Result<GpuHandle, String> {
        if let Some(indirect) = self.rep.indirect {
            return Ok(indirect);
        }
        direct_draw_args_buffer(viewer, label, [self.vertex_count, 1, 0, 0])
    }
}

struct ArtifactPlan<'a> {
    color_lut: GpuHandle,
    sphere_reps: Vec<SphereArtifactRep<'a>>,
    cylinder_reps: Vec<CylinderArtifactRep<'a>>,
    capsule_reps: Vec<CapsuleArtifactRep<'a>>,
    triangle_reps: Vec<TriangleArtifactRep<'a>>,
    sphere_count: u32,
    cylinder_count: u32,
    capsule_count: u32,
    triangle_count: u32,
}

#[derive(Clone, Copy)]
struct PrimitiveBuffers {
    spheres: GpuHandle,
    cylinders: GpuHandle,
    capsules: GpuHandle,
    triangles: GpuHandle,
}

#[derive(Clone, Copy)]
struct PrimitiveCounts {
    spheres: u32,
    cylinders: u32,
    capsules: u32,
    triangles: u32,
}

impl ArtifactPlan<'_> {
    fn primitive_counts(&self) -> PrimitiveCounts {
        PrimitiveCounts {
            spheres: self.sphere_count,
            cylinders: self.cylinder_count,
            capsules: self.capsule_count,
            triangles: self.triangle_count,
        }
    }
}

#[derive(Clone, Copy)]
struct BvhBuffers {
    nodes: GpuHandle,
    indices: GpuHandle,
}

struct BvhBuildInput<'a> {
    primitives: PrimitiveBuffers,
    bvh: BvhBuffers,
    counts: PrimitiveCounts,
    primitive_count: u32,
    shape: &'a BvhShape,
    max_dispatch_dimension: u32,
}

struct RaytraceDispatchInput<'a> {
    params: &'a RaytraceParams,
    primitives: PrimitiveBuffers,
    bvh: BvhBuffers,
    counts: PrimitiveCounts,
    bvh_node_count: u32,
    max_dispatch_dimension: u32,
}

struct DownsampleInput {
    source_buffer: GpuHandle,
    src_width: u32,
    dst_width: u32,
    dst_height: u32,
    factor: u32,
    max_dispatch_dimension: u32,
}

impl ArtifactPlan<'_> {
    fn primitive_count(&self) -> Result<u32, String> {
        self.sphere_count
            .checked_add(self.cylinder_count)
            .and_then(|count| count.checked_add(self.capsule_count))
            .and_then(|count| count.checked_add(self.triangle_count))
            .ok_or_else(|| "artifact primitive count overflow".to_string())
    }
}

struct BvhShape {
    leaf_slots: u32,
    leaf_start: u32,
    node_count: u32,
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

    let plan = plan::plan_artifact_primitives(snapshot)?;
    let primitive_count = plan.primitive_count()?;
    if primitive_count == 0 {
        return Err("GPU artifact snapshot contains no raytraceable GPU primitives".to_string());
    }

    let setup_start = std::time::Instant::now();
    let mut batch_commands = Vec::new();
    let sphere_buffer = create_buffer(
        viewer,
        "ray.artifact.spheres",
        storage_bytes_for::<GpuSphere>(plan.sphere_count, "ray artifact spheres")?,
        storage_usage(),
        None,
    )?;
    let cylinder_buffer = create_buffer(
        viewer,
        "ray.artifact.cylinders",
        storage_bytes_for::<GpuCylinder>(plan.cylinder_count, "ray artifact cylinders")?,
        storage_usage(),
        None,
    )?;
    let capsule_buffer = create_buffer(
        viewer,
        "ray.artifact.capsules",
        storage_bytes_for::<GpuCapsule>(plan.capsule_count, "ray artifact capsules")?,
        storage_usage(),
        None,
    )?;
    let triangle_buffer = create_buffer(
        viewer,
        "ray.artifact.triangles",
        storage_bytes_for::<GpuTriangle>(plan.triangle_count, "ray artifact triangles")?,
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
        max_dispatch_dimension,
        &mut batch_commands,
    )?;

    let bvh_shape = bvh::bvh_shape_for(primitive_count)?;
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
        u64::from(primitive_count) * std::mem::size_of::<u32>() as u64,
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
