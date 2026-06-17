use patinae_plugin::prelude::ViewerLike;
use patinae_scene::{GpuDeviceLimits, GpuHandle};

use super::layout::ArtifactVisibleTriangleParams;
use super::reps::{ArtifactPlan, TriangleArtifactRep};
use super::resources::{
    checked_storage_buffer_size, checked_storage_bytes, create_buffer, storage_usage,
};
use crate::gpu::RaytraceParams;

pub(super) struct SurfaceVisibility {
    pub(super) counter_buffer: GpuHandle,
    pub(super) counter_bytes: u64,
    pub(super) rep_count: usize,
}

#[derive(Clone, Copy)]
pub(super) struct VisibleTriangleParamWindow {
    pub(super) source_triangle_start: u32,
    pub(super) source_triangle_count: u32,
    pub(super) triangle_offset: u32,
    pub(super) output_triangle_capacity: u32,
    pub(super) counter_index: u32,
    pub(super) dispatch_width: u32,
}

pub(super) fn create_surface_visibility_counters(
    viewer: &mut dyn ViewerLike,
    plan: &ArtifactPlan<'_>,
    device_limits: &GpuDeviceLimits,
) -> Result<Option<SurfaceVisibility>, String> {
    let surface_rep_count = plan
        .triangle_reps
        .iter()
        .filter_map(|rep| rep.visibility_counter_index)
        .max()
        .map(|index| index as usize + 1)
        .unwrap_or(0);
    if surface_rep_count == 0 {
        return Ok(None);
    }

    let counter_bytes = checked_storage_buffer_size(
        checked_storage_bytes(
            surface_rep_count as u64,
            std::mem::size_of::<u32>() as u64,
            "surface visibility counters",
        )?,
        device_limits,
        "surface visibility counters",
    )?;
    let counter_buffer = create_buffer(
        viewer,
        "ray.artifact.surface_visible_counts",
        counter_bytes,
        storage_usage(),
        None,
    )?;

    Ok(Some(SurfaceVisibility {
        counter_buffer,
        counter_bytes,
        rep_count: surface_rep_count,
    }))
}

pub(super) fn visible_triangle_params(
    params: &RaytraceParams,
    rep: &TriangleArtifactRep<'_>,
    window: VisibleTriangleParamWindow,
) -> ArtifactVisibleTriangleParams {
    ArtifactVisibleTriangleParams {
        view_matrix: params.view_matrix,
        proj_matrix: params.proj_matrix,
        source_vertex_count: rep.vertex_count,
        source_triangle_start: window.source_triangle_start,
        source_triangle_count: window.source_triangle_count,
        triangle_offset: window.triangle_offset,
        output_triangle_capacity: window.output_triangle_capacity,
        atom_offset: rep.rep.atom_offset,
        rep_slot: rep.rep_slot,
        transparency: rep.rep.transparency,
        dispatch_width: window.dispatch_width,
        counter_index: window.counter_index,
        _pad0: 0,
        _pad1: 0,
    }
}
