use patinae_plugin::prelude::ViewerLike;

use super::indirect::{active_triangle_vertex_count, read_indirect_draw_args};
use super::layout::STD_VERTEX_STRIDE;
use super::reps::ArtifactPlan;
use super::resources::checked_storage_bytes;

pub(super) fn resolve_indirect_triangle_counts(
    viewer: &mut dyn ViewerLike,
    plan: &mut ArtifactPlan<'_>,
) -> Result<(), String> {
    let mut triangle_reps = Vec::with_capacity(plan.triangle_reps.len());
    let mut triangle_count = 0_u32;

    for mut rep in std::mem::take(&mut plan.triangle_reps) {
        if let Some(indirect) = rep.rep.indirect {
            let draw_args = read_indirect_draw_args(viewer, indirect, rep.rep.rep_kind)?;
            rep.vertex_count = active_triangle_vertex_count(draw_args[0], rep.vertex_count);
            rep.triangle_count = rep.vertex_count / 3;
            rep.geometry_binding_size = checked_storage_bytes(
                u64::from(rep.vertex_count),
                STD_VERTEX_STRIDE,
                "triangle geometry",
            )?;
        }
        if rep.triangle_count == 0 {
            continue;
        }
        rep.triangle_offset = triangle_count;
        triangle_count = triangle_count
            .checked_add(rep.triangle_count)
            .ok_or_else(|| "artifact triangle count overflow".to_string())?;
        triangle_reps.push(rep);
    }

    plan.triangle_reps = triangle_reps;
    plan.triangle_count = triangle_count;
    Ok(())
}

pub(super) fn apply_visible_surface_triangle_counts(
    plan: &mut ArtifactPlan<'_>,
    counts: &[u32],
) -> Result<(), String> {
    let mut triangle_reps = Vec::with_capacity(plan.triangle_reps.len());
    let mut triangle_count = 0_u32;
    let mut counter_index = 0_u32;

    for mut rep in std::mem::take(&mut plan.triangle_reps) {
        if rep.uses_surface_visibility_culling() {
            let visible_count = counts
                .get(counter_index as usize)
                .copied()
                .ok_or_else(|| "surface visibility count readback is missing a rep".to_string())?;
            if visible_count > rep.source_triangle_count() {
                return Err(format!(
                    "surface visibility count {visible_count} exceeds source triangle count {}",
                    rep.source_triangle_count()
                ));
            }
            rep.triangle_count = visible_count;
            rep.visibility_counter_index = Some(counter_index);
            counter_index = counter_index
                .checked_add(1)
                .ok_or_else(|| "surface visibility counter index overflow".to_string())?;
        }

        if rep.triangle_count == 0 {
            continue;
        }
        rep.triangle_offset = triangle_count;
        triangle_count = triangle_count
            .checked_add(rep.triangle_count)
            .ok_or_else(|| "artifact triangle count overflow".to_string())?;
        triangle_reps.push(rep);
    }

    if counter_index as usize != counts.len() {
        return Err(format!(
            "surface visibility readback returned {} counters for {counter_index} surface reps",
            counts.len()
        ));
    }

    plan.triangle_reps = triangle_reps;
    plan.triangle_count = triangle_count;
    Ok(())
}
