use super::reps::ArtifactPlan;

pub(super) fn prepare_triangle_gpu_metadata(plan: &mut ArtifactPlan<'_>) -> Result<usize, String> {
    let mut counter_index = 0_u32;

    for rep in &mut plan.triangle_reps {
        if rep.uses_surface_visibility_culling() {
            rep.visibility_counter_index = Some(counter_index);
            counter_index = counter_index
                .checked_add(1)
                .ok_or_else(|| "surface visibility counter index overflow".to_string())?;
        }
    }

    Ok(counter_index as usize)
}
