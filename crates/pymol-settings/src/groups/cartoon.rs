//! Cartoon representation settings (object-overridable).

use crate::define_settings_group;

define_settings_group! {
    /// Cartoon representation parameters.
    group CartoonSettings / CartoonOverrides {
        sampling: i32 = -1,
            name = "cartoon_sampling",
            side_effects = [RepresentationRebuild];
        power: f32 = 2.0,
            name = "cartoon_power",
            min = 0.0, max = 10.0,
            side_effects = [RepresentationRebuild];
        power_b: f32 = 0.52,
            name = "cartoon_power_b",
            min = 0.0, max = 10.0,
            side_effects = [RepresentationRebuild];
        throw: f32 = 1.35,
            name = "cartoon_throw",
            min = 0.0, max = 10.0,
            side_effects = [RepresentationRebuild];
        smooth_first: i32 = 1,
            name = "cartoon_smooth_first",
            side_effects = [RepresentationRebuild];
        smooth_last: i32 = 1,
            name = "cartoon_smooth_last",
            side_effects = [RepresentationRebuild];
        smooth_cycles: i32 = 2,
            name = "cartoon_smooth_cycles",
            side_effects = [RepresentationRebuild];
        smooth_loops: bool = false,
            name = "cartoon_smooth_loops",
            side_effects = [RepresentationRebuild];
        flat_cycles: i32 = 4,
            name = "cartoon_flat_cycles",
            side_effects = [RepresentationRebuild];
        gap_cutoff: i32 = 10,
            name = "cartoon_gap_cutoff",
            side_effects = [RepresentationRebuild];
        refine: i32 = 5,
            name = "cartoon_refine",
            side_effects = [RepresentationRebuild];
        refine_normals: i32 = -1,
            name = "cartoon_refine_normals",
            side_effects = [RepresentationRebuild];
        round_helices: bool = true,
            name = "cartoon_round_helices",
            side_effects = [RepresentationRebuild];
        color: i32 = -1,
            name = "cartoon_color",
            side_effects = [ColorRebuild];
        nucleic_acid_color: i32 = -1,
            name = "cartoon_nucleic_acid_color",
            side_effects = [ColorRebuild];
    }
}
