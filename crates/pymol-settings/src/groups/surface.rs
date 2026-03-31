//! Surface representation settings (object-overridable).

use crate::define_settings_group;

define_settings_group! {
    /// Surface representation parameters.
    group SurfaceSettings / SurfaceOverrides {
        surface_type: i32 = 0,
            name = "surface_type",
            side_effects = [RepresentationRebuild];
        solvent: bool = false,
            name = "surface_solvent",
            side_effects = [RepresentationRebuild];
        quality: i32 = 0,
            name = "surface_quality",
            side_effects = [RepresentationRebuild];
        color: i32 = -1,
            name = "surface_color",
            side_effects = [ColorRebuild];
        individual_chains: bool = false,
            name = "surface_individual_chains",
            side_effects = [RepresentationRebuild];
        transparency: f32 = 0.0,
            name = "transparency",
            min = 0.0, max = 1.0,
            side_effects = [SceneInvalidate];
        transparency_mode: i32 = 2,
            name = "transparency_mode",
            side_effects = [SceneInvalidate];
    }
}
