//! Sphere representation settings (object-overridable).

use crate::define_settings_group;

define_settings_group! {
    /// Sphere representation parameters.
    group SphereSettings / SphereOverrides {
        scale: f32 = 1.0,
            name = "sphere_scale",
            min = 0.01, max = 10.0,
            side_effects = [RepresentationRebuild];
        quality: i32 = 1,
            name = "sphere_quality",
            min = 0.0, max = 4.0,
            side_effects = [RepresentationRebuild];
        color: i32 = -1,
            name = "sphere_color",
            side_effects = [ColorRebuild];
    }
}
