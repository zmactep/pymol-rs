//! Sphere representation settings (object-overridable).

use crate::define_settings_group;
use crate::Color;

define_settings_group! {
    /// Sphere representation parameters.
    group SphereSettings / SphereOverrides {
        scale: f32 = 1.0,
            name = "sphere_scale",
            min = 0.01, max = 10.0,
            side_effects = [RepresentationRebuild];
        color: Color = Color::UNSET,
            name = "sphere_color",
            side_effects = [ColorRebuild];
        transparency: f32 = 0.0,
            name = "sphere_transparency",
            min = 0.0, max = 1.0,
            side_effects = [ColorRebuild];
    }
}
