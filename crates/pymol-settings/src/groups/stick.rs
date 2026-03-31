//! Stick representation settings (object-overridable).

use crate::define_settings_group;

define_settings_group! {
    /// Stick representation parameters.
    group StickSettings / StickOverrides {
        radius: f32 = 0.25,
            name = "stick_radius",
            min = 0.01, max = 5.0,
            side_effects = [RepresentationRebuild];
        color: i32 = -1,
            name = "stick_color",
            side_effects = [ColorRebuild];
        valence_scale: f32 = 1.0,
            name = "stick_valence_scale",
            min = 0.0, max = 5.0,
            side_effects = [RepresentationRebuild];
        valence: bool = true,
            name = "valence",
            side_effects = [RepresentationRebuild];
        valence_size: f32 = 0.06,
            name = "valence_size",
            min = 0.0, max = 1.0,
            side_effects = [RepresentationRebuild];
    }
}
