//! Line representation settings (object-overridable).

use crate::define_settings_group;

define_settings_group! {
    /// Line representation parameters.
    group LineSettings / LineOverrides {
        color: i32 = -1,
            name = "line_color",
            side_effects = [ColorRebuild];
    }
}
