//! Line representation settings (object-overridable).

use crate::define_settings_group;
use crate::Color;

define_settings_group! {
    /// Line representation parameters.
    group LineSettings / LineOverrides {
        color: Color = Color::UNSET,
            name = "line_color",
            side_effects = [ColorRebuild];
    }
}
