//! General object settings.

use crate::define_settings_group;

define_settings_group! {
    /// Generic object settings that do not belong to a representation group.
    group ObjectSettings / ObjectSettingOverrides {
        state: i32 = 1,
            name = "state",
            side_effects = [ViewportUpdate];
    }
}
