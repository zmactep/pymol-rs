//! Mesh representation settings (object-overridable).

use crate::define_settings_group;

define_settings_group! {
    /// Mesh representation parameters.
    group MeshSettings / MeshOverrides {
        color: i32 = -1,
            name = "mesh_color",
            side_effects = [ColorRebuild];
        lighting: bool = false,
            name = "mesh_lighting",
            side_effects = [SceneInvalidate];
    }
}
