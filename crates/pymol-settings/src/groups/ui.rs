//! UI / interface settings (global-only).

use crate::define_settings_group;
use crate::enums::MouseSelectionMode;

define_settings_group! {
    /// User interface and interaction settings.
    group_global UiSettings {
        mouse_selection_mode: MouseSelectionMode = MouseSelectionMode::Residues,
            name = "mouse_selection_mode",
            hints = MouseSelectionMode,
            side_effects = [OrthoDirty];
        mouse_wheel_scale: f32 = 0.5,
            name = "mouse_wheel_scale",
            min = 0.01, max = 10.0;
        antialias: i32 = 1,
            name = "antialias",
            side_effects = [SceneInvalidate];
        selection_width: f32 = 3.0,
            name = "selection_width",
            min = 0.5, max = 20.0,
            side_effects = [SceneInvalidate];
        transparent_panels: bool = false,
            name = "transparent_panels",
            side_effects = [OrthoDirty];
        opaque_background: bool = false,
            name = "opaque_background",
            side_effects = [SceneInvalidate];
        orthoscopic: bool = false,
            name = "orthoscopic",
            side_effects = [SceneInvalidate];
        bg_rgb: i32 = 0,
            name = "bg_rgb",
            side_effects = [ViewportUpdate];
        bg_rgb_top: i32 = 0x00004D,
            name = "bg_rgb_top",
            side_effects = [ViewportUpdate];
        bg_rgb_bottom: i32 = 0x333380,
            name = "bg_rgb_bottom",
            side_effects = [ViewportUpdate];
    }
}
