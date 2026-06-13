//! Fast approximate anti-aliasing (global-only).
//!
//! Single full-screen postprocess (Lottes 2009 / FXAA 3.11) that smooths
//! impostor-representation
//! silhouettes. Default-enabled because it is cheap (~0.3 ms at 1080p
//! on M-series Apple GPUs) and noticeably improves edge quality.

use crate::define_settings_group;

define_settings_group! {
    /// Fast approximate anti-aliasing (FXAA 3.11).
    group_global FxaaSettings {
        enabled: bool = true,
            name = "fxaa_enabled",
            side_effects = [SceneInvalidate];
    }
}
