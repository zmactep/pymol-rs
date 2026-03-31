//! Movie / animation settings (global-only).

use crate::define_settings_group;

define_settings_group! {
    /// Movie playback and animation settings.
    group_global MovieSettings {
        movie_loop: bool = true,
            name = "movie_loop";
        movie_fps: f32 = 30.0,
            name = "movie_fps",
            min = 1.0, max = 120.0;
        movie_auto_interpolate: bool = true,
            name = "movie_auto_interpolate";
    }
}
