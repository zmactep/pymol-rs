//! Ray tracing settings (global-only).

use crate::define_settings_group;

define_settings_group! {
    /// Ray tracing parameters.
    group_global RayTraceSettings {
        max_passes: i32 = 25,
            name = "ray_max_passes",
            min = 1.0, max = 100.0;
        shadow: bool = true,
            name = "ray_shadow";
        trace_mode: i32 = 0,
            name = "ray_trace_mode";
        trace_fog: f32 = -1.0,
            name = "ray_trace_fog";
        trace_color: i32 = -6,
            name = "ray_trace_color";
        opaque_background: i32 = -1,
            name = "ray_opaque_background";
        transparency_shadows: bool = true,
            name = "ray_transparency_shadows";
        trace_depth_factor: f32 = 0.1,
            name = "ray_trace_depth_factor",
            min = 0.0, max = 1.0;
        trace_slope_factor: f32 = 0.6,
            name = "ray_trace_slope_factor",
            min = 0.0, max = 2.0;
        trace_disco_factor: f32 = 0.05,
            name = "ray_trace_disco_factor",
            min = 0.0, max = 1.0;
        trace_gain: f32 = 0.12,
            name = "ray_trace_gain",
            min = 0.0, max = 1.0;
    }
}
