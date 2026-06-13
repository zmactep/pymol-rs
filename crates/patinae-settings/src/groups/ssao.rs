//! Screen-space ambient occlusion (global-only).
//!
//! Reconstructs view-space position and normal from the depth buffer, samples a
//! 32-point hemisphere
//! per pixel, writes a single-channel occlusion factor that the compose
//! pass multiplies into the opaque-pass colour. Translucent surfaces do
//! not participate (no depth write under WBOIT).

use crate::define_settings_group;

define_settings_group! {
    /// Screen-space ambient occlusion (depth-only Crytek-style SSAO).
    group_global SsaoSettings {
        /// Master toggle. Off by default — large assemblies pay a measurable
        /// per-frame cost (compute + bilateral blur + compose).
        enabled: bool = false,
            name = "ssao_enabled";
        /// Sampling sphere radius in world units. Smaller → tighter
        /// crevice darkening; larger → broader ambient shadow.
        radius: f32 = 0.5,
            name = "ssao_radius",
            min = 0.05, max = 5.0;
        /// AO factor multiplier on the final compose. `1.0` = full
        /// darkening at fully-occluded pixels; `0.5` = half-strength.
        intensity: f32 = 1.0,
            name = "ssao_intensity",
            min = 0.0, max = 4.0;
        /// Depth bias to suppress self-occlusion on near-flat surfaces.
        /// In view-space units (Å for molecular scenes).
        bias: f32 = 0.025,
            name = "ssao_bias",
            min = 0.0, max = 1.0;
    }
}
