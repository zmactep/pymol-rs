/// Smooth density falloff used by production SAS surface generation.
pub(super) const SAS_ALPHA: f32 = 4.0;

/// Iso level at which an isolated atom's SAS surface lies exactly at
/// `r + probe`. `SAS_ISO = exp(-SAS_ALPHA)`.
pub(super) const SAS_ISO: f32 = 0.018_315_64;
