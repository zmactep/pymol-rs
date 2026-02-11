//! Shared utility functions for cartoon and ribbon representations

use lin_alg::f32::Vec3;
use pymol_mol::SecondaryStructure;

use super::frame::FrameWithMetadata;

/// Normalize a vector safely, returning a default for zero-length vectors
#[inline]
pub fn normalize_safe(v: Vec3) -> Vec3 {
    let len_sq = v.magnitude_squared();
    if len_sq > 1e-10 {
        v / len_sq.sqrt()
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    }
}

/// PyMOL's smooth easing function for tapering and interpolation
///
/// Creates an ease-in/ease-out curve for smooth transitions.
/// Matches PyMOL's `smooth(x, power)` function from Vector.cpp.
///
/// # Arguments
/// * `x` - Input value in range [0, 1]
/// * `power` - Power factor (typically 2.0 for quadratic easing)
///
/// # Returns
/// Smoothed value in range [0, 1]
#[inline]
pub fn smooth(x: f32, power: f32) -> f32 {
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }
    if x <= 0.5 {
        0.5 * (2.0 * x).powf(power)
    } else {
        1.0 - 0.5 * (2.0 * (1.0 - x)).powf(power)
    }
}

/// Check if a secondary structure type is a helix
#[inline]
pub fn is_helix(ss: SecondaryStructure) -> bool {
    matches!(
        ss,
        SecondaryStructure::Helix | SecondaryStructure::Helix310 | SecondaryStructure::HelixPi
    )
}

/// Find contiguous regions of frames matching a secondary structure predicate
///
/// Returns a vector of `(start_idx, end_idx)` tuples for each contiguous run
/// where the predicate returns true.
pub fn find_ss_regions(
    frames: &[FrameWithMetadata],
    predicate: impl Fn(SecondaryStructure) -> bool,
) -> Vec<(usize, usize)> {
    let mut regions = Vec::new();
    let mut in_region = false;
    let mut region_start = 0;

    for (i, frame) in frames.iter().enumerate() {
        let matches = predicate(frame.ss_type);

        if matches && !in_region {
            region_start = i;
            in_region = true;
        } else if !matches && in_region {
            regions.push((region_start, i - 1));
            in_region = false;
        }
    }

    // Handle region at end
    if in_region {
        regions.push((region_start, frames.len() - 1));
    }

    regions
}
