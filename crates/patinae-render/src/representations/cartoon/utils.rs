//! Cartoon utility helpers.
//!
//! Used by the run-aware cartoon pipeline. Functions here are CPU-only;
//! equivalents needed inside WGSL kernels are duplicated in those shaders.

use lin_alg::f32::Vec3;
use patinae_mol::SecondaryStructure;

/// Piecewise sigmoid blend on `[0, 1]`, used as the displacement envelope
/// during cartoon sampling.
pub fn smooth(x: f32, power: f32) -> f32 {
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }
    if x <= 0.5 {
        return 0.5 * (2.0 * x).powf(power);
    }
    1.0 - 0.5 * (2.0 * (1.0 - x)).powf(power)
}

/// Normalize `v`. If magnitude is below `1e-6`, return the zero vector.
pub fn normalize_safe(v: Vec3) -> Vec3 {
    let l = v.magnitude();
    if l > 1e-6 {
        v / l
    } else {
        Vec3::new(0.0, 0.0, 0.0)
    }
}

/// True if `ss` is any alpha-helix variant: alpha, 3-10, or pi.
pub fn is_helix(ss: SecondaryStructure) -> bool {
    matches!(
        ss,
        SecondaryStructure::Helix | SecondaryStructure::Helix310 | SecondaryStructure::HelixPi
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn smooth_endpoints() {
        assert_eq!(smooth(0.0, 2.0), 0.0);
        assert_eq!(smooth(1.0, 2.0), 1.0);
        assert!((smooth(0.5, 2.0) - 0.5).abs() < 1e-6);
    }

    #[test]
    fn smooth_power_2_is_quadratic_at_quarter() {
        // smooth(0.25, 2) = 0.5 * (0.5)^2 = 0.125
        assert!((smooth(0.25, 2.0) - 0.125).abs() < 1e-6);
        // smooth(0.75, 2) = 1 - 0.5 * (0.5)^2 = 0.875
        assert!((smooth(0.75, 2.0) - 0.875).abs() < 1e-6);
    }

    #[test]
    fn normalize_safe_zero_returns_zero() {
        let z = Vec3::new(1e-9, 0.0, 0.0);
        assert_eq!(normalize_safe(z), Vec3::new(0.0, 0.0, 0.0));
    }

    #[test]
    fn normalize_safe_returns_unit() {
        let v = Vec3::new(3.0, 4.0, 0.0);
        let n = normalize_safe(v);
        assert!((n.magnitude() - 1.0).abs() < 1e-6);
    }

    #[test]
    fn is_helix_classification() {
        assert!(is_helix(SecondaryStructure::Helix));
        assert!(is_helix(SecondaryStructure::Helix310));
        assert!(is_helix(SecondaryStructure::HelixPi));
        assert!(!is_helix(SecondaryStructure::Sheet));
        assert!(!is_helix(SecondaryStructure::Loop));
        assert!(!is_helix(SecondaryStructure::Turn));
    }
}
