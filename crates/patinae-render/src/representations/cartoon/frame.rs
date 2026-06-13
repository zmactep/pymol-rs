//! Reference frame helper.
//!
//! Used by CPU sampling/refinement. Per-vertex frames at extrude time are
//! computed on the GPU, but the CPU sample stage occasionally needs
//! `ReferenceFrame` for picking-vertex placement and unit tests.

use lin_alg::f32::Vec3;
use patinae_mol::SecondaryStructure;

use super::utils::normalize_safe;

/// Right-handed orthonormal frame at one point along the cartoon spline:
///   T = tangent (along chain direction, unit length)
///   B = T × O (perpendicular to ribbon plane)
///   N = B × T (in ribbon plane, perpendicular to T)
///
/// For helix residues with `O = axis × tangent` (set by
/// `compute_round_helix_orientations`), this puts B mostly along the
/// helix axis and N radially outward — the orientation that makes the
/// oval cross-section read as a smooth ribbon wrapped around the cylinder.
///
/// For sheet residues with `O = CA→O` (after enforce-consistency +
/// flatten), O lies in-plane perpendicular to the strand, B = T × O is
/// out of strand plane, and N = B × T is in-plane perpendicular to chain.
/// The flat rectangle ribbon's wide axis runs along N, lying flat in the
/// strand plane (`generate_explicit_sheet` convention).
#[derive(Debug, Clone, Copy)]
pub struct ReferenceFrame {
    pub position: Vec3,
    pub tangent: Vec3,
    pub normal: Vec3,
    pub binormal: Vec3,
}

impl ReferenceFrame {
    /// Build a frame from `(position, tangent, orientation)`:
    ///   T = normalize(tangent)
    ///   B = normalize(T × orientation)
    ///   N = B × T
    pub fn new(position: Vec3, tangent: Vec3, orientation: Vec3) -> Self {
        let t = normalize_safe(tangent);
        let b = normalize_safe(t.cross(orientation));
        let n = b.cross(t);
        Self {
            position,
            tangent: t,
            normal: n,
            binormal: b,
        }
    }

    /// Lift a 2D `(x, y)` profile point in the (normal, binormal) plane
    /// into world-space.
    #[inline]
    pub fn transform_local(&self, local: (f32, f32)) -> Vec3 {
        self.position + self.normal * local.0 + self.binormal * local.1
    }

    /// Lift a 2D `(nx, ny)` profile normal into world-space.
    #[inline]
    pub fn local_normal(&self, local: (f32, f32)) -> Vec3 {
        normalize_safe(self.normal * local.0 + self.binormal * local.1)
    }
}

/// Reference frame plus colour, segment, and secondary-structure metadata.
#[derive(Debug, Clone, Copy)]
pub struct FrameWithMetadata {
    pub frame: ReferenceFrame,
    pub color: [f32; 4],
    pub ss_type: SecondaryStructure,
    /// Source segment index used by ribbon mesh generation.
    pub segment_idx: usize,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn frame_is_orthonormal() {
        let pos = Vec3::new(0.0, 0.0, 0.0);
        let t = Vec3::new(1.0, 0.0, 0.0);
        let o = Vec3::new(0.0, 1.0, 0.0);
        let f = ReferenceFrame::new(pos, t, o);
        // T = (1, 0, 0), O = (0, 1, 0), B = T × O = (0, 0, 1), N = B × T = (0, 1, 0).
        assert!((f.tangent - Vec3::new(1.0, 0.0, 0.0)).magnitude() < 1e-6);
        assert!((f.binormal - Vec3::new(0.0, 0.0, 1.0)).magnitude() < 1e-6);
        assert!((f.normal - Vec3::new(0.0, 1.0, 0.0)).magnitude() < 1e-6);
        // Orthogonal pairs.
        assert!(f.tangent.dot(f.binormal).abs() < 1e-6);
        assert!(f.tangent.dot(f.normal).abs() < 1e-6);
        assert!(f.binormal.dot(f.normal).abs() < 1e-6);
    }

    #[test]
    fn transform_local_lifts_into_n_b_plane() {
        let f = ReferenceFrame::new(
            Vec3::new(1.0, 2.0, 3.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        );
        // (x=2, y=3) in (n, b) → 1 + 2*(0,1,0) + 3*(0,0,1) = (1, 4, 6) when
        // pos = (1, 2, 3). N = (0,1,0), B = (0,0,1).
        let p = f.transform_local((2.0, 3.0));
        assert!((p - Vec3::new(1.0, 4.0, 6.0)).magnitude() < 1e-6);
    }
}
