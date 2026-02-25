//! Crystallographic unit cell math
//!
//! Provides conversion between Cartesian (real-space) and fractional coordinates
//! using unit cell parameters. The math follows PyMOL's `Crystal.cpp` implementation.

use lin_alg::f32::{Mat4, Vec3};

use crate::linalg::{mat3x3_to_mat4, invert_3x3, transform_3x3};

/// Crystallographic unit cell with precomputed transformation matrices
pub struct CrystalCell {
    /// Cell edge lengths (a, b, c) in Angstroms
    pub lengths: [f32; 3],
    /// Cell angles (alpha, beta, gamma) in degrees
    pub angles: [f32; 3],
    /// 3x3 fractional-to-real matrix (row-major)
    frac_to_real: [f32; 9],
    /// 3x3 real-to-fractional matrix (row-major)
    real_to_frac: [f32; 9],
}

impl CrystalCell {
    /// Create a new crystal cell from unit cell parameters
    ///
    /// Computes the fractional↔real transformation matrices on construction.
    pub fn new(lengths: [f32; 3], angles: [f32; 3]) -> Self {
        let frac_to_real = compute_frac_to_real(&lengths, &angles);
        let real_to_frac = invert_3x3(&frac_to_real);
        CrystalCell {
            lengths,
            angles,
            frac_to_real,
            real_to_frac,
        }
    }

    /// Get the 3x3 fractional-to-real matrix (row-major)
    pub fn frac_to_real(&self) -> &[f32; 9] {
        &self.frac_to_real
    }

    /// Get the 3x3 real-to-fractional matrix (row-major)
    pub fn real_to_frac(&self) -> &[f32; 9] {
        &self.real_to_frac
    }

    /// Get fractional-to-real as a 4x4 homogeneous matrix
    pub fn frac_to_real_4x4(&self) -> Mat4 {
        mat3x3_to_mat4(&self.frac_to_real)
    }

    /// Get real-to-fractional as a 4x4 homogeneous matrix
    pub fn real_to_frac_4x4(&self) -> Mat4 {
        mat3x3_to_mat4(&self.real_to_frac)
    }

    /// Transform a real-space vector to fractional coordinates
    pub fn to_fractional(&self, v: Vec3) -> Vec3 {
        transform_3x3(&self.real_to_frac, v)
    }

    /// Transform fractional coordinates to real-space
    pub fn to_real(&self, v: Vec3) -> Vec3 {
        transform_3x3(&self.frac_to_real, v)
    }
}

/// Compute the frac-to-real 3x3 matrix from cell parameters
///
/// Follows PyMOL's `Crystal.cpp:fracToReal()` (lines 189-221).
/// The matrix is row-major: `[row0_col0, row0_col1, row0_col2, row1_col0, ...]`
fn compute_frac_to_real(lengths: &[f32; 3], angles: &[f32; 3]) -> [f32; 9] {
    let [a, b, c] = *lengths;
    let [alpha_deg, beta_deg, gamma_deg] = *angles;

    // Handle degenerate case
    if a == 0.0 || b == 0.0 || c == 0.0 || alpha_deg == 0.0 || beta_deg == 0.0 || gamma_deg == 0.0
    {
        return [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
    }

    let alpha = alpha_deg.to_radians();
    let beta = beta_deg.to_radians();
    let gamma = gamma_deg.to_radians();

    let ca = alpha.cos();
    let cb = beta.cos();
    let cg = gamma.cos();
    let sb = beta.sin();
    let sg = gamma.sin();

    // cos(alpha*) = (cos(beta)*cos(gamma) - cos(alpha)) / (sin(beta)*sin(gamma))
    let cabgs = (cb * cg - ca) / (sb * sg);
    let sabgs = (1.0 - cabgs * cabgs).max(0.0).sqrt();

    // Row-major 3x3 matching PyMOL's FracToReal layout:
    // f2r[0] = a,       f2r[1] = cg*b,      f2r[2] = cb*c
    // f2r[3] = 0,       f2r[4] = sg*b,      f2r[5] = -sb*cabgs*c
    // f2r[6] = 0,       f2r[7] = 0,          f2r[8] = sb*sabgs*c
    [
        a,
        cg * b,
        cb * c,
        0.0,
        sg * b,
        -sb * cabgs * c,
        0.0,
        0.0,
        sb * sabgs * c,
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_orthogonal_cell() {
        // Orthogonal cell (90/90/90) should produce near-identity scaling
        let cell = CrystalCell::new([10.0, 20.0, 30.0], [90.0, 90.0, 90.0]);

        // frac (0.5, 0.5, 0.5) should map to real (5, 10, 15)
        let frac = Vec3::new(0.5, 0.5, 0.5);
        let real = cell.to_real(frac);
        assert!((real.x - 5.0).abs() < 1e-4);
        assert!((real.y - 10.0).abs() < 1e-4);
        assert!((real.z - 15.0).abs() < 1e-4);

        // Round-trip
        let back = cell.to_fractional(real);
        assert!((back.x - 0.5).abs() < 1e-4);
        assert!((back.y - 0.5).abs() < 1e-4);
        assert!((back.z - 0.5).abs() < 1e-4);
    }

    #[test]
    fn test_monoclinic_cell() {
        // Monoclinic cell: alpha=90, beta=100, gamma=90
        let cell = CrystalCell::new([10.0, 20.0, 30.0], [90.0, 100.0, 90.0]);

        // Round-trip test
        let v = Vec3::new(5.0, 10.0, 15.0);
        let frac = cell.to_fractional(v);
        let back = cell.to_real(frac);
        assert!((back.x - v.x).abs() < 1e-3);
        assert!((back.y - v.y).abs() < 1e-3);
        assert!((back.z - v.z).abs() < 1e-3);
    }
}
