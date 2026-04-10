//! Crystallographic unit cell math
//!
//! Provides conversion between Cartesian (real-space) and fractional coordinates
//! using unit cell parameters.
//!
//! Reference: Giacovazzo et al., "Fundamentals of Crystallography", 3rd ed.,
//! IUCr/Oxford University Press, Ch. 2 — defines the fractional↔Cartesian
//! transformation from cell lengths (a, b, c) and angles (α, β, γ).

use lin_alg::f32::{Mat3, Mat4, Vec3};

/// Crystallographic unit cell with precomputed transformation matrices
pub struct CrystalCell {
    /// Cell edge lengths (a, b, c) in Angstroms
    pub lengths: Vec3,
    /// Cell angles (alpha, beta, gamma) in degrees
    pub angles: Vec3,
    /// 3x3 fractional-to-real matrix (column-major, lin_alg convention)
    frac_to_real: Mat3,
    /// 3x3 real-to-fractional matrix (column-major, lin_alg convention)
    real_to_frac: Mat3,
}

impl CrystalCell {
    /// Create a new crystal cell from unit cell parameters
    ///
    /// Computes the fractional↔real transformation matrices on construction.
    pub fn new(lengths: Vec3, angles: Vec3) -> Self {
        let frac_to_real = compute_frac_to_real(lengths, angles);
        let real_to_frac = frac_to_real.inverse().unwrap_or(Mat3::new_identity());
        CrystalCell {
            lengths,
            angles,
            frac_to_real,
            real_to_frac,
        }
    }

    /// Get the 3x3 fractional-to-real matrix
    pub fn frac_to_real(&self) -> &Mat3 {
        &self.frac_to_real
    }

    /// Get the 3x3 real-to-fractional matrix
    pub fn real_to_frac(&self) -> &Mat3 {
        &self.real_to_frac
    }

    /// Get fractional-to-real as a 4x4 homogeneous matrix (column-major)
    pub fn frac_to_real_4x4(&self) -> Mat4 {
        mat3_to_mat4(&self.frac_to_real)
    }

    /// Get real-to-fractional as a 4x4 homogeneous matrix (column-major)
    pub fn real_to_frac_4x4(&self) -> Mat4 {
        mat3_to_mat4(&self.real_to_frac)
    }

    /// Transform a real-space vector to fractional coordinates
    pub fn to_fractional(&self, v: Vec3) -> Vec3 {
        self.real_to_frac.clone() * v
    }

    /// Transform fractional coordinates to real-space
    pub fn to_real(&self, v: Vec3) -> Vec3 {
        self.frac_to_real.clone() * v
    }
}

/// Expand a Mat3 into a Mat4.
fn mat3_to_mat4(m: &Mat3) -> Mat4 {
    let d = &m.data;
    Mat4::new([
        d[0], d[1], d[2], 0.0, // col 0
        d[3], d[4], d[5], 0.0, // col 1
        d[6], d[7], d[8], 0.0, // col 2
        0.0,  0.0,  0.0,  1.0, // col 3
    ])
}

/// Compute the frac-to-real 3x3 matrix from cell parameters.
///
/// The matrix M satisfies r = M · f where r is Cartesian and f is fractional.
/// Derived from Giacovazzo et al., "Fundamentals of Crystallography", eq. 2.18.
fn compute_frac_to_real(lengths: Vec3, angles: Vec3) -> Mat3 {
    let (a, b, c) = (lengths.x, lengths.y, lengths.z);
    let (alpha_deg, beta_deg, gamma_deg) = (angles.x, angles.y, angles.z);

    // Handle degenerate case
    if a == 0.0 || b == 0.0 || c == 0.0 || alpha_deg == 0.0 || beta_deg == 0.0 || gamma_deg == 0.0
    {
        return Mat3::new_identity();
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

    // The 3x3 matrix in row-major form:
    // row 0: a,    cg*b,        cb*c
    // row 1: 0,    sg*b,        -sb*cabgs*c
    // row 2: 0,    0,           sb*sabgs*c
    //
    // In column-major (lin_alg): data[col*3 + row]
    Mat3::new([
        a,   0.0, 0.0,           // col 0: (a, 0, 0)
        cg * b, sg * b, 0.0,     // col 1: (cg*b, sg*b, 0)
        cb * c, -sb * cabgs * c, sb * sabgs * c, // col 2
    ])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_orthogonal_cell() {
        // Orthogonal cell (90/90/90) should produce near-identity scaling
        let cell = CrystalCell::new(
            Vec3::new(10.0, 20.0, 30.0),
            Vec3::new(90.0, 90.0, 90.0),
        );

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
        let cell = CrystalCell::new(
            Vec3::new(10.0, 20.0, 30.0),
            Vec3::new(90.0, 100.0, 90.0),
        );

        // Round-trip test
        let v = Vec3::new(5.0, 10.0, 15.0);
        let frac = cell.to_fractional(v);
        let back = cell.to_real(frac);
        assert!((back.x - v.x).abs() < 1e-3);
        assert!((back.y - v.y).abs() < 1e-3);
        assert!((back.z - v.z).abs() < 1e-3);
    }
}
