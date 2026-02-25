//! 4×4 and 3×3 row-major matrix utilities
//!
//! Homogeneous transformation helpers matching PyMOL's row-major matrix
//! conventions (e.g. `copy33f44f`, `left_multiply44f44f`).

use lin_alg::f32::{Mat4, Vec3};

/// Expand a 3×3 row-major matrix into a 4×4 homogeneous Mat4 (row-major data)
///
/// Equivalent to PyMOL's `copy33f44f`: top-left 3×3 from src,
/// column 3 = 0, row 3 = [0,0,0,1].
pub fn mat3x3_to_mat4(m: &[f32; 9]) -> Mat4 {
    Mat4 {
        data: [
            m[0], m[1], m[2], 0.0, // row 0
            m[3], m[4], m[5], 0.0, // row 1
            m[6], m[7], m[8], 0.0, // row 2
            0.0, 0.0, 0.0, 1.0, // row 3
        ],
    }
}

/// Left-multiply: result = left * right (row-major 4×4)
///
/// Equivalent to PyMOL's `left_multiply44f44f` but returns a new matrix
/// instead of mutating `right` in place.
pub fn left_multiply_mat4(left: &Mat4, right: &Mat4) -> Mat4 {
    let l = &left.data;
    let r = &right.data;
    let mut out = [0.0f32; 16];
    for row in 0..4 {
        for col in 0..4 {
            out[row * 4 + col] = l[row * 4] * r[col]
                + l[row * 4 + 1] * r[4 + col]
                + l[row * 4 + 2] * r[8 + col]
                + l[row * 4 + 3] * r[12 + col];
        }
    }
    Mat4 { data: out }
}

/// Transform a Vec3 by a 4×4 row-major matrix (homogeneous, w=1)
pub fn transform_mat4(m: &Mat4, v: Vec3) -> Vec3 {
    Vec3::new(
        m.data[0] * v.x + m.data[1] * v.y + m.data[2] * v.z + m.data[3],
        m.data[4] * v.x + m.data[5] * v.y + m.data[6] * v.z + m.data[7],
        m.data[8] * v.x + m.data[9] * v.y + m.data[10] * v.z + m.data[11],
    )
}

/// Check if a 4×4 row-major matrix is approximately identity
pub fn is_identity_mat4(m: &Mat4) -> bool {
    let id = [
        1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    ];
    m.data
        .iter()
        .zip(id.iter())
        .all(|(a, b)| (a - b).abs() < 1e-4)
}

/// Transform a Vec3 by a 3×3 row-major matrix
pub fn transform_3x3(m: &[f32; 9], v: Vec3) -> Vec3 {
    Vec3::new(
        m[0] * v.x + m[1] * v.y + m[2] * v.z,
        m[3] * v.x + m[4] * v.y + m[5] * v.z,
        m[6] * v.x + m[7] * v.y + m[8] * v.z,
    )
}

/// Invert a 3×3 row-major matrix using Cramer's rule
///
/// Returns identity if the matrix is singular (determinant ≈ 0).
pub fn invert_3x3(m: &[f32; 9]) -> [f32; 9] {
    let a = m[0];
    let b = m[1];
    let c = m[2];
    let d = m[3];
    let e = m[4];
    let f = m[5];
    let g = m[6];
    let h = m[7];
    let i = m[8];

    let det = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);

    if det.abs() < 1e-30 {
        return [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
    }

    let inv_det = 1.0 / det;

    [
        (e * i - f * h) * inv_det,
        (c * h - b * i) * inv_det,
        (b * f - c * e) * inv_det,
        (f * g - d * i) * inv_det,
        (a * i - c * g) * inv_det,
        (c * d - a * f) * inv_det,
        (d * h - e * g) * inv_det,
        (b * g - a * h) * inv_det,
        (a * e - b * d) * inv_det,
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mat3x3_to_mat4() {
        let m3 = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
        let m4 = mat3x3_to_mat4(&m3);
        assert_eq!(m4.data[0], 1.0);
        assert_eq!(m4.data[3], 0.0); // col 3
        assert_eq!(m4.data[12], 0.0); // row 3
        assert_eq!(m4.data[15], 1.0); // [3][3]
    }

    #[test]
    fn test_left_multiply_identity() {
        let id = Mat4::new_identity();
        let m = Mat4 {
            data: [
                1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 0.0, 0.0, 0.0,
                1.0,
            ],
        };
        let result = left_multiply_mat4(&id, &m);
        for i in 0..16 {
            assert!(
                (result.data[i] - m.data[i]).abs() < 1e-6,
                "mismatch at index {i}"
            );
        }
    }

    #[test]
    fn test_is_identity() {
        assert!(is_identity_mat4(&Mat4::new_identity()));
        let mut m = Mat4::new_identity();
        m.data[3] = 0.1;
        assert!(!is_identity_mat4(&m));
    }

    #[test]
    fn test_transform_mat4() {
        let m = Mat4 {
            data: [
                1.0, 0.0, 0.0, 5.0, // translate x+5
                0.0, 1.0, 0.0, 0.0, //
                0.0, 0.0, 1.0, 0.0, //
                0.0, 0.0, 0.0, 1.0, //
            ],
        };
        let v = Vec3::new(1.0, 2.0, 3.0);
        let result = transform_mat4(&m, v);
        assert!((result.x - 6.0).abs() < 1e-6);
        assert!((result.y - 2.0).abs() < 1e-6);
        assert!((result.z - 3.0).abs() < 1e-6);
    }

    #[test]
    fn test_transform_3x3() {
        // Identity
        let id = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let v = Vec3::new(1.0, 2.0, 3.0);
        let result = transform_3x3(&id, v);
        assert!((result.x - 1.0).abs() < 1e-6);
        assert!((result.y - 2.0).abs() < 1e-6);
        assert!((result.z - 3.0).abs() < 1e-6);
    }

    #[test]
    fn test_invert_3x3_identity() {
        let id = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let inv = invert_3x3(&id);
        for i in 0..9 {
            assert!(
                (inv[i] - id[i]).abs() < 1e-6,
                "mismatch at index {i}: {} vs {}",
                inv[i],
                id[i]
            );
        }
    }

    #[test]
    fn test_invert_3x3_roundtrip() {
        // Scale matrix
        let m = [2.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 4.0];
        let inv = invert_3x3(&m);
        assert!((inv[0] - 0.5).abs() < 1e-6);
        assert!((inv[4] - 1.0 / 3.0).abs() < 1e-6);
        assert!((inv[8] - 0.25).abs() < 1e-6);
    }
}
