//! Quaternion math for smooth rotation interpolation
//!
//! Provides quaternion-matrix conversion and SLERP for use in
//! movie keyframe interpolation.

use lin_alg::f32::{Mat4, Vec3};

/// A unit quaternion representing a rotation.
///
/// Stored as (w, x, y, z) where w is the scalar part.
#[derive(Debug, Clone, Copy)]
pub struct Quat {
    pub w: f32,
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

impl Quat {
    /// Create a new quaternion
    pub fn new(w: f32, x: f32, y: f32, z: f32) -> Self {
        Self { w, x, y, z }
    }

    /// Identity quaternion (no rotation)
    pub fn identity() -> Self {
        Self {
            w: 1.0,
            x: 0.0,
            y: 0.0,
            z: 0.0,
        }
    }

    /// Extract a quaternion from a rotation matrix.
    ///
    /// The matrix is assumed to be a proper rotation (orthonormal, det=+1).
    /// Uses Shepperd's method for numerical stability.
    ///
    /// lin_alg Mat4 uses column-major storage:
    /// data[0..4] = column 0, data[4..8] = column 1, etc.
    /// So data[col*4 + row] = M[row][col].
    pub fn from_mat4(m: &Mat4) -> Self {
        // m.data[col*4 + row]
        // m00 = data[0], m11 = data[5], m22 = data[10]
        // m01 = data[4], m10 = data[1]
        // m02 = data[8], m20 = data[2]
        // m12 = data[9], m21 = data[6]
        let m00 = m.data[0];
        let m11 = m.data[5];
        let m22 = m.data[10];
        let trace = m00 + m11 + m22;

        if trace > 0.0 {
            let s = (trace + 1.0).sqrt() * 2.0; // s = 4*w
            Self {
                w: 0.25 * s,
                x: (m.data[6] - m.data[9]) / s, // (m21 - m12) / s
                y: (m.data[8] - m.data[2]) / s, // (m02 - m20) / s
                z: (m.data[1] - m.data[4]) / s, // (m10 - m01) / s
            }
        } else if m00 > m11 && m00 > m22 {
            let s = (1.0 + m00 - m11 - m22).sqrt() * 2.0; // s = 4*x
            Self {
                w: (m.data[6] - m.data[9]) / s,
                x: 0.25 * s,
                y: (m.data[4] + m.data[1]) / s, // (m01 + m10) / s
                z: (m.data[8] + m.data[2]) / s, // (m02 + m20) / s
            }
        } else if m11 > m22 {
            let s = (1.0 + m11 - m00 - m22).sqrt() * 2.0; // s = 4*y
            Self {
                w: (m.data[8] - m.data[2]) / s,
                x: (m.data[4] + m.data[1]) / s,
                y: 0.25 * s,
                z: (m.data[9] + m.data[6]) / s, // (m12 + m21) / s
            }
        } else {
            let s = (1.0 + m22 - m00 - m11).sqrt() * 2.0; // s = 4*z
            Self {
                w: (m.data[1] - m.data[4]) / s,
                x: (m.data[8] + m.data[2]) / s,
                y: (m.data[9] + m.data[6]) / s,
                z: 0.25 * s,
            }
        }
    }

    /// Convert this quaternion to a 4x4 rotation matrix.
    pub fn to_mat4(&self) -> Mat4 {
        let Quat { w, x, y, z } = *self;

        let x2 = x + x;
        let y2 = y + y;
        let z2 = z + z;
        let xx = x * x2;
        let xy = x * y2;
        let xz = x * z2;
        let yy = y * y2;
        let yz = y * z2;
        let zz = z * z2;
        let wx = w * x2;
        let wy = w * y2;
        let wz = w * z2;

        // Column-major: data[col*4 + row]
        let mut m = Mat4::new_identity();
        m.data[0] = 1.0 - (yy + zz);
        m.data[1] = xy + wz;
        m.data[2] = xz - wy;
        m.data[3] = 0.0;

        m.data[4] = xy - wz;
        m.data[5] = 1.0 - (xx + zz);
        m.data[6] = yz + wx;
        m.data[7] = 0.0;

        m.data[8] = xz + wy;
        m.data[9] = yz - wx;
        m.data[10] = 1.0 - (xx + yy);
        m.data[11] = 0.0;

        m.data[12] = 0.0;
        m.data[13] = 0.0;
        m.data[14] = 0.0;
        m.data[15] = 1.0;

        m
    }

    /// Dot product of two quaternions.
    pub fn dot(&self, other: &Self) -> f32 {
        self.w * other.w + self.x * other.x + self.y * other.y + self.z * other.z
    }

    /// Normalize this quaternion.
    pub fn normalized(&self) -> Self {
        let len = (self.w * self.w + self.x * self.x + self.y * self.y + self.z * self.z).sqrt();
        if len < 1e-10 {
            Self::identity()
        } else {
            let inv = 1.0 / len;
            Self {
                w: self.w * inv,
                x: self.x * inv,
                y: self.y * inv,
                z: self.z * inv,
            }
        }
    }

    /// Negate quaternion (represents same rotation).
    pub fn negated(&self) -> Self {
        Self {
            w: -self.w,
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }

    /// Spherical linear interpolation between two quaternions.
    ///
    /// Always takes the shortest path (flips if dot product is negative).
    pub fn slerp(a: &Self, b: &Self, t: f32) -> Self {
        let mut dot = a.dot(b);

        // If dot is negative, negate one quaternion to take the shorter path
        let b = if dot < 0.0 {
            dot = -dot;
            b.negated()
        } else {
            *b
        };

        // Clamp dot to valid range for acos
        let dot = dot.min(1.0);

        if dot > 0.9995 {
            // Very close: use linear interpolation for stability
            Self {
                w: a.w + (b.w - a.w) * t,
                x: a.x + (b.x - a.x) * t,
                y: a.y + (b.y - a.y) * t,
                z: a.z + (b.z - a.z) * t,
            }
            .normalized()
        } else {
            let theta = dot.acos();
            let sin_theta = theta.sin();
            let s0 = ((1.0 - t) * theta).sin() / sin_theta;
            let s1 = (t * theta).sin() / sin_theta;

            Self {
                w: a.w * s0 + b.w * s1,
                x: a.x * s0 + b.x * s1,
                y: a.y * s0 + b.y * s1,
                z: a.z * s0 + b.z * s1,
            }
        }
    }
}

/// Linearly interpolate two f32 values.
pub fn lerp(a: f32, b: f32, t: f32) -> f32 {
    a + (b - a) * t
}

/// Linearly interpolate two Vec3 values.
pub fn lerp_vec3(a: &Vec3, b: &Vec3, t: f32) -> Vec3 {
    Vec3::new(
        lerp(a.x, b.x, t),
        lerp(a.y, b.y, t),
        lerp(a.z, b.z, t),
    )
}

/// Cubic Hermite interpolation for f32 values.
///
/// Given values and tangents at t=0 and t=1, compute the interpolated
/// value at parameter `t`.
pub fn hermite(p0: f32, m0: f32, p1: f32, m1: f32, t: f32) -> f32 {
    let t2 = t * t;
    let t3 = t2 * t;
    let h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
    let h10 = t3 - 2.0 * t2 + t;
    let h01 = -2.0 * t3 + 3.0 * t2;
    let h11 = t3 - t2;
    h00 * p0 + h10 * m0 + h01 * p1 + h11 * m1
}

/// Cubic Hermite interpolation for Vec3 values.
pub fn hermite_vec3(p0: &Vec3, m0: &Vec3, p1: &Vec3, m1: &Vec3, t: f32) -> Vec3 {
    Vec3::new(
        hermite(p0.x, m0.x, p1.x, m1.x, t),
        hermite(p0.y, m0.y, p1.y, m1.y, t),
        hermite(p0.z, m0.z, p1.z, m1.z, t),
    )
}

/// Compute Catmull-Rom tangent for a point in a sequence.
///
/// Given the previous value, current value, and next value,
/// compute the tangent at the current point.
pub fn catmull_rom_tangent(prev: f32, next: f32) -> f32 {
    (next - prev) * 0.5
}

/// Compute Catmull-Rom tangent for Vec3.
pub fn catmull_rom_tangent_vec3(prev: &Vec3, next: &Vec3) -> Vec3 {
    Vec3::new(
        catmull_rom_tangent(prev.x, next.x),
        catmull_rom_tangent(prev.y, next.y),
        catmull_rom_tangent(prev.z, next.z),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f32::consts::PI;

    #[test]
    fn test_identity_roundtrip() {
        let m = Mat4::new_identity();
        let q = Quat::from_mat4(&m);
        let m2 = q.to_mat4();

        for i in 0..16 {
            assert!(
                (m.data[i] - m2.data[i]).abs() < 1e-5,
                "index {}: {} vs {}",
                i,
                m.data[i],
                m2.data[i]
            );
        }
    }

    #[test]
    fn test_rotation_roundtrip() {
        // Create a 90-degree rotation around Y
        let angle = PI / 2.0;
        let c = angle.cos();
        let s = angle.sin();

        let mut m = Mat4::new_identity();
        // Rotation around Y: column-major
        m.data[0] = c;
        m.data[2] = -s;
        m.data[8] = s;
        m.data[10] = c;

        let q = Quat::from_mat4(&m);
        let m2 = q.to_mat4();

        for i in 0..16 {
            assert!(
                (m.data[i] - m2.data[i]).abs() < 1e-5,
                "index {}: {} vs {}",
                i,
                m.data[i],
                m2.data[i]
            );
        }
    }

    #[test]
    fn test_slerp_endpoints() {
        let q1 = Quat::identity();
        let angle = PI / 4.0;
        let q2 = Quat::new((angle / 2.0).cos(), 0.0, (angle / 2.0).sin(), 0.0);

        // t=0 should give q1
        let r0 = Quat::slerp(&q1, &q2, 0.0);
        assert!((r0.w - q1.w).abs() < 1e-5);
        assert!((r0.x - q1.x).abs() < 1e-5);
        assert!((r0.y - q1.y).abs() < 1e-5);
        assert!((r0.z - q1.z).abs() < 1e-5);

        // t=1 should give q2
        let r1 = Quat::slerp(&q1, &q2, 1.0);
        assert!((r1.w - q2.w).abs() < 1e-5);
        assert!((r1.x - q2.x).abs() < 1e-5);
        assert!((r1.y - q2.y).abs() < 1e-5);
        assert!((r1.z - q2.z).abs() < 1e-5);
    }

    #[test]
    fn test_slerp_midpoint() {
        let q1 = Quat::identity();
        // 90 degrees around Y
        let angle = PI / 2.0;
        let q2 = Quat::new((angle / 2.0).cos(), 0.0, (angle / 2.0).sin(), 0.0);

        let mid = Quat::slerp(&q1, &q2, 0.5);
        // Should be 45 degrees around Y
        let expected_angle = PI / 4.0;
        let expected = Quat::new(
            (expected_angle / 2.0).cos(),
            0.0,
            (expected_angle / 2.0).sin(),
            0.0,
        );

        assert!((mid.w - expected.w).abs() < 1e-4);
        assert!((mid.y - expected.y).abs() < 1e-4);
    }

    #[test]
    fn test_hermite_endpoints() {
        assert!((hermite(0.0, 1.0, 1.0, 1.0, 0.0) - 0.0).abs() < 1e-5);
        assert!((hermite(0.0, 1.0, 1.0, 1.0, 1.0) - 1.0).abs() < 1e-5);
    }

    #[test]
    fn test_hermite_linear() {
        // Zero tangents should give smooth but still monotonic interpolation
        let v = hermite(0.0, 0.0, 1.0, 0.0, 0.5);
        assert!(v > 0.0 && v < 1.0);
    }
}
