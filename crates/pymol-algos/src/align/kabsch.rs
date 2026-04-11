//! Kabsch algorithm for optimal rigid-body superposition
//!
//! Given two sets of corresponding 3D points, finds the rotation and
//! translation that minimizes RMSD.

use lin_alg::f32::{Mat3, Mat4, Vec3, Vec4};

use crate::AlignError;
use crate::linalg::svd3::svd3;

/// Result of Kabsch superposition
#[derive(Debug, Clone)]
pub struct KabschResult {
    /// 3×3 rotation stored in the upper-left of a Mat4
    pub rotation: Mat4,
    /// Translation vector (applied after rotation)
    pub translation: Vec3,
    /// RMSD after superposition
    pub rmsd: f32,
    /// Number of atom pairs used
    pub n_atoms: usize,
}

/// Compute optimal superposition of source onto target.
///
/// Returns the transformation that maps source → target.
/// Both slices must have the same length (≥ 3).
pub fn kabsch(
    source: &[Vec3],
    target: &[Vec3],
    weights: Option<&[f32]>,
) -> Result<KabschResult, AlignError> {
    let n = source.len();
    if n != target.len() {
        return Err(AlignError::LengthMismatch(n, target.len()));
    }
    if n < 3 {
        return Err(AlignError::TooFewAtoms(n));
    }

    // 1. Compute (weighted) centroids
    let total_weight: f32;
    let centroid_src: Vec3;
    let centroid_tgt: Vec3;

    if let Some(w) = weights {
        total_weight = w.iter().sum();
        let mut cs = Vec3::new(0.0, 0.0, 0.0);
        let mut ct = Vec3::new(0.0, 0.0, 0.0);
        for i in 0..n {
            cs += source[i] * w[i];
            ct += target[i] * w[i];
        }
        centroid_src = cs * (1.0 / total_weight);
        centroid_tgt = ct * (1.0 / total_weight);
    } else {
        total_weight = n as f32;
        let mut cs = Vec3::new(0.0, 0.0, 0.0);
        let mut ct = Vec3::new(0.0, 0.0, 0.0);
        for i in 0..n {
            cs += source[i];
            ct += target[i];
        }
        centroid_src = cs * (1.0 / total_weight);
        centroid_tgt = ct * (1.0 / total_weight);
    }

    // 2. Center both sets
    let centered_src: Vec<Vec3> = source.iter().map(|s| *s - centroid_src).collect();
    let centered_tgt: Vec<Vec3> = target.iter().map(|t| *t - centroid_tgt).collect();

    // 3. Compute cross-covariance matrix H = Σ (w_i * p_i ⊗ q_i)
    // H = Pᵀ · W · Q, column-major
    let mut h_cols = [Vec3::new_zero(); 3];
    for i in 0..n {
        let w = weights.map_or(1.0, |ws| ws[i]);
        let s = centered_src[i] * w;
        let t = centered_tgt[i];
        h_cols[0] += s * t.x;
        h_cols[1] += s * t.y;
        h_cols[2] += s * t.z;
    }
    let h = Mat3::from_cols(h_cols[0], h_cols[1], h_cols[2]);

    // 4. SVD of H
    let h_arr = [
        [h.data[0], h.data[1], h.data[2]],
        [h.data[3], h.data[4], h.data[5]],
        [h.data[6], h.data[7], h.data[8]],
    ];
    let svd = svd3(&h_arr);

    // 5. Rotation R = V · diag(1, 1, d) · Uᵀ where d = det(V·Uᵀ)
    let u_mat = Mat3::from(svd.u);
    let vt_mat = Mat3::from(svd.vt);

    // Compute det(V · Uᵀ) to handle reflections
    let det_u = u_mat.determinant();
    let det_v = vt_mat.transpose().determinant();

    let d = if det_u * det_v < 0.0 { -1.0f32 } else { 1.0f32 };

    // R = V · diag(1,1,d) · Uᵀ
    let diag = Mat3::new([1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, d]);
    let rot_mat = vt_mat.transpose() * diag * u_mat.transpose();

    // 6. Translation t = centroid_target - R · centroid_source
    let translation = centroid_tgt - rot_mat.clone() * centroid_src;

    // 7. Compute RMSD
    let rotated_src: Vec<Vec3> = centered_src.iter().map(|s| rot_mat.clone() * *s).collect();
    let rmsd = if let Some(w) = weights {
        let (sum_sq, tw) = rotated_src
            .iter()
            .zip(centered_tgt.iter())
            .zip(w.iter())
            .fold((0.0f32, 0.0f32), |(ss, tw), ((r, t), &wi)| {
                (ss + wi * (*r - *t).magnitude_squared(), tw + wi)
            });
        (sum_sq / tw).sqrt()
    } else {
        self::rmsd(&rotated_src, &centered_tgt)
    };

    // Build Mat4 from rotation Mat3
    let rd = &rot_mat.data;
    let rotation = Mat4::new([
        rd[0], rd[1], rd[2], 0.0, // col 0
        rd[3], rd[4], rd[5], 0.0, // col 1
        rd[6], rd[7], rd[8], 0.0, // col 2
        0.0,   0.0,   0.0,  1.0, // col 3
    ]);

    Ok(KabschResult {
        rotation,
        translation,
        rmsd,
        n_atoms: n,
    })
}

/// Compute RMSD between two equal-length coordinate sets (no superposition)
pub fn rmsd(coords_a: &[Vec3], coords_b: &[Vec3]) -> f32 {
    assert_eq!(coords_a.len(), coords_b.len());
    let n = coords_a.len();
    if n == 0 {
        return 0.0;
    }
    let sum: f32 = coords_a
        .iter()
        .zip(coords_b.iter())
        .map(|(a, b)| (*a - *b).magnitude_squared())
        .sum();
    (sum / n as f32).sqrt()
}

/// Apply a Kabsch transformation to coordinates in-place
pub fn apply_transform(coords: &mut [Vec3], result: &KabschResult) {
    let m = &result.rotation;
    let t = result.translation;
    for coord in coords.iter_mut() {
        let v = m.clone() * Vec4::new(coord.x, coord.y, coord.z, 1.0);
        *coord = Vec3::new(v.x, v.y, v.z) + t;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn v(x: f32, y: f32, z: f32) -> Vec3 {
        Vec3::new(x, y, z)
    }

    #[test]
    fn test_identity_case() {
        let points = vec![v(0.0, 0.0, 0.0), v(1.0, 0.0, 0.0), v(0.0, 1.0, 0.0), v(0.0, 0.0, 1.0)];
        let result = kabsch(&points, &points, None).unwrap();
        assert!(result.rmsd < 1e-5, "RMSD should be ~0, got {}", result.rmsd);
        assert_eq!(result.n_atoms, 4);
    }

    #[test]
    fn test_pure_translation() {
        let source = vec![v(0.0, 0.0, 0.0), v(1.0, 0.0, 0.0), v(0.0, 1.0, 0.0), v(0.0, 0.0, 1.0)];
        let target: Vec<Vec3> = source.iter().map(|p| *p + Vec3::new(5.0, 3.0, 1.0)).collect();
        let result = kabsch(&source, &target, None).unwrap();
        assert!(result.rmsd < 1e-4, "RMSD should be ~0, got {}", result.rmsd);
        assert!((result.translation.x - 5.0).abs() < 0.1);
        assert!((result.translation.y - 3.0).abs() < 0.1);
        assert!((result.translation.z - 1.0).abs() < 0.1);
    }

    #[test]
    fn test_known_rotation() {
        // Rotate 90° around Z axis
        let source = vec![v(1.0, 0.0, 0.0), v(0.0, 1.0, 0.0), v(-1.0, 0.0, 0.0), v(0.0, -1.0, 0.0), v(0.0, 0.0, 1.0)];
        let target: Vec<Vec3> = source.iter().map(|p| v(-p.y, p.x, p.z)).collect();
        let result = kabsch(&source, &target, None).unwrap();
        assert!(result.rmsd < 1e-3, "RMSD should be ~0, got {}", result.rmsd);
    }

    #[test]
    fn test_apply_transform() {
        let source = vec![v(0.0, 0.0, 0.0), v(1.0, 0.0, 0.0), v(0.0, 1.0, 0.0), v(0.0, 0.0, 1.0)];
        let target: Vec<Vec3> = source.iter().map(|p| *p + Vec3::new(5.0, 3.0, 1.0)).collect();
        let result = kabsch(&source, &target, None).unwrap();

        let mut transformed = source.clone();
        apply_transform(&mut transformed, &result);

        for (t, tgt) in transformed.iter().zip(target.iter()) {
            assert!((t.x - tgt.x).abs() < 0.1, "Transform mismatch: {:?} vs {:?}", t, tgt);
            assert!((t.y - tgt.y).abs() < 0.1, "Transform mismatch: {:?} vs {:?}", t, tgt);
            assert!((t.z - tgt.z).abs() < 0.1, "Transform mismatch: {:?} vs {:?}", t, tgt);
        }
    }

    #[test]
    fn test_length_mismatch() {
        let a = vec![Vec3::new(0.0, 0.0, 0.0); 5];
        let b = vec![Vec3::new(0.0, 0.0, 0.0); 4];
        assert!(kabsch(&a, &b, None).is_err());
    }

    #[test]
    fn test_too_few_atoms() {
        let a = vec![Vec3::new(0.0, 0.0, 0.0); 2];
        let b = vec![Vec3::new(0.0, 0.0, 0.0); 2];
        assert!(kabsch(&a, &b, None).is_err());
    }

    #[test]
    fn test_rmsd_function() {
        let a = vec![v(0.0, 0.0, 0.0), v(1.0, 0.0, 0.0)];
        let b = vec![v(0.0, 0.0, 0.0), v(0.0, 0.0, 0.0)];
        let r = rmsd(&a, &b);
        // sqrt((0 + 1) / 2) = sqrt(0.5) ≈ 0.707
        assert!((r - std::f32::consts::FRAC_1_SQRT_2).abs() < 0.01);
    }

    #[test]
    fn test_reflection_handling() {
        let source = vec![v(1.0, 0.0, 0.0), v(0.0, 1.0, 0.0), v(0.0, 0.0, 1.0), v(1.0, 1.0, 1.0)];
        // Mirror through XY plane
        let target: Vec<Vec3> = source.iter().map(|p| v(p.x, p.y, -p.z)).collect();
        let result = kabsch(&source, &target, None).unwrap();

        // Verify rotation matrix has det = +1 (not a reflection)
        // Column-major: M[row,col] = data[col*4 + row]
        let r = &result.rotation.data;
        let det = r[0] * (r[5] * r[10] - r[9] * r[6])
            - r[1] * (r[4] * r[10] - r[8] * r[6])
            + r[2] * (r[4] * r[9] - r[8] * r[5]);
        assert!((det - 1.0).abs() < 0.1, "det(R) should be +1, got {}", det);
    }
}
