//! Kabsch algorithm for optimal rigid-body superposition
//!
//! Given two sets of corresponding 3D points, finds the rotation and
//! translation that minimizes RMSD.

use lin_alg::f32::{Mat4, Vec3, Vec4};

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
    // H[col][row] (column-major) where H = Pᵀ · W · Q
    // H_jk = Σ_i w_i * src_i[j] * tgt_i[k]
    // In column-major: H[col][row] = Σ_i w_i * src_i[row] * tgt_i[col]
    let mut h = [[0.0f32; 3]; 3];
    for i in 0..n {
        let w = weights.map_or(1.0, |ws| ws[i]);
        let s = centered_src[i];
        let t = centered_tgt[i];
        let sv = [s.x, s.y, s.z];
        let tv = [t.x, t.y, t.z];
        for col in 0..3 {
            for row in 0..3 {
                h[col][row] += w * sv[row] * tv[col];
            }
        }
    }

    // 4. SVD of H
    let svd = svd3(&h);

    // 5. Rotation R = V · diag(1, 1, d) · Uᵀ where d = det(V·Uᵀ)
    // Reconstruct V from Vᵀ: V[col][row] = Vt[row][col]
    let vt = &svd.vt;
    let u = &svd.u;

    // Compute det(V · Uᵀ) to handle reflections
    let det_u = u[0][0] * (u[1][1] * u[2][2] - u[1][2] * u[2][1])
        - u[1][0] * (u[0][1] * u[2][2] - u[0][2] * u[2][1])
        + u[2][0] * (u[0][1] * u[1][2] - u[0][2] * u[1][1]);

    // V[col][row] = Vt[row][col], so det(V):
    let v00 = vt[0][0]; let v10 = vt[1][0]; let v20 = vt[2][0];
    let v01 = vt[0][1]; let v11 = vt[1][1]; let v21 = vt[2][1];
    let v02 = vt[0][2]; let v12 = vt[1][2]; let v22 = vt[2][2];
    let det_v = v00 * (v11 * v22 - v12 * v21)
        - v01 * (v10 * v22 - v12 * v20)
        + v02 * (v10 * v21 - v11 * v20);

    let d = if det_u * det_v < 0.0 { -1.0f32 } else { 1.0f32 };

    // R = V · diag(1,1,d) · Uᵀ
    // Build rotation in column-major
    let mut rot = [[0.0f32; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let diag = [1.0, 1.0, d];
            let mut sum = 0.0f32;
            for k in 0..3 {
                sum += vt[i][k] * diag[k] * u[k][j];
            }
            rot[j][i] = sum;
        }
    }

    // 6. Translation t = centroid_target - R · centroid_source
    let r_cs = Vec3::new(
        rot[0][0] * centroid_src.x + rot[1][0] * centroid_src.y + rot[2][0] * centroid_src.z,
        rot[0][1] * centroid_src.x + rot[1][1] * centroid_src.y + rot[2][1] * centroid_src.z,
        rot[0][2] * centroid_src.x + rot[1][2] * centroid_src.y + rot[2][2] * centroid_src.z,
    );
    let translation = centroid_tgt - r_cs;

    // 7. Compute RMSD
    let mut sum_sq = 0.0f32;
    let mut tw = 0.0f32;
    for i in 0..n {
        let w = weights.map_or(1.0, |ws| ws[i]);
        let s = centered_src[i];
        let rotated = Vec3::new(
            rot[0][0] * s.x + rot[1][0] * s.y + rot[2][0] * s.z,
            rot[0][1] * s.x + rot[1][1] * s.y + rot[2][1] * s.z,
            rot[0][2] * s.x + rot[1][2] * s.y + rot[2][2] * s.z,
        );
        let diff = rotated - centered_tgt[i];
        sum_sq += w * (diff.x * diff.x + diff.y * diff.y + diff.z * diff.z);
        tw += w;
    }
    let rmsd = (sum_sq / tw).sqrt();

    // Build column-major Mat4: data[col*4 + row]
    // rot[col][row] is already column-major 3×3 — place directly
    let rotation = Mat4::new([
        rot[0][0], rot[0][1], rot[0][2], 0.0, // col 0
        rot[1][0], rot[1][1], rot[1][2], 0.0, // col 1
        rot[2][0], rot[2][1], rot[2][2], 0.0, // col 2
        0.0,       0.0,       0.0,       1.0, // col 3
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
        .map(|(a, b)| {
            let d = *a - *b;
            d.x * d.x + d.y * d.y + d.z * d.z
        })
        .sum();
    (sum / n as f32).sqrt()
}

/// Apply a Kabsch transformation to coordinates in-place
pub fn apply_transform(coords: &mut [Vec3], result: &KabschResult) {
    let m = &result.rotation;
    let t = &result.translation;
    for coord in coords.iter_mut() {
        let v = m.clone() * Vec4::new(coord.x, coord.y, coord.z, 1.0);
        coord.x = v.x + t.x;
        coord.y = v.y + t.y;
        coord.z = v.z + t.z;
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
