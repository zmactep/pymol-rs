//! Kabsch algorithm for optimal rigid-body superposition
//!
//! Given two sets of corresponding 3D points, finds the rotation and
//! translation that minimizes RMSD.

use lin_alg::f32::{Mat4, Vec3};

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
    source: &[[f32; 3]],
    target: &[[f32; 3]],
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
    let centroid_src: [f32; 3];
    let centroid_tgt: [f32; 3];

    if let Some(w) = weights {
        total_weight = w.iter().sum();
        let mut cs = [0.0f32; 3];
        let mut ct = [0.0f32; 3];
        for i in 0..n {
            for k in 0..3 {
                cs[k] += w[i] * source[i][k];
                ct[k] += w[i] * target[i][k];
            }
        }
        centroid_src = [cs[0] / total_weight, cs[1] / total_weight, cs[2] / total_weight];
        centroid_tgt = [ct[0] / total_weight, ct[1] / total_weight, ct[2] / total_weight];
    } else {
        total_weight = n as f32;
        let mut cs = [0.0f32; 3];
        let mut ct = [0.0f32; 3];
        for i in 0..n {
            for k in 0..3 {
                cs[k] += source[i][k];
                ct[k] += target[i][k];
            }
        }
        centroid_src = [cs[0] / total_weight, cs[1] / total_weight, cs[2] / total_weight];
        centroid_tgt = [ct[0] / total_weight, ct[1] / total_weight, ct[2] / total_weight];
    }

    // 2. Center both sets
    let mut centered_src = vec![[0.0f32; 3]; n];
    let mut centered_tgt = vec![[0.0f32; 3]; n];
    for i in 0..n {
        for k in 0..3 {
            centered_src[i][k] = source[i][k] - centroid_src[k];
            centered_tgt[i][k] = target[i][k] - centroid_tgt[k];
        }
    }

    // 3. Compute cross-covariance matrix H = Σ (w_i * p_i ⊗ q_i)
    // H[col][row] (column-major) where H = Pᵀ · W · Q
    // H_jk = Σ_i w_i * src_i[j] * tgt_i[k]
    // In column-major: H[col][row] = Σ_i w_i * src_i[row] * tgt_i[col]
    let mut h = [[0.0f32; 3]; 3];
    for i in 0..n {
        let w = weights.map_or(1.0, |ws| ws[i]);
        for col in 0..3 {
            for row in 0..3 {
                h[col][row] += w * centered_src[i][row] * centered_tgt[i][col];
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
    // V·Uᵀ in column-major: for each output col j, output row i:
    // (V·Uᵀ)[j][i] = Σ_k V[k][i] * Uᵀ[j][k] = Σ_k Vt[i][k] * U[j][k]
    // But simpler: det(V·Uᵀ) = det(V) * det(U)
    // Since U and V are orthogonal, det = ±1, so det(V·Uᵀ) = det_v * det_u
    let det_u = u[0][0] * (u[1][1] * u[2][2] - u[1][2] * u[2][1])
        - u[1][0] * (u[0][1] * u[2][2] - u[0][2] * u[2][1])
        + u[2][0] * (u[0][1] * u[1][2] - u[0][2] * u[1][1]);

    // V[col][row] = Vt[row][col], so det(V):
    // V columns are: col0 = [vt[0][0], vt[1][0], vt[2][0]], etc.
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
            // R[col=j][row=i] = Σ_k V[k][i] * diag[k] * Uᵀ[j][k]
            // V[k][i] = Vt[i][k]
            // Uᵀ[j][k] = U[k][j] ... no wait.
            // Uᵀ[col][row] = U[row][col]
            // So Uᵀ[j][k] = U[k][j]
            let diag = [1.0, 1.0, d];
            let mut sum = 0.0f32;
            for k in 0..3 {
                sum += vt[i][k] * diag[k] * u[k][j];
            }
            rot[j][i] = sum;
        }
    }

    // 6. Translation t = centroid_target - R · centroid_source
    let r_cs = [
        rot[0][0] * centroid_src[0] + rot[1][0] * centroid_src[1] + rot[2][0] * centroid_src[2],
        rot[0][1] * centroid_src[0] + rot[1][1] * centroid_src[1] + rot[2][1] * centroid_src[2],
        rot[0][2] * centroid_src[0] + rot[1][2] * centroid_src[1] + rot[2][2] * centroid_src[2],
    ];
    let translation = Vec3::new(
        centroid_tgt[0] - r_cs[0],
        centroid_tgt[1] - r_cs[1],
        centroid_tgt[2] - r_cs[2],
    );

    // 7. Compute RMSD
    // Apply rotation to centered source, compare to centered target
    let mut sum_sq = 0.0f32;
    let mut tw = 0.0f32;
    for i in 0..n {
        let w = weights.map_or(1.0, |ws| ws[i]);
        let rotated = [
            rot[0][0] * centered_src[i][0] + rot[1][0] * centered_src[i][1] + rot[2][0] * centered_src[i][2],
            rot[0][1] * centered_src[i][0] + rot[1][1] * centered_src[i][1] + rot[2][1] * centered_src[i][2],
            rot[0][2] * centered_src[i][0] + rot[1][2] * centered_src[i][1] + rot[2][2] * centered_src[i][2],
        ];
        for k in 0..3 {
            let d = rotated[k] - centered_tgt[i][k];
            sum_sq += w * d * d;
        }
        tw += w;
    }
    let rmsd = (sum_sq / tw).sqrt();

    // Build Mat4 rotation (row-major: data[row*4 + col])
    // rot[col][row] is our column-major 3×3, so element(i,j) = rot[j][i]
    let rotation = Mat4::new([
        rot[0][0], rot[1][0], rot[2][0], 0.0,  // row 0
        rot[0][1], rot[1][1], rot[2][1], 0.0,  // row 1
        rot[0][2], rot[1][2], rot[2][2], 0.0,  // row 2
        0.0,       0.0,       0.0,       1.0,  // row 3
    ]);

    Ok(KabschResult {
        rotation,
        translation,
        rmsd,
        n_atoms: n,
    })
}

/// Compute RMSD between two equal-length coordinate sets (no superposition)
pub fn rmsd(coords_a: &[[f32; 3]], coords_b: &[[f32; 3]]) -> f32 {
    assert_eq!(coords_a.len(), coords_b.len());
    let n = coords_a.len();
    if n == 0 {
        return 0.0;
    }
    let sum: f32 = coords_a
        .iter()
        .zip(coords_b.iter())
        .map(|(a, b)| {
            let dx = a[0] - b[0];
            let dy = a[1] - b[1];
            let dz = a[2] - b[2];
            dx * dx + dy * dy + dz * dz
        })
        .sum();
    (sum / n as f32).sqrt()
}

/// Apply a Kabsch transformation to coordinates in-place
pub fn apply_transform(coords: &mut [[f32; 3]], result: &KabschResult) {
    let r = &result.rotation.data;
    let t = &result.translation;
    // Mat4 is row-major: data[row*4 + col]
    // data[0]=r00, data[1]=r01, data[2]=r02
    // data[4]=r10, data[5]=r11, data[6]=r12
    // data[8]=r20, data[9]=r21, data[10]=r22
    for coord in coords.iter_mut() {
        let x = coord[0];
        let y = coord[1];
        let z = coord[2];
        coord[0] = r[0] * x + r[1] * y + r[2] * z + t.x;
        coord[1] = r[4] * x + r[5] * y + r[6] * z + t.y;
        coord[2] = r[8] * x + r[9] * y + r[10] * z + t.z;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_identity_case() {
        let points = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let result = kabsch(&points, &points, None).unwrap();
        assert!(result.rmsd < 1e-5, "RMSD should be ~0, got {}", result.rmsd);
        assert_eq!(result.n_atoms, 4);
    }

    #[test]
    fn test_pure_translation() {
        let source = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let target: Vec<[f32; 3]> = source.iter().map(|p| [p[0] + 5.0, p[1] + 3.0, p[2] + 1.0]).collect();
        let result = kabsch(&source, &target, None).unwrap();
        assert!(result.rmsd < 1e-4, "RMSD should be ~0, got {}", result.rmsd);
        assert!((result.translation.x - 5.0).abs() < 0.1);
        assert!((result.translation.y - 3.0).abs() < 0.1);
        assert!((result.translation.z - 1.0).abs() < 0.1);
    }

    #[test]
    fn test_known_rotation() {
        // Rotate 90° around Z axis
        let source = vec![
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let target: Vec<[f32; 3]> = source
            .iter()
            .map(|p| [-p[1], p[0], p[2]]) // 90° Z rotation
            .collect();
        let result = kabsch(&source, &target, None).unwrap();
        assert!(result.rmsd < 1e-3, "RMSD should be ~0, got {}", result.rmsd);
    }

    #[test]
    fn test_apply_transform() {
        let source = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let target: Vec<[f32; 3]> = source.iter().map(|p| [p[0] + 5.0, p[1] + 3.0, p[2] + 1.0]).collect();
        let result = kabsch(&source, &target, None).unwrap();

        let mut transformed = source.clone();
        apply_transform(&mut transformed, &result);

        for (t, tgt) in transformed.iter().zip(target.iter()) {
            for k in 0..3 {
                assert!(
                    (t[k] - tgt[k]).abs() < 0.1,
                    "Transform mismatch: {:?} vs {:?}",
                    t, tgt
                );
            }
        }
    }

    #[test]
    fn test_length_mismatch() {
        let a = vec![[0.0; 3]; 5];
        let b = vec![[0.0; 3]; 4];
        assert!(kabsch(&a, &b, None).is_err());
    }

    #[test]
    fn test_too_few_atoms() {
        let a = vec![[0.0; 3]; 2];
        let b = vec![[0.0; 3]; 2];
        assert!(kabsch(&a, &b, None).is_err());
    }

    #[test]
    fn test_rmsd_function() {
        let a = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let b = vec![[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        let r = rmsd(&a, &b);
        // sqrt((0 + 1) / 2) = sqrt(0.5) ≈ 0.707
        assert!((r - 0.7071).abs() < 0.01);
    }

    #[test]
    fn test_reflection_handling() {
        // Points that would require a reflection — Kabsch should give rotation only
        let source = vec![
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
        ];
        // Mirror through XY plane
        let target: Vec<[f32; 3]> = source.iter().map(|p| [p[0], p[1], -p[2]]).collect();
        let result = kabsch(&source, &target, None).unwrap();

        // Verify rotation matrix has det = +1 (not a reflection)
        let r = &result.rotation.data;
        let det = r[0] * (r[5] * r[10] - r[6] * r[9])
            - r[4] * (r[1] * r[10] - r[2] * r[9])
            + r[8] * (r[1] * r[6] - r[2] * r[5]);
        assert!((det - 1.0).abs() < 0.1, "det(R) should be +1, got {}", det);
    }
}
