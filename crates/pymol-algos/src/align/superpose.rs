//! Iterative superposition with outlier rejection
//!
//! Combines Kabsch superposition (Acta Crystallogr. A32:922-923, 1976) with
//! iterative outlier rejection for robust structural alignment.

use lin_alg::f32::Vec3;

use super::kabsch::{self, KabschResult};
use crate::AlignError;

/// Parameters for iterative superposition
#[derive(Debug, Clone)]
pub struct SuperposeParams {
    /// Number of outlier rejection cycles (0 = no rejection)
    pub cycles: u32,
    /// Outlier rejection cutoff — an atom pair is rejected when
    /// `distance / current_RMSD > cutoff` (PyMOL convention).
    /// Default: 2.0
    pub cutoff: f32,
}

impl Default for SuperposeParams {
    fn default() -> Self {
        Self {
            cycles: 5,
            cutoff: 2.0,
        }
    }
}

/// Result of iterative superposition
#[derive(Debug, Clone)]
pub struct SuperposeResult {
    /// Final Kabsch result (rotation + translation)
    pub transform: KabschResult,
    /// Initial RMSD (before outlier rejection)
    pub initial_rmsd: f32,
    /// Final RMSD (after outlier rejection)
    pub final_rmsd: f32,
    /// Number of cycles performed
    pub cycles_performed: u32,
    /// Number of atom pairs rejected as outliers
    pub n_rejected: usize,
    /// Number of atom pairs used in final fit
    pub n_aligned: usize,
    /// Indices of pairs that were NOT rejected (into original pairs array)
    pub aligned_pairs: Vec<(usize, usize)>,
}

/// Iterative superposition with outlier rejection.
///
/// `pairs` is a list of (source_index, target_index) correspondences.
/// Coordinates are looked up in `source_coords` and `target_coords`.
pub fn superpose(
    source_coords: &[Vec3],
    target_coords: &[Vec3],
    pairs: &[(usize, usize)],
    params: &SuperposeParams,
) -> Result<SuperposeResult, AlignError> {
    if pairs.len() < 3 {
        return Err(AlignError::TooFewAtoms(pairs.len()));
    }

    // Extract initial paired coordinates
    let mut active_pairs: Vec<(usize, usize)> = pairs.to_vec();

    let extract = |pairs: &[(usize, usize)]| -> (Vec<Vec3>, Vec<Vec3>) {
        let src: Vec<Vec3> = pairs.iter().map(|&(si, _)| source_coords[si]).collect();
        let tgt: Vec<Vec3> = pairs.iter().map(|&(_, ti)| target_coords[ti]).collect();
        (src, tgt)
    };

    // Initial Kabsch
    let (src, tgt) = extract(&active_pairs);
    let initial_result = kabsch::kabsch(&src, &tgt, None)?;
    let initial_rmsd = initial_result.rmsd;

    let mut current_result = initial_result;
    let mut cycles_performed = 0;

    // Iterative outlier rejection
    for cycle in 0..params.cycles {
        if active_pairs.len() < 3 {
            break;
        }

        // Apply current transform to source pairs and compute per-pair distances
        let (src, tgt) = extract(&active_pairs);
        let mut transformed_src = src.clone();
        kabsch::apply_transform(&mut transformed_src, &current_result);

        // Reject pairs where distance/RMSD > cutoff (PyMOL convention).
        let rms = current_result.rmsd;
        let mut keep = Vec::new();
        if rms > 1e-6 {
            for (i, (ts, t)) in transformed_src.iter().zip(tgt.iter()).enumerate() {
                let d = *ts - *t;
                let dist = d.magnitude();
                if (dist / rms) <= params.cutoff {
                    keep.push(i);
                }
            }
        } else {
            // RMSD ≈ 0 → perfect fit, keep everything
            keep.extend(0..active_pairs.len());
        }

        // If no pairs were rejected, we've converged
        if keep.len() == active_pairs.len() {
            cycles_performed = cycle + 1;
            break;
        }

        // If all pairs would be rejected, stop
        if keep.len() < 3 {
            if keep.is_empty() {
                return Err(AlignError::AllRejected);
            }
            // Keep what we have but don't re-fit
            cycles_performed = cycle + 1;
            break;
        }

        // Filter to kept pairs
        active_pairs = keep.iter().map(|&i| active_pairs[i]).collect();

        // Re-run Kabsch on remaining pairs
        let (src, tgt) = extract(&active_pairs);
        current_result = kabsch::kabsch(&src, &tgt, None)?;
        cycles_performed = cycle + 1;
    }

    let n_aligned = active_pairs.len();
    let n_rejected = pairs.len() - n_aligned;

    Ok(SuperposeResult {
        final_rmsd: current_result.rmsd,
        transform: current_result,
        initial_rmsd,
        cycles_performed,
        n_rejected,
        n_aligned,
        aligned_pairs: active_pairs,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn v(x: f32, y: f32, z: f32) -> Vec3 {
        Vec3::new(x, y, z)
    }

    #[test]
    fn test_no_outliers() {
        let source = vec![v(0.0, 0.0, 0.0), v(1.0, 0.0, 0.0), v(0.0, 1.0, 0.0), v(0.0, 0.0, 1.0)];
        let target = source.clone();
        let pairs: Vec<(usize, usize)> = (0..4).map(|i| (i, i)).collect();

        let result = superpose(&source, &target, &pairs, &SuperposeParams::default()).unwrap();
        assert!(result.final_rmsd < 1e-4);
        assert_eq!(result.n_rejected, 0);
        assert_eq!(result.n_aligned, 4);
    }

    #[test]
    fn test_with_outliers() {
        let source: Vec<Vec3> = vec![
            v(0.0, 0.0, 0.0), v(1.0, 0.0, 0.0), v(0.0, 1.0, 0.0), v(0.0, 0.0, 1.0),
            v(1.0, 1.0, 0.0), v(1.0, 0.0, 1.0), v(0.0, 1.0, 1.0), v(1.0, 1.0, 1.0),
            v(0.5, 0.5, 0.5), v(0.5, 0.5, 0.0),
        ];
        let target: Vec<Vec3> = vec![
            v(0.0, 0.0, 0.0), v(1.0, 0.0, 0.0), v(0.0, 1.0, 0.0), v(0.0, 0.0, 1.0),
            v(1.0, 1.0, 0.0), v(1.0, 0.0, 1.0), v(0.0, 1.0, 1.0), v(1.0, 1.0, 1.0),
            v(20.5, 20.5, 20.5), v(20.5, 20.5, 20.0),
        ];

        let pairs: Vec<(usize, usize)> = (0..10).map(|i| (i, i)).collect();

        let result = superpose(
            &source, &target, &pairs,
            &SuperposeParams { cycles: 5, cutoff: 2.0 },
        ).unwrap();

        assert!(result.n_rejected >= 2, "Expected ≥2 rejected, got {}", result.n_rejected);
        assert!(result.final_rmsd < result.initial_rmsd);
    }

    #[test]
    fn test_zero_cycles() {
        let source = vec![v(0.0, 0.0, 0.0), v(1.0, 0.0, 0.0), v(0.0, 1.0, 0.0), v(0.0, 0.0, 1.0)];
        let target = source.clone();
        let pairs: Vec<(usize, usize)> = (0..4).map(|i| (i, i)).collect();

        let result = superpose(
            &source, &target, &pairs,
            &SuperposeParams { cycles: 0, cutoff: 2.0 },
        ).unwrap();

        assert_eq!(result.cycles_performed, 0);
        assert_eq!(result.n_rejected, 0);
    }

    #[test]
    fn test_too_few_pairs() {
        let source = vec![v(0.0, 0.0, 0.0), v(1.0, 0.0, 0.0)];
        let target = source.clone();
        let pairs: Vec<(usize, usize)> = (0..2).map(|i| (i, i)).collect();

        assert!(superpose(&source, &target, &pairs, &SuperposeParams::default()).is_err());
    }
}
