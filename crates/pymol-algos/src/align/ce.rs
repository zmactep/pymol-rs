//! Combinatorial Extension (CE) structural alignment
//!
//! Aligns protein structures by Cα distance-matrix fingerprints and
//! combinatorial path extension.
//!
//! # References
//!
//! - Shindyalov IN, Bourne PE (1998). "Protein structure alignment by
//!   incremental combinatorial extension (CE) of the optimal path."
//!   Protein Engineering 11(9):739-747.
//! - Zhang Y, Skolnick J (2004). "Scoring function for automated assessment
//!   of protein structure template quality." Proteins 57(4):702-710.
//!   (Z-score normalization)

use std::cmp::Ordering;

use lin_alg::f32::Vec3;

use super::kabsch;
use crate::AlignError;

/// Maximum candidate paths retained during search.
const MAX_KEPT: usize = 20;

/// Parameters for the CE structural alignment algorithm.
///
/// Default values match those recommended by Shindyalov & Bourne (1998).
#[derive(Debug, Clone)]
pub struct CeParams {
    /// Fragment window size (number of residues per aligned fragment pair)
    pub win_size: usize,
    /// Maximum gap allowed between consecutive AFPs in the path
    pub gap_max: usize,
    /// Distance cutoff D_0 for single-AFP similarity (Å) — Eq. 9
    pub d0: f32,
    /// Distance cutoff D_1 for path extension and whole-path scoring (Å) — Eq. 10–11
    pub d1: f32,
    /// Maximum number of best paths to keep
    pub max_paths: usize,
}

impl Default for CeParams {
    fn default() -> Self {
        Self {
            win_size: 8,
            gap_max: 30,
            d0: 3.0,
            d1: 4.0,
            max_paths: MAX_KEPT,
        }
    }
}

/// Result of CE structural alignment
#[derive(Debug, Clone)]
pub struct CeResult {
    /// Aligned residue pairs: (source_index, target_index) into the
    /// original Cα coordinate arrays
    pub pairs: Vec<(usize, usize)>,
    /// Number of aligned residues
    pub n_aligned: usize,
    /// RMSD of the aligned pairs (after Kabsch superposition)
    pub rmsd: f32,
    /// Z-score of the alignment (quality metric)
    pub z_score: f32,
}

/// Perform Combinatorial Extension structural alignment.
///
/// Takes two sets of Cα coordinates and finds the best structural alignment
/// using fragment-based comparison and combinatorial extension.
///
/// Returns aligned residue pairs that can be fed to `superpose()` for
/// the final rigid-body fit with outlier rejection.
pub fn ce_align(
    source_ca: &[Vec3],
    target_ca: &[Vec3],
    params: &CeParams,
) -> Result<CeResult, AlignError> {
    let len_a = source_ca.len();
    let len_b = target_ca.len();
    let w = params.win_size;

    if len_a < w || len_b < w {
        return Err(AlignError::TooFewAtoms(len_a.min(len_b)));
    }

    // Step 1: Compute intra-molecular Cα distance matrices
    let dm_a = distance_matrix(source_ca);
    let dm_b = distance_matrix(target_ca);

    // Step 2: Compute single-AFP similarity matrix D_ii (Eq. 7 variant —
    // full non-neighboring intra-fragment distances)
    let sim = afp_similarity_matrix(&dm_a, &dm_b, len_a, len_b, w);

    // Step 3: Find optimal AFP paths by combinatorial extension (Eq. 9–11)
    let paths = extend_paths(&sim, &dm_a, &dm_b, len_a, len_b, params);

    if paths.is_empty() {
        return Err(AlignError::NoMatches);
    }

    // Step 4: Select path with best Kabsch RMSD
    // NOTE: The paper's final optimization (gap relocation ±m/2 and
    // iterative Needleman-Wunsch DP refinement) is not implemented.
    let (pairs, rmsd) = select_best_path(source_ca, target_ca, &paths, w)
        .ok_or(AlignError::NoMatches)?;

    let n_aligned = pairs.len();
    if n_aligned < 3 {
        return Err(AlignError::TooFewAtoms(n_aligned));
    }

    let z = z_score(n_aligned, rmsd, len_a, len_b);

    Ok(CeResult {
        pairs,
        n_aligned,
        rmsd,
        z_score: z,
    })
}

// ============================================================================
// Intra-molecular distance matrix
// ============================================================================

/// Compute symmetric NxN pairwise Cα distance matrix.
fn distance_matrix(coords: &[Vec3]) -> Vec<Vec<f32>> {
    let n = coords.len();
    let mut dm = vec![vec![0.0f32; n]; n];
    for i in 0..n {
        for j in (i + 1)..n {
            let dist = (coords[i] - coords[j]).magnitude();
            dm[i][j] = dist;
            dm[j][i] = dist;
        }
    }
    dm
}

// ============================================================================
// AFP similarity matrix — single-AFP distance D_ii (Eq. 7 variant)
// ============================================================================

/// Compute single-AFP similarity scores D_ii using the full non-neighboring
/// distance set (Eq. 7 variant, Shindyalov & Bourne 1998).
///
/// For each pair of starting positions (i_a, i_b) — one fragment from each
/// protein — computes the average absolute difference of all unique
/// non-neighboring intra-fragment Cα distance pairs. This uses
/// (m-1)(m-2)/2 pairs per fragment (all |k-l| >= 2), unlike the
/// "independent" set (Eq. 6) which uses only m anti-diagonal distances.
///
/// Lower scores indicate more similar local structure.
/// A value of -1.0 indicates the fragment extends beyond the sequence end.
fn afp_similarity_matrix(
    dm_a: &[Vec<f32>],
    dm_b: &[Vec<f32>],
    len_a: usize,
    len_b: usize,
    win_size: usize,
) -> Vec<Vec<f32>> {
    let w = win_size;
    let pair_count = ((w - 1) * (w - 2)) as f32 / 2.0;

    let mut sim = vec![vec![-1.0f32; len_b]; len_a];

    for i_a in 0..=len_a.saturating_sub(w) {
        for i_b in 0..=len_b.saturating_sub(w) {
            let score: f32 = (0..(w - 2))
                .flat_map(|row| ((row + 2)..w).map(move |col| (row, col)))
                .map(|(row, col)| {
                    (dm_a[i_a + row][i_a + col] - dm_b[i_b + row][i_b + col]).abs()
                })
                .sum();

            sim[i_a][i_b] = score / pair_count;
        }
    }

    sim
}

// ============================================================================
// Path helpers
// ============================================================================

/// Resolve a gap index to AFP anchor positions (Eq. 1–5).
///
/// Implements path continuity conditions: no gap (Eq. 1), gap in A only
/// (Eq. 2), or gap in B only (Eq. 3), bounded by the gap limit G (Eq. 4–5).
/// Gap indices alternate insertion between structures:
/// g=0 → +0 in B, g=1 → +1 in A, g=2 → +1 in B, g=3 → +2 in A, ...
fn resolve_gap(last: (usize, usize), win_size: usize, gap_idx: usize) -> (usize, usize) {
    let offset = (gap_idx + 1) / 2;
    let (mut j_a, mut j_b) = (last.0 + win_size, last.1 + win_size);
    if (gap_idx + 1) % 2 == 0 {
        j_a += offset;
    } else {
        j_b += offset;
    }
    (j_a, j_b)
}

/// Compute inter-fragment distance D_ij between a candidate AFP at (j_a, j_b)
/// and all existing path anchors, using the independent distance set (Eq. 6).
///
/// For each pair of AFPs, computes m anti-diagonal distance comparisons
/// (one per residue position), averaged over m × n_path. This corresponds
/// to the "Ind." (independent) method in Table I of the paper.
fn fragment_d_score(
    path: &[(usize, usize)],
    j_a: usize,
    j_b: usize,
    dm_a: &[Vec<f32>],
    dm_b: &[Vec<f32>],
    w: usize,
) -> f64 {
    let total: f64 = path
        .iter()
        .map(|&(pa, pb)| {
            let first = (dm_a[pa][j_a] as f64 - dm_b[pb][j_b] as f64).abs();
            let last = (dm_a[pa + w - 1][j_a + w - 1] as f64
                - dm_b[pb + w - 1][j_b + w - 1] as f64)
                .abs();
            let mid: f64 = (1..(w - 1))
                .map(|k| {
                    (dm_a[pa + k][j_a + w - 1 - k] as f64
                        - dm_b[pb + k][j_b + w - 1 - k] as f64)
                        .abs()
                })
                .sum();
            first + last + mid
        })
        .sum();
    total / (w * path.len()) as f64
}

// ============================================================================
// Combinatorial extension of the alignment path (Eq. 9–11)
// ============================================================================

/// Grow an AFP path greedily from a seed.
///
/// Returns all intermediate path snapshots (at each extension step) so that
/// `select_best_path` can evaluate truncation points at different lengths.
/// This matches the original CE algorithm's behavior of considering paths
/// at every length, not just the final greedy result.
fn grow_path(
    sim: &[Vec<f32>],
    dm_a: &[Vec<f32>],
    dm_b: &[Vec<f32>],
    len_a: usize,
    len_b: usize,
    seed: (usize, usize),
    seed_score: f64,
    params: &CeParams,
    norm: &[f64],
) -> Vec<(Vec<(usize, usize)>, f64)> {
    let w = params.win_size;
    let d0 = params.d0 as f64;
    let d1 = params.d1 as f64;
    let win_sum = ((w - 1) * (w - 2)) / 2;
    let n_gaps = params.gap_max * 2 + 1;

    let mut path = vec![seed];
    let mut cumulative_scores: Vec<f64> = Vec::new();
    let mut snapshots: Vec<(Vec<(usize, usize)>, f64)> = Vec::new();

    loop {
        let last = *path.last().unwrap();

        // Find the best gap extension: single AFP must satisfy D_an < D_0
        // (Eq. 9), then inter-fragment D-score must satisfy Eq. 10.
        let best = (0..n_gaps)
            .filter_map(|g| {
                let (j_a, j_b) = resolve_gap(last, w, g);
                if j_a + w > len_a || j_b + w > len_b {
                    return None;
                }
                let afp_sim = sim[j_a][j_b] as f64;
                if afp_sim >= d0 || afp_sim < 0.0 {
                    return None; // Eq. 9: D_an < D_0
                }
                let d = fragment_d_score(&path, j_a, j_b, dm_a, dm_b, w);
                (d < d1).then_some(((j_a, j_b), afp_sim, d)) // Eq. 10
            })
            .min_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(Ordering::Equal));

        let Some(((j_a, j_b), afp_sim, d_score)) = best else {
            break;
        };

        // Incrementally update whole-path score (approximates Eq. 11:
        // (1/n²) ΣΣ D_ij < D_1), combining inter-AFP D-scores (Eq. 6)
        // with intra-AFP similarities (Eq. 7).
        let score1 = (d_score * (w * path.len()) as f64 + afp_sim * win_sum as f64)
            / (w * path.len() + win_sum) as f64;

        let prev = cumulative_scores.last().copied().unwrap_or(seed_score);
        let step = path.len();
        let cumulative =
            (prev * norm[step - 1] + score1 * (norm[step] - norm[step - 1])) / norm[step];

        if cumulative > d1 {
            break; // Eq. 11: whole-path score exceeds D_1
        }

        path.push((j_a, j_b));
        cumulative_scores.push(cumulative);

        // Snapshot the path at this length for later evaluation
        snapshots.push((path.clone(), cumulative));
    }

    snapshots
}

/// Extend AFP paths through the similarity matrix.
///
/// Starting from each seed AFP with score < d0, greedily extends the path
/// by appending the best-scoring neighboring AFP within the gap limit.
/// The cumulative D-score (average inter-fragment distance deviation)
/// must remain below d1 for the path to continue growing.
///
/// Returns the set of best candidate paths found, each represented as
/// a list of (i_a, i_b) AFP anchor positions.
fn extend_paths(
    sim: &[Vec<f32>],
    dm_a: &[Vec<f32>],
    dm_b: &[Vec<f32>],
    len_a: usize,
    len_b: usize,
    params: &CeParams,
) -> Vec<Vec<(usize, usize)>> {
    let w = params.win_size;
    let d0 = params.d0 as f64;
    let max_kept = params.max_paths.min(MAX_KEPT);
    let smaller = len_a.min(len_b);
    let win_sum = ((w - 1) * (w - 2)) / 2;

    // Precompute cumulative normalization factors for D-score averaging
    let norm: Vec<f64> = (0..smaller)
        .map(|i| ((i + 1) * i * w / 2 + (i + 1) * win_sum) as f64)
        .collect();

    let mut best_len = 0usize;
    let mut candidates: Vec<(Vec<(usize, usize)>, f64)> = Vec::new();

    for i_a in 0..len_a {
        // Early termination: can't form a longer path than the best found
        if best_len > 1 && i_a + w * (best_len - 1) > len_a {
            break;
        }

        for i_b in 0..len_b {
            let seed_score = sim[i_a][i_b] as f64;
            if seed_score >= d0 || seed_score < 0.0 {
                continue;
            }

            if best_len > 1 && i_b + w * (best_len - 1) > len_b {
                break;
            }

            let snapshots = grow_path(
                sim, dm_a, dm_b, len_a, len_b, (i_a, i_b), seed_score, params, &norm,
            );
            for (path, score) in snapshots {
                best_len = best_len.max(path.len());
                candidates.push((path, score));
            }

            // Periodically trim to bound memory
            if candidates.len() >= max_kept * 4 {
                sort_and_trim(&mut candidates, max_kept);
            }
        }
    }

    sort_and_trim(&mut candidates, max_kept);
    candidates.into_iter().map(|(path, _)| path).collect()
}

/// Sort candidates by quality (longer paths first, lower scores first)
/// and truncate to `max`.
fn sort_and_trim(candidates: &mut Vec<(Vec<(usize, usize)>, f64)>, max: usize) {
    candidates.sort_by(|a, b| {
        b.0.len()
            .cmp(&a.0.len())
            .then_with(|| a.1.partial_cmp(&b.1).unwrap_or(Ordering::Equal))
    });
    candidates.truncate(max);
}

// ============================================================================
// Path selection by Kabsch RMSD
// ============================================================================

/// Expand an AFP anchor path into residue-level pairs and evaluate by
/// Kabsch RMSD.
fn expand_and_evaluate(
    source_ca: &[Vec3],
    target_ca: &[Vec3],
    path: &[(usize, usize)],
    win_size: usize,
) -> Option<(Vec<(usize, usize)>, f32)> {
    let pairs: Vec<(usize, usize)> = path
        .iter()
        .flat_map(|&(first, second)| {
            (0..win_size).filter_map(move |k| {
                let (si, ti) = (first + k, second + k);
                (si < source_ca.len() && ti < target_ca.len()).then_some((si, ti))
            })
        })
        .collect();

    if pairs.len() < 3 {
        return None;
    }

    let src: Vec<Vec3> = pairs.iter().map(|&(si, _)| source_ca[si]).collect();
    let tgt: Vec<Vec3> = pairs.iter().map(|&(_, ti)| target_ca[ti]).collect();

    kabsch::kabsch(&src, &tgt, None)
        .ok()
        .map(|r| (pairs, r.rmsd))
}

/// Select the best candidate path by computing Kabsch RMSD for each.
///
/// Expands AFP anchor paths into full residue-level correspondences,
/// then evaluates each via optimal rigid-body superposition (Kabsch).
/// Prefers more aligned residues first, then lower RMSD as tiebreaker.
fn select_best_path(
    source_ca: &[Vec3],
    target_ca: &[Vec3],
    paths: &[Vec<(usize, usize)>],
    win_size: usize,
) -> Option<(Vec<(usize, usize)>, f32)> {
    paths
        .iter()
        .filter(|path| !path.is_empty())
        .filter_map(|path| expand_and_evaluate(source_ca, target_ca, path, win_size))
        .max_by(|(pairs_a, rmsd_a), (pairs_b, rmsd_b)| {
            pairs_a
                .len()
                .cmp(&pairs_b.len())
                .then_with(|| rmsd_b.partial_cmp(rmsd_a).unwrap_or(Ordering::Equal))
        })
}

// ============================================================================
// Z-score (Zhang & Skolnick, 2004)
// ============================================================================

/// Compute alignment quality Z-score.
///
/// Uses the length-dependent d0 normalization from Zhang & Skolnick (2004):
/// `d0 = 1.24 * (L_min - 15)^(1/3) - 1.8`
fn z_score(n_aligned: usize, rmsd: f32, len_source: usize, len_target: usize) -> f32 {
    let l_min = len_source.min(len_target) as f32;
    if l_min < 1.0 || rmsd < 1e-6 {
        return n_aligned as f32;
    }
    let coverage = n_aligned as f32 / l_min;
    let n = n_aligned as f32;
    let d0 = 1.24 * (l_min - 15.0).max(1.0).cbrt() - 1.8;
    let d0 = d0.max(0.5);
    n * coverage / (1.0 + rmsd / d0)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn v(x: f32, y: f32, z: f32) -> Vec3 {
        Vec3::new(x, y, z)
    }

    /// Generate an alpha-helix-like Cα trace with realistic ~3.8Å CA-CA spacing.
    fn make_helix(n: usize, offset: Vec3) -> Vec<Vec3> {
        (0..n)
            .map(|i| {
                let t = i as f32;
                v(
                    2.3 * (t * 1.745).cos() + offset.x,
                    2.3 * (t * 1.745).sin() + offset.y,
                    1.5 * t + offset.z,
                )
            })
            .collect()
    }

    /// Generate a beta-strand-like Cα trace (extended, 3.3Å spacing)
    fn make_strand(n: usize, offset: Vec3) -> Vec<Vec3> {
        (0..n)
            .map(|i| {
                let t = i as f32;
                v(
                    3.3 * t + offset.x,
                    offset.y + if i % 2 == 0 { 0.5 } else { -0.5 },
                    offset.z,
                )
            })
            .collect()
    }

    #[test]
    fn test_distance_matrix_symmetric() {
        let coords = make_helix(5, v(0.0, 0.0, 0.0));
        let dm = distance_matrix(&coords);
        assert_eq!(dm.len(), 5);
        for (i, row) in dm.iter().enumerate().take(5) {
            assert_eq!(row[i], 0.0);
            for j in 0..5 {
                assert!((row[j] - dm[j][i]).abs() < 1e-6);
            }
        }
    }

    #[test]
    fn test_afp_similarity_dimensions() {
        let a = make_helix(15, v(0.0, 0.0, 0.0));
        let b = make_helix(20, v(5.0, 5.0, 5.0));
        let dm_a = distance_matrix(&a);
        let dm_b = distance_matrix(&b);
        let s = afp_similarity_matrix(&dm_a, &dm_b, a.len(), b.len(), 8);
        assert_eq!(s.len(), 15);
        assert_eq!(s[0].len(), 20);
        assert!(s[0][0] >= 0.0);
        assert_eq!(s[14][0], -1.0);
    }

    #[test]
    fn test_afp_similarity_identical_fragments() {
        let coords = make_helix(20, v(0.0, 0.0, 0.0));
        let dm = distance_matrix(&coords);
        let s = afp_similarity_matrix(&dm, &dm, coords.len(), coords.len(), 8);
        assert!(s[0][0].abs() < 1e-6);
        assert!(s[5][5].abs() < 1e-6);
    }

    #[test]
    fn test_identical_structures() {
        let coords = make_helix(20, v(0.0, 0.0, 0.0));
        let result = ce_align(&coords, &coords, &CeParams::default()).unwrap();
        assert!(result.n_aligned >= 15, "Expected >=15 aligned, got {}", result.n_aligned);
        assert!(result.rmsd < 0.1, "Expected low RMSD for identical, got {}", result.rmsd);
    }

    #[test]
    fn test_translated_structure() {
        let source = make_helix(25, v(0.0, 0.0, 0.0));
        let target = make_helix(25, v(10.0, 20.0, 30.0));
        let result = ce_align(&source, &target, &CeParams::default()).unwrap();
        assert!(result.n_aligned >= 18, "Expected >=18 aligned, got {}", result.n_aligned);
        assert!(result.rmsd < 0.5, "Expected low RMSD for translated, got {}", result.rmsd);
    }

    #[test]
    fn test_different_lengths() {
        let source = make_helix(25, v(0.0, 0.0, 0.0));
        let target = make_helix(40, v(5.0, 5.0, 5.0));
        let result = ce_align(&source, &target, &CeParams::default()).unwrap();
        assert!(result.n_aligned >= 8, "Expected >=8 aligned, got {}", result.n_aligned);
    }

    #[test]
    fn test_dissimilar_structures() {
        let helix = make_helix(20, v(0.0, 0.0, 0.0));
        let strand = make_strand(20, v(0.0, 0.0, 0.0));
        let result = ce_align(&helix, &strand, &CeParams::default());
        if let Ok(r) = result {
            assert!(r.z_score < 5.0, "Z-score should be low for dissimilar, got {}", r.z_score);
        }
    }

    #[test]
    fn test_too_short() {
        let short = make_helix(5, v(0.0, 0.0, 0.0));
        let result = ce_align(&short, &short, &CeParams::default());
        assert!(result.is_err());
    }

    #[test]
    fn test_custom_params() {
        let source = make_helix(20, v(0.0, 0.0, 0.0));
        let target = make_helix(20, v(5.0, 0.0, 0.0));
        let params = CeParams {
            win_size: 6,
            gap_max: 10,
            d0: 2.5,
            d1: 3.5,
            max_paths: 10,
        };
        let result = ce_align(&source, &target, &params).unwrap();
        assert!(result.n_aligned >= 8, "Expected >=8 aligned, got {}", result.n_aligned);
    }

    #[test]
    fn test_z_score_positive() {
        let z = z_score(50, 2.0, 100, 120);
        assert!(z > 0.0, "Z-score should be positive, got {}", z);
    }
}
