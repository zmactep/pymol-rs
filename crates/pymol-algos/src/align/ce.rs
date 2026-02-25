//! Combinatorial Extension (CE) structural alignment
//!
//! Implements the CE algorithm (Shindyalov & Bourne, 1998) for structure-based
//! protein alignment. Works on Cα coordinates, finding residue correspondences
//! purely from 3D structure by comparing local distance-matrix fingerprints
//! and chaining compatible fragment matches.
//!
//! This is a faithful port of PyMOL's `ccealignmodule.cpp`.

use super::kabsch;
use crate::AlignError;

const MAX_KEPT: usize = 20;
const SENTINEL_SCORE: f64 = 1e6;

/// Parameters for the CE structural alignment algorithm
#[derive(Debug, Clone)]
pub struct CeParams {
    /// Fragment window size (number of residues per fragment)
    pub win_size: usize,
    /// Maximum gap allowed between consecutive AFPs in the path
    pub gap_max: usize,
    /// Distance cutoff for fragment similarity scoring (Angstroms)
    pub d0: f32,
    /// Distance cutoff for path extension (Angstroms)
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
    source_ca: &[[f32; 3]],
    target_ca: &[[f32; 3]],
    params: &CeParams,
) -> Result<CeResult, AlignError> {
    let len_a = source_ca.len();
    let len_b = target_ca.len();
    let w = params.win_size;

    if len_a < w || len_b < w {
        return Err(AlignError::TooFewAtoms(len_a.min(len_b)));
    }

    // Step 1: Compute full NxN distance matrices
    let dm_a = calc_dm(source_ca);
    let dm_b = calc_dm(target_ca);

    // Step 2: Compute CE similarity matrix
    let s = calc_s(&dm_a, &dm_b, len_a, len_b, w);

    // Step 3: Find best AFP paths
    let paths = find_path(&s, &dm_a, &dm_b, len_a, len_b, params);

    if paths.is_empty() {
        return Err(AlignError::NoMatches);
    }

    // Step 4: Select path with best Kabsch RMSD
    let (pairs, rmsd) = find_best(source_ca, target_ca, &paths, w)
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

/// Compute full NxN pairwise distance matrix.
/// Equivalent to PyMOL's `calcDM`.
fn calc_dm(coords: &[[f32; 3]]) -> Vec<Vec<f32>> {
    let n = coords.len();
    let mut dm = vec![vec![0.0f32; n]; n];
    for i in 0..n {
        for j in (i + 1)..n {
            let dx = coords[i][0] - coords[j][0];
            let dy = coords[i][1] - coords[j][1];
            let dz = coords[i][2] - coords[j][2];
            let d = (dx * dx + dy * dy + dz * dz).sqrt();
            dm[i][j] = d;
            dm[j][i] = d;
        }
    }
    dm
}

/// Compute CE similarity matrix.
/// Equivalent to PyMOL's `calcS`.
///
/// For each fragment pair (iA, iB), compares intra-fragment distance patterns.
/// Skips adjacent residue distances (constant ~3.8Å CA-CA bond).
/// S[iA][iB] = mean |d_A[iA+row][iA+col] - d_B[iB+row][iB+col]|
/// for row=0..w-3, col=row+2..w-1.
/// S[iA][iB] = -1.0 if the fragment extends past the protein end.
fn calc_s(
    dm_a: &[Vec<f32>],
    dm_b: &[Vec<f32>],
    len_a: usize,
    len_b: usize,
    win_size: usize,
) -> Vec<Vec<f32>> {
    let w = win_size;
    // Number of distance pairs compared (skipping adjacent)
    let sum_size = ((w - 1) * (w - 2)) as f32 / 2.0;

    let mut s = vec![vec![-1.0f32; len_b]; len_a];

    for i_a in 0..len_a {
        for i_b in 0..len_b {
            if i_a > len_a - w || i_b > len_b - w {
                continue; // S[i_a][i_b] stays -1.0
            }

            let mut score = 0.0f32;
            // Skip distance from residue to its immediate neighbor (row, row+1)
            // PyMOL: for (row = 0; row < wSize - 2; row++)
            //          for (col = row + 2; col < wSize; col++)
            for row in 0..(w - 2) {
                for col in (row + 2)..w {
                    score += (dm_a[i_a + row][i_a + col] - dm_b[i_b + row][i_b + col]).abs();
                }
            }

            s[i_a][i_b] = score / sum_size;
        }
    }

    s
}

/// Find optimal AFP paths through the similarity matrix.
/// Equivalent to PyMOL's `findPath`.
///
/// Returns up to `max_paths` best paths, each as a list of (first, second) AFP pairs.
fn find_path(
    s: &[Vec<f32>],
    dm_a: &[Vec<f32>],
    dm_b: &[Vec<f32>],
    len_a: usize,
    len_b: usize,
    params: &CeParams,
) -> Vec<Vec<(usize, usize)>> {
    let win_size = params.win_size;
    let gap_max = params.gap_max;
    let d0 = params.d0 as f64;
    let d1 = params.d1 as f64;
    let max_kept = params.max_paths.min(MAX_KEPT);

    let smaller = len_a.min(len_b);
    let win_sum = ((win_size - 1) * (win_size - 2)) / 2;

    // Best path seen so far (across all seeds)
    let mut best_path: Vec<(usize, usize)> = Vec::new();
    let mut best_path_score: f64 = SENTINEL_SCORE;
    let mut best_path_length: usize = 0;

    // Ring buffer of top paths
    let mut path_buffer: Vec<Option<Vec<(usize, usize)>>> = vec![None; max_kept];
    let mut score_buffer = vec![SENTINEL_SCORE; max_kept];
    let mut len_buffer = vec![0usize; max_kept];
    let mut buffer_index: usize = 0;
    let mut buffer_size: usize = 0;

    // winCache: precomputed cumulative weight denominators
    let win_cache: Vec<i64> = (0..smaller)
        .map(|i| ((i + 1) * i * win_size / 2 + (i + 1) * win_sum) as i64)
        .collect();

    // allScoreBuffer: partial gapped scores per path position per gap index
    let mut all_score_buffer = vec![vec![SENTINEL_SCORE; gap_max * 2 + 1]; smaller];

    // tIndex: tracks which gap index was chosen at each path position
    let mut t_index = vec![0usize; smaller];

    for i_a in 0..len_a {
        // Pruning: can't possibly build a longer path than best_path_length from here
        if best_path_length > 1 && i_a > len_a - win_size * (best_path_length - 1) {
            break;
        }

        for i_b in 0..len_b {
            let s_val = s[i_a][i_b] as f64;
            if s_val >= d0 || s_val == -1.0 {
                continue;
            }

            if best_path_length > 1 && i_b > len_b - win_size * (best_path_length - 1) {
                break;
            }

            // Start a new path from (iA, iB)
            let mut cur_path = vec![(0usize, 0usize); smaller];
            cur_path[0] = (i_a, i_b);
            let mut cur_path_length: usize = 1;
            t_index[0] = 0;

            // Reset allScoreBuffer for this path
            for row in all_score_buffer.iter_mut().take(smaller) {
                for val in row.iter_mut().take(gap_max * 2 + 1) {
                    *val = SENTINEL_SCORE;
                }
            }

            let mut done = false;
            while !done {
                let mut gap_best_score: f64 = SENTINEL_SCORE;
                let mut gap_best_index: Option<usize> = None;

                // Try all gaps (g used for both indexing and arithmetic)
                #[allow(clippy::needless_range_loop)]
                for g in 0..(gap_max * 2 + 1) {
                    let last = cur_path[cur_path_length - 1];
                    let mut j_a = last.0 + win_size;
                    let mut j_b = last.1 + win_size;

                    // Alternating gap offset: even g+1 gaps A, odd g+1 gaps B
                    #[allow(clippy::manual_div_ceil)]
                    let gap_offset = (g + 1) / 2;
                    if (g + 1) % 2 == 0 {
                        j_a += gap_offset;
                    } else {
                        j_b += gap_offset;
                    }

                    // Boundary checks
                    if j_a > len_a - win_size || j_b > len_b - win_size {
                        continue;
                    }

                    let s_jab = s[j_a][j_b] as f64;
                    if s_jab >= d0 || s_jab == -1.0 {
                        continue;
                    }

                    // Score extension using cross-fragment distances
                    let mut cur_score: f64 = 0.0;
                    for &prev in cur_path.iter().take(cur_path_length) {
                        // First atom of prev fragment vs first atom of new fragment
                        cur_score += (dm_a[prev.0][j_a] as f64
                            - dm_b[prev.1][j_b] as f64)
                            .abs();
                        // Last atom of prev fragment vs last atom of new fragment
                        cur_score += (dm_a[prev.0 + win_size - 1][j_a + win_size - 1] as f64
                            - dm_b[prev.1 + win_size - 1][j_b + win_size - 1] as f64)
                            .abs();
                        // Anti-diagonal: k=1..winSize-2
                        for k in 1..(win_size - 1) {
                            cur_score += (dm_a[prev.0 + k][j_a + win_size - 1 - k] as f64
                                - dm_b[prev.1 + k][j_b + win_size - 1 - k] as f64)
                                .abs();
                        }
                    }

                    cur_score /= (win_size * cur_path_length) as f64;

                    if cur_score >= d1 {
                        continue;
                    }

                    // Store gapped best
                    if cur_score < gap_best_score {
                        cur_path[cur_path_length] = (j_a, j_b);
                        gap_best_score = cur_score;
                        gap_best_index = Some(g);
                        all_score_buffer[cur_path_length - 1][g] = cur_score;
                    }
                }

                // Calculate cumulative total score
                if let Some(gbi) = gap_best_index {
                    #[allow(clippy::manual_div_ceil)]
                    let j_gap = (gbi + 1) / 2;
                    let (g_a, g_b) = if (gbi + 1) % 2 == 0 {
                        let prev = cur_path[cur_path_length - 1];
                        (prev.0 + win_size + j_gap, prev.1 + win_size)
                    } else {
                        let prev = cur_path[cur_path_length - 1];
                        (prev.0 + win_size, prev.1 + win_size + j_gap)
                    };

                    // score1: weighted combination of cross-fragment score with local S score
                    let score1 = (all_score_buffer[cur_path_length - 1][gbi]
                        * (win_size * cur_path_length) as f64
                        + s[g_a][g_b] as f64 * win_sum as f64)
                        / (win_size * cur_path_length + win_sum) as f64;

                    // score2: cumulative weighted average with previous path score
                    let prev_score = if cur_path_length > 1 {
                        all_score_buffer[cur_path_length - 2][t_index[cur_path_length - 1]]
                    } else {
                        s[i_a][i_b] as f64
                    };

                    let score2 = (prev_score * win_cache[cur_path_length - 1] as f64
                        + score1
                            * (win_cache[cur_path_length] - win_cache[cur_path_length - 1]) as f64)
                        / win_cache[cur_path_length] as f64;

                    let cur_total_score = score2;

                    // Heuristic: path getting sloppy, stop
                    if cur_total_score > d1 {
                        done = true;
                        // gap_best_index effectively -1
                    } else {
                        all_score_buffer[cur_path_length - 1][gbi] = cur_total_score;
                        t_index[cur_path_length] = gbi;
                        cur_path_length += 1;
                    }

                    // Update best path if this one is better
                    if !done
                        && (cur_path_length > best_path_length
                            || (cur_path_length == best_path_length
                                && cur_total_score < best_path_score))
                    {
                        best_path_length = cur_path_length;
                        best_path_score = cur_total_score;
                        best_path = cur_path[..cur_path_length].to_vec();
                    }
                } else {
                    // No good extension found
                    done = true;
                    cur_path_length = cur_path_length.saturating_sub(1);
                }
            }

            // Add to ring buffer if best path is good enough
            if best_path_length > len_buffer[buffer_index]
                || (best_path_length == len_buffer[buffer_index]
                    && best_path_score < score_buffer[buffer_index])
            {
                buffer_index = if buffer_index == max_kept - 1 {
                    0
                } else {
                    buffer_index + 1
                };
                buffer_size = if buffer_size < max_kept {
                    buffer_size + 1
                } else {
                    max_kept
                };

                let path_copy = best_path.clone();

                let idx = if buffer_index == 0 && buffer_size == max_kept {
                    max_kept - 1
                } else {
                    buffer_index - 1
                };

                path_buffer[idx] = Some(path_copy);
                score_buffer[idx] = best_path_score;
                len_buffer[idx] = best_path_length;
            }
        }
    }

    // Collect non-empty paths from buffer
    path_buffer.into_iter().flatten().collect()
}

/// Select the best path by computing Kabsch RMSD for each.
/// Equivalent to PyMOL's `findBest`.
///
/// Expands each AFP path into individual residue pairs, computes RMSD via Kabsch,
/// returns the path with lowest RMSD.
fn find_best(
    source_ca: &[[f32; 3]],
    target_ca: &[[f32; 3]],
    paths: &[Vec<(usize, usize)>],
    win_size: usize,
) -> Option<(Vec<(usize, usize)>, f32)> {
    let mut best_rmsd = f32::MAX;
    let mut best_pairs: Option<Vec<(usize, usize)>> = None;

    for path in paths {
        if path.is_empty() {
            continue;
        }

        // Expand AFPs into individual residue pairs
        let mut src_coords = Vec::new();
        let mut tgt_coords = Vec::new();
        let mut pairs = Vec::new();

        for &(first, second) in path {
            for k in 0..win_size {
                let si = first + k;
                let ti = second + k;
                if si < source_ca.len() && ti < target_ca.len() {
                    src_coords.push(source_ca[si]);
                    tgt_coords.push(target_ca[ti]);
                    pairs.push((si, ti));
                }
            }
        }

        if src_coords.len() < 3 {
            continue;
        }

        // Compute Kabsch RMSD
        match kabsch::kabsch(&src_coords, &tgt_coords, None) {
            Ok(result) => {
                if result.rmsd < best_rmsd
                    || (result.rmsd == best_rmsd
                        && pairs.len() > best_pairs.as_ref().map_or(0, |p| p.len()))
                {
                    best_rmsd = result.rmsd;
                    best_pairs = Some(pairs);
                }
            }
            Err(_) => continue,
        }
    }

    best_pairs.map(|p| (p, best_rmsd))
}

/// Compute Z-score for alignment quality assessment.
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

    /// Generate an alpha-helix-like Cα trace with realistic ~3.8Å CA-CA spacing.
    /// Parameters match a real alpha helix: radius ~2.3Å, rise ~1.5Å/residue,
    /// turn ~100°/residue.
    fn make_helix(n: usize, offset: [f32; 3]) -> Vec<[f32; 3]> {
        (0..n)
            .map(|i| {
                let t = i as f32;
                [
                    2.3 * (t * 1.745).cos() + offset[0],
                    2.3 * (t * 1.745).sin() + offset[1],
                    1.5 * t + offset[2],
                ]
            })
            .collect()
    }

    /// Generate a beta-strand-like Cα trace (extended, 3.3Å spacing)
    fn make_strand(n: usize, offset: [f32; 3]) -> Vec<[f32; 3]> {
        (0..n)
            .map(|i| {
                let t = i as f32;
                [
                    3.3 * t + offset[0],
                    offset[1] + if i % 2 == 0 { 0.5 } else { -0.5 },
                    offset[2],
                ]
            })
            .collect()
    }

    #[test]
    fn test_calc_dm_symmetric() {
        let coords = make_helix(5, [0.0, 0.0, 0.0]);
        let dm = calc_dm(&coords);
        assert_eq!(dm.len(), 5);
        for i in 0..5 {
            assert_eq!(dm[i][i], 0.0);
            for j in 0..5 {
                assert!((dm[i][j] - dm[j][i]).abs() < 1e-6);
            }
        }
    }

    #[test]
    fn test_calc_s_dimensions() {
        let a = make_helix(15, [0.0, 0.0, 0.0]);
        let b = make_helix(20, [5.0, 5.0, 5.0]);
        let dm_a = calc_dm(&a);
        let dm_b = calc_dm(&b);
        let s = calc_s(&dm_a, &dm_b, a.len(), b.len(), 8);
        assert_eq!(s.len(), 15);
        assert_eq!(s[0].len(), 20);
        // Valid positions should not be -1.0
        assert!(s[0][0] >= 0.0);
        // Past-end positions should be -1.0
        assert_eq!(s[14][0], -1.0);
    }

    #[test]
    fn test_calc_s_identical_fragments() {
        let coords = make_helix(20, [0.0, 0.0, 0.0]);
        let dm = calc_dm(&coords);
        let s = calc_s(&dm, &dm, coords.len(), coords.len(), 8);
        // Same fragment at same position should score 0
        assert!(s[0][0].abs() < 1e-6);
        assert!(s[5][5].abs() < 1e-6);
    }

    #[test]
    fn test_identical_structures() {
        let coords = make_helix(20, [0.0, 0.0, 0.0]);
        let result = ce_align(&coords, &coords, &CeParams::default()).unwrap();
        assert!(
            result.n_aligned >= 15,
            "Expected >=15 aligned, got {}",
            result.n_aligned
        );
        assert!(
            result.rmsd < 0.1,
            "Expected low RMSD for identical, got {}",
            result.rmsd
        );
    }

    #[test]
    fn test_translated_structure() {
        let source = make_helix(25, [0.0, 0.0, 0.0]);
        let target = make_helix(25, [10.0, 20.0, 30.0]);
        let result = ce_align(&source, &target, &CeParams::default()).unwrap();
        // Distance matrices are translation-invariant
        assert!(
            result.n_aligned >= 18,
            "Expected >=18 aligned, got {}",
            result.n_aligned
        );
        assert!(
            result.rmsd < 0.5,
            "Expected low RMSD for translated, got {}",
            result.rmsd
        );
    }

    #[test]
    fn test_different_lengths() {
        // Use longer helices so there are enough fragment positions for CE
        let source = make_helix(25, [0.0, 0.0, 0.0]);
        let target = make_helix(40, [5.0, 5.0, 5.0]);
        let result = ce_align(&source, &target, &CeParams::default()).unwrap();
        assert!(
            result.n_aligned >= 8,
            "Expected >=8 aligned, got {}",
            result.n_aligned
        );
    }

    #[test]
    fn test_dissimilar_structures() {
        let helix = make_helix(20, [0.0, 0.0, 0.0]);
        let strand = make_strand(20, [0.0, 0.0, 0.0]);
        let result = ce_align(&helix, &strand, &CeParams::default());
        match result {
            Ok(r) => assert!(
                r.z_score < 5.0,
                "Z-score should be low for dissimilar, got {}",
                r.z_score
            ),
            Err(_) => {} // Also acceptable
        }
    }

    #[test]
    fn test_too_short() {
        let short = make_helix(5, [0.0, 0.0, 0.0]);
        let result = ce_align(&short, &short, &CeParams::default());
        assert!(result.is_err());
    }

    #[test]
    fn test_custom_params() {
        let source = make_helix(20, [0.0, 0.0, 0.0]);
        let target = make_helix(20, [5.0, 0.0, 0.0]);
        let params = CeParams {
            win_size: 6,
            gap_max: 10,
            d0: 2.5,
            d1: 3.5,
            max_paths: 10,
        };
        let result = ce_align(&source, &target, &params).unwrap();
        assert!(
            result.n_aligned >= 8,
            "Expected >=8 aligned, got {}",
            result.n_aligned
        );
    }

    #[test]
    fn test_z_score_positive() {
        let z = z_score(50, 2.0, 100, 120);
        assert!(z > 0.0, "Z-score should be positive, got {}", z);
    }
}
