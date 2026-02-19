//! PyMOL-compatible cartoon pipeline
//!
//! Implements the exact sequence from RepCartoon.cpp:
//! 1. Compute differences, normals, tangents from guide points
//! 2. Compute round helix orientations (using tangent vectors)
//! 3. Refine normals (orthogonalize, fix flips)
//! 4. Flatten sheets, smooth loops
//! 5. Re-compute tangents after position modifications
//! 6. Per-run interpolation (CartoonGenerateSample) and mesh generation

use lin_alg::f32::Vec3;
use pymol_mol::SecondaryStructure;

use super::backbone::{BackboneSegment, GuidePoint};
use super::frame::{FrameWithMetadata, ReferenceFrame};
use super::geometry::{
    connect_rings, find_helix_regions, find_sheet_termini,
    generate_explicit_sheet, generate_face_strips,
    CartoonGeometrySettings, Profile,
};
use super::utils::{is_helix, normalize_safe, smooth};
use crate::vertex::MeshVertex;

// ============================================================================
// Pipeline Settings
// ============================================================================

/// Settings for the cartoon pipeline
#[derive(Debug, Clone)]
pub struct PipelineSettings {
    /// Interpolation samples per residue pair (cartoon_sampling, default 7)
    pub sampling: u32,
    /// Power for sigmoid blend (cartoon_power, default 2.0)
    pub power_a: f32,
    /// Power for displacement envelope (cartoon_power_b, default 0.52)
    pub power_b: f32,
    /// Throw factor for displacement magnitude (cartoon_throw, default 1.35)
    pub throw_factor: f32,
    /// Number of sheet flattening cycles (cartoon_flat_cycles, default 4)
    pub flat_cycles: u32,
    /// First window size for loop smoothing (cartoon_smooth_first, default 1)
    pub smooth_first: u32,
    /// Last window size for loop smoothing (cartoon_smooth_last, default 1)
    pub smooth_last: u32,
    /// Number of smoothing cycles (cartoon_smooth_cycles, default 2)
    pub smooth_cycles: u32,
    /// Whether to refine normals (cartoon_refine_normals, default true)
    pub refine_normals: bool,
    /// Whether to compute round helix orientations (cartoon_round_helices, default true)
    pub round_helices: bool,
    /// Position refinement cycles after each interpolation (cartoon_refine, default 5)
    pub refine: u32,
    /// Whether to smooth loop regions (cartoon_smooth_loops, default false)
    pub smooth_loops: bool,
}

impl Default for PipelineSettings {
    fn default() -> Self {
        Self {
            sampling: 7,
            power_a: 2.0,
            power_b: 0.52,
            throw_factor: 1.35,
            flat_cycles: 4,
            smooth_first: 1,
            smooth_last: 1,
            smooth_cycles: 2,
            refine_normals: true,
            round_helices: true,
            refine: 5,
            smooth_loops: false,
        }
    }
}

// ============================================================================
// CartoonType — maps SS type + settings to rendering style
// ============================================================================

/// Cartoon rendering type (matches PyMOL's cur_car values)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CartoonType {
    Loop,
    Oval,
    Rect,
    Arrow,
    Dumbbell,
    /// Nucleic acid flat ribbon (same geometry as Rect but no arrow terminus)
    NucleicRect,
}

/// Map from SS type + settings to CartoonType
fn cartoon_type_for(ss: SecondaryStructure, settings: &CartoonGeometrySettings) -> CartoonType {
    if settings.uniform_tube {
        return CartoonType::Loop;
    }
    if is_helix(ss) {
        if settings.fancy_helices {
            CartoonType::Dumbbell
        } else {
            CartoonType::Oval
        }
    } else if ss == SecondaryStructure::Sheet {
        if settings.fancy_sheets {
            CartoonType::Arrow
        } else {
            CartoonType::Rect
        }
    } else if ss == SecondaryStructure::NucleicRibbon {
        // Nucleic ribbon: flat rectangle without arrow terminus
        CartoonType::NucleicRect
    } else {
        CartoonType::Loop
    }
}

// ============================================================================
// Phase 2: Compute Differences, Normals, Tangents
// ============================================================================

/// Compute difference vectors and their lengths between consecutive guide points.
///
/// Returns (distances, normalized_directions) — both Vec<f32>/Vec<Vec3> of length n-1.
/// Matches RepCartoonComputeDifferencesAndNormals.
fn compute_differences_and_normals(gps: &[GuidePoint]) -> (Vec<f32>, Vec<Vec3>) {
    let n = gps.len();
    if n < 2 {
        return (Vec::new(), Vec::new());
    }

    let mut dl = Vec::with_capacity(n - 1);
    let mut nv = Vec::with_capacity(n - 1);

    for a in 0..n - 1 {
        let diff = gps[a + 1].position - gps[a].position;
        let d = diff.magnitude();
        dl.push(d);
        if d > 1e-4 {
            nv.push(diff / d);
        } else if a > 0 {
            // Copy previous if zero length
            nv.push(*nv.last().unwrap_or(&Vec3::new(0.0, 0.0, 0.0)));
        } else {
            nv.push(Vec3::new(0.0, 0.0, 0.0));
        }
    }

    (dl, nv)
}

/// Compute tangent vectors from normalized direction vectors.
///
/// Matches RepCartoonComputeTangents exactly:
/// - First: tv[0] = nv[0]
/// - Interior: tv[i] = normalize(nv[i] + nv[i-1])
/// - Last: tv[n-1] = nv[n-2]
fn compute_tangents(nv: &[Vec3], n: usize) -> Vec<Vec3> {
    if n < 2 || nv.is_empty() {
        return vec![Vec3::new(0.0, 0.0, 1.0); n];
    }

    let mut tv = Vec::with_capacity(n);

    // First point
    tv.push(nv[0]);

    // Interior points
    for a in 1..n - 1 {
        if a < nv.len() && a >= 1 {
            let sum = nv[a] + nv[a - 1];
            tv.push(normalize_safe(sum));
        } else if a >= 1 && a - 1 < nv.len() {
            tv.push(nv[a - 1]);
        } else if a < nv.len() {
            tv.push(nv[a]);
        } else {
            tv.push(Vec3::new(0.0, 0.0, 1.0));
        }
    }

    // Last point
    tv.push(nv[nv.len() - 1]);

    tv
}

// ============================================================================
// Phase 3: Round Helix Orientations
// ============================================================================

/// Compute round helix orientations using sliding window of 4-5 CA positions.
///
/// Matches RepCartoonComputeRoundHelices exactly:
/// - Maintains sliding window v1..v5 of helix CA positions
/// - Computes weighted center from 4 CAs: 0.2130*(v1+v4) + 0.2870*(v2+v3)
/// - Helix axis = normalize(prev_center - center)
/// - Orientation = normalize(cross(axis, tangent))
#[allow(unused_assignments)]
fn compute_round_helix_orientations(gps: &mut [GuidePoint], tv: &[Vec3]) {
    let n = gps.len();
    if n < 2 {
        return;
    }

    // Sliding window of helix CA positions (indices into gps)
    let mut v1: Option<usize> = None;
    let mut v2: Option<usize> = None;
    let mut v3: Option<usize> = None;
    let mut v4: Option<usize> = None;
    let mut v5: Option<usize> = None;
    let mut last = 0i32;
    let mut prev_center = Vec3::new(0.0, 0.0, 0.0);

    for a in 0..n {
        // Shift window
        v5 = v4;
        v4 = v3;
        v3 = v2;
        v2 = v1;

        if is_helix(gps[a].ss_type) {
            v1 = Some(a);
        } else {
            // Early termination: helix ended with < 2 centers computed
            if last < 2 {
                // Compute axis from chain direction vectors
                if let (Some(i2), Some(i3)) = (v2, v3) {
                    let mut t0 = gps[i2].position - gps[a].position;
                    t0 = normalize_safe(t0);
                    let mut t1 = gps[i3].position - gps[i2].position;
                    t1 = normalize_safe(t1);
                    t0 = t0 + t1;

                    if let Some(i4) = v4 {
                        t1 = gps[i4].position - gps[i3].position;
                        t1 = normalize_safe(t1);
                        t0 = t0 + t1;
                    }
                    if let Some(i5) = v5 {
                        if let Some(i4) = v4 {
                            t1 = gps[i5].position - gps[i4].position;
                            t1 = normalize_safe(t1);
                            t0 = t0 + t1;
                        }
                    }
                    t0 = normalize_safe(t0);

                    // Set orientations using cross(axis, tangent)
                    // a-1 corresponds to v2 (if exists)
                    if a >= 1 {
                        gps[a - 1].orientation = normalize_safe(t0.cross(tv[a - 1]));
                    }
                    if a >= 2 {
                        gps[a - 2].orientation = normalize_safe(t0.cross(tv[a - 2]));
                    }
                    if v4.is_some() && a >= 3 {
                        gps[a - 3].orientation = normalize_safe(t0.cross(tv[a - 3]));
                    }
                    if v5.is_some() && a >= 4 {
                        gps[a - 4].orientation = normalize_safe(t0.cross(tv[a - 4]));
                    }

                    // Check for goofy flip on short tight helices
                    if v4.is_some() && v5.is_some() && a >= 4 {
                        if gps[a - 3].orientation.dot(gps[a - 4].orientation) < -0.8 {
                            gps[a - 4].orientation = gps[a - 4].orientation * -1.0;
                        }
                    }
                }
            }

            v1 = None;
            v2 = None;
            v3 = None;
            v4 = None;
            v5 = None;
            last = 0;
        }

        // When we have 4 helix CAs, compute center and axis
        if let (Some(i1), Some(i2), Some(i3), Some(i4)) = (v1, v2, v3, v4) {
            let t0_raw = gps[i1].position + gps[i4].position;
            let t1_raw = gps[i2].position + gps[i3].position;
            let center = t0_raw * 0.2130 + t1_raw * 0.2870;

            if last > 0 {
                // We have a previous center — compute axis
                let axis = normalize_safe(prev_center - center);

                // Set orientations = cross(axis, tangent)
                gps[a].orientation = normalize_safe(axis.cross(tv[a]));
                if a >= 1 {
                    gps[a - 1].orientation = normalize_safe(axis.cross(tv[a - 1]));
                }
                if a >= 2 {
                    gps[a - 2].orientation = normalize_safe(axis.cross(tv[a - 2]));
                }

                if last == 1 {
                    // 5th CA — also fix first residues of helix
                    if a >= 3 {
                        gps[a - 3].orientation = normalize_safe(axis.cross(tv[a - 3]));
                    }
                    if a >= 4 {
                        gps[a - 4].orientation = normalize_safe(axis.cross(tv[a - 4]));
                    }
                }
            }

            last += 1;
            prev_center = center;
        }
    }
}

// ============================================================================
// Phase 4: Refine Normals
// ============================================================================

/// Refine orientation vectors to prevent flips and kinks.
///
/// Matches RepCartoonRefineNormals exactly:
/// 1. Make orientations orthogonal to tangents (interior residues only)
/// 2. Generate alternative inverted orientations (NOT for helices)
/// 3. Forward iterate: pick orientation with highest dot product with previous
/// 4. Soften kinks where dot(prev,curr) * dot(curr,next) < -0.10
fn refine_normals(gps: &mut [GuidePoint], tv: &[Vec3], nv: &[Vec3]) {
    let n = gps.len();
    if n < 3 {
        return;
    }

    // Step 1: Make orientations orthogonal to tangents (interior only)
    for a in 1..n - 1 {
        let orient = gps[a].orientation;
        let tangent = tv[a];
        let removed = orient - tangent * orient.dot(tangent);
        let mag = removed.magnitude();
        if mag > 1e-6 {
            gps[a].orientation = removed / mag;
        }
    }

    // Step 2: Generate alternative inverted orientations
    // For helices, keep same (don't allow inversion)
    let mut alternatives: Vec<[Vec3; 2]> = Vec::with_capacity(n);
    for a in 0..n {
        let orient = gps[a].orientation;
        let alt = if is_helix(gps[a].ss_type) {
            orient // No inversion for helices
        } else {
            orient * -1.0
        };
        alternatives.push([orient, alt]);
    }

    // Step 3: Forward iterate through pairs to select optimal orientation
    for a in 1..n - 1 {
        // nv is the chain direction vectors (normalized differences)
        // PyMOL removes chain component from both previous orientation and candidates
        let chain_dir = if a > 0 && a - 1 < nv.len() {
            nv[a - 1]
        } else {
            tv[a]
        };

        // Previous orientation, projected perpendicular to chain direction
        let prev_orient = gps[a - 1].orientation;
        let o0 = prev_orient - chain_dir * prev_orient.dot(chain_dir);
        let o0_len = o0.magnitude();
        let o0 = if o0_len > 1e-6 { o0 / o0_len } else { prev_orient };

        // Candidate 0 (original), projected
        let c0 = alternatives[a][0];
        let o1_0 = c0 - chain_dir * c0.dot(chain_dir);
        let o1_0_len = o1_0.magnitude();
        let o1_0 = if o1_0_len > 1e-6 { o1_0 / o1_0_len } else { c0 };

        // Candidate 1 (inverted), projected
        let c1 = alternatives[a][1];
        let o1_1 = c1 - chain_dir * c1.dot(chain_dir);
        let o1_1_len = o1_1.magnitude();
        let o1_1 = if o1_1_len > 1e-6 { o1_1 / o1_1_len } else { c1 };

        let dot0 = o0.dot(o1_0);
        let dot1 = o0.dot(o1_1);

        if dot1 > dot0 {
            gps[a].orientation = c1;
        } else {
            gps[a].orientation = c0;
        }
    }

    // Step 4: Soften kinks
    // Store results in alternatives buffer, then apply
    for a in 1..n - 1 {
        let curr = gps[a].orientation;
        let prev = gps[a - 1].orientation;
        let next = gps[a + 1].orientation;

        let dp = curr.dot(next) * curr.dot(prev);
        if dp < -0.10 {
            // Kink detected
            let tangent = tv[a];

            // Average neighbors + tiny bit of current
            let mut t0 = next + prev;
            let t1 = curr * 0.001;
            t0 = t0 + t1;

            // Remove tangent component
            t0 = t0 - tangent * t0.dot(tangent);
            let t0_len = t0.magnitude();

            if t0_len > 1e-6 {
                let t0_norm = t0 / t0_len;

                let t2 = if curr.dot(t0_norm) < 0.0 {
                    curr - t0_norm
                } else {
                    curr + t0_norm
                };
                let t2 = normalize_safe(t2);

                // Blend factor based on kink severity
                let blend = (2.0 * (-0.10 - dp)).clamp(0.0, 1.0);

                // mix3f: result = curr*(1-blend) + t2*blend
                let result = curr * (1.0 - blend) + t2 * blend;
                alternatives[a][0] = result;
            } else {
                alternatives[a][0] = curr;
            }
        } else {
            alternatives[a][0] = curr;
        }
    }

    // Apply kink-softened orientations
    for a in 1..n - 1 {
        gps[a].orientation = alternatives[a][0];
    }
}

// ============================================================================
// Phase 5: Flatten Sheets, Smooth Loops
// ============================================================================

/// Flatten sheet regions to make them planar.
///
/// Matches RepCartoonFlattenSheets: iterative window averaging of positions
/// and orientations, then removing tangent component from orientations.
fn flatten_sheets(gps: &mut [GuidePoint], flat_cycles: u32) {
    let n = gps.len();
    if n < 3 {
        return;
    }

    // Find sheet runs
    let runs = find_sheet_runs(gps);

    for (first, last) in runs {
        let f = 1usize; // PyMOL uses fixed window f=1

        for _ in 0..flat_cycles {
            // Temporary buffers
            let mut tmp_pos = vec![Vec3::new(0.0, 0.0, 0.0); n];
            let mut tmp_orient = vec![Vec3::new(0.0, 0.0, 0.0); n];

            // Average positions
            for b in (first + f)..=(last.saturating_sub(f)) {
                let mut sum = Vec3::new(0.0, 0.0, 0.0);
                for e in (b - f)..=(b + f) {
                    sum = sum + gps[e].position;
                }
                tmp_pos[b] = sum / (2 * f + 1) as f32;
            }
            for b in (first + f)..=(last.saturating_sub(f)) {
                gps[b].position = tmp_pos[b];
            }

            // Average orientations
            for b in (first + f)..=(last.saturating_sub(f)) {
                let mut sum = Vec3::new(0.0, 0.0, 0.0);
                for e in (b - f)..=(b + f) {
                    sum = sum + gps[e].orientation;
                }
                tmp_orient[b] = sum / (2 * f + 1) as f32;
            }
            for b in (first + f)..=(last.saturating_sub(f)) {
                gps[b].orientation = tmp_orient[b];
            }

            // Remove tangent component from orientations
            // PyMOL computes tangent from smoothed positions: tmp[b] = normalize(pv[b+1] - pv[b-1])
            for b in (first + f)..=(last.saturating_sub(f)) {
                let prev = if b > 0 { b - 1 } else { 0 };
                let next = (b + 1).min(n - 1);
                let tangent = normalize_safe(gps[next].position - gps[prev].position);

                let orient = gps[b].orientation;
                let removed = orient - tangent * orient.dot(tangent);
                gps[b].orientation = normalize_safe(removed);
            }
        }
    }
}

/// Smooth loop regions with progressive window sizes.
///
/// Matches RepCartoonSmoothLoops: identifies loop regions (ss == NONE/Loop),
/// extends them by 1 residue into adjacent SS regions if within same segment,
/// then applies windowed averaging.
fn smooth_loops(
    gps: &mut [GuidePoint],
    smooth_first: u32,
    smooth_last: u32,
    smooth_cycles: u32,
) {
    let n = gps.len();
    if n < 3 {
        return;
    }

    // Find loop runs (non-helix, non-sheet)
    let runs = find_loop_runs(gps);
    let mut tmp = vec![Vec3::new(0.0, 0.0, 0.0); n];

    for (mut first, mut last) in runs {
        // PyMOL extends loop regions by 1 into adjacent segments
        if first > 0 {
            first -= 1;
        }
        if last < n - 1 {
            last += 1;
        }

        for f in smooth_first..=smooth_last {
            let f = f as usize;

            for _ in 0..smooth_cycles {
                // Average positions
                for b in (first + f)..=(last.saturating_sub(f)) {
                    let mut sum = Vec3::new(0.0, 0.0, 0.0);
                    for e in (b - f)..=(b + f) {
                        sum = sum + gps[e].position;
                    }
                    tmp[b] = sum / (2 * f + 1) as f32;
                }
                for b in (first + f)..=(last.saturating_sub(f)) {
                    gps[b].position = tmp[b];
                }

                // Average orientations
                for b in (first + f)..=(last.saturating_sub(f)) {
                    let mut sum = Vec3::new(0.0, 0.0, 0.0);
                    for e in (b - f)..=(b + f) {
                        sum = sum + gps[e].orientation;
                    }
                    tmp[b] = sum / (2 * f + 1) as f32;
                }
                for b in (first + f)..=(last.saturating_sub(f)) {
                    gps[b].orientation = normalize_safe(tmp[b]);
                }
            }
        }
    }
}

/// Find contiguous sheet/nucleic-ribbon runs. Returns (first, last) index pairs.
fn find_sheet_runs(gps: &[GuidePoint]) -> Vec<(usize, usize)> {
    let mut runs = Vec::new();
    let mut first: Option<usize> = None;

    for (a, gp) in gps.iter().enumerate() {
        if gp.ss_type.is_flat_ribbon() {
            if first.is_none() {
                first = Some(a);
            }
        } else {
            if let Some(f) = first {
                runs.push((f, a - 1));
                first = None;
            }
        }
    }
    if let Some(f) = first {
        runs.push((f, gps.len() - 1));
    }
    runs
}

/// Find contiguous loop runs (not helix, not sheet, not nucleic ribbon). Returns (first, last) index pairs.
fn find_loop_runs(gps: &[GuidePoint]) -> Vec<(usize, usize)> {
    let mut runs = Vec::new();
    let mut first: Option<usize> = None;

    for (a, gp) in gps.iter().enumerate() {
        let is_loop = !is_helix(gp.ss_type) && !gp.ss_type.is_flat_ribbon();
        if is_loop {
            if first.is_none() {
                first = Some(a);
            }
        } else {
            if let Some(f) = first {
                runs.push((f, a - 1));
                first = None;
            }
        }
    }
    if let Some(f) = first {
        runs.push((f, gps.len() - 1));
    }
    runs
}

// ============================================================================
// Phase 6: Per-Run Interpolation and Mesh Generation
// ============================================================================

/// Interpolated point in the extrude buffer
struct ExtrudePoint {
    position: Vec3,
    orientation: Vec3,
    color: [f32; 4],
}

/// CartoonGenerateSample — interpolate between two guide points.
///
/// For each pair of guide atoms, generates `sampling` interpolated points.
/// The first invocation (n_p==0) also generates a starting point (total: sampling+1
/// for the first pair, sampling for subsequent pairs).
fn cartoon_generate_sample(
    buffer: &mut Vec<ExtrudePoint>,
    n_p: &mut usize,
    gp1: &GuidePoint,
    gp2: &GuidePoint,
    tv1: &Vec3,
    tv2: &Vec3,
    dev: f32,
    sampling: u32,
    power_a: f32,
    power_b: f32,
) {
    let pos1 = gp1.position;
    let pos2 = gp2.position;
    let orient1 = gp1.orientation;
    let orient2 = gp2.orientation;

    for b in 0..sampling {
        if *n_p == 0 {
            // First point of first segment: include starting position
            let f0_raw = b as f32 / sampling as f32;

            // Color: discrete switch at 0.5
            let color = if f0_raw <= 0.5 { gp1.color } else { gp2.color };

            let f0 = smooth(f0_raw, power_a);
            let f1 = 1.0 - f0;
            let f2 = smooth(f0, power_b);
            let f3 = smooth(f1, power_b);
            let f4 = dev * f2 * f3;

            let pos = pos1 * f1 + pos2 * f0 + (*tv1 * f3 - *tv2 * f2) * f4;

            // Store orientation (will be used later for frame construction)
            // PyMOL: copy3f(vo, vn - 6) — copies the orientation of guide point 1
            // This is only for the very first point, it gets the raw orientation
            let orient = orient1;

            buffer.push(ExtrudePoint {
                position: pos,
                orientation: orient,
                color,
            });
            *n_p += 1;
        }

        // Generate the next interpolated point
        let f0_raw = (b as f32 + 1.0) / sampling as f32;

        // Color: discrete switch at 0.5
        let color = if f0_raw <= 0.5 { gp1.color } else { gp2.color };

        let f0 = smooth(f0_raw, power_a);
        let f1 = 1.0 - f0;
        let f2 = smooth(f0, power_b);
        let f3 = smooth(f1, power_b);
        let f4 = dev * f2 * f3;

        let pos = pos1 * f1 + pos2 * f0 + (*tv1 * f3 - *tv2 * f2) * f4;
        let orient = normalize_safe(orient1 * (f1 * f2) + orient2 * (f0 * f3));

        buffer.push(ExtrudePoint {
            position: pos,
            orientation: orient,
            color,
        });

        // Last sample: override orientation with gp2's raw orientation
        // PyMOL: if(b == sampling - 1) copy3f(vo + 3, vn - 6)
        if b == sampling - 1 {
            buffer.last_mut().unwrap().orientation = orient2;
        }

        *n_p += 1;
    }
}

/// Refine interpolated positions after each CartoonGenerateSample call.
///
/// Matches CartoonGenerateRefine (RepCartoon.cpp line 2165):
/// Smooths positions along the axis perpendicular to the ribbon plane
/// (cross product of consecutive guide point orientations).
/// This removes the "wobble" from the helical Cα trace.
fn cartoon_generate_refine(
    buffer: &mut [ExtrudePoint],
    n_p: usize,
    sampling: u32,
    refine_cycles: u32,
    orient_a: &Vec3,
    orient_b: &Vec3,
) {
    let sampling = sampling as usize;
    if sampling <= 1 || refine_cycles == 0 || n_p < sampling {
        return;
    }

    // t0 = cross(orient_a, orient_b) — axis perpendicular to ribbon plane
    let t0 = orient_a.cross(*orient_b);
    if t0.magnitude() < 1e-4 {
        return;
    }
    let t0 = t0 / t0.magnitude();

    // The last `sampling` points in the buffer start at n_p - sampling
    // But we also need the point BEFORE them (at n_p - sampling - 1)
    let start = n_p - sampling;
    if start == 0 {
        return; // Need at least one point before the window
    }

    let mut tmp = vec![Vec3::new(0.0, 0.0, 0.0); sampling];

    for _ in 0..refine_cycles {
        // For b in 0..sampling-1: smooth positions along t0
        // p0 = buffer[start - 1 + b], p1 = buffer[start + b], p2 = buffer[start + 1 + b]
        for b in 0..sampling - 1 {
            let idx0 = start - 1 + b;
            let idx1 = start + b;
            let idx2 = start + 1 + b;

            let f0 = t0.dot(buffer[idx0].position);
            let f1 = t0.dot(buffer[idx1].position);
            let f2 = t0.dot(buffer[idx2].position);

            let f3 = (f2 + f0) / 2.0;
            // Move p1 along t0 by (f3 - f1)
            tmp[b] = buffer[idx1].position + t0 * (f3 - f1);
        }

        // Apply
        for b in 0..sampling - 1 {
            buffer[start + b].position = tmp[b];
        }
    }
}

/// Compute tangent vectors from interpolated positions.
///
/// Matches ExtrudeComputeTangents: uses adjacent differences, averaged at interior.
fn compute_tangents_from_positions(positions: &[Vec3]) -> Vec<Vec3> {
    let n = positions.len();
    if n < 2 {
        return vec![Vec3::new(0.0, 0.0, 1.0); n];
    }

    // First compute normalized differences
    let mut nv = Vec::with_capacity(n - 1);
    for a in 0..n - 1 {
        let diff = positions[a + 1] - positions[a];
        nv.push(normalize_safe(diff));
    }

    let mut tangents = Vec::with_capacity(n);

    // First: copy first difference
    tangents.push(nv[0]);

    // Interior: average adjacent differences
    for a in 1..n - 1 {
        tangents.push(normalize_safe(nv[a] + nv[a - 1]));
    }

    // Last: copy last difference
    tangents.push(nv[n - 2]);

    tangents
}

/// Extrude a run of interpolated points into mesh geometry.
///
/// Builds reference frames from (position, tangent, orientation), selects
/// the appropriate profile, and generates mesh vertices/indices.
fn extrude_run(
    buffer: &[ExtrudePoint],
    car_type: CartoonType,
    geom: &CartoonGeometrySettings,
    all_vertices: &mut Vec<MeshVertex>,
    all_indices: &mut Vec<u32>,
) {
    let n_p = buffer.len();
    if n_p < 2 {
        return;
    }

    // Collect positions for tangent computation
    let positions: Vec<Vec3> = buffer.iter().map(|p| p.position).collect();
    let tangents = compute_tangents_from_positions(&positions);

    // Build reference frames: ExtrudeBuildNormals2f
    let frames: Vec<FrameWithMetadata> = buffer
        .iter()
        .enumerate()
        .map(|(i, ep)| {
            let frame = ReferenceFrame::new(positions[i], tangents[i], ep.orientation);
            FrameWithMetadata {
                frame,
                color: ep.color,
                ss_type: car_type_to_ss(car_type),
                segment_idx: 0,
            }
        })
        .collect();

    // Generate mesh based on cartoon type
    match car_type {
        CartoonType::Oval => {
            // PyMOL: wide axis (oval_length=1.35) along binormal, thin (oval_width=0.25) along normal
            // Profile::ellipse(width=normal_axis, height=binormal_axis)
            let profile = Profile::ellipse(geom.helix_width, geom.helix_height, geom.quality);
            extrude_tube(all_vertices, all_indices, &frames, &profile);
            // Caps using the same corrected profile
            if frames.len() >= 2 {
                generate_cap_from_profile(all_vertices, all_indices, &frames[0], &profile, true);
                generate_cap_from_profile(all_vertices, all_indices, frames.last().unwrap(), &profile, false);
            }
        }
        CartoonType::Loop => {
            let profile = Profile::circle(geom.loop_radius, geom.quality);
            extrude_tube(all_vertices, all_indices, &frames, &profile);
            if frames.len() >= 2 {
                generate_cap_from_profile(all_vertices, all_indices, &frames[0], &profile, true);
                generate_cap_from_profile(all_vertices, all_indices, frames.last().unwrap(), &profile, false);
            }
        }
        CartoonType::Rect => {
            // Sheet body — use explicit sheet rendering
            let sheet_termini = find_sheet_termini(&frames);
            if sheet_termini.is_empty() {
                // Fallback: entire run is sheet
                generate_explicit_sheet(
                    all_vertices,
                    all_indices,
                    &frames,
                    0,
                    frames.len() - 1,
                    geom,
                );
            } else {
                for &(start, end) in &sheet_termini {
                    generate_explicit_sheet(all_vertices, all_indices, &frames, start, end, geom);
                }
            }
        }
        CartoonType::Arrow => {
            // Sheet with arrow — same as Rect (generate_explicit_sheet handles arrows
            // when fancy_sheets is enabled)
            let sheet_termini = find_sheet_termini(&frames);
            if sheet_termini.is_empty() {
                generate_explicit_sheet(
                    all_vertices,
                    all_indices,
                    &frames,
                    0,
                    frames.len() - 1,
                    geom,
                );
            } else {
                for &(start, end) in &sheet_termini {
                    generate_explicit_sheet(all_vertices, all_indices, &frames, start, end, geom);
                }
            }
        }
        CartoonType::NucleicRect => {
            // Nucleic acid flat ribbon — same as sheet body but NO arrow terminus.
            // Render the entire run as one explicit sheet with fancy_sheets disabled.
            let mut no_arrow_geom = geom.clone();
            no_arrow_geom.fancy_sheets = false;
            generate_explicit_sheet(
                all_vertices,
                all_indices,
                &frames,
                0,
                frames.len() - 1,
                &no_arrow_geom,
            );
        }
        CartoonType::Dumbbell => {
            // Dumbbell helix: flat ribbon + edge tubes
            // Rotated dumbbell: wide along binormal, face normals along ±normal (radial)
            let profile = make_rotated_dumbbell(geom.dumbbell_width, geom.dumbbell_length);
            let helix_regions = find_helix_regions(&frames);
            let taper_frames = 8usize;

            // Generate flat ribbon with taper
            let mut ring_starts: Vec<usize> = Vec::new();

            for (i, fm) in frames.iter().enumerate() {
                let ring_start = all_vertices.len();
                ring_starts.push(ring_start);

                let taper = calculate_dumbbell_taper_inline(i, &helix_regions, taper_frames);

                for j in 0..profile.len() {
                    let local_pos = profile.points[j];
                    // Taper the wide dimension (binormal = second coord)
                    let tapered_pos = (local_pos.0, local_pos.1 * taper);

                    let world_pos = fm.frame.transform_local(tapered_pos);
                    let world_normal = fm.frame.local_normal(profile.normals[j]);

                    all_vertices.push(MeshVertex {
                        position: [world_pos.x, world_pos.y, world_pos.z],
                        normal: [world_normal.x, world_normal.y, world_normal.z],
                        color: fm.color,
                    });
                }
            }

            // Connect rings with face strips (2 faces for dumbbell)
            if ring_starts.len() >= 2 {
                generate_face_strips(all_indices, &ring_starts, 2, profile.len());
            }

            // Generate edge tubes along binormal (wide axis)
            if !helix_regions.is_empty() {
                let sin45 = std::f32::consts::FRAC_1_SQRT_2;
                let base_offset = sin45 * geom.dumbbell_length;
                let tube_radius = geom.dumbbell_radius;
                let tube_profile = Profile::circle(tube_radius, 16);

                for &(start, end) in &helix_regions {
                    let region_len = end - start + 1;
                    let sub_n = region_len.saturating_sub(taper_frames);

                    for sign in [-1.0_f32, 1.0_f32] {
                        let mut edge_ring_starts: Vec<usize> = Vec::new();

                        for i in start..=end {
                            let fm = &frames[i];
                            let ring_start = all_vertices.len();
                            edge_ring_starts.push(ring_start);

                            let frame_in_region = i - start;
                            let taper_factor = if frame_in_region < taper_frames {
                                smooth(frame_in_region as f32 / taper_frames as f32, 2.0)
                            } else if frame_in_region > sub_n {
                                smooth((end - i) as f32 / taper_frames as f32, 2.0)
                            } else {
                                1.0
                            };

                            let offset = base_offset * sign * taper_factor;
                            let offset_pos = fm.frame.position + fm.frame.binormal * offset;

                            for j in 0..tube_profile.len() {
                                let local = tube_profile.points[j];
                                let world_pos = offset_pos
                                    + fm.frame.normal * local.0
                                    + fm.frame.binormal * local.1;
                                let world_normal = fm.frame.local_normal(tube_profile.normals[j]);
                                all_vertices.push(MeshVertex {
                                    position: [world_pos.x, world_pos.y, world_pos.z],
                                    normal: [world_normal.x, world_normal.y, world_normal.z],
                                    color: fm.color,
                                });
                            }
                        }

                        for i in 0..edge_ring_starts.len().saturating_sub(1) {
                            connect_rings(
                                all_indices,
                                edge_ring_starts[i] as u32,
                                edge_ring_starts[i + 1] as u32,
                                tube_profile.len() as u32,
                            );
                        }
                    }
                }
            }
        }
    }
}

/// Map CartoonType back to SecondaryStructure for frame metadata.
///
/// Note: for NucleicRibbon, the source SS type is preserved by passing it
/// through the guide points, not through this mapping. This mapping is only
/// used when creating FrameWithMetadata from extrude points, where the
/// original SS information has already been mapped to CartoonType.
fn car_type_to_ss(car: CartoonType) -> SecondaryStructure {
    match car {
        CartoonType::Oval | CartoonType::Dumbbell => SecondaryStructure::Helix,
        CartoonType::Rect | CartoonType::Arrow => SecondaryStructure::Sheet,
        CartoonType::NucleicRect => SecondaryStructure::NucleicRibbon,
        CartoonType::Loop => SecondaryStructure::Loop,
    }
}

/// Inline dumbbell taper calculation (avoids importing from geometry)
fn calculate_dumbbell_taper_inline(
    frame_idx: usize,
    helix_regions: &[(usize, usize)],
    taper_frames: usize,
) -> f32 {
    if taper_frames == 0 {
        return 1.0;
    }
    for &(start, end) in helix_regions {
        if frame_idx >= start && frame_idx <= end {
            let region_len = end - start + 1;
            let sub_n = region_len.saturating_sub(taper_frames);
            if frame_idx < start + taper_frames {
                let f = (frame_idx - start) as f32 / taper_frames as f32;
                return smooth(f, 2.0);
            } else if frame_idx > start + sub_n {
                let f = (end - frame_idx) as f32 / taper_frames as f32;
                return smooth(f, 2.0);
            }
            return 1.0;
        }
    }
    1.0
}

/// Create a dumbbell (flat ribbon) cross-section profile rotated 90°.
///
/// The resulting profile has:
/// - **Wide dimension** (`length`) along the binormal (second local coord)
/// - **Thin dimension** (`width`) along the normal (first local coord)
/// - **Face normals** pointing along ±normal (radially outward)
///
/// This matches PyMOL's convention when used with round helix orientations,
/// producing a flat ribbon with two visible faces (top and bottom).
///
/// The profile has 4 vertices forming 2 faces, connected by `generate_face_strips`.
fn make_rotated_dumbbell(width: f32, length: f32) -> Profile {
    let w = std::f32::consts::FRAC_1_SQRT_2 * width;
    let l = std::f32::consts::FRAC_1_SQRT_2 * length;

    // Rotated: width along normal (first coord), length along binormal (second coord)
    // Face normals along ±normal
    let points = vec![
        (w, -l),  // top face, left edge
        (w, l),   // top face, right edge
        (-w, -l), // bottom face, left edge
        (-w, l),  // bottom face, right edge
    ];

    let normals = vec![
        (1.0, 0.0),  // top face normal → +normal (radially outward)
        (1.0, 0.0),
        (-1.0, 0.0), // bottom face normal → -normal (radially inward)
        (-1.0, 0.0),
    ];

    Profile {
        points,
        normals,
        profile_type: super::geometry::ProfileType::Flat2Face,
    }
}

/// Generate a flat cap disc from the same profile used for the tube body.
///
/// Creates a triangle fan centered at the frame position, with edge vertices
/// placed at the profile points transformed into world space. The cap normal
/// points along the tangent direction (negated for start caps, positive for end caps).
///
/// This ensures cap dimensions match the tube exactly, regardless of profile shape
/// (circle for loops, ellipse for helices).
fn generate_cap_from_profile(
    vertices: &mut Vec<MeshVertex>,
    indices: &mut Vec<u32>,
    frame_meta: &FrameWithMetadata,
    profile: &Profile,
    is_start: bool,
) {
    let normal_dir = if is_start { -1.0 } else { 1.0 };
    let cap_normal = frame_meta.frame.tangent * normal_dir;

    // Center vertex
    let center_idx = vertices.len() as u32;
    vertices.push(MeshVertex {
        position: [
            frame_meta.frame.position.x,
            frame_meta.frame.position.y,
            frame_meta.frame.position.z,
        ],
        normal: [cap_normal.x, cap_normal.y, cap_normal.z],
        color: frame_meta.color,
    });

    // Edge vertices from profile
    let edge_start = vertices.len() as u32;
    for point in &profile.points {
        let world_pos = frame_meta.frame.transform_local(*point);
        vertices.push(MeshVertex {
            position: [world_pos.x, world_pos.y, world_pos.z],
            normal: [cap_normal.x, cap_normal.y, cap_normal.z],
            color: frame_meta.color,
        });
    }

    // Fan triangles
    let n = profile.len() as u32;
    for j in 0..n {
        let j_next = (j + 1) % n;
        if is_start {
            indices.push(center_idx);
            indices.push(edge_start + j_next);
            indices.push(edge_start + j);
        } else {
            indices.push(center_idx);
            indices.push(edge_start + j);
            indices.push(edge_start + j_next);
        }
    }
}

/// Generate tube mesh from frames + round profile (circle/ellipse)
fn extrude_tube(
    vertices: &mut Vec<MeshVertex>,
    indices: &mut Vec<u32>,
    frames: &[FrameWithMetadata],
    profile: &Profile,
) {
    let mut ring_starts: Vec<usize> = Vec::new();

    for fm in frames.iter() {
        let ring_start = vertices.len();
        ring_starts.push(ring_start);

        for j in 0..profile.len() {
            let world_pos = fm.frame.transform_local(profile.points[j]);
            let world_normal = fm.frame.local_normal(profile.normals[j]);

            vertices.push(MeshVertex {
                position: [world_pos.x, world_pos.y, world_pos.z],
                normal: [world_normal.x, world_normal.y, world_normal.z],
                color: fm.color,
            });
        }
    }

    // Connect rings
    for i in 0..ring_starts.len().saturating_sub(1) {
        connect_rings(
            indices,
            ring_starts[i] as u32,
            ring_starts[i + 1] as u32,
            profile.len() as u32,
        );
    }
}

// ============================================================================
// Top-level Pipeline
// ============================================================================

/// Complete cartoon pipeline for one backbone segment.
///
/// Implements the full PyMOL cartoon rendering pipeline:
/// 1. Compute differences, normals, tangents
/// 2. Compute round helix orientations
/// 3. Refine normals
/// 4. Flatten sheets, smooth loops
/// 5. Recompute tangents after position modifications
/// 6. Per-run interpolation and mesh generation
pub fn generate_segment_cartoon(
    segment: &mut BackboneSegment,
    settings: &PipelineSettings,
    geom_settings: &CartoonGeometrySettings,
) -> (Vec<MeshVertex>, Vec<u32>) {
    let n = segment.guide_points.len();
    if n < 2 {
        return (Vec::new(), Vec::new());
    }

    let gps = &mut segment.guide_points;

    // Phase 2: differences, normals, tangents
    let (_dl, nv) = compute_differences_and_normals(gps);
    let tv = compute_tangents(&nv, n);

    // Phase 3: round helix orientations
    if settings.round_helices {
        compute_round_helix_orientations(gps, &tv);
    }

    // Phase 4: refine normals
    if settings.refine_normals {
        refine_normals(gps, &tv, &nv);
    }

    // Phase 5: flatten sheets, smooth loops
    flatten_sheets(gps, settings.flat_cycles);
    if settings.smooth_loops {
        smooth_loops(gps, settings.smooth_first, settings.smooth_last, settings.smooth_cycles);
    }

    // Recompute after position changes
    let (dl, nv) = compute_differences_and_normals(gps);
    let tv = compute_tangents(&nv, n);

    // Phase 6: per-run interpolation and mesh generation
    generate_per_run_mesh(gps, &tv, &dl, geom_settings, settings)
}

/// Generate mesh from guide points using per-run extrusion.
///
/// Implements GenerateRepCartoonCGO: iterates through guide points,
/// detects cartoon type changes, interpolates each run, then extrudes.
fn generate_per_run_mesh(
    gps: &[GuidePoint],
    tv: &[Vec3],
    dl: &[f32],
    geom: &CartoonGeometrySettings,
    settings: &PipelineSettings,
) -> (Vec<MeshVertex>, Vec<u32>) {
    let mut all_vertices = Vec::new();
    let mut all_indices = Vec::new();

    let n = gps.len();
    if n < 2 {
        return (all_vertices, all_indices);
    }

    // Extrude buffer for current run
    let mut buffer: Vec<ExtrudePoint> = Vec::new();
    let mut n_p: usize = 0;
    let mut cur_car: Option<CartoonType> = None;

    let sampling = settings.sampling;
    let power_a = settings.power_a;
    let power_b = settings.power_b;

    // PyMOL-style loop: determine cartoon type from the PAIR (a, a+1).
    //
    // **Boundary pair rule**: when two adjacent guide points have different SS types
    // (e.g. Helix→Loop, Sheet→Loop), the pair is rendered as Loop type. This creates
    // smooth tube transitions between SS elements. Only pairs where BOTH endpoints
    // share the same SS type get that type's rendering (Oval, Rect, etc.).
    //
    // When the cartoon type changes between consecutive pairs, the current run is
    // flushed (extruded) and a new run begins. The boundary pair always starts
    // a new Loop run.
    let mut a = 0;
    while a < n {
        // Determine pair type: boundary pairs (different SS types) → always Loop
        let pair_car = if a < n - 1 {
            let car_a = cartoon_type_for(gps[a].ss_type, geom);
            let car_b = cartoon_type_for(gps[a + 1].ss_type, geom);
            if car_a == car_b { car_a } else { CartoonType::Loop }
        } else {
            cartoon_type_for(gps[a].ss_type, geom)
        };

        // Type changed from previous pair — flush the accumulated run
        let type_changed = cur_car.map_or(false, |c| c != pair_car);

        if type_changed && n_p > 0 {
            extrude_run(
                &buffer[..n_p],
                cur_car.unwrap(),
                geom,
                &mut all_vertices,
                &mut all_indices,
            );
            buffer.clear();
            n_p = 0;
        }

        cur_car = Some(pair_car);

        // Interpolate between a and a+1
        if a < n - 1 {
            let dev = settings.throw_factor * dl[a];

            cartoon_generate_sample(
                &mut buffer,
                &mut n_p,
                &gps[a],
                &gps[a + 1],
                &tv[a],
                &tv[a + 1],
                dev,
                sampling,
                power_a,
                power_b,
            );

            // Refine positions (CartoonGenerateRefine) — only for same-type pairs
            let car_a = cartoon_type_for(gps[a].ss_type, geom);
            let car_b = cartoon_type_for(gps[a + 1].ss_type, geom);
            if settings.refine > 0 && car_a == car_b {
                cartoon_generate_refine(
                    &mut buffer,
                    n_p,
                    sampling,
                    settings.refine,
                    &gps[a].orientation,
                    &gps[a + 1].orientation,
                );
            }
        }

        a += 1;
    }

    // Extrude final run
    if n_p > 0 {
        if let Some(car) = cur_car {
            extrude_run(&buffer[..n_p], car, geom, &mut all_vertices, &mut all_indices);
        }
    }

    (all_vertices, all_indices)
}

#[cfg(test)]
mod tests {
    use super::*;
    use lin_alg::f32::Vec3;
    use pymol_mol::{AtomIndex, SecondaryStructure};

    fn make_gp(x: f32, ss: SecondaryStructure) -> GuidePoint {
        GuidePoint::new(
            Vec3::new(x, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            [1.0, 1.0, 1.0, 1.0],
            ss,
            AtomIndex(0),
            1,
        )
    }

    #[test]
    fn test_compute_differences_and_normals() {
        let gps = vec![
            make_gp(0.0, SecondaryStructure::Loop),
            make_gp(3.8, SecondaryStructure::Loop),
            make_gp(7.6, SecondaryStructure::Loop),
        ];

        let (dl, nv) = compute_differences_and_normals(&gps);
        assert_eq!(dl.len(), 2);
        assert!((dl[0] - 3.8).abs() < 0.01);
        assert!((dl[1] - 3.8).abs() < 0.01);
        assert!((nv[0].x - 1.0).abs() < 0.01);
        assert!((nv[1].x - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_compute_tangents() {
        let nv = vec![
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        ];
        let tv = compute_tangents(&nv, 3);
        assert_eq!(tv.len(), 3);
        // First = nv[0]
        assert!((tv[0].x - 1.0).abs() < 0.01);
        // Middle = normalize(nv[0] + nv[1]) = normalize((1,1,0))
        assert!((tv[1].x - 0.707).abs() < 0.01);
        assert!((tv[1].y - 0.707).abs() < 0.01);
        // Last = nv[1]
        assert!((tv[2].y - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_cartoon_generate_sample_first_pair() {
        let gp1 = make_gp(0.0, SecondaryStructure::Loop);
        let gp2 = make_gp(3.8, SecondaryStructure::Loop);
        let tv1 = Vec3::new(1.0, 0.0, 0.0);
        let tv2 = Vec3::new(1.0, 0.0, 0.0);

        let mut buffer = Vec::new();
        let mut n_p = 0usize;

        cartoon_generate_sample(
            &mut buffer,
            &mut n_p,
            &gp1,
            &gp2,
            &tv1,
            &tv2,
            1.35 * 3.8,
            7,
            2.0,
            0.52,
        );

        // First pair: sampling + 1 points (starts with b=0 point, plus sampling regular points)
        assert_eq!(n_p, 8); // 7 + 1 extra start point
        assert_eq!(buffer.len(), 8);

        // First point should be near pos1
        assert!((buffer[0].position.x - 0.0).abs() < 0.1);
        // Last point should be near pos2
        assert!((buffer[7].position.x - 3.8).abs() < 0.1);
    }

    #[test]
    fn test_cartoon_generate_sample_subsequent_pair() {
        let gp1 = make_gp(0.0, SecondaryStructure::Loop);
        let gp2 = make_gp(3.8, SecondaryStructure::Loop);
        let gp3 = make_gp(7.6, SecondaryStructure::Loop);
        let tv1 = Vec3::new(1.0, 0.0, 0.0);
        let tv2 = Vec3::new(1.0, 0.0, 0.0);
        let tv3 = Vec3::new(1.0, 0.0, 0.0);

        let mut buffer = Vec::new();
        let mut n_p = 0usize;

        // First pair
        cartoon_generate_sample(&mut buffer, &mut n_p, &gp1, &gp2, &tv1, &tv2, 1.35 * 3.8, 7, 2.0, 0.52);
        assert_eq!(n_p, 8);

        // Second pair: n_p > 0, so no extra start point
        cartoon_generate_sample(&mut buffer, &mut n_p, &gp2, &gp3, &tv2, &tv3, 1.35 * 3.8, 7, 2.0, 0.52);
        assert_eq!(n_p, 15); // 8 + 7
    }

    #[test]
    fn test_generate_segment_cartoon_loop() {
        let mut segment = super::super::backbone::BackboneSegment::new("A");
        for i in 0..5 {
            segment.push(make_gp(i as f32 * 3.8, SecondaryStructure::Loop));
        }

        let settings = PipelineSettings::default();
        let geom = CartoonGeometrySettings::default();
        let (vertices, indices) = generate_segment_cartoon(&mut segment, &settings, &geom);

        assert!(!vertices.is_empty(), "Should have vertices");
        assert!(!indices.is_empty(), "Should have indices");
    }

    #[test]
    fn test_generate_segment_cartoon_mixed_ss() {
        let mut segment = super::super::backbone::BackboneSegment::new("A");

        // Helix region
        for i in 0..4 {
            segment.push(GuidePoint::new(
                Vec3::new(i as f32 * 3.8, 2.3 * (i as f32 * 1.7).cos(), 2.3 * (i as f32 * 1.7).sin()),
                Vec3::new(0.0, 1.0, 0.0),
                [1.0, 0.0, 0.0, 1.0],
                SecondaryStructure::Helix,
                AtomIndex(i),
                i as i32 + 1,
            ));
        }
        // Loop region
        for i in 4..7 {
            segment.push(make_gp(i as f32 * 3.8, SecondaryStructure::Loop));
        }
        // Sheet region
        for i in 7..10 {
            segment.push(make_gp(i as f32 * 3.8, SecondaryStructure::Sheet));
        }

        let settings = PipelineSettings::default();
        let geom = CartoonGeometrySettings::default();
        let (vertices, indices) = generate_segment_cartoon(&mut segment, &settings, &geom);

        assert!(!vertices.is_empty(), "Mixed SS should have vertices");
        assert!(!indices.is_empty(), "Mixed SS should have indices");
    }
}
