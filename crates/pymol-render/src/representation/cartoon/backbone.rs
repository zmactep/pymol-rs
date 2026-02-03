//! Backbone trace extraction for cartoon representation
//!
//! Extracts C-alpha (CA) positions and carbonyl oxygen (O) directions from
//! protein residues to create guide points for cartoon rendering.
//!
//! This module also implements PyMOL-style smoothing algorithms:
//! - Helix axis centering (ExtrudeShiftToAxis)
//! - Sheet flattening (RepCartoonFlattenSheets)
//! - Loop smoothing (RepCartoonSmoothLoops)

use lin_alg::f32::Vec3;
use pymol_mol::{AtomIndex, CoordSet, ObjectMolecule, RepMask, SecondaryStructure};

/// Settings for cartoon smoothing operations
#[derive(Debug, Clone)]
pub struct CartoonSmoothSettings {
    /// Number of orientation smoothing cycles (cartoon_smooth_cycles)
    pub smooth_cycles: u32,
    /// Number of sheet flattening cycles (cartoon_flat_cycles)
    pub flat_cycles: u32,
    /// First window size for loop smoothing (cartoon_smooth_first)
    pub smooth_first: u32,
    /// Last window size for loop smoothing (cartoon_smooth_last)
    pub smooth_last: u32,
    /// Number of spline refinement cycles (cartoon_refine)
    #[allow(dead_code)]
    pub refine_cycles: u32,
    /// Whether to refine normals (cartoon_refine_normals)
    pub refine_normals: bool,
    /// Power for the smooth sigmoid function (cartoon_power)
    #[allow(dead_code)]
    pub power_a: f32,
    /// Power for displacement calculation (cartoon_power_b)
    #[allow(dead_code)]
    pub power_b: f32,
    /// Displacement throw factor (cartoon_throw)
    #[allow(dead_code)]
    pub throw_factor: f32,
}

impl Default for CartoonSmoothSettings {
    fn default() -> Self {
        Self {
            smooth_cycles: 2,   // PyMOL default
            flat_cycles: 4,     // PyMOL default
            smooth_first: 1,    // PyMOL default
            smooth_last: 1,     // PyMOL default
            refine_cycles: 5,   // PyMOL default
            refine_normals: true,
            power_a: 2.0,       // PyMOL default (cartoon_power)
            power_b: 0.52,      // PyMOL default (cartoon_power_b)
            throw_factor: 1.35, // PyMOL default (cartoon_throw)
        }
    }
}

/// A guide point for cartoon rendering
///
/// Each guide point represents a residue's contribution to the cartoon trace.
/// It contains the CA position, O direction for ribbon orientation, and color.
#[derive(Debug, Clone)]
pub struct GuidePoint {
    /// Position (typically CA or interpolated backbone position)
    pub position: Vec3,
    /// Orientation vector (typically CA -> O direction for ribbon normal)
    pub orientation: Vec3,
    /// RGBA color for this guide point
    pub color: [f32; 4],
    /// Secondary structure type
    pub ss_type: SecondaryStructure,
    /// Reference to the CA atom index (for selection/picking)
    #[allow(dead_code)]
    pub atom_idx: AtomIndex,
    /// Residue sequence number (for gap detection)
    #[allow(dead_code)]
    pub resv: i32,
    /// B-factor (for putty representation)
    pub b_factor: f32,
}

impl GuidePoint {
    /// Create a new guide point
    pub fn new(
        position: Vec3,
        orientation: Vec3,
        color: [f32; 4],
        ss_type: SecondaryStructure,
        atom_idx: AtomIndex,
        resv: i32,
        b_factor: f32,
    ) -> Self {
        Self {
            position,
            orientation,
            color,
            ss_type,
            atom_idx,
            resv,
            b_factor,
        }
    }
}

/// A continuous segment of backbone for cartoon rendering
///
/// Chain breaks or missing residues split the backbone into multiple segments.
#[derive(Debug, Clone)]
pub struct BackboneSegment {
    /// Chain identifier
    #[allow(dead_code)]
    pub chain_id: String,
    /// Guide points along this segment
    pub guide_points: Vec<GuidePoint>,
}

impl BackboneSegment {
    /// Create a new backbone segment
    pub fn new(chain_id: impl Into<String>) -> Self {
        Self {
            chain_id: chain_id.into(),
            guide_points: Vec::new(),
        }
    }

    /// Add a guide point to the segment
    pub fn push(&mut self, point: GuidePoint) {
        self.guide_points.push(point);
    }

    /// Check if the segment is empty
    pub fn is_empty(&self) -> bool {
        self.guide_points.is_empty()
    }

    /// Get the number of guide points
    pub fn len(&self) -> usize {
        self.guide_points.len()
    }
}

/// Extract backbone segments from a molecule
///
/// This function iterates through chains and residues, extracting CA positions
/// and O directions to create guide points. It detects chain breaks based on
/// residue numbering gaps.
///
/// The `rep_mask` parameter specifies which representation visibility to check
/// (e.g., `RepMask::CARTOON` for cartoon or `RepMask::RIBBON` for ribbon).
pub fn extract_backbone_segments(
    molecule: &ObjectMolecule,
    coord_set: &CoordSet,
    colors: &crate::color_resolver::ColorResolver,
    gap_cutoff: i32,
    rep_mask: u32,
) -> Vec<BackboneSegment> {
    let mut segments = Vec::new();

    for chain in molecule.chains() {
        let mut current_segment = BackboneSegment::new(chain.id());
        let mut prev_resv: Option<i32> = None;

        for residue in chain.residues() {
            // Only process protein residues
            if !residue.is_protein() {
                continue;
            }

            // Find CA atom
            let ca = match residue.ca() {
                Some((idx, atom)) => (idx, atom),
                None => continue,
            };

            let (ca_idx, ca_atom) = ca;

            // Check if the specified representation is visible for this atom
            if !ca_atom.visible_reps.is_visible(rep_mask) {
                continue;
            }

            // Get CA coordinates
            let ca_pos = match coord_set.get_atom_coord(ca_idx) {
                Some(pos) => pos,
                None => continue,
            };

            // Check for gap (chain break)
            if let Some(prev) = prev_resv {
                let gap = (residue.resv() - prev).abs();
                if gap > gap_cutoff {
                    // Start a new segment
                    if !current_segment.is_empty() {
                        segments.push(current_segment);
                    }
                    current_segment = BackboneSegment::new(chain.id());
                }
            }

            // Find O atom for orientation
            let orientation = find_orientation_vector(&residue, coord_set, ca_pos);

            // Resolve color based on representation type
            let color = if rep_mask == RepMask::RIBBON {
                colors.resolve_ribbon(ca_atom, molecule)
            } else {
                colors.resolve_cartoon(ca_atom, molecule)
            };

            // Create guide point
            let guide_point = GuidePoint::new(
                ca_pos,
                orientation,
                color,
                ca_atom.ss_type,
                ca_idx,
                residue.resv(),
                ca_atom.b_factor,
            );

            current_segment.push(guide_point);
            prev_resv = Some(residue.resv());
        }

        // Push the final segment
        if !current_segment.is_empty() {
            segments.push(current_segment);
        }
    }

    segments
}

/// Find the orientation vector for a residue
///
/// The orientation is typically the CA -> O direction, which defines the
/// ribbon normal direction. Falls back to a default if O is not found.
fn find_orientation_vector(
    residue: &pymol_mol::ResidueView,
    coord_set: &CoordSet,
    ca_pos: Vec3,
) -> Vec3 {
    // Try to find O (carbonyl oxygen) atom
    if let Some((o_idx, _)) = residue.find_by_name("O") {
        if let Some(o_pos) = coord_set.get_atom_coord(o_idx) {
            let dir = o_pos - ca_pos;
            let len_sq = dir.magnitude_squared();
            if len_sq > 1e-6 {
                return dir / len_sq.sqrt();
            }
        }
    }

    // Try C (carbonyl carbon) as fallback
    if let Some((c_idx, _)) = residue.find_by_name("C") {
        if let Some(c_pos) = coord_set.get_atom_coord(c_idx) {
            let dir = c_pos - ca_pos;
            let len_sq = dir.magnitude_squared();
            if len_sq > 1e-6 {
                return dir / len_sq.sqrt();
            }
        }
    }

    // Default orientation if nothing found
    Vec3::new(0.0, 1.0, 0.0)
}

/// Smooth the orientation vectors along a segment
///
/// This helps prevent sudden flips in ribbon orientation.
pub fn smooth_orientations(segment: &mut BackboneSegment, cycles: u32) {
    if segment.len() < 3 {
        return;
    }

    for _ in 0..cycles {
        let len = segment.guide_points.len();
        let mut new_orientations = Vec::with_capacity(len);

        for i in 0..len {
            let prev = if i > 0 { i - 1 } else { 0 };
            let next = if i < len - 1 { i + 1 } else { len - 1 };

            // Average neighboring orientations
            let o_prev = segment.guide_points[prev].orientation;
            let o_curr = segment.guide_points[i].orientation;
            let o_next = segment.guide_points[next].orientation;

            let mut avg = o_prev + o_curr * 2.0 + o_next;
            let len_sq = avg.magnitude_squared();
            if len_sq > 1e-6 {
                avg = avg / len_sq.sqrt();
            }

            new_orientations.push(avg);
        }

        // Apply smoothed orientations
        for (i, orientation) in new_orientations.into_iter().enumerate() {
            segment.guide_points[i].orientation = orientation;
        }
    }
}

/// Apply all PyMOL-style smoothing operations to a backbone segment
///
/// This applies smoothing in the correct order:
/// 1. Helix smoothing (heavy orientation smoothing to prevent ribbon twist)
/// 2. Sheet flattening (make sheets planar) - RepCartoonFlattenSheets
/// 3. Loop smoothing (progressive window smoothing) - RepCartoonSmoothLoops  
/// 4. Normal refinement (orthogonalize and fix kinks) - RepCartoonRefineNormals
pub fn apply_pymol_smoothing(segment: &mut BackboneSegment, settings: &CartoonSmoothSettings) {
    if segment.len() < 2 {
        return;
    }

    // 1. Smooth helices (heavy orientation smoothing to prevent ribbon twist)
    smooth_helices(segment, settings.smooth_cycles);

    // 2. Flatten sheets (RepCartoonFlattenSheets)
    flatten_sheets(segment, settings.flat_cycles);

    // 3. Smooth loops (RepCartoonSmoothLoops)
    smooth_loops(segment, settings.smooth_first, settings.smooth_last, settings.smooth_cycles);

    // 4. Refine normals (RepCartoonRefineNormals)
    if settings.refine_normals {
        refine_normals(segment);
    }
}

/// Smooth positions across the entire segment (including SS boundaries)
#[allow(dead_code)]
fn global_smooth_positions(segment: &mut BackboneSegment, cycles: u32) {
    let len = segment.guide_points.len();
    if len < 2 {
        return;
    }

    for _ in 0..cycles {
        let mut new_positions = Vec::with_capacity(len);

        for i in 0..len {
            let prev = i.saturating_sub(1);
            let next = (i + 1).min(len - 1);

            let p_prev = segment.guide_points[prev].position;
            let p_curr = segment.guide_points[i].position;
            let p_next = segment.guide_points[next].position;

            // Weighted average: center gets more weight
            let avg = (p_prev + p_curr * 4.0 + p_next) / 6.0;
            new_positions.push(avg);
        }

        for (i, pos) in new_positions.into_iter().enumerate() {
            segment.guide_points[i].position = pos;
        }
    }
}

/// Smooth orientations across the entire segment
#[allow(dead_code)]
fn global_smooth_orientations(segment: &mut BackboneSegment, cycles: u32) {
    let len = segment.guide_points.len();
    if len < 2 {
        return;
    }

    for _ in 0..cycles {
        // First, ensure consistent direction
        for i in 1..len {
            let prev = segment.guide_points[i - 1].orientation;
            let curr = segment.guide_points[i].orientation;
            if prev.dot(curr) < 0.0 {
                segment.guide_points[i].orientation = curr * -1.0;
            }
        }

        // Then smooth
        let mut new_orientations = Vec::with_capacity(len);

        for i in 0..len {
            let prev = i.saturating_sub(1);
            let next = (i + 1).min(len - 1);

            let o_prev = segment.guide_points[prev].orientation;
            let o_curr = segment.guide_points[i].orientation;
            let o_next = segment.guide_points[next].orientation;

            // Weighted average
            let mut avg = o_prev + o_curr * 4.0 + o_next;
            let mag = avg.magnitude();
            if mag > 1e-6 {
                avg = avg / mag;
            }
            new_orientations.push(avg);
        }

        for (i, orient) in new_orientations.into_iter().enumerate() {
            segment.guide_points[i].orientation = orient;
        }
    }
}

// ============================================================================
// Helix Axis Centering (ExtrudeShiftToAxis equivalent)
// ============================================================================

/// Shift helix CA positions toward the helix axis
///
/// Alpha helices have CA atoms positioned ~2.3Å from the central axis.
/// This function detects helix regions and shifts their positions toward
/// the computed axis, creating a smoother, more cylindrical appearance.
///
/// The algorithm (based on PyMOL's ExtrudeShiftToAxis):
/// 1. For each helix residue, shift position toward axis along the NORMAL direction
/// 2. The normal direction (CA->O) points roughly outward from the helix axis
/// 3. Apply multiple passes of heavy smoothing to create smooth cylindrical trace
///
/// Note: This is for cylindrical helix mode only, not used for ribbon helices.
#[allow(dead_code)]
fn shift_helix_to_axis(segment: &mut BackboneSegment, shift_distance: f32) {
    let runs = detect_ss_runs(segment);

    for (start, end, ss_type) in runs {
        if !is_helix(ss_type) {
            continue;
        }

        let helix_len = end - start + 1;
        if helix_len < 2 {
            continue; // Need at least 2 residues for any smoothing
        }

        // Step 1: Shift each position toward the axis along the normal direction
        // In a helix, the normal (CA->O direction) points outward from the axis
        // So we shift OPPOSITE to the normal to move toward the axis
        for i in start..=end {
            let normal = segment.guide_points[i].orientation;
            let normal_len = normal.magnitude();
            
            if normal_len > 0.1 {
                let shift_dir = normal / normal_len;
                // Shift OPPOSITE to normal (toward the helix axis)
                segment.guide_points[i].position = 
                    segment.guide_points[i].position - shift_dir * shift_distance;
            }
        }

        // Step 2: Heavy window smoothing - multiple passes with increasing window
        // This is critical for creating smooth helices
        for window in 1..=4 {
            window_smooth_positions(segment, start, end, 3, window);
        }

        // Step 3: Smooth orientations within the helix to prevent twisting
        window_smooth_orientations_range(segment, start, end, 3, 2);

        // Step 4: Make orientations consistent along the helix
        // Propagate from start, flipping if necessary
        for i in (start + 1)..=end.min(segment.guide_points.len() - 1) {
            let prev_orient = segment.guide_points[i - 1].orientation;
            let curr_orient = segment.guide_points[i].orientation;
            
            if prev_orient.dot(curr_orient) < 0.0 {
                segment.guide_points[i].orientation = curr_orient * -1.0;
            }
        }
    }
}

// ============================================================================
// Sheet Flattening (RepCartoonFlattenSheets equivalent)
// ============================================================================

/// Flatten sheet regions to make them planar (RepCartoonFlattenSheets)
///
/// This implements PyMOL's RepCartoonFlattenSheets algorithm exactly:
/// - Uses a fixed window size of 1 for iterative averaging
/// - For each flat_cycles iteration:
///   1. Average positions with window [-1, 1]
///   2. Average orientations with window [-1, 1]
///   3. Remove tangent component from orientations (make perpendicular)
fn flatten_sheets(segment: &mut BackboneSegment, flat_cycles: u32) {
    let len = segment.guide_points.len();
    if len < 3 {
        return;
    }

    let runs = detect_ss_runs(segment);

    for (first, last, ss_type) in runs {
        if ss_type != SecondaryStructure::Sheet {
            continue;
        }

        let sheet_len = last - first + 1;
        if sheet_len < 3 {
            continue;
        }

        // PyMOL uses fixed window size f=1
        let f = 1usize;

        // Apply flat_cycles iterations of smoothing
        for _ in 0..flat_cycles {
            // Temporary storage for smoothed values
            let mut tmp_positions = vec![Vec3::new(0.0, 0.0, 0.0); len];
            let mut tmp_orientations = vec![Vec3::new(0.0, 0.0, 0.0); len];

            // Step 1: Average positions with window [-f, f]
            for b in (first + f)..=(last.saturating_sub(f)).min(len - 1) {
                let mut sum = Vec3::new(0.0, 0.0, 0.0);
                for e in (b.saturating_sub(f))..=(b + f).min(len - 1) {
                    sum = sum + segment.guide_points[e].position;
                }
                tmp_positions[b] = sum / (2 * f + 1) as f32;
            }

            // Apply smoothed positions
            for b in (first + f)..=(last.saturating_sub(f)).min(len - 1) {
                segment.guide_points[b].position = tmp_positions[b];
            }

            // Step 2: Average orientations with window [-f, f]
            for b in (first + f)..=(last.saturating_sub(f)).min(len - 1) {
                let mut sum = Vec3::new(0.0, 0.0, 0.0);
                for e in (b.saturating_sub(f))..=(b + f).min(len - 1) {
                    sum = sum + segment.guide_points[e].orientation;
                }
                tmp_orientations[b] = sum / (2 * f + 1) as f32;
            }

            // Apply smoothed orientations (without normalizing yet, like PyMOL)
            for b in (first + f)..=(last.saturating_sub(f)).min(len - 1) {
                segment.guide_points[b].orientation = tmp_orientations[b];
            }

            // Step 3: Remove tangent component from orientations (make perpendicular)
            for b in (first + f)..=(last.saturating_sub(f)).min(len - 1) {
                // Compute tangent from neighboring positions
                let prev_idx = (b.saturating_sub(1)).max(first);
                let next_idx = (b + 1).min(last).min(len - 1);
                
                let tangent = segment.guide_points[next_idx].position 
                            - segment.guide_points[prev_idx].position;
                let tangent_len = tangent.magnitude();
                
                if tangent_len > 1e-6 {
                    let tangent_norm = tangent / tangent_len;
                    
                    // Remove tangent component: orient = orient - tangent*(orient·tangent)
                    let orient = segment.guide_points[b].orientation;
                    let tangent_component = tangent_norm * orient.dot(tangent_norm);
                    let new_orient = orient - tangent_component;
                    
                    // Normalize
                    let new_len = new_orient.magnitude();
                    if new_len > 1e-6 {
                        segment.guide_points[b].orientation = new_orient / new_len;
                    }
                }
            }
        }
    }
}

/// Compute average orientation vector for a range
#[allow(dead_code)]
fn compute_average_orientation(segment: &BackboneSegment, start: usize, end: usize) -> Vec3 {
    let mut sum = Vec3::new(0.0, 0.0, 0.0);
    let mut count = 0;

    // Use first orientation as reference for consistent direction
    let reference = if start < segment.guide_points.len() {
        segment.guide_points[start].orientation
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };

    for i in start..=end.min(segment.guide_points.len() - 1) {
        let orient = segment.guide_points[i].orientation;
        // Add with consistent direction relative to reference
        if orient.dot(reference) < 0.0 {
            sum = sum - orient;
        } else {
            sum = sum + orient;
        }
        count += 1;
    }

    if count > 0 {
        let avg = sum / count as f32;
        let len = avg.magnitude();
        if len > 1e-6 {
            return avg / len;
        }
    }

    reference
}

/// Make orientations perpendicular to the local tangent direction
#[allow(dead_code)]
fn orthogonalize_orientations_to_tangent(segment: &mut BackboneSegment, start: usize, end: usize) {
    for i in start..=end.min(segment.guide_points.len() - 1) {
        // Compute local tangent
        let tangent = compute_local_tangent(segment, i);

        // Remove tangent component from orientation
        let orient = segment.guide_points[i].orientation;
        let tangent_component = tangent * orient.dot(tangent);
        let new_orient = orient - tangent_component;

        // Normalize
        let len = new_orient.magnitude();
        if len > 1e-6 {
            segment.guide_points[i].orientation = new_orient / len;
        }
    }
}

/// Compute local tangent direction at a guide point
fn compute_local_tangent(segment: &BackboneSegment, idx: usize) -> Vec3 {
    let len = segment.guide_points.len();
    if len < 2 {
        return Vec3::new(0.0, 0.0, 1.0);
    }

    let prev = idx.saturating_sub(1);
    let next = (idx + 1).min(len - 1);

    let p_prev = segment.guide_points[prev].position;
    let p_next = segment.guide_points[next].position;

    let tangent = p_next - p_prev;
    let mag = tangent.magnitude();

    if mag > 1e-6 {
        tangent / mag
    } else {
        Vec3::new(0.0, 0.0, 1.0)
    }
}

/// Refine orientations at sheet termini for arrow heads
#[allow(dead_code)]
fn refine_sheet_tips(segment: &mut BackboneSegment, start: usize, end: usize) {
    if end <= start || start >= segment.guide_points.len() {
        return;
    }

    // At the C-terminus (arrow tip), align orientation with sheet plane
    // by averaging with neighbors
    let last_idx = end.min(segment.guide_points.len() - 1);
    if last_idx > start + 1 {
        let prev_orient = segment.guide_points[last_idx - 1].orientation;
        let curr_orient = segment.guide_points[last_idx].orientation;

        let avg = (prev_orient + curr_orient) * 0.5;
        let len = avg.magnitude();
        if len > 1e-6 {
            segment.guide_points[last_idx].orientation = avg / len;
        }
    }
}

// ============================================================================
// Loop Smoothing (RepCartoonSmoothLoops equivalent)
// ============================================================================

/// Smooth loop regions with progressive window sizes (RepCartoonSmoothLoops)
///
/// This implements PyMOL's RepCartoonSmoothLoops algorithm exactly:
/// - Outer loop: window sizes from smooth_first to smooth_last
/// - Inner loop: smooth_cycles iterations per window size
/// - Averages positions and orientations within the window
/// - Normalizes orientations after smoothing
fn smooth_loops(segment: &mut BackboneSegment, smooth_first: u32, smooth_last: u32, smooth_cycles: u32) {
    let len = segment.guide_points.len();
    if len < 3 {
        return;
    }

    let runs = detect_ss_runs(segment);

    for (first, last, ss_type) in runs {
        // Only process loop/coil regions (not helices or sheets)
        if is_helix(ss_type) || ss_type == SecondaryStructure::Sheet {
            continue;
        }

        let loop_len = last - first + 1;
        if loop_len < 3 {
            continue;
        }

        // Temporary storage for smoothed values
        let mut tmp = vec![Vec3::new(0.0, 0.0, 0.0); len];

        // Progressive window smoothing: from smooth_first to smooth_last
        for f in smooth_first..=smooth_last {
            let f = f as usize;
            
            // Apply smooth_cycles iterations for this window size
            for _ in 0..smooth_cycles {
                // Step 1: Average positions with window [-f, f]
                for b in (first + f)..=(last.saturating_sub(f)).min(len - 1) {
                    let mut sum = Vec3::new(0.0, 0.0, 0.0);
                    for e in (b.saturating_sub(f))..=(b + f).min(len - 1) {
                        sum = sum + segment.guide_points[e].position;
                    }
                    tmp[b] = sum / (2 * f + 1) as f32;
                }

                // Apply smoothed positions (respecting no_smooth flag would go here)
                for b in (first + f)..=(last.saturating_sub(f)).min(len - 1) {
                    segment.guide_points[b].position = tmp[b];
                }

                // Step 2: Average orientations with window [-f, f]
                for b in (first + f)..=(last.saturating_sub(f)).min(len - 1) {
                    let mut sum = Vec3::new(0.0, 0.0, 0.0);
                    for e in (b.saturating_sub(f))..=(b + f).min(len - 1) {
                        sum = sum + segment.guide_points[e].orientation;
                    }
                    tmp[b] = sum / (2 * f + 1) as f32;
                }

                // Apply smoothed orientations and normalize
                for b in (first + f)..=(last.saturating_sub(f)).min(len - 1) {
                    let smoothed = tmp[b];
                    let mag = smoothed.magnitude();
                    if mag > 1e-6 {
                        segment.guide_points[b].orientation = smoothed / mag;
                    }
                }
            }
        }
    }
}

// ============================================================================
// Helix Smoothing
// ============================================================================

/// Ensure helix orientations are consistent (no sudden flips)
///
/// Smooth helices and handle helix-loop transitions
///
/// This ensures orientation consistency within helices and creates smooth
/// transitions at helix termini to prevent pinching/twisting artifacts.
fn smooth_helices(segment: &mut BackboneSegment, _smooth_cycles: u32) {
    let len = segment.guide_points.len();
    if len < 2 {
        return;
    }

    let runs = detect_ss_runs(segment);

    for (first, last, ss_type) in runs {
        // Only process helix regions
        if !is_helix(ss_type) {
            continue;
        }

        // Ensure orientation consistency within the helix
        // Propagate from the first residue, flipping any that point opposite
        for i in (first + 1)..=last.min(len - 1) {
            let prev_orient = segment.guide_points[i - 1].orientation;
            let curr_orient = segment.guide_points[i].orientation;

            if prev_orient.dot(curr_orient) < 0.0 {
                segment.guide_points[i].orientation = curr_orient * -1.0;
            }
        }

        // Handle helix-loop transition at the C-terminus (end of helix)
        // Extend the helix orientation pattern into 1-2 adjacent loop residues
        if last + 1 < len {
            let helix_end_orient = segment.guide_points[last].orientation;
            let next_idx = last + 1;
            let next_orient = segment.guide_points[next_idx].orientation;
            
            // Check if the next residue is a loop (not helix or sheet)
            let next_ss = segment.guide_points[next_idx].ss_type;
            if !is_helix(next_ss) && next_ss != SecondaryStructure::Sheet {
                // Ensure consistent direction with helix end
                if helix_end_orient.dot(next_orient) < 0.0 {
                    segment.guide_points[next_idx].orientation = next_orient * -1.0;
                }
                
                // Blend the first loop residue's orientation toward the helix end
                let blended = normalize_safe(
                    helix_end_orient * 0.6 + segment.guide_points[next_idx].orientation * 0.4
                );
                segment.guide_points[next_idx].orientation = blended;
            }
        }

        // Handle loop-helix transition at the N-terminus (start of helix)
        if first > 0 {
            let helix_start_orient = segment.guide_points[first].orientation;
            let prev_idx = first - 1;
            let prev_orient = segment.guide_points[prev_idx].orientation;
            
            // Check if the previous residue is a loop
            let prev_ss = segment.guide_points[prev_idx].ss_type;
            if !is_helix(prev_ss) && prev_ss != SecondaryStructure::Sheet {
                // Ensure consistent direction with helix start
                if helix_start_orient.dot(prev_orient) < 0.0 {
                    segment.guide_points[prev_idx].orientation = prev_orient * -1.0;
                }
                
                // Blend the last loop residue's orientation toward the helix start
                let blended = normalize_safe(
                    helix_start_orient * 0.6 + segment.guide_points[prev_idx].orientation * 0.4
                );
                segment.guide_points[prev_idx].orientation = blended;
            }
        }
    }
}

/// Normalize a vector safely
fn normalize_safe(v: Vec3) -> Vec3 {
    let len_sq = v.magnitude_squared();
    if len_sq > 1e-10 {
        v / len_sq.sqrt()
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    }
}

// ============================================================================
// Normal Refinement (RepCartoonRefineNormals equivalent)
// ============================================================================

/// Refine orientation vectors to prevent flips and kinks (RepCartoonRefineNormals)
///
/// This implements PyMOL's RepCartoonRefineNormals algorithm exactly:
/// 1. Make orientations orthogonal to tangents (for interior residues only)
/// 2. Generate alternative inverted orientations (but NOT for helices)
/// 3. Forward iterate through pairs to select optimal orientation
/// 4. Detect and soften kinks where dot(prev,curr) * dot(curr,next) < -0.10
fn refine_normals(segment: &mut BackboneSegment) {
    let len = segment.guide_points.len();
    if len < 3 {
        return;
    }

    // Compute tangent vectors for each position (direction along the chain)
    let tangents: Vec<Vec3> = (0..len)
        .map(|i| compute_local_tangent(segment, i))
        .collect();

    // Step 1: Make orientations orthogonal to tangents
    // Only operate on interior residues (not the end vectors)
    for i in 1..len - 1 {
        let tangent = tangents[i];
        let orient = segment.guide_points[i].orientation;

        // Remove tangent component: orient = orient - tangent*(orient·tangent)
        let tangent_component = tangent * orient.dot(tangent);
        let new_orient = orient - tangent_component;

        let mag = new_orient.magnitude();
        if mag > 1e-6 {
            segment.guide_points[i].orientation = new_orient / mag;
        }
    }

    // Step 2: Generate alternative orientations (original and inverted)
    // For helices, don't allow inversion (would confuse inside/outside)
    let mut alternatives: Vec<[Vec3; 2]> = vec![[Vec3::new(0.0, 1.0, 0.0); 2]; len];
    for i in 0..len {
        let orient = segment.guide_points[i].orientation;
        let ss_type = segment.guide_points[i].ss_type;
        
        // Original orientation
        alternatives[i][0] = orient;
        
        // Inverted orientation (only for non-helices)
        if is_helix(ss_type) {
            alternatives[i][1] = orient; // Keep same for helices
        } else {
            alternatives[i][1] = orient * -1.0;
        }
    }

    // Step 3: Forward iterate through pairs to select optimal orientation
    // Compare with previous orientation in the chain direction
    for i in 1..len - 1 {
        let tangent = tangents[i];
        
        // Get previous orientation (already selected, made perpendicular to chain)
        let prev_orient = segment.guide_points[i - 1].orientation;
        let mut prev_perp = prev_orient - tangent * prev_orient.dot(tangent);
        let prev_len = prev_perp.magnitude();
        if prev_len > 1e-6 {
            prev_perp = prev_perp / prev_len;
        }

        // Try both candidate orientations and select the one most aligned with previous
        let cand0 = alternatives[i][0];
        let cand1 = alternatives[i][1];
        
        // Make candidates perpendicular to tangent
        let mut cand0_perp = cand0 - tangent * cand0.dot(tangent);
        let mut cand1_perp = cand1 - tangent * cand1.dot(tangent);
        
        let len0 = cand0_perp.magnitude();
        let len1 = cand1_perp.magnitude();
        
        if len0 > 1e-6 {
            cand0_perp = cand0_perp / len0;
        }
        if len1 > 1e-6 {
            cand1_perp = cand1_perp / len1;
        }

        // Select the candidate with higher dot product with previous
        let dot0 = prev_perp.dot(cand0_perp);
        let dot1 = prev_perp.dot(cand1_perp);

        let best = if dot1 > dot0 { cand1 } else { cand0 };
        segment.guide_points[i].orientation = best;
    }

    // Step 4: Detect and soften kinks
    // A kink occurs when dot(prev, curr) * dot(curr, next) < -0.10
    // This is PyMOL's threshold value
    for i in 1..len - 1 {
        let tangent = tangents[i];
        let prev = segment.guide_points[i - 1].orientation;
        let curr = segment.guide_points[i].orientation;
        let next = segment.guide_points[i + 1].orientation;

        let dot_prev_curr = prev.dot(curr);
        let dot_curr_next = curr.dot(next);

        if dot_prev_curr * dot_curr_next < -0.10 {
            // Kink detected - soften by averaging with neighbors
            let mut avg = prev + next;
            // Add small amount of current to avoid division issues
            avg = avg + curr * 0.001;
            
            // Remove tangent component
            avg = avg - tangent * avg.dot(tangent);
            let avg_len = avg.magnitude();
            
            if avg_len > 1e-6 {
                let avg_norm = avg / avg_len;
                
                // If the averaged vector points opposite to current, adjust
                let mut smoothed = if curr.dot(avg_norm) < 0.0 {
                    curr - avg_norm
                } else {
                    curr + avg_norm
                };
                
                let smooth_len = smoothed.magnitude();
                if smooth_len > 1e-6 {
                    smoothed = smoothed / smooth_len;
                }

                // Calculate blend factor based on kink severity
                let dp = 2.0 * (-0.10 - dot_prev_curr * dot_curr_next);
                let blend = dp.min(1.0).max(0.0);

                // Mix current with smoothed based on kink severity
                let result = curr * (1.0 - blend) + smoothed * blend;
                let result_len = result.magnitude();
                if result_len > 1e-6 {
                    segment.guide_points[i].orientation = result / result_len;
                }
            }
        }
    }
}

// ============================================================================
// Helper functions for smoothing
// ============================================================================

/// Window-based position smoothing for a range of guide points
#[allow(dead_code)]
fn window_smooth_positions(
    segment: &mut BackboneSegment,
    start: usize,
    end: usize,
    cycles: u32,
    window_size: usize,
) {
    let len = segment.guide_points.len();
    if len < 3 || window_size == 0 {
        return;
    }

    for _ in 0..cycles {
        let mut new_positions = Vec::with_capacity(end - start + 1);

        for i in start..=end.min(len - 1) {
            // Calculate window boundaries (clamped to range)
            let win_start = i.saturating_sub(window_size).max(start);
            let win_end = (i + window_size).min(end).min(len - 1);
            let win_count = (win_end - win_start + 1) as f32;

            // Average positions in window
            let mut avg = Vec3::new(0.0, 0.0, 0.0);
            for j in win_start..=win_end {
                avg = avg + segment.guide_points[j].position;
            }
            avg = avg / win_count;

            // Weight center position more heavily
            let center_weight = 2.0;
            let total_weight = win_count - 1.0 + center_weight;
            let weighted_avg = (avg * (win_count - 1.0)
                + segment.guide_points[i].position * center_weight)
                / total_weight;

            new_positions.push(weighted_avg);
        }

        // Apply smoothed positions
        for (idx, pos) in new_positions.into_iter().enumerate() {
            segment.guide_points[start + idx].position = pos;
        }
    }
}

/// Window-based orientation smoothing for a range of guide points
#[allow(dead_code)]
fn window_smooth_orientations_range(
    segment: &mut BackboneSegment,
    start: usize,
    end: usize,
    cycles: u32,
    window_size: usize,
) {
    let len = segment.guide_points.len();
    if len < 3 || window_size == 0 {
        return;
    }

    for _ in 0..cycles {
        let mut new_orientations = Vec::with_capacity(end - start + 1);

        for i in start..=end.min(len - 1) {
            // Calculate window boundaries (clamped to range)
            let win_start = i.saturating_sub(window_size).max(start);
            let win_end = (i + window_size).min(end).min(len - 1);

            // Average orientations in window (with flip detection)
            let mut avg = Vec3::new(0.0, 0.0, 0.0);
            let reference = segment.guide_points[i].orientation;

            for j in win_start..=win_end {
                let orient = segment.guide_points[j].orientation;
                // Flip if pointing opposite to reference
                if orient.dot(reference) < 0.0 {
                    avg = avg - orient;
                } else {
                    avg = avg + orient;
                }
            }

            // Weight center orientation more heavily
            avg = avg + reference;

            let mag = avg.magnitude();
            if mag > 1e-6 {
                new_orientations.push(avg / mag);
            } else {
                new_orientations.push(reference);
            }
        }

        // Apply smoothed orientations
        for (idx, orient) in new_orientations.into_iter().enumerate() {
            segment.guide_points[start + idx].orientation = orient;
        }
    }
}

/// Check if a secondary structure type is a helix
fn is_helix(ss: SecondaryStructure) -> bool {
    matches!(
        ss,
        SecondaryStructure::Helix | SecondaryStructure::Helix310 | SecondaryStructure::HelixPi
    )
}

/// Detect secondary structure runs within a segment
///
/// Returns tuples of (start_index, end_index, ss_type) for each contiguous run
/// of secondary structure.
#[allow(dead_code)]
pub fn detect_ss_runs(segment: &BackboneSegment) -> Vec<(usize, usize, SecondaryStructure)> {
    if segment.is_empty() {
        return Vec::new();
    }

    let mut runs = Vec::new();
    let mut run_start = 0;
    let mut current_ss = segment.guide_points[0].ss_type;

    for (i, point) in segment.guide_points.iter().enumerate().skip(1) {
        if point.ss_type != current_ss {
            runs.push((run_start, i - 1, current_ss));
            run_start = i;
            current_ss = point.ss_type;
        }
    }

    // Push the last run
    runs.push((run_start, segment.guide_points.len() - 1, current_ss));

    runs
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_guide_point_creation() {
        let point = GuidePoint::new(
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            [1.0, 1.0, 1.0, 1.0],
            SecondaryStructure::Helix,
            AtomIndex(0),
            1,
            20.0,
        );

        assert_eq!(point.ss_type, SecondaryStructure::Helix);
        assert_eq!(point.resv, 1);
    }

    #[test]
    fn test_backbone_segment() {
        let mut segment = BackboneSegment::new("A");
        assert!(segment.is_empty());

        segment.push(GuidePoint::new(
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            [1.0, 1.0, 1.0, 1.0],
            SecondaryStructure::Loop,
            AtomIndex(0),
            1,
            20.0,
        ));

        assert_eq!(segment.len(), 1);
        assert!(!segment.is_empty());
    }

    #[test]
    fn test_detect_ss_runs() {
        let mut segment = BackboneSegment::new("A");

        // Add some guide points with different secondary structures
        for (i, ss) in [
            SecondaryStructure::Loop,
            SecondaryStructure::Helix,
            SecondaryStructure::Helix,
            SecondaryStructure::Helix,
            SecondaryStructure::Sheet,
            SecondaryStructure::Sheet,
            SecondaryStructure::Loop,
        ]
        .iter()
        .enumerate()
        {
            segment.push(GuidePoint::new(
                Vec3::new(i as f32, 0.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
                [1.0, 1.0, 1.0, 1.0],
                *ss,
                AtomIndex(i as u32),
                i as i32 + 1,
                20.0,
            ));
        }

        let runs = detect_ss_runs(&segment);

        assert_eq!(runs.len(), 4);
        assert_eq!(runs[0], (0, 0, SecondaryStructure::Loop));
        assert_eq!(runs[1], (1, 3, SecondaryStructure::Helix));
        assert_eq!(runs[2], (4, 5, SecondaryStructure::Sheet));
        assert_eq!(runs[3], (6, 6, SecondaryStructure::Loop));
    }

    #[test]
    fn test_cartoon_smooth_settings_default() {
        let settings = CartoonSmoothSettings::default();
        // These match PyMOL's defaults exactly
        assert_eq!(settings.smooth_cycles, 2);
        assert_eq!(settings.flat_cycles, 4);
        assert_eq!(settings.smooth_first, 1);
        assert_eq!(settings.smooth_last, 1);  // PyMOL default
        assert_eq!(settings.refine_cycles, 5);
        assert!(settings.refine_normals);
        assert!((settings.power_a - 2.0).abs() < 0.01);       // cartoon_power
        assert!((settings.power_b - 0.52).abs() < 0.01);      // cartoon_power_b
        assert!((settings.throw_factor - 1.35).abs() < 0.01); // cartoon_throw
    }

    #[test]
    fn test_is_helix() {
        assert!(is_helix(SecondaryStructure::Helix));
        assert!(is_helix(SecondaryStructure::Helix310));
        assert!(is_helix(SecondaryStructure::HelixPi));
        assert!(!is_helix(SecondaryStructure::Sheet));
        assert!(!is_helix(SecondaryStructure::Loop));
        assert!(!is_helix(SecondaryStructure::Turn));
    }

    #[test]
    fn test_compute_local_tangent() {
        let mut segment = BackboneSegment::new("A");
        
        // Create a simple linear segment along X axis
        for i in 0..5 {
            segment.push(GuidePoint::new(
                Vec3::new(i as f32 * 3.8, 0.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
                [1.0, 1.0, 1.0, 1.0],
                SecondaryStructure::Loop,
                AtomIndex(i),
                i as i32 + 1,
                20.0,
            ));
        }

        // Tangent at middle point should be along X
        let tangent = compute_local_tangent(&segment, 2);
        assert!(tangent.x > 0.9);
        assert!(tangent.y.abs() < 0.1);
        assert!(tangent.z.abs() < 0.1);
    }

    #[test]
    fn test_refine_normals() {
        let mut segment = BackboneSegment::new("A");
        
        // Create a segment with alternating flipped orientations
        for i in 0..5 {
            let flip = if i % 2 == 0 { 1.0 } else { -1.0 };
            segment.push(GuidePoint::new(
                Vec3::new(i as f32 * 3.8, 0.0, 0.0),
                Vec3::new(0.0, flip, 0.0),
                [1.0, 1.0, 1.0, 1.0],
                SecondaryStructure::Loop,
                AtomIndex(i),
                i as i32 + 1,
                20.0,
            ));
        }

        refine_normals(&mut segment);

        // After refinement, all normals should point in roughly the same direction
        let first_orient = segment.guide_points[0].orientation;
        for gp in &segment.guide_points {
            let dot = first_orient.dot(gp.orientation);
            assert!(dot > 0.0, "Orientations should not be flipped after refinement");
        }
    }

    #[test]
    fn test_apply_pymol_smoothing() {
        let mut segment = BackboneSegment::new("A");
        
        // Create a helix-like segment
        for i in 0..8 {
            // Simulate helical positions (slight twist)
            let angle = i as f32 * 100.0_f32.to_radians();
            let x = i as f32 * 1.5;
            let y = 2.3 * angle.cos();
            let z = 2.3 * angle.sin();
            
            segment.push(GuidePoint::new(
                Vec3::new(x, y, z),
                Vec3::new(0.0, angle.sin(), angle.cos()),
                [1.0, 0.0, 0.0, 1.0],
                SecondaryStructure::Helix,
                AtomIndex(i),
                i as i32 + 1,
                20.0,
            ));
        }

        let settings = CartoonSmoothSettings::default();
        apply_pymol_smoothing(&mut segment, &settings);

        // After smoothing, positions should be closer to the axis
        // Verify that the segment is still valid
        assert_eq!(segment.len(), 8);
        for gp in &segment.guide_points {
            // Orientations should be unit vectors
            let mag = gp.orientation.magnitude();
            assert!((mag - 1.0).abs() < 0.1, "Orientation should be normalized");
        }
    }
}
