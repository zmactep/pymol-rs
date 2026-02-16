//! Reference frame calculation for cartoon representation
//!
//! Calculates local coordinate frames along the backbone curve for
//! orienting ribbon cross-sections.

use lin_alg::f32::Vec3;

use super::backbone::GuidePoint;
use super::spline::{compute_tangents, interpolate_backbone, InterpolationSettings};
use super::utils::normalize_safe;

/// A reference frame at a point along the cartoon backbone
///
/// The frame consists of three orthogonal vectors:
/// - Tangent (T): Along the direction of the curve
/// - Normal (N): Perpendicular to tangent, in the ribbon plane (toward O)
/// - Binormal (B): T × N, perpendicular to the ribbon
#[derive(Debug, Clone, Copy)]
pub struct ReferenceFrame {
    /// Position in world space
    pub position: Vec3,
    /// Tangent vector (along curve direction, normalized)
    pub tangent: Vec3,
    /// Normal vector (ribbon width direction, normalized)
    pub normal: Vec3,
    /// Binormal vector (ribbon thickness direction, normalized)
    pub binormal: Vec3,
}

impl ReferenceFrame {
    /// Create a new reference frame from position, tangent, and orientation
    ///
    /// Computes an orthonormal frame where:
    /// - tangent: points along the backbone direction
    /// - normal: points in the ribbon width direction
    /// - binormal: perpendicular to both (ribbon thickness direction)
    pub fn new(position: Vec3, tangent: Vec3, orientation: Vec3) -> Self {
        let t = normalize_safe(tangent);

        // Binormal = tangent × orientation (perpendicular to both)
        let b = t.cross(orientation);
        let b = normalize_safe(b);

        // Normal = binormal × tangent (ensures orthogonality)
        let n = b.cross(t);

        Self {
            position,
            tangent: t,
            normal: n,
            binormal: b,
        }
    }

    /// Transform a local 2D point to world space
    ///
    /// The local point is in the normal-binormal plane, centered at position.
    #[inline]
    pub fn transform_local(&self, local: (f32, f32)) -> Vec3 {
        self.position + self.normal * local.0 + self.binormal * local.1
    }

    /// Get the normal for a vertex given local coordinates
    ///
    /// For lighting calculations on extruded geometry.
    #[inline]
    pub fn local_normal(&self, local: (f32, f32)) -> Vec3 {
        let n = self.normal * local.0 + self.binormal * local.1;
        normalize_safe(n)
    }
}

/// Reference frame with associated metadata
#[derive(Debug, Clone)]
pub struct FrameWithMetadata {
    /// The reference frame
    pub frame: ReferenceFrame,
    /// Interpolated color
    pub color: [f32; 4],
    /// Secondary structure type
    pub ss_type: pymol_mol::SecondaryStructure,
    /// Source segment index (which guide point pair this frame was interpolated from)
    pub segment_idx: usize,
}

/// Generate reference frames along a backbone segment
///
/// Uses displacement-based interpolation to create smooth curves through
/// the guide points, then builds orthonormal frames for geometry extrusion.
pub fn generate_frames(
    guide_points: &[GuidePoint],
    settings: &InterpolationSettings,
    subdivisions: u32,
    smooth_cycles: u32,
) -> Vec<FrameWithMetadata> {
    if guide_points.is_empty() {
        return Vec::new();
    }

    if guide_points.len() == 1 {
        let gp = &guide_points[0];
        return vec![FrameWithMetadata {
            frame: ReferenceFrame::new(gp.position, Vec3::new(0.0, 0.0, 1.0), gp.orientation),
            color: gp.color,
            ss_type: gp.ss_type,
            segment_idx: 0,
        }];
    }

    // Extract data from guide points
    let positions: Vec<Vec3> = guide_points.iter().map(|gp| gp.position).collect();
    let orientations: Vec<Vec3> = guide_points.iter().map(|gp| gp.orientation).collect();

    // Ensure orientation consistency (no sudden flips)
    let consistent_orientations = ensure_orientation_consistency(&orientations);

    // Compute tangent vectors at each guide point
    let tangents = compute_tangents(&positions);

    // Interpolate using displacement-based algorithm
    let interp_points = interpolate_backbone(
        &positions,
        &tangents,
        &consistent_orientations,
        subdivisions,
        settings,
    );

    // Build frames from interpolated points
    let mut frames: Vec<FrameWithMetadata> = interp_points
        .iter()
        .map(|ip| {
            // Get metadata from source guide points
            let seg = ip.segment.min(guide_points.len() - 2);
            let gp1 = &guide_points[seg];
            let gp2 = &guide_points[seg + 1];

            // Interpolate color
            let color = interpolate_color_linear(&gp1.color, &gp2.color, ip.local_t);

            // Use SS type from the nearest guide point (discrete assignment)
            let ss_type = if ip.local_t < 0.5 { gp1.ss_type } else { gp2.ss_type };

            // Create frame from interpolated point
            // The interpolation already provides position, tangent, and normal
            let frame = ReferenceFrame::new(ip.position, ip.tangent, ip.normal);

            FrameWithMetadata {
                frame,
                color,
                ss_type,
                segment_idx: ip.segment,
            }
        })
        .collect();

    // Override orientations for sheet frames to keep sheets flat.
    // Compute a single averaged orientation per contiguous sheet run,
    // then apply it to all frames in that run. This prevents twist caused
    // by per-residue orientation variation.
    {
        // Find contiguous sheet runs in guide points
        let mut sheet_runs: Vec<(usize, usize)> = Vec::new();
        let mut in_sheet = false;
        let mut run_start = 0;
        for (i, gp) in guide_points.iter().enumerate() {
            if gp.ss_type == pymol_mol::SecondaryStructure::Sheet {
                if !in_sheet {
                    run_start = i;
                    in_sheet = true;
                }
            } else if in_sheet {
                sheet_runs.push((run_start, i - 1));
                in_sheet = false;
            }
        }
        if in_sheet {
            sheet_runs.push((run_start, guide_points.len() - 1));
        }

        // Compute averaged orientation for each sheet run
        let mut avg_orientations: Vec<(usize, usize, Vec3)> = Vec::new();
        for &(start, end) in &sheet_runs {
            let mut sum = Vec3::new(0.0, 0.0, 0.0);
            for i in start..=end {
                sum = sum + consistent_orientations[i];
            }
            let avg = normalize_safe(sum);
            avg_orientations.push((start, end, avg));
        }

        // Apply averaged orientation to all sheet frames
        for frame in &mut frames {
            if frame.ss_type == pymol_mol::SecondaryStructure::Sheet {
                let seg = frame.segment_idx.min(guide_points.len() - 2);
                // Find which sheet run this frame belongs to
                if let Some(&(_, _, avg_orient)) = avg_orientations.iter().find(|&&(s, e, _)| {
                    seg >= s && seg <= e
                }) {
                    frame.frame = ReferenceFrame::new(
                        frame.frame.position,
                        frame.frame.tangent,
                        avg_orient,
                    );
                }
            }
        }
    }

    // Apply smoothing to prevent normal flips (sheets are already skipped)
    smooth_frames(&mut frames, smooth_cycles);

    frames
}

/// Ensure orientation vectors are consistent (no sudden flips)
fn ensure_orientation_consistency(orientations: &[Vec3]) -> Vec<Vec3> {
    if orientations.is_empty() {
        return Vec::new();
    }

    let mut consistent = Vec::with_capacity(orientations.len());
    consistent.push(orientations[0]);

    for i in 1..orientations.len() {
        let prev = consistent[i - 1];
        let curr = orientations[i];

        // If current orientation points opposite to previous, flip it
        if prev.dot(curr) < 0.0 {
            consistent.push(curr * -1.0);
        } else {
            consistent.push(curr);
        }
    }

    consistent
}

/// Smooth reference frames while preserving consistent ribbon orientation
///
/// Uses parallel transport with blending to ensure smooth, continuous frames.
/// IMPORTANT: Sheets are excluded from smoothing - PyMOL computes sheet frames
/// independently at each point from the orientation vector to prevent twist.
fn smooth_frames(frames: &mut [FrameWithMetadata], cycles: u32) {
    use pymol_mol::SecondaryStructure;

    if frames.len() < 2 {
        return;
    }

    // Store original binormals (from orientation-based construction)
    let original_binormals: Vec<Vec3> = frames.iter().map(|f| f.frame.binormal).collect();

    // First pass: use parallel transport to establish continuous binormals
    // SKIP sheets - they should keep their orientation-based frames (PyMOL behavior)
    for i in 1..frames.len() {
        // Skip sheets - they need independent frame computation to stay flat
        if frames[i].ss_type == SecondaryStructure::Sheet {
            continue;
        }

        let tangent = frames[i].frame.tangent;
        let prev_binormal = frames[i - 1].frame.binormal;
        let orig_binormal = original_binormals[i];

        // Parallel transport: project previous binormal to current tangent's perpendicular plane
        let transported = prev_binormal - tangent * prev_binormal.dot(tangent);
        let transport_len = transported.magnitude();

        if transport_len > 1e-6 {
            let transported_norm = transported / transport_len;

            // Blend between transported and original based on agreement
            let agreement = transported_norm.dot(orig_binormal);

            let new_binormal = if agreement > 0.5 {
                // Good agreement - blend towards original
                let blend = agreement;
                normalize_safe(transported_norm * (1.0 - blend * 0.5) + orig_binormal * (blend * 0.5))
            } else if agreement < -0.5 {
                // Opposite - use transported for continuity
                normalize_safe(transported_norm)
            } else {
                // Poor agreement - use transported
                transported_norm
            };

            // Compute normal from binormal: N = B × T
            let new_normal = new_binormal.cross(tangent);

            frames[i].frame.binormal = new_binormal;
            frames[i].frame.normal = new_normal;
        }
    }

    // Apply smoothing cycles
    // SKIP sheets - they should not be smoothed (PyMOL behavior)
    for _ in 0..cycles {
        let len = frames.len();
        let mut new_binormals = Vec::with_capacity(len);

        // Smooth binormals by averaging with neighbors
        for i in 0..len {
            // Skip sheets - preserve their orientation-based binormals
            if frames[i].ss_type == SecondaryStructure::Sheet {
                new_binormals.push(frames[i].frame.binormal);
                continue;
            }

            let prev = if i > 0 { i - 1 } else { 0 };
            let next = if i < len - 1 { i + 1 } else { len - 1 };

            let b_prev = frames[prev].frame.binormal;
            let b_curr = frames[i].frame.binormal;
            let b_next = frames[next].frame.binormal;

            // Weighted average (center has more weight)
            let avg = b_prev + b_curr * 2.0 + b_next;
            new_binormals.push(normalize_safe(avg));
        }

        // Apply smoothed binormals and recompute frames
        for (i, new_binormal) in new_binormals.into_iter().enumerate() {
            // Skip sheets
            if frames[i].ss_type == SecondaryStructure::Sheet {
                continue;
            }

            let tangent = frames[i].frame.tangent;

            // Make binormal perpendicular to tangent
            let projected = new_binormal - tangent * new_binormal.dot(tangent);
            let binormal = normalize_safe(projected);

            // Compute normal: N = B × T
            let normal = binormal.cross(tangent);

            frames[i].frame.binormal = binormal;
            frames[i].frame.normal = normal;
        }
    }
}

/// Linear interpolation between two colors
fn interpolate_color_linear(c1: &[f32; 4], c2: &[f32; 4], t: f32) -> [f32; 4] {
    [
        c1[0] + (c2[0] - c1[0]) * t,
        c1[1] + (c2[1] - c1[1]) * t,
        c1[2] + (c2[2] - c1[2]) * t,
        c1[3] + (c2[3] - c1[3]) * t,
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use pymol_mol::{AtomIndex, SecondaryStructure};

    #[test]
    fn test_reference_frame() {
        let frame = ReferenceFrame::new(
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        );

        // Check orthogonality
        assert!(frame.tangent.dot(frame.normal).abs() < 1e-6);
        assert!(frame.tangent.dot(frame.binormal).abs() < 1e-6);
        assert!(frame.normal.dot(frame.binormal).abs() < 1e-6);
    }

    #[test]
    fn test_transform_local() {
        let frame = ReferenceFrame::new(
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        );

        let local = (1.0, 0.0);
        let world = frame.transform_local(local);

        assert!((world.y - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_generate_frames() {
        let guide_points = vec![
            GuidePoint::new(
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
                [1.0, 0.0, 0.0, 1.0],
                SecondaryStructure::Helix,
                AtomIndex(0),
                1
            ),
            GuidePoint::new(
                Vec3::new(3.8, 0.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
                [0.0, 1.0, 0.0, 1.0],
                SecondaryStructure::Helix,
                AtomIndex(1),
                2
            ),
            GuidePoint::new(
                Vec3::new(7.6, 0.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
                [0.0, 0.0, 1.0, 1.0],
                SecondaryStructure::Helix,
                AtomIndex(2),
                3
            ),
        ];

        let settings = InterpolationSettings::default();
        let frames = generate_frames(&guide_points, &settings, 7, 2);

        // Should have interpolated frames
        assert!(frames.len() > guide_points.len());

        // First position should be near first guide point
        assert!((frames[0].frame.position - guide_points[0].position).magnitude() < 0.1);
    }

    #[test]
    fn test_ensure_orientation_consistency() {
        let orientations = vec![
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, -1.0, 0.0), // Flipped
            Vec3::new(0.0, 1.0, 0.0),
        ];

        let consistent = ensure_orientation_consistency(&orientations);

        // All should now point in same general direction
        for i in 1..consistent.len() {
            assert!(consistent[i - 1].dot(consistent[i]) > 0.0);
        }
    }
}
