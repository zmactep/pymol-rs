//! Reference frame calculation for cartoon representation
//!
//! Calculates local coordinate frames along the backbone spline for
//! orienting ribbon cross-sections.

use lin_alg::f32::Vec3;

use super::backbone::GuidePoint;
use super::spline::{interpolate_with_tangents, CatmullRomSpline, InterpolatedPoint};

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
    /// Create a new reference frame using PyMOL's get_system2f3f algorithm
    ///
    /// This is equivalent to PyMOL's `get_system2f3f` function:
    /// 1. Takes tangent (along curve) and orientation (from CA->O) as input
    /// 2. Computes binormal = tangent × orientation (then normalize)
    /// 3. Recomputes normal = binormal × tangent (ensures orthogonality)
    ///
    /// The result is an orthonormal frame where:
    /// - tangent: points along the backbone direction
    /// - normal: points in the ribbon width direction (perpendicular to tangent, in the plane defined by tangent and orientation)
    /// - binormal: perpendicular to both (ribbon thickness direction)
    pub fn new(position: Vec3, tangent: Vec3, orientation: Vec3) -> Self {
        let t = normalize_safe(tangent);
        
        // Binormal = tangent × orientation (perpendicular to both)
        let b = t.cross(orientation);
        let b = normalize_safe(b);
        
        // Normal = binormal × tangent (ensures orthogonality, lies in tangent-orientation plane)
        let n = b.cross(t);
        // n is already normalized since b and t are orthonormal

        Self {
            position,
            tangent: t,
            normal: n,
            binormal: b,
        }
    }

    /// Create a frame with explicit binormal (for pre-computed frames)
    pub fn from_tnb(position: Vec3, tangent: Vec3, normal: Vec3, binormal: Vec3) -> Self {
        Self {
            position,
            tangent: normalize_safe(tangent),
            normal: normalize_safe(normal),
            binormal: normalize_safe(binormal),
        }
    }
    
    /// Build frame using PyMOL's get_system1f approach (for tubes/loops)
    ///
    /// This is for circular cross-sections where orientation doesn't matter
    /// as much - it uses parallel transport to minimize twist.
    #[allow(dead_code)]
    pub fn new_tube(position: Vec3, tangent: Vec3, prev_binormal: Option<Vec3>) -> Self {
        let t = normalize_safe(tangent);
        
        // Use previous binormal if available (parallel transport)
        // Otherwise, find an arbitrary perpendicular
        let b = match prev_binormal {
            Some(prev_b) => {
                // Project previous binormal to be perpendicular to new tangent
                let proj = prev_b - t * prev_b.dot(t);
                let len = proj.magnitude();
                if len > 1e-6 {
                    proj / len
                } else {
                    find_perpendicular(t)
                }
            }
            None => find_perpendicular(t),
        };
        
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

/// Normalize a vector, returning a default if zero-length
#[inline]
fn normalize_safe(v: Vec3) -> Vec3 {
    let len_sq = v.magnitude_squared();
    if len_sq > 1e-10 {
        v / len_sq.sqrt()
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    }
}

/// Generate reference frames along a backbone segment
///
/// This interpolates positions and orientations, then calculates proper
/// orthonormal frames for geometry extrusion.
pub fn generate_frames(
    guide_points: &[GuidePoint],
    spline: &CatmullRomSpline,
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
            b_factor: gp.b_factor,
            segment_idx: 0,
            local_t: 0.0,
        }];
    }

    // Extract positions for spline interpolation
    let positions: Vec<Vec3> = guide_points.iter().map(|gp| gp.position).collect();

    // Interpolate positions with tangents
    let interp_points = interpolate_with_tangents(spline, &positions, subdivisions);

    // Interpolate orientations
    let orientations: Vec<Vec3> = guide_points.iter().map(|gp| gp.orientation).collect();
    let interp_orientations = interpolate_orientations(&orientations, &interp_points);

    // Build initial frames
    let mut frames: Vec<FrameWithMetadata> = interp_points
        .iter()
        .zip(interp_orientations.iter())
        .map(|(ip, orientation)| {
            // Get metadata from source guide points
            let seg = ip.segment.min(guide_points.len() - 2);
            let gp1 = &guide_points[seg];
            let gp2 = &guide_points[seg + 1];

            // Interpolate color
            let color = interpolate_color_linear(&gp1.color, &gp2.color, ip.local_t);

            // Interpolate B-factor
            let b_factor = gp1.b_factor + (gp2.b_factor - gp1.b_factor) * ip.local_t;

            // Use SS type from the nearest guide point (discrete assignment)
            // This ensures proper sheet termini detection and smooth visual transitions
            let ss_type = if ip.local_t < 0.5 { gp1.ss_type } else { gp2.ss_type };

            // Create frame
            let frame = ReferenceFrame::new(ip.position, ip.tangent, *orientation);

            FrameWithMetadata {
                frame,
                color,
                ss_type,
                b_factor,
                segment_idx: ip.segment,
                local_t: ip.local_t,
            }
        })
        .collect();

    // Apply smoothing to prevent normal flips
    smooth_frames(&mut frames, smooth_cycles);

    frames
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
    /// Interpolated B-factor (for putty representation)
    #[allow(dead_code)]
    pub b_factor: f32,
    /// Source segment index
    #[allow(dead_code)]
    pub segment_idx: usize,
    /// Local parameter within segment
    #[allow(dead_code)]
    pub local_t: f32,
}

/// Interpolate orientation vectors along the spline
/// 
/// First ensures orientation continuity (no sudden flips), then interpolates.
fn interpolate_orientations(orientations: &[Vec3], interp_points: &[InterpolatedPoint]) -> Vec<Vec3> {
    if orientations.is_empty() {
        return vec![Vec3::new(0.0, 1.0, 0.0); interp_points.len()];
    }

    if orientations.len() == 1 {
        return vec![orientations[0]; interp_points.len()];
    }

    // First pass: ensure orientation continuity by flipping vectors that point
    // opposite to their predecessor. This prevents sudden 180° rotations.
    let mut consistent_orientations = Vec::with_capacity(orientations.len());
    consistent_orientations.push(orientations[0]);
    
    for i in 1..orientations.len() {
        let prev = consistent_orientations[i - 1];
        let curr = orientations[i];
        
        // If current orientation points opposite to previous, flip it
        if prev.dot(curr) < 0.0 {
            consistent_orientations.push(curr * -1.0);
        } else {
            consistent_orientations.push(curr);
        }
    }

    // Now interpolate using the consistent orientations
    interp_points
        .iter()
        .map(|ip| {
            let seg = ip.segment.min(consistent_orientations.len() - 2);
            let o1 = consistent_orientations[seg];
            let o2 = consistent_orientations[seg + 1];

            // Spherical linear interpolation (slerp) for better orientation blending
            slerp_vec3(o1, o2, ip.local_t)
        })
        .collect()
}

/// Spherical linear interpolation between two direction vectors
/// 
/// Takes the shortest path between vectors - if they point in opposite directions,
/// v2 is negated to ensure continuity.
fn slerp_vec3(v1: Vec3, v2: Vec3, t: f32) -> Vec3 {
    let v1 = normalize_safe(v1);
    let mut v2 = normalize_safe(v2);

    let mut dot = v1.dot(v2);

    // If vectors point in opposite directions, negate v2 to take the shortest path
    // This is crucial for maintaining orientation continuity along the backbone
    if dot < 0.0 {
        v2 = v2 * -1.0;
        dot = -dot;
    }

    // Clamp to avoid numerical issues
    dot = dot.min(1.0);

    // If vectors are nearly parallel, use linear interpolation
    if dot > 0.9995 {
        let result = v1 + (v2 - v1) * t;
        return normalize_safe(result);
    }

    let theta = dot.acos();
    let sin_theta = theta.sin();

    if sin_theta.abs() < 1e-6 {
        return v1;
    }

    let a = ((1.0 - t) * theta).sin() / sin_theta;
    let b = (t * theta).sin() / sin_theta;

    normalize_safe(v1 * a + v2 * b)
}

/// Smooth reference frames while preserving consistent ribbon orientation
///
/// This function ensures smooth, continuous frames by:
/// 1. Propagating a consistent binormal direction using parallel transport
/// 2. Blending with the original orientation-based binormal to preserve ribbon direction
/// 3. Smoothing to reduce local variations
///
/// The approach uses parallel transport (projecting previous binormal forward)
/// but blends it with the original orientation-based frame to maintain the
/// intended ribbon direction from the CA→O vectors.
fn smooth_frames(frames: &mut [FrameWithMetadata], cycles: u32) {
    if frames.len() < 2 {
        return;
    }

    // Store original binormals (from orientation-based construction)
    let original_binormals: Vec<Vec3> = frames.iter().map(|f| f.frame.binormal).collect();

    // First pass: use parallel transport to establish continuous binormals
    for i in 1..frames.len() {
        let tangent = frames[i].frame.tangent;
        let prev_binormal = frames[i - 1].frame.binormal;
        let orig_binormal = original_binormals[i];
        
        // Parallel transport: project previous binormal to current tangent's perpendicular plane
        let transported = prev_binormal - tangent * prev_binormal.dot(tangent);
        let transport_len = transported.magnitude();
        
        if transport_len > 1e-6 {
            let transported_norm = transported / transport_len;
            
            // Blend between transported and original based on how much they agree
            // If they point the same way, use more of the original (preserves orientation)
            // If they point opposite, use the transported (ensures continuity)
            let agreement = transported_norm.dot(orig_binormal);
            
            let new_binormal = if agreement > 0.5 {
                // Good agreement - blend towards original
                let blend = agreement; // 0.5 to 1.0 → more original
                normalize_safe(transported_norm * (1.0 - blend * 0.5) + orig_binormal * (blend * 0.5))
            } else if agreement < -0.5 {
                // Opposite - transported is flipped, use negated original
                normalize_safe(transported_norm)
            } else {
                // Poor agreement - use transported for continuity
                transported_norm
            };
            
            // Compute normal from binormal: N = B × T
            let new_normal = new_binormal.cross(tangent);
            
            frames[i].frame.binormal = new_binormal;
            frames[i].frame.normal = new_normal;
        }
    }

    // Apply smoothing cycles
    for _ in 0..cycles {
        let len = frames.len();
        let mut new_binormals = Vec::with_capacity(len);

        // Smooth binormals by averaging with neighbors
        for i in 0..len {
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

/// Calculate frames using the parallel transport method
///
/// This is an alternative to using O-atom orientations that produces
/// frames with minimal twist along the curve.
#[allow(dead_code)]
pub fn parallel_transport_frames(
    positions: &[Vec3],
    subdivisions: u32,
    spline: &CatmullRomSpline,
) -> Vec<ReferenceFrame> {
    let interp_points = interpolate_with_tangents(spline, positions, subdivisions);

    if interp_points.is_empty() {
        return Vec::new();
    }

    let mut frames = Vec::with_capacity(interp_points.len());

    // Initialize first frame
    let first_tangent = interp_points[0].tangent;
    let initial_normal = find_perpendicular(first_tangent);
    let initial_binormal = first_tangent.cross(initial_normal);

    frames.push(ReferenceFrame::from_tnb(
        interp_points[0].position,
        first_tangent,
        initial_normal,
        initial_binormal,
    ));

    // Propagate frames using parallel transport
    for i in 1..interp_points.len() {
        let prev_frame = &frames[i - 1];
        let curr_tangent = interp_points[i].tangent;

        // Rotate previous normal to be perpendicular to new tangent
        let axis = prev_frame.tangent.cross(curr_tangent);
        let axis_len_sq = axis.magnitude_squared();

        let new_normal = if axis_len_sq > 1e-10 {
            let axis = axis / axis_len_sq.sqrt();
            let cos_angle = prev_frame.tangent.dot(curr_tangent).max(-1.0).min(1.0);
            let angle = cos_angle.acos();

            rotate_around_axis(prev_frame.normal, axis, angle)
        } else {
            prev_frame.normal
        };

        let new_binormal = curr_tangent.cross(new_normal);

        frames.push(ReferenceFrame::from_tnb(
            interp_points[i].position,
            curr_tangent,
            new_normal,
            new_binormal,
        ));
    }

    frames
}

/// Find a vector perpendicular to the given vector
#[allow(dead_code)]
fn find_perpendicular(v: Vec3) -> Vec3 {
    let abs_x = v.x.abs();
    let abs_y = v.y.abs();
    let abs_z = v.z.abs();

    let other = if abs_x <= abs_y && abs_x <= abs_z {
        Vec3::new(1.0, 0.0, 0.0)
    } else if abs_y <= abs_z {
        Vec3::new(0.0, 1.0, 0.0)
    } else {
        Vec3::new(0.0, 0.0, 1.0)
    };

    let perp = v.cross(other);
    normalize_safe(perp)
}

/// Rotate a vector around an axis by an angle (radians)
#[allow(dead_code)]
fn rotate_around_axis(v: Vec3, axis: Vec3, angle: f32) -> Vec3 {
    let cos_a = angle.cos();
    let sin_a = angle.sin();

    // Rodrigues' rotation formula
    let term1 = v * cos_a;
    let term2 = axis.cross(v) * sin_a;
    let term3 = axis * axis.dot(v) * (1.0 - cos_a);

    normalize_safe(term1 + term2 + term3)
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
    fn test_slerp() {
        let v1 = Vec3::new(1.0, 0.0, 0.0);
        let v2 = Vec3::new(0.0, 1.0, 0.0);

        // At t=0, should be v1
        let at_0 = slerp_vec3(v1, v2, 0.0);
        assert!((at_0 - v1).magnitude_squared() < 1e-6);

        // At t=1, should be v2
        let at_1 = slerp_vec3(v1, v2, 1.0);
        assert!((at_1 - v2).magnitude_squared() < 1e-6);

        // At t=0.5, should be normalized average direction
        let at_05 = slerp_vec3(v1, v2, 0.5);
        assert!((at_05.magnitude_squared() - 1.0).abs() < 1e-6);
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
                1,
                20.0,
            ),
            GuidePoint::new(
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
                [0.0, 1.0, 0.0, 1.0],
                SecondaryStructure::Helix,
                AtomIndex(1),
                2,
                25.0,
            ),
            GuidePoint::new(
                Vec3::new(2.0, 0.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
                [0.0, 0.0, 1.0, 1.0],
                SecondaryStructure::Helix,
                AtomIndex(2),
                3,
                30.0,
            ),
        ];

        let spline = CatmullRomSpline::new();
        let frames = generate_frames(&guide_points, &spline, 2, 1);

        // Should have interpolated frames
        assert!(frames.len() > guide_points.len());

        // First and last positions should match guide points
        assert!((frames[0].frame.position - guide_points[0].position).magnitude_squared() < 1e-6);
    }
}
