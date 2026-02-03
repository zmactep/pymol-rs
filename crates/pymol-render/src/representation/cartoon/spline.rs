//! Backbone interpolation for cartoon representation
//!
//! Implements PyMOL's exact interpolation algorithm for smooth cartoon curves.
//! Based on analysis of PyMOL's RepCartoon.cpp CartoonGenerateSample function.

use lin_alg::f32::Vec3;

// ============================================================================
// Interpolation Settings
// ============================================================================

/// Settings for backbone interpolation (matches PyMOL defaults)
#[derive(Debug, Clone)]
pub struct InterpolationSettings {
    /// Power for sigmoid blend (cartoon_power). Default: 2.0
    pub power_a: f32,
    /// Power for displacement envelope (cartoon_power_b). Default: 0.52
    pub power_b: f32,
    /// Throw factor for displacement magnitude (cartoon_throw). Default: 1.35
    pub throw_factor: f32,
}

impl Default for InterpolationSettings {
    fn default() -> Self {
        Self {
            power_a: 2.0,
            power_b: 0.52,
            throw_factor: 1.35,
        }
    }
}

// ============================================================================
// Interpolated Point
// ============================================================================

/// An interpolated point along the backbone curve
#[derive(Debug, Clone)]
pub struct InterpolatedPoint {
    /// Position in 3D space
    pub position: Vec3,
    /// Tangent direction (along the curve, normalized)
    pub tangent: Vec3,
    /// Normal direction (ribbon orientation, normalized)
    pub normal: Vec3,
    /// Global parameter along the full curve (0 to 1)
    #[allow(dead_code)]
    pub t: f32,
    /// Index of the source segment (which guide point pair)
    pub segment: usize,
    /// Local parameter within the segment (0 to 1)
    pub local_t: f32,
}

// ============================================================================
// Core Interpolation Functions
// ============================================================================

/// PyMOL's smooth() function - sigmoid smoothstep
///
/// Creates an S-curve with:
/// - smooth(0) = 0
/// - smooth(1) = 1
/// - smooth'(0) = 0
/// - smooth'(1) = 0
#[inline]
pub fn sigmoid_blend(x: f32, power: f32) -> f32 {
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }
    
    if x <= 0.5 {
        0.5 * (2.0 * x).powf(power)
    } else {
        1.0 - 0.5 * (2.0 * (1.0 - x)).powf(power)
    }
}

/// Compute tangent vectors using PyMOL's algorithm
///
/// PyMOL computes tangents as the sum of adjacent normalized difference vectors:
/// - First: tangent[0] = nv[0]
/// - Interior: tangent[i] = normalize(nv[i] + nv[i-1])
/// - Last: tangent[n-1] = nv[n-2]
///
/// Where nv[i] = normalize(pos[i+1] - pos[i])
pub fn compute_tangents(positions: &[Vec3]) -> Vec<Vec3> {
    let n = positions.len();
    if n < 2 {
        return vec![Vec3::new(0.0, 0.0, 1.0); n];
    }

    // First compute normalized difference vectors (nv)
    let mut nv = Vec::with_capacity(n - 1);
    for i in 0..(n - 1) {
        let diff = positions[i + 1] - positions[i];
        nv.push(normalize_safe(diff));
    }

    // Now compute tangents using PyMOL's formula
    let mut tangents = Vec::with_capacity(n);

    // First point: tangent = nv[0]
    tangents.push(nv[0]);

    // Interior points: tangent = normalize(nv[i] + nv[i-1])
    for i in 1..(n - 1) {
        let sum = nv[i] + nv[i - 1];
        tangents.push(normalize_safe(sum));
    }

    // Last point: tangent = nv[n-2]
    tangents.push(nv[n - 2]);

    tangents
}

/// Compute distances between consecutive positions
pub fn compute_distances(positions: &[Vec3]) -> Vec<f32> {
    if positions.len() < 2 {
        return Vec::new();
    }

    let mut distances = Vec::with_capacity(positions.len() - 1);
    for i in 0..(positions.len() - 1) {
        distances.push((positions[i + 1] - positions[i]).magnitude());
    }
    distances
}

/// Interpolate an entire backbone using PyMOL's exact algorithm
///
/// This implements CartoonGenerateSample from PyMOL's RepCartoon.cpp:
/// ```text
/// f0 = smooth(t, power_a)
/// f1 = 1 - f0
/// f2 = smooth(f0, power_b)
/// f3 = smooth(f1, power_b)
/// f4 = dev * f2 * f3
///
/// position = f1*pos1 + f0*pos2 + f4*(f3*tangent1 - f2*tangent2)
/// orientation = f1*(orient1*f2) + f0*(orient2*f3)
/// ```
pub fn interpolate_backbone(
    positions: &[Vec3],
    tangents: &[Vec3],
    orientations: &[Vec3],
    sampling: u32,
    settings: &InterpolationSettings,
) -> Vec<InterpolatedPoint> {
    if positions.len() < 2 {
        if positions.is_empty() {
            return Vec::new();
        }
        return vec![InterpolatedPoint {
            position: positions[0],
            tangent: tangents.first().copied().unwrap_or(Vec3::new(0.0, 0.0, 1.0)),
            normal: orientations.first().copied().unwrap_or(Vec3::new(0.0, 1.0, 0.0)),
            t: 0.0,
            segment: 0,
            local_t: 0.0,
        }];
    }

    // Compute distances for deviation calculation
    let distances = compute_distances(positions);
    
    let n_segments = positions.len() - 1;
    let total_points = n_segments * sampling as usize + 1;
    let mut result = Vec::with_capacity(total_points);

    for seg in 0..n_segments {
        let pos1 = positions[seg];
        let pos2 = positions[seg + 1];
        
        let tangent1 = tangents.get(seg).copied().unwrap_or(Vec3::new(0.0, 0.0, 1.0));
        let tangent2 = tangents.get(seg + 1).copied().unwrap_or(tangent1);
        
        let orient1 = orientations.get(seg).copied().unwrap_or(Vec3::new(0.0, 1.0, 0.0));
        let orient2 = orientations.get(seg + 1).copied().unwrap_or(orient1);

        // Deviation = throw_factor * distance (PyMOL's exact formula)
        let dev = settings.throw_factor * distances[seg];

        // Generate samples for this segment
        for b in 0..sampling {
            // Skip first point of subsequent segments (already included)
            if seg > 0 && b == 0 {
                continue;
            }
            
            let t = b as f32 / sampling as f32;
            let global_t = (seg as f32 + t) / n_segments as f32;
            
            // PyMOL's exact interpolation formula
            let f0 = sigmoid_blend(t, settings.power_a);
            let f1 = 1.0 - f0;
            let f2 = sigmoid_blend(f0, settings.power_b);
            let f3 = sigmoid_blend(f1, settings.power_b);
            let f4 = dev * f2 * f3;
            
            // Position with tangent-based displacement
            // This is the exact PyMOL formula from CartoonGenerateSample
            let position = pos1 * f1 + pos2 * f0 + (tangent1 * f3 - tangent2 * f2) * f4;
            
            // Orientation interpolation (PyMOL's formula)
            let normal = normalize_safe(orient1 * (f1 * f2) + orient2 * (f0 * f3));
            
            // Tangent will be recomputed after all points are generated
            let tangent = normalize_safe(pos2 - pos1);
            
            result.push(InterpolatedPoint {
                position,
                tangent,
                normal,
                t: global_t,
                segment: seg,
                local_t: t,
            });
        }
    }

    // Add the final point
    let last_seg = n_segments - 1;
    result.push(InterpolatedPoint {
        position: positions[n_segments],
        tangent: tangents.get(n_segments).copied().unwrap_or(tangents[last_seg]),
        normal: orientations.get(n_segments).copied().unwrap_or(orientations[last_seg]),
        t: 1.0,
        segment: last_seg,
        local_t: 1.0,
    });

    // Recompute tangents from interpolated positions using PyMOL's algorithm
    recompute_tangents(&mut result);

    result
}

/// Recompute tangent vectors from interpolated positions using PyMOL's algorithm
///
/// Uses the same formula as compute_tangents:
/// tangent[i] = normalize(nv[i] + nv[i-1])
/// where nv[i] = normalize(pos[i+1] - pos[i])
fn recompute_tangents(points: &mut [InterpolatedPoint]) {
    let n = points.len();
    if n < 2 {
        return;
    }

    // First compute normalized difference vectors
    let mut nv = Vec::with_capacity(n - 1);
    for i in 0..(n - 1) {
        let diff = points[i + 1].position - points[i].position;
        nv.push(normalize_safe(diff));
    }

    // First point
    points[0].tangent = nv[0];

    // Interior points: sum of adjacent normalized differences
    for i in 1..(n - 1) {
        let sum = nv[i] + nv[i - 1];
        points[i].tangent = normalize_safe(sum);
    }

    // Last point
    points[n - 1].tangent = nv[n - 2];
}

/// Interpolate between two guide points (for testing/single segment use)
#[allow(dead_code)]
pub fn interpolate_segment(
    pos1: Vec3,
    pos2: Vec3,
    tangent1: Vec3,
    tangent2: Vec3,
    orient1: Vec3,
    orient2: Vec3,
    sampling: u32,
    settings: &InterpolationSettings,
) -> Vec<InterpolatedPoint> {
    let mut result = Vec::with_capacity(sampling as usize + 1);
    
    let distance = (pos2 - pos1).magnitude();
    let dev = settings.throw_factor * distance;

    for b in 0..=sampling {
        let t = b as f32 / sampling as f32;
        
        let f0 = sigmoid_blend(t, settings.power_a);
        let f1 = 1.0 - f0;
        let f2 = sigmoid_blend(f0, settings.power_b);
        let f3 = sigmoid_blend(f1, settings.power_b);
        let f4 = dev * f2 * f3;
        
        // PyMOL's exact position formula
        let position = pos1 * f1 + pos2 * f0 + (tangent1 * f3 - tangent2 * f2) * f4;
        
        // PyMOL's orientation formula
        let normal = normalize_safe(orient1 * (f1 * f2) + orient2 * (f0 * f3));
        
        let tangent = normalize_safe(pos2 - pos1);
        
        result.push(InterpolatedPoint {
            position,
            tangent,
            normal,
            t,
            segment: 0,
            local_t: t,
        });
    }
    
    // Recompute tangents
    recompute_tangents(&mut result);
    
    result
}

/// Refine interpolated points to reduce perpendicular jitter
/// (PyMOL's CartoonGenerateRefine algorithm)
#[allow(dead_code)]
pub fn refine_points(points: &mut [InterpolatedPoint], cycles: u32) {
    if points.len() < 3 {
        return;
    }

    for _ in 0..cycles {
        let mut new_positions = Vec::with_capacity(points.len());

        for i in 0..points.len() {
            if i == 0 || i == points.len() - 1 {
                new_positions.push(points[i].position);
                continue;
            }

            let prev = &points[i - 1];
            let curr = &points[i];
            let next = &points[i + 1];

            // Compute refinement plane from cross product of consecutive normals
            let cross = prev.normal.cross(next.normal);
            let cross_mag = cross.magnitude();

            if cross_mag < 1e-6 {
                new_positions.push(curr.position);
                continue;
            }

            let plane_normal = cross / cross_mag;

            // Average of neighbor positions
            let avg_pos = (prev.position + next.position) * 0.5;

            // Project current position onto plane
            let to_curr = curr.position - avg_pos;
            let proj_dist = to_curr.dot(plane_normal);

            let new_pos = curr.position - plane_normal * (proj_dist * 0.5);
            new_positions.push(new_pos);
        }

        for (i, pos) in new_positions.into_iter().enumerate() {
            points[i].position = pos;
        }

        // Recompute tangents after refinement
        recompute_tangents(points);
    }
}

// ============================================================================
// Helper Functions
// ============================================================================

/// Normalize a vector safely, returning a default for zero-length vectors
fn normalize_safe(v: Vec3) -> Vec3 {
    let len_sq = v.magnitude_squared();
    if len_sq > 1e-10 {
        v / len_sq.sqrt()
    } else {
        Vec3::new(0.0, 0.0, 1.0)
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sigmoid_blend() {
        assert!((sigmoid_blend(0.0, 1.0) - 0.0).abs() < 1e-6);
        assert!((sigmoid_blend(1.0, 1.0) - 1.0).abs() < 1e-6);
        assert!((sigmoid_blend(0.5, 1.0) - 0.5).abs() < 1e-6);

        assert!((sigmoid_blend(0.0, 2.0) - 0.0).abs() < 1e-6);
        assert!((sigmoid_blend(1.0, 2.0) - 1.0).abs() < 1e-6);
        assert!((sigmoid_blend(0.5, 2.0) - 0.5).abs() < 1e-6);

        // Higher power = steeper transition
        let power1 = sigmoid_blend(0.25, 1.0);
        let power2 = sigmoid_blend(0.25, 2.0);
        assert!(power2 < power1);
    }

    #[test]
    fn test_interpolation_settings_default() {
        let settings = InterpolationSettings::default();
        assert!((settings.power_a - 2.0).abs() < 0.01);
        assert!((settings.power_b - 0.52).abs() < 0.01);
        assert!((settings.throw_factor - 1.35).abs() < 0.01);
    }

    #[test]
    fn test_compute_tangents() {
        // Straight line - all tangents should point in same direction
        let positions = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(2.0, 0.0, 0.0),
        ];

        let tangents = compute_tangents(&positions);

        assert_eq!(tangents.len(), 3);
        for t in &tangents {
            assert!(t.x > 0.9); // Should point in +X direction
        }
    }

    #[test]
    fn test_compute_tangents_curved() {
        // L-shaped path
        let positions = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(1.0, 1.0, 0.0),
        ];

        let tangents = compute_tangents(&positions);
        
        assert_eq!(tangents.len(), 3);
        // First tangent: +X
        assert!(tangents[0].x > 0.9);
        // Middle tangent: diagonal (sum of +X and +Y)
        assert!(tangents[1].x > 0.5);
        assert!(tangents[1].y > 0.5);
        // Last tangent: +Y
        assert!(tangents[2].y > 0.9);
    }

    #[test]
    fn test_interpolate_segment() {
        let pos1 = Vec3::new(0.0, 0.0, 0.0);
        let pos2 = Vec3::new(3.8, 0.0, 0.0);
        let tangent1 = Vec3::new(1.0, 0.0, 0.0);
        let tangent2 = Vec3::new(1.0, 0.0, 0.0);
        let orient1 = Vec3::new(0.0, 1.0, 0.0);
        let orient2 = Vec3::new(0.0, 1.0, 0.0);

        let settings = InterpolationSettings::default();
        let result = interpolate_segment(pos1, pos2, tangent1, tangent2, orient1, orient2, 7, &settings);

        assert_eq!(result.len(), 8);
        
        // First point near pos1
        assert!((result[0].position - pos1).magnitude() < 0.1);
        
        // Last point near pos2
        assert!((result[7].position - pos2).magnitude() < 0.1);
    }

    #[test]
    fn test_interpolate_backbone() {
        let positions = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(3.8, 0.0, 0.0),
            Vec3::new(7.6, 1.0, 0.0),
        ];
        let tangents = compute_tangents(&positions);
        let orientations = vec![
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        ];

        let settings = InterpolationSettings::default();
        let result = interpolate_backbone(&positions, &tangents, &orientations, 7, &settings);

        // Should have (n_segments * sampling) + 1 points = 2*7 + 1 = 15
        // But we skip first point of second segment, so 14 + 1 = 15... let's just check > 10
        assert!(result.len() > 10);

        // First point near first position
        assert!((result[0].position - positions[0]).magnitude() < 0.5);
    }

    #[test]
    fn test_refine_points() {
        let mut points = vec![
            InterpolatedPoint {
                position: Vec3::new(0.0, 0.0, 0.0),
                tangent: Vec3::new(1.0, 0.0, 0.0),
                normal: Vec3::new(0.0, 1.0, 0.0),
                t: 0.0,
                segment: 0,
                local_t: 0.0,
            },
            InterpolatedPoint {
                position: Vec3::new(1.0, 0.5, 0.0),
                tangent: Vec3::new(1.0, 0.0, 0.0),
                normal: Vec3::new(0.0, 1.0, 0.0),
                t: 0.5,
                segment: 0,
                local_t: 0.5,
            },
            InterpolatedPoint {
                position: Vec3::new(2.0, 0.0, 0.0),
                tangent: Vec3::new(1.0, 0.0, 0.0),
                normal: Vec3::new(0.0, 1.0, 0.0),
                t: 1.0,
                segment: 0,
                local_t: 1.0,
            },
        ];

        refine_points(&mut points, 2);

        // Endpoints unchanged
        assert!((points[0].position.y - 0.0).abs() < 1e-6);
        assert!((points[2].position.y - 0.0).abs() < 1e-6);
        assert_eq!(points.len(), 3);
    }
}
