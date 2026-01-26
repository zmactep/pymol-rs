//! Geometry generation for cartoon representation
//!
//! Provides cross-section profiles and mesh extrusion for different
//! secondary structure types (helices, sheets, loops).

use pymol_mol::SecondaryStructure;
use pymol_settings::SettingResolver;

use super::frame::FrameWithMetadata;
use crate::vertex::MeshVertex;

/// Cross-section profile for cartoon extrusion
///
/// A profile is a 2D shape that gets extruded along the backbone spline.
/// The coordinates are in the normal-binormal plane of the reference frame.
#[derive(Debug, Clone)]
pub struct Profile {
    /// Points defining the profile (in local 2D coordinates)
    /// The profile is assumed to be closed (last point connects to first)
    pub points: Vec<(f32, f32)>,
    /// Outward-facing normals for each point (for lighting)
    pub normals: Vec<(f32, f32)>,
}

impl Profile {
    /// Create an elliptical profile (for helices)
    ///
    /// - `width`: Half-width in the normal direction (ribbon width)
    /// - `height`: Half-height in the binormal direction (ribbon thickness)
    /// - `quality`: Number of points around the ellipse
    pub fn ellipse(width: f32, height: f32, quality: u32) -> Self {
        let n = quality.max(32) as usize; // Minimum 32 for smooth, round appearance
        let mut points = Vec::with_capacity(n);
        let mut normals = Vec::with_capacity(n);

        for i in 0..n {
            let angle = 2.0 * std::f32::consts::PI * (i as f32) / (n as f32);
            let cos_a = angle.cos();
            let sin_a = angle.sin();

            // Point on ellipse
            let x = width * cos_a;
            let y = height * sin_a;
            points.push((x, y));

            // Normal to ellipse (perpendicular to tangent)
            // For ellipse: normal = (x/a², y/b²) normalized
            let nx = cos_a / width;
            let ny = sin_a / height;
            let len = (nx * nx + ny * ny).sqrt();
            normals.push((nx / len, ny / len));
        }

        Profile { points, normals }
    }

    /// Create a circular profile (for loops/tubes)
    ///
    /// - `radius`: Tube radius
    /// - `quality`: Number of points around the circle
    pub fn circle(radius: f32, quality: u32) -> Self {
        Self::ellipse(radius, radius, quality)
    }

    /// Create a rectangular profile (for sheets) with specified quality
    ///
    /// - `width`: Half-width in the normal direction (sheet width)
    /// - `height`: Half-height in the binormal direction (sheet thickness)
    /// - `quality`: Number of points around the rectangle (must match circle/ellipse for smooth transitions)
    ///
    /// Uses a superellipse formula to create a rounded rectangle shape that
    /// transitions smoothly to/from circular profiles.
    pub fn rectangle(width: f32, height: f32, quality: u32) -> Self {
        let n = quality.max(32) as usize; // Minimum 32 for smooth appearance
        let mut points = Vec::with_capacity(n);
        let mut normals = Vec::with_capacity(n);

        // Use superellipse with high exponent for rectangular shape
        // Higher exponent = more rectangular, lower = more elliptical
        let exponent = 4.0_f32; // Gives a nice rounded rectangle

        for i in 0..n {
            let angle = 2.0 * std::f32::consts::PI * (i as f32) / (n as f32);
            let cos_a = angle.cos();
            let sin_a = angle.sin();

            // Superellipse: |x/a|^n + |y/b|^n = 1
            // Parametric form with consistent vertex distribution
            let sign_x = cos_a.signum();
            let sign_y = sin_a.signum();
            
            let abs_cos = cos_a.abs();
            let abs_sin = sin_a.abs();
            
            // Superellipse parametric coordinates
            let x = width * sign_x * abs_cos.powf(2.0 / exponent);
            let y = height * sign_y * abs_sin.powf(2.0 / exponent);
            points.push((x, y));

            // Normal to superellipse (gradient of implicit function)
            // For |x/a|^n + |y/b|^n = 1, gradient is:
            // (n * sign(x) * |x/a|^(n-1) / a, n * sign(y) * |y/b|^(n-1) / b)
            let nx = if abs_cos > 1e-6 {
                sign_x * abs_cos.powf((exponent - 2.0) / exponent) / width
            } else {
                0.0
            };
            let ny = if abs_sin > 1e-6 {
                sign_y * abs_sin.powf((exponent - 2.0) / exponent) / height
            } else {
                0.0
            };
            
            let len = (nx * nx + ny * ny).sqrt();
            if len > 1e-6 {
                normals.push((nx / len, ny / len));
            } else {
                // Fallback: use direction from center
                let len2 = (x * x + y * y).sqrt();
                if len2 > 1e-6 {
                    normals.push((x / len2, y / len2));
                } else {
                    normals.push((1.0, 0.0));
                }
            }
        }

        Profile { points, normals }
    }

    /// Scale the profile uniformly
    #[allow(dead_code)]
    pub fn scale(&self, factor: f32) -> Self {
        Profile {
            points: self.points.iter().map(|(x, y)| (x * factor, y * factor)).collect(),
            normals: self.normals.clone(),
        }
    }

    /// Get the number of points in the profile
    pub fn len(&self) -> usize {
        self.points.len()
    }

    /// Check if the profile is empty
    #[allow(dead_code)]
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }
    
    /// Blend two profiles together
    ///
    /// Creates a new profile by linearly interpolating between `self` and `other`.
    /// Both profiles must have the same number of points.
    ///
    /// - `t`: Blend factor (0.0 = self, 1.0 = other)
    pub fn blend(&self, other: &Profile, t: f32) -> Self {
        debug_assert_eq!(self.points.len(), other.points.len(), "Profiles must have same vertex count");
        
        let t = t.clamp(0.0, 1.0);
        let one_minus_t = 1.0 - t;
        
        let points: Vec<(f32, f32)> = self.points.iter()
            .zip(other.points.iter())
            .map(|((x1, y1), (x2, y2))| {
                (x1 * one_minus_t + x2 * t, y1 * one_minus_t + y2 * t)
            })
            .collect();
        
        let normals: Vec<(f32, f32)> = self.normals.iter()
            .zip(other.normals.iter())
            .map(|((nx1, ny1), (nx2, ny2))| {
                let nx = nx1 * one_minus_t + nx2 * t;
                let ny = ny1 * one_minus_t + ny2 * t;
                let len = (nx * nx + ny * ny).sqrt();
                if len > 1e-6 {
                    (nx / len, ny / len)
                } else {
                    (*nx1, *ny1)
                }
            })
            .collect();
        
        Profile { points, normals }
    }
}

/// Settings for cartoon geometry generation
#[derive(Debug, Clone)]
pub struct CartoonGeometrySettings {
    /// Helix ribbon width (cartoon_oval_width)
    pub helix_width: f32,
    /// Helix ribbon length/height (cartoon_oval_length)
    pub helix_height: f32,
    /// Sheet arrow width (cartoon_rect_width)
    pub sheet_width: f32,
    /// Sheet arrow height/thickness (cartoon_rect_length)
    pub sheet_height: f32,
    /// Loop tube radius (cartoon_loop_radius)
    pub loop_radius: f32,
    /// Quality for curved profiles
    pub quality: u32,
    /// Whether to use round helices
    pub round_helices: bool,
    /// Whether to use flat sheets (affects sheet geometry)
    #[allow(dead_code)]
    pub flat_sheets: bool,
    /// Whether to draw arrows on beta sheets (cartoon_fancy_sheets)
    pub fancy_sheets: bool,
    /// Arrow tip scale factor (how much wider the arrow gets before tapering)
    pub arrow_tip_scale: f32,
    /// Number of frames for arrow taper region (calculated from subdivisions)
    pub arrow_length: usize,
    /// Number of residues the arrow should span
    pub arrow_residues: usize,
    /// If true, use uniform circular tube profile for ALL secondary structures (ribbon mode)
    pub uniform_tube: bool,
}

impl CartoonGeometrySettings {
    /// Create settings from the settings resolver
    pub fn from_resolver(settings: &SettingResolver) -> Self {
        // Setting IDs from pymol-settings definitions
        const CARTOON_OVAL_WIDTH: u16 = 101;
        const CARTOON_OVAL_LENGTH: u16 = 100;
        const CARTOON_RECT_WIDTH: u16 = 97;
        const CARTOON_RECT_LENGTH: u16 = 96;
        const CARTOON_LOOP_RADIUS: u16 = 92;
        const CARTOON_OVAL_QUALITY: u16 = 102;
        const CARTOON_ROUND_HELICES: u16 = 111;
        const CARTOON_FLAT_SHEETS: u16 = 113;
        const CARTOON_FANCY_SHEETS: u16 = 119;

        // Quality: PyMOL uses -1 for "auto", so we need to handle negative values
        // Higher quality = smoother curves. Using 32 for smooth, round helices
        let quality_raw = settings.get_int_if_defined(CARTOON_OVAL_QUALITY).unwrap_or(32);
        let quality = if quality_raw < 0 { 32u32 } else { (quality_raw as u32).max(32) };

        Self {
            helix_width: settings.get_float_if_defined(CARTOON_OVAL_WIDTH).unwrap_or(0.25),
            helix_height: settings.get_float_if_defined(CARTOON_OVAL_LENGTH).unwrap_or(1.35),
            sheet_width: settings.get_float_if_defined(CARTOON_RECT_WIDTH).unwrap_or(0.4),
            sheet_height: settings.get_float_if_defined(CARTOON_RECT_LENGTH).unwrap_or(1.4),
            loop_radius: settings.get_float_if_defined(CARTOON_LOOP_RADIUS).unwrap_or(0.2),
            quality,
            round_helices: settings.get_bool_if_defined(CARTOON_ROUND_HELICES).unwrap_or(true),
            flat_sheets: settings.get_bool_if_defined(CARTOON_FLAT_SHEETS).unwrap_or(true),
            fancy_sheets: settings.get_bool_if_defined(CARTOON_FANCY_SHEETS).unwrap_or(true),
            arrow_tip_scale: 1.5,
            arrow_length: 0, // Will be calculated from subdivisions
            arrow_residues: 2, // Arrow spans ~2 residues (like PyMOL's sampling parameter)
            uniform_tube: false, // Default: different profiles for different SS types
        }
    }
    
    /// Set the arrow length based on spline subdivisions
    ///
    /// The arrow should span approximately `arrow_residues` residues.
    /// Each residue has `subdivisions + 1` frames, so:
    /// arrow_length = arrow_residues * (subdivisions + 1)
    pub fn with_subdivisions(mut self, subdivisions: u32) -> Self {
        self.arrow_length = self.arrow_residues * (subdivisions as usize + 1);
        self
    }
}

impl Default for CartoonGeometrySettings {
    fn default() -> Self {
        Self {
            helix_width: 0.25,
            helix_height: 1.35,
            sheet_width: 0.4,
            sheet_height: 1.4,
            loop_radius: 0.2,
            quality: 32, // High quality for smooth, round helices
            round_helices: true,
            flat_sheets: true,
            fancy_sheets: true,
            arrow_tip_scale: 1.5,
            arrow_length: 16, // Default for subdivisions=7: 2 * (7+1) = 16
            arrow_residues: 2,
            uniform_tube: false,
        }
    }
}

/// Generate cartoon mesh from reference frames
pub fn generate_cartoon_mesh(
    frames: &[FrameWithMetadata],
    settings: &CartoonGeometrySettings,
    sheet_termini: &[(usize, usize)], // (start, end) indices of sheet runs for arrow heads
) -> (Vec<MeshVertex>, Vec<u32>) {
    if frames.is_empty() {
        return (Vec::new(), Vec::new());
    }

    let mut vertices = Vec::new();
    let mut indices = Vec::new();

    // Track previous profile for connecting segments
    let mut prev_ring_start: Option<usize> = None;
    let mut prev_profile_len: usize = 0;
    
    // Detect SS transition boundaries for blending
    let transition_frames = detect_ss_transitions(frames);
    // Number of frames to blend across at each transition
    // Higher value = smoother transitions between secondary structure types
    let blend_range = 8usize;

    for (i, frame_meta) in frames.iter().enumerate() {
        // Get the base profile for this frame's SS type
        let base_profile = select_profile(frame_meta, settings, i, frames.len(), sheet_termini);
        
        // Check if we need to blend with a different profile
        let profile = apply_transition_blending(
            i,
            frames,
            &base_profile,
            settings,
            &transition_frames,
            blend_range,
            sheet_termini,
        );
        
        let ring_start = vertices.len();

        // Generate vertices for this profile ring
        for j in 0..profile.len() {
            let local_pos = profile.points[j];
            let local_normal = profile.normals[j];

            let world_pos = frame_meta.frame.transform_local(local_pos);
            
            // Transform the 2D profile normal to world space using the frame
            // This is the correct approach for ellipses: the profile stores the
            // mathematically correct normal (gradient of ellipse equation), and
            // we transform it through the frame's basis vectors.
            // 
            // For an ellipse (x/a)² + (y/b)² = 1, the normal at point (a*cos(θ), b*sin(θ))
            // is proportional to (cos(θ)/a, sin(θ)/b), NOT the radial direction.
            let world_normal = frame_meta.frame.local_normal(local_normal);

            vertices.push(MeshVertex {
                position: [world_pos.x, world_pos.y, world_pos.z],
                normal: [world_normal.x, world_normal.y, world_normal.z],
                color: frame_meta.color,
            });
        }

        // Connect to previous ring if profiles have same vertex count
        if let Some(prev_start) = prev_ring_start {
            if profile.len() == prev_profile_len {
                connect_rings(
                    &mut indices,
                    prev_start as u32,
                    ring_start as u32,
                    profile.len() as u32,
                );
            }
        }

        prev_ring_start = Some(ring_start);
        prev_profile_len = profile.len();
    }

    // Cap the ends
    if frames.len() >= 2 {
        // Cap the start (never an arrow tip)
        cap_tube_end(&mut vertices, &mut indices, &frames[0], settings, true, false);
        
        // Check if the end is an arrow tip (at a sheet terminus)
        let last_idx = frames.len() - 1;
        let is_arrow_tip = is_at_arrow_tip(last_idx, sheet_termini);
        cap_tube_end(
            &mut vertices,
            &mut indices,
            &frames[last_idx],
            settings,
            false,
            is_arrow_tip,
        );
    }

    (vertices, indices)
}

/// Select the appropriate profile for a frame based on secondary structure
fn select_profile(
    frame_meta: &FrameWithMetadata,
    settings: &CartoonGeometrySettings,
    frame_idx: usize,
    _total_frames: usize,
    sheet_termini: &[(usize, usize)],
) -> Profile {
    // Ribbon mode: use uniform circular tube for ALL secondary structures
    if settings.uniform_tube {
        return Profile::circle(settings.loop_radius, settings.quality);
    }
    
    match frame_meta.ss_type {
        SecondaryStructure::Helix
        | SecondaryStructure::Helix310
        | SecondaryStructure::HelixPi => {
            // Note: helix_height is the ribbon WIDTH (in normal direction)
            //       helix_width is the ribbon THICKNESS (in binormal direction)
            // The ellipse first param goes to normal (x), second to binormal (y)
            if settings.round_helices {
                Profile::ellipse(settings.helix_height, settings.helix_width, settings.quality)
            } else {
                Profile::rectangle(settings.helix_height, settings.helix_width, settings.quality)
            }
        }
        SecondaryStructure::Sheet => {
            // Check if this is near a sheet terminus for arrow tapering
            let scale = calculate_arrow_scale(
                frame_idx,
                sheet_termini,
                settings.arrow_length,
                settings.arrow_tip_scale,
            );
            
            // sheet_height is the ribbon WIDTH (in normal direction)
            // sheet_width is the ribbon THICKNESS (in binormal direction)
            // Scale affects the width (normal direction) for arrow tapering
            let ribbon_width = settings.sheet_height * scale;
            let ribbon_thickness = settings.sheet_width * 0.5; // Sheets are flatter
            
            if ribbon_width < 0.01 {
                // At or near the arrow tip - create minimal profile to avoid artifacts
                Profile::rectangle(0.01, ribbon_thickness, settings.quality)
            } else {
                Profile::rectangle(ribbon_width, ribbon_thickness, settings.quality)
            }
        }
        _ => {
            // Loop/coil/turn/bend - use tube
            Profile::circle(settings.loop_radius, settings.quality)
        }
    }
}

/// Calculate arrow width scale factor for sheet termini (PyMOL-style)
///
/// Returns:
/// - 1.0 for normal sheet regions (no arrow)
/// - 1.5 at the start of the arrow (widest point)
/// - 0.0 at the arrow tip (point)
/// - Values between 1.5 and 0.0 for the taper region
fn calculate_arrow_scale(
    frame_idx: usize,
    sheet_termini: &[(usize, usize)],
    arrow_length: usize,
    tip_scale: f32,
) -> f32 {
    if arrow_length == 0 {
        return 1.0;
    }
    
    for &(_start, end) in sheet_termini {
        // Check if near the end of a sheet (C-terminus gets arrow)
        if frame_idx >= end.saturating_sub(arrow_length) && frame_idx <= end {
            let distance_from_end = end - frame_idx;
            // PyMOL formula: tip_scale * distance_from_end / arrow_length
            // At distance_from_end = arrow_length: scale = tip_scale (widest)
            // At distance_from_end = 0 (tip): scale = 0 (point)
            let scale = tip_scale * (distance_from_end as f32) / (arrow_length as f32);
            return scale;
        }
    }
    1.0 // Not in an arrow region, use normal width
}

/// Check if a frame index is at an arrow tip (end of a sheet segment)
fn is_at_arrow_tip(frame_idx: usize, sheet_termini: &[(usize, usize)]) -> bool {
    sheet_termini.iter().any(|&(_start, end)| frame_idx == end)
}

/// Detect frame indices where SS type transitions occur
fn detect_ss_transitions(frames: &[FrameWithMetadata]) -> Vec<usize> {
    let mut transitions = Vec::new();
    
    for i in 1..frames.len() {
        if frames[i].ss_type != frames[i - 1].ss_type {
            transitions.push(i);
        }
    }
    
    transitions
}

/// Apply smooth profile blending at SS transitions
///
/// When near a transition boundary, blend between the current profile and
/// the adjacent SS type's profile for smooth visual transitions.
fn apply_transition_blending(
    frame_idx: usize,
    frames: &[FrameWithMetadata],
    base_profile: &Profile,
    settings: &CartoonGeometrySettings,
    transitions: &[usize],
    blend_range: usize,
    sheet_termini: &[(usize, usize)],
) -> Profile {
    if blend_range == 0 || transitions.is_empty() {
        return base_profile.clone();
    }
    
    // Find the nearest transition
    let mut nearest_transition: Option<usize> = None;
    let mut nearest_distance = usize::MAX;
    
    for &trans_idx in transitions {
        let distance = if frame_idx >= trans_idx {
            frame_idx - trans_idx
        } else {
            trans_idx - frame_idx
        };
        
        if distance < nearest_distance && distance <= blend_range {
            nearest_distance = distance;
            nearest_transition = Some(trans_idx);
        }
    }
    
    let trans_idx = match nearest_transition {
        Some(t) => t,
        None => return base_profile.clone(),
    };
    
    // Determine which SS types are involved in the transition
    let (before_ss, after_ss) = if trans_idx > 0 && trans_idx < frames.len() {
        (frames[trans_idx - 1].ss_type, frames[trans_idx].ss_type)
    } else {
        return base_profile.clone();
    };
    
    // Get the "other" profile (the one we're transitioning to/from)
    let other_ss = if frame_idx < trans_idx { after_ss } else { before_ss };
    
    // Create profile for the other SS type
    let other_profile = create_profile_for_ss(other_ss, settings, frame_idx, frames.len(), sheet_termini);
    
    // Calculate blend factor based on distance from transition
    // At transition: blend = 0.5
    // At blend_range distance: blend = 0.0 (before) or 1.0 (after)
    let blend_factor = if frame_idx < trans_idx {
        // Before transition: blend from base (0) towards other (1)
        let dist = trans_idx - frame_idx;
        1.0 - (dist as f32 / (blend_range as f32 + 1.0))
    } else {
        // After transition: blend from other towards base
        let dist = frame_idx - trans_idx;
        dist as f32 / (blend_range as f32 + 1.0)
    };
    
    // Clamp blend factor
    let blend_factor = blend_factor.clamp(0.0, 1.0);
    
    // If blend factor is very close to 0 or 1, skip blending
    if blend_factor < 0.01 {
        return base_profile.clone();
    }
    if blend_factor > 0.99 {
        return base_profile.clone();
    }
    
    // Blend the profiles
    if frame_idx < trans_idx {
        // Before transition: base -> other
        base_profile.blend(&other_profile, blend_factor)
    } else {
        // After transition: other -> base
        other_profile.blend(base_profile, blend_factor)
    }
}

/// Create a profile for a given SS type (helper for blending)
fn create_profile_for_ss(
    ss_type: SecondaryStructure,
    settings: &CartoonGeometrySettings,
    frame_idx: usize,
    _total_frames: usize,
    sheet_termini: &[(usize, usize)],
) -> Profile {
    // Ribbon mode: use uniform circular tube for ALL secondary structures
    if settings.uniform_tube {
        return Profile::circle(settings.loop_radius, settings.quality);
    }
    
    match ss_type {
        SecondaryStructure::Helix
        | SecondaryStructure::Helix310
        | SecondaryStructure::HelixPi => {
            // helix_height = ribbon width (normal), helix_width = ribbon thickness (binormal)
            if settings.round_helices {
                Profile::ellipse(settings.helix_height, settings.helix_width, settings.quality)
            } else {
                Profile::rectangle(settings.helix_height, settings.helix_width, settings.quality)
            }
        }
        SecondaryStructure::Sheet => {
            let scale = calculate_arrow_scale(
                frame_idx,
                sheet_termini,
                settings.arrow_length,
                settings.arrow_tip_scale,
            );
            // sheet_height = ribbon width (normal), sheet_width = ribbon thickness (binormal)
            let ribbon_width = settings.sheet_height * scale;
            let ribbon_thickness = settings.sheet_width * 0.5;
            if ribbon_width < 0.01 {
                Profile::rectangle(0.01, ribbon_thickness, settings.quality)
            } else {
                Profile::rectangle(ribbon_width, ribbon_thickness, settings.quality)
            }
        }
        _ => Profile::circle(settings.loop_radius, settings.quality),
    }
}

/// Connect two profile rings with triangles
fn connect_rings(indices: &mut Vec<u32>, ring1_start: u32, ring2_start: u32, profile_len: u32) {
    for j in 0..profile_len {
        let j_next = (j + 1) % profile_len;

        let v00 = ring1_start + j;
        let v01 = ring1_start + j_next;
        let v10 = ring2_start + j;
        let v11 = ring2_start + j_next;

        // Two triangles per quad - CCW winding when viewed from outside
        // (reversed from previous CW winding to match outward-pointing vertex normals)
        indices.push(v00);
        indices.push(v01);
        indices.push(v10);

        indices.push(v01);
        indices.push(v11);
        indices.push(v10);
    }
}

/// Add end cap to the tube
///
/// - `is_start`: true for the start of a segment, false for the end
/// - `is_arrow_tip`: true if this is at an arrow tip (skip cap since geometry tapers to a point)
fn cap_tube_end(
    vertices: &mut Vec<MeshVertex>,
    indices: &mut Vec<u32>,
    frame_meta: &FrameWithMetadata,
    settings: &CartoonGeometrySettings,
    is_start: bool,
    is_arrow_tip: bool,
) {
    // Skip cap for arrow tips - the geometry tapers to a point
    if is_arrow_tip && !is_start {
        return;
    }

    // Ribbon mode: use uniform circular tube for ALL secondary structures
    let profile = if settings.uniform_tube {
        Profile::circle(settings.loop_radius, settings.quality)
    } else {
        // helix_height/sheet_height = ribbon width (normal direction)
        // helix_width/sheet_width = ribbon thickness (binormal direction)
        match frame_meta.ss_type {
            SecondaryStructure::Helix
            | SecondaryStructure::Helix310
            | SecondaryStructure::HelixPi => {
                if settings.round_helices {
                    Profile::ellipse(settings.helix_height, settings.helix_width, settings.quality)
                } else {
                    Profile::rectangle(settings.helix_height, settings.helix_width, settings.quality)
                }
            }
            SecondaryStructure::Sheet => {
                Profile::rectangle(settings.sheet_height, settings.sheet_width * 0.5, settings.quality)
            }
            _ => Profile::circle(settings.loop_radius, settings.quality),
        }
    };

    let normal_dir = if is_start { -1.0 } else { 1.0 };
    let cap_normal = frame_meta.frame.tangent * normal_dir;

    // Add center vertex
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

    // Add edge vertices
    let edge_start = vertices.len() as u32;
    for point in &profile.points {
        let world_pos = frame_meta.frame.transform_local(*point);
        vertices.push(MeshVertex {
            position: [world_pos.x, world_pos.y, world_pos.z],
            normal: [cap_normal.x, cap_normal.y, cap_normal.z],
            color: frame_meta.color,
        });
    }

    // Create fan triangles
    let n = profile.len() as u32;
    for j in 0..n {
        let j_next = (j + 1) % n;
        if is_start {
            // Reverse winding for start cap
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

/// Find sheet termini indices in a sequence of frames
pub fn find_sheet_termini(frames: &[FrameWithMetadata]) -> Vec<(usize, usize)> {
    let mut termini = Vec::new();
    let mut in_sheet = false;
    let mut sheet_start = 0;

    for (i, frame) in frames.iter().enumerate() {
        let is_sheet = frame.ss_type == SecondaryStructure::Sheet;

        if is_sheet && !in_sheet {
            // Sheet starts
            sheet_start = i;
            in_sheet = true;
        } else if !is_sheet && in_sheet {
            // Sheet ends
            termini.push((sheet_start, i - 1));
            in_sheet = false;
        }
    }

    // Handle sheet at end
    if in_sheet {
        termini.push((sheet_start, frames.len() - 1));
    }

    termini
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::frame::ReferenceFrame;
    use lin_alg::f32::Vec3;

    #[test]
    fn test_ellipse_profile() {
        let profile = Profile::ellipse(1.0, 0.5, 32);
        assert_eq!(profile.len(), 32);

        // Check that points lie on ellipse
        for (x, y) in &profile.points {
            let dist = (x / 1.0).powi(2) + (y / 0.5).powi(2);
            assert!((dist - 1.0).abs() < 0.01);
        }
    }

    #[test]
    fn test_circle_profile() {
        let profile = Profile::circle(0.5, 32);
        assert_eq!(profile.len(), 32);

        // Check that all points are equidistant from center
        for (x, y) in &profile.points {
            let dist = (x * x + y * y).sqrt();
            assert!((dist - 0.5).abs() < 0.01);
        }
    }

    #[test]
    fn test_rectangle_profile() {
        let profile = Profile::rectangle(1.0, 0.5, 32);
        assert_eq!(profile.len(), 32); // Same vertex count as ellipse

        // Check normals are unit vectors
        for (nx, ny) in &profile.normals {
            let len = (nx * nx + ny * ny).sqrt();
            assert!((len - 1.0).abs() < 0.01);
        }
        
        // Verify rectangle and ellipse have same vertex count for smooth transitions
        let ellipse = Profile::ellipse(1.0, 0.5, 32);
        assert_eq!(profile.len(), ellipse.len());
    }

    #[test]
    fn test_profile_scale() {
        let profile = Profile::circle(1.0, 8);
        let scaled = profile.scale(2.0);

        assert_eq!(profile.len(), scaled.len());

        for i in 0..profile.len() {
            assert!((scaled.points[i].0 - profile.points[i].0 * 2.0).abs() < 1e-6);
            assert!((scaled.points[i].1 - profile.points[i].1 * 2.0).abs() < 1e-6);
        }
    }

    #[test]
    fn test_find_sheet_termini() {
        let make_frame = |ss| FrameWithMetadata {
            frame: ReferenceFrame::new(
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
            ),
            color: [1.0, 1.0, 1.0, 1.0],
            ss_type: ss,
            b_factor: 20.0,
            segment_idx: 0,
            local_t: 0.0,
        };

        let frames = vec![
            make_frame(SecondaryStructure::Loop),
            make_frame(SecondaryStructure::Sheet),
            make_frame(SecondaryStructure::Sheet),
            make_frame(SecondaryStructure::Sheet),
            make_frame(SecondaryStructure::Loop),
            make_frame(SecondaryStructure::Sheet),
            make_frame(SecondaryStructure::Sheet),
            make_frame(SecondaryStructure::Loop),
        ];

        let termini = find_sheet_termini(&frames);
        assert_eq!(termini.len(), 2);
        assert_eq!(termini[0], (1, 3));
        assert_eq!(termini[1], (5, 6));
    }

    #[test]
    fn test_arrow_scale_calculation() {
        // Sheet termini: frames 5-10 (end at 10)
        // With arrow_length=3, arrow region is frames 7-10
        let sheet_termini = vec![(5, 10)];
        let arrow_length = 3;
        let tip_scale = 1.5;

        // Outside arrow region - should return 1.0
        assert!((calculate_arrow_scale(0, &sheet_termini, arrow_length, tip_scale) - 1.0).abs() < 0.01);
        assert!((calculate_arrow_scale(5, &sheet_termini, arrow_length, tip_scale) - 1.0).abs() < 0.01);
        assert!((calculate_arrow_scale(6, &sheet_termini, arrow_length, tip_scale) - 1.0).abs() < 0.01);
        assert!((calculate_arrow_scale(11, &sheet_termini, arrow_length, tip_scale) - 1.0).abs() < 0.01);

        // At the arrow tip (frame_idx == end) - should return 0.0
        // distance_from_end = 10 - 10 = 0, scale = 1.5 * 0/3 = 0.0
        assert!((calculate_arrow_scale(10, &sheet_termini, arrow_length, tip_scale) - 0.0).abs() < 0.01);

        // At the start of arrow region (distance_from_end == arrow_length)
        // Frame 7: distance_from_end = 10 - 7 = 3, scale = 1.5 * 3/3 = 1.5
        assert!((calculate_arrow_scale(7, &sheet_termini, arrow_length, tip_scale) - 1.5).abs() < 0.01);
        
        // Frame 8: distance_from_end = 10 - 8 = 2, scale = 1.5 * 2/3 = 1.0
        assert!((calculate_arrow_scale(8, &sheet_termini, arrow_length, tip_scale) - 1.0).abs() < 0.01);
        
        // Frame 9: distance_from_end = 10 - 9 = 1, scale = 1.5 * 1/3 = 0.5
        assert!((calculate_arrow_scale(9, &sheet_termini, arrow_length, tip_scale) - 0.5).abs() < 0.01);
    }

    #[test]
    fn test_is_at_arrow_tip() {
        let sheet_termini = vec![(1, 3), (5, 6)];

        // At arrow tips
        assert!(is_at_arrow_tip(3, &sheet_termini));
        assert!(is_at_arrow_tip(6, &sheet_termini));

        // Not at arrow tips
        assert!(!is_at_arrow_tip(0, &sheet_termini));
        assert!(!is_at_arrow_tip(1, &sheet_termini));
        assert!(!is_at_arrow_tip(2, &sheet_termini));
        assert!(!is_at_arrow_tip(4, &sheet_termini));
        assert!(!is_at_arrow_tip(5, &sheet_termini));
        assert!(!is_at_arrow_tip(7, &sheet_termini));
    }
    
    #[test]
    fn test_profile_blend() {
        // Create two profiles with the same vertex count
        let circle = Profile::circle(1.0, 12);
        let rect = Profile::rectangle(1.0, 0.5, 12);
        
        assert_eq!(circle.len(), rect.len());
        
        // Blend at t=0 should give circle
        let blend_0 = circle.blend(&rect, 0.0);
        assert_eq!(blend_0.len(), circle.len());
        for i in 0..circle.len() {
            assert!((blend_0.points[i].0 - circle.points[i].0).abs() < 1e-5);
            assert!((blend_0.points[i].1 - circle.points[i].1).abs() < 1e-5);
        }
        
        // Blend at t=1 should give rect
        let blend_1 = circle.blend(&rect, 1.0);
        for i in 0..rect.len() {
            assert!((blend_1.points[i].0 - rect.points[i].0).abs() < 1e-5);
            assert!((blend_1.points[i].1 - rect.points[i].1).abs() < 1e-5);
        }
        
        // Blend at t=0.5 should be midway
        let blend_05 = circle.blend(&rect, 0.5);
        for i in 0..circle.len() {
            let expected_x = (circle.points[i].0 + rect.points[i].0) / 2.0;
            let expected_y = (circle.points[i].1 + rect.points[i].1) / 2.0;
            assert!((blend_05.points[i].0 - expected_x).abs() < 1e-5);
            assert!((blend_05.points[i].1 - expected_y).abs() < 1e-5);
        }
    }
    
    #[test]
    fn test_detect_ss_transitions() {
        let make_frame = |ss| FrameWithMetadata {
            frame: ReferenceFrame::new(
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
            ),
            color: [1.0, 1.0, 1.0, 1.0],
            ss_type: ss,
            b_factor: 20.0,
            segment_idx: 0,
            local_t: 0.0,
        };

        let frames = vec![
            make_frame(SecondaryStructure::Loop),
            make_frame(SecondaryStructure::Loop),
            make_frame(SecondaryStructure::Helix), // Transition at index 2
            make_frame(SecondaryStructure::Helix),
            make_frame(SecondaryStructure::Sheet), // Transition at index 4
            make_frame(SecondaryStructure::Sheet),
            make_frame(SecondaryStructure::Loop),  // Transition at index 6
        ];

        let transitions = detect_ss_transitions(&frames);
        assert_eq!(transitions.len(), 3);
        assert_eq!(transitions[0], 2);
        assert_eq!(transitions[1], 4);
        assert_eq!(transitions[2], 6);
    }
}
