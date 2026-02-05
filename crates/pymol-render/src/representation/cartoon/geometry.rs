//! Geometry generation for cartoon representation
//!
//! Provides cross-section profiles and mesh extrusion for different
//! secondary structure types (helices, sheets, loops).

use lin_alg::f32::Vec3;
use pymol_mol::SecondaryStructure;
use pymol_settings::SettingResolver;

use super::frame::FrameWithMetadata;
use crate::vertex::MeshVertex;

/// Profile rendering type - determines how mesh is generated
/// 
/// Round profiles use `connect_rings()` which wraps around all vertices.
/// Flat profiles use `generate_face_strips()` which renders each face pair separately.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProfileType {
    /// Round/closed profiles (circle, ellipse, rectangle) - use connect_rings()
    Round,
    /// Flat profiles with 4 faces (sheets) - 8 vertices, face-pair rendering
    Flat4Face,
    /// Flat profiles with 2 faces (dumbbells) - 4 vertices, top/bottom only
    Flat2Face,
}

/// Normalize a vector safely, returning a default for zero-length vectors
#[inline]
fn normalize_safe(v: Vec3) -> Vec3 {
    let len_sq = v.magnitude_squared();
    if len_sq > 1e-10 {
        v / len_sq.sqrt()
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    }
}

/// PyMOL's smooth easing function for tapering
///
/// Creates an ease-in/ease-out curve for smooth transitions at helix ends.
/// Matches PyMOL's `smooth(x, power)` function from Vector.cpp.
///
/// # Arguments
/// * `x` - Input value in range [0, 1]
/// * `power` - Power factor (typically 2.0 for quadratic easing)
///
/// # Returns
/// Smoothed value in range [0, 1]
#[inline]
fn smooth(x: f32, power: f32) -> f32 {
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

/// Cross-section profile for cartoon extrusion
///
/// A profile is a 2D shape that gets extruded along the backbone spline.
/// The coordinates are in the normal-binormal plane of the reference frame.
#[derive(Debug, Clone)]
pub struct Profile {
    /// Points defining the profile (in local 2D coordinates)
    pub points: Vec<(f32, f32)>,
    /// Outward-facing normals for each point (for lighting)
    pub normals: Vec<(f32, f32)>,
    /// Rendering type - determines mesh generation method
    pub profile_type: ProfileType,
}

impl Profile {
    /// Create an elliptical profile (for helices) - PyMOL-compatible
    ///
    /// Matches PyMOL's ExtrudeOval function:
    /// - Position: (binormal = cos(θ) * width, normal = sin(θ) * length)
    /// - Normal: (binormal = cos(θ) * length, normal = sin(θ) * width)
    ///
    /// In Rust's (x=normal, y=binormal) coordinate system:
    /// - `width`: oval_length - Half-width in the normal direction (ribbon width, ~1.35)
    /// - `height`: oval_width - Half-height in the binormal direction (ribbon thickness, ~0.25)
    /// - `quality`: Number of points around the ellipse
    pub fn ellipse(width: f32, height: f32, quality: u32) -> Self {
        let n = quality.max(32) as usize; // Minimum 32 for smooth, round appearance
        let mut points = Vec::with_capacity(n);
        let mut normals = Vec::with_capacity(n);

        for i in 0..n {
            let angle = 2.0 * std::f32::consts::PI * (i as f32) / (n as f32);
            let cos_a = angle.cos();
            let sin_a = angle.sin();

            // Point on ellipse - PyMOL mapping:
            // PyMOL: (binormal = cos * width, normal = sin * length)
            // Rust (x=normal, y=binormal): x = sin * width, y = cos * height
            let x = width * sin_a;   // normal direction (wide, ~1.35)
            let y = height * cos_a;  // binormal direction (thin, ~0.25)
            points.push((x, y));

            // Normal to ellipse - PyMOL's formula (swapped width/length for shading):
            // PyMOL: (binormal = cos * length, normal = sin * width)
            // Rust (x=normal, y=binormal): nx = sin * height, ny = cos * width
            let nx = height * sin_a;
            let ny = width * cos_a;
            let len = (nx * nx + ny * ny).sqrt();
            normals.push((nx / len, ny / len));
        }

        Profile { points, normals, profile_type: ProfileType::Round }
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

        Profile { points, normals, profile_type: ProfileType::Round }
    }

    /// Create a flat-faced rectangular profile (for sheets) matching PyMOL's ExtrudeRectangle
    ///
    /// This creates a flat ribbon with 4 faces (top, right, bottom, left) organized
    /// as face pairs for PyMOL-style rendering. Each face has 2 vertices with the same normal.
    ///
    /// - `thickness`: Half-extent in the normal direction (perpendicular to surface, small ~0.2)
    /// - `width`: Half-extent in the binormal direction (lateral ribbon width, large ~1.4)
    ///
    /// PyMOL coordinate mapping:
    /// - Normal direction (first coord): perpendicular to ribbon surface (+X = top, -X = bottom)
    /// - Binormal direction (second coord): lateral ribbon extent (+Y = right, -Y = left)
    ///
    /// 8 vertices organized as 4 face pairs (like PyMOL's ExtrudeRectangle):
    /// - Face 0-1: Top face (normal +X)
    /// - Face 2-3: Right side (normal +Y)  
    /// - Face 4-5: Bottom face (normal -X)
    /// - Face 6-7: Left side (normal -Y)
    ///
    /// IMPORTANT: All horizontal faces (top, bottom) must have vertices in the same
    /// order (left-to-right) for consistent triangle winding in generate_face_strips.
    /// Vertical faces (left, right) must go in the same direction (top-to-bottom).
    pub fn flat_rectangle(thickness: f32, width: f32) -> Self {
        // PyMOL uses cos(π/4) scaling for the rectangular shape
        let cos45 = std::f32::consts::FRAC_1_SQRT_2;
        let t = cos45 * thickness;  // half-thickness scaled
        let w = cos45 * width;      // half-width scaled
        
        // 8 vertices: 2 per face, organized as face pairs
        // CRITICAL: All faces must have consistent vertex ordering for correct winding
        // Horizontal faces go left-to-right, vertical faces go top-to-bottom
        let points = vec![
            // Top face (vertices 0-1): at +t, going left to right (-w to +w)
            (t, -w), (t, w),
            // Right side (vertices 2-3): at +w, going top to bottom (+t to -t)
            (t, w), (-t, w),
            // Bottom face (vertices 4-5): at -t, going left to right (-w to +w) - SAME as top
            (-t, -w), (-t, w),
            // Left side (vertices 6-7): at -w, going top to bottom (+t to -t) - SAME as right
            (t, -w), (-t, -w),
        ];
        
        // Flat normals for each face - all vertices in a face pair have same normal
        let normals = vec![
            // Top face normals (+X direction = perpendicular up from ribbon)
            (1.0, 0.0), (1.0, 0.0),
            // Right side normals (+Y direction = along ribbon width)
            (0.0, 1.0), (0.0, 1.0),
            // Bottom face normals (-X direction = perpendicular down from ribbon)
            (-1.0, 0.0), (-1.0, 0.0),
            // Left side normals (-Y direction = along ribbon width)
            (0.0, -1.0), (0.0, -1.0),
        ];
        
        Profile { points, normals, profile_type: ProfileType::Flat4Face }
    }

    /// Create a dumbbell profile (for fancy helices) matching PyMOL's ExtrudeDumbbell1
    ///
    /// This creates a flat ribbon with top/bottom faces only (no side faces).
    /// Used when cartoon_fancy_helices is enabled.
    ///
    /// - `width`: Half-thickness perpendicular to ribbon surface (dumbbell_width, default 0.17)
    /// - `length`: Half-width of the ribbon laterally (dumbbell_length, default 1.6)
    ///
    /// IMPORTANT: To match ellipse orientation (wide in NORMAL direction):
    /// - First coord (NORMAL direction): ribbon WIDTH, spanning ±l
    /// - Second coord (BINORMAL direction): face position at ±w (ribbon thickness)
    /// - Face normals point in ±binormal direction (perpendicular to flat surfaces)
    ///
    /// PyMOL renders dumbbells as:
    /// 1. A flat ribbon (this profile) - top and bottom faces only
    /// 2. Two separate edge tubes (handled by generate_dumbbell_edge_tubes)
    ///
    /// 4 vertices organized as 2 face pairs:
    /// - Face 0-1: Top face (at +w in binormal, spanning ±l in normal)
    /// - Face 2-3: Bottom face (at -w in binormal, spanning ±l in normal)
    ///
    /// No side faces - the edge tubes provide the lateral boundary at ±l in normal.
    /// IMPORTANT: Both faces must have vertices in the same order (left-to-right)
    /// for consistent triangle winding in generate_face_strips.
    pub fn dumbbell(width: f32, length: f32) -> Self {
        // PyMOL's dumbbell uses cos(π/4) and sin(π/4) for positioning
        let w = std::f32::consts::FRAC_1_SQRT_2 * width;
        let l = std::f32::consts::FRAC_1_SQRT_2 * length;
        
        // 4 vertices: 2 per face, organized as face pairs
        // To match ellipse orientation (wide in normal direction):
        // - First coord (normal): ribbon width, spanning ±l
        // - Second coord (binormal): face position at ±w
        let points = vec![
            // Top face (vertices 0-1): at +w in binormal, spanning -l to +l in normal
            (-l, w), (l, w),
            // Bottom face (vertices 2-3): at -w in binormal, spanning -l to +l in normal
            (-l, -w), (l, -w),
        ];
        
        // Face normals point perpendicular to the flat surfaces (in ±binormal direction)
        let normals = vec![
            // Top face normals point up (+binormal, +Y direction)
            (0.0, 1.0), (0.0, 1.0),
            // Bottom face normals point down (-binormal, -Y direction)
            (0.0, -1.0), (0.0, -1.0),
        ];
        
        Profile { points, normals, profile_type: ProfileType::Flat2Face }
    }

    /// Scale the profile uniformly
    #[allow(dead_code)]
    pub fn scale(&self, factor: f32) -> Self {
        Profile {
            points: self.points.iter().map(|(x, y)| (x * factor, y * factor)).collect(),
            normals: self.normals.clone(),
            profile_type: self.profile_type,
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
        
        // Preserve the profile type from the source profile
        Profile { points, normals, profile_type: self.profile_type }
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
    /// Whether to use dumbbell helices (cartoon_fancy_helices)
    pub fancy_helices: bool,
    /// Dumbbell helix length/height in binormal direction (cartoon_dumbbell_length)
    pub dumbbell_length: f32,
    /// Dumbbell helix width in normal direction (cartoon_dumbbell_width)
    pub dumbbell_width: f32,
    /// Dumbbell edge tube radius (cartoon_dumbbell_radius)
    pub dumbbell_radius: f32,
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
        const CARTOON_FANCY_HELICES: u16 = 118;
        const CARTOON_DUMBBELL_LENGTH: u16 = 115;
        const CARTOON_DUMBBELL_WIDTH: u16 = 116;
        const CARTOON_DUMBBELL_RADIUS: u16 = 117;

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
            fancy_helices: settings.get_bool_if_defined(CARTOON_FANCY_HELICES).unwrap_or(false),
            dumbbell_length: settings.get_float_if_defined(CARTOON_DUMBBELL_LENGTH).unwrap_or(1.6),
            dumbbell_width: settings.get_float_if_defined(CARTOON_DUMBBELL_WIDTH).unwrap_or(0.17),
            dumbbell_radius: settings.get_float_if_defined(CARTOON_DUMBBELL_RADIUS).unwrap_or(0.16),
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
            fancy_helices: false,
            dumbbell_length: 1.6,
            dumbbell_width: 0.17,
            dumbbell_radius: 0.16,
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

    // Track segments of same profile type for batched rendering
    let mut segment_ring_starts: Vec<usize> = Vec::new();
    let mut segment_frame_indices: Vec<usize> = Vec::new(); // Track frame indices for capping
    let mut current_profile_type: Option<ProfileType> = None;
    let mut current_profile_len: usize = 0;
    
    // Detect SS transition boundaries for blending
    let transition_frames = detect_ss_transitions(frames);
    // Number of frames to blend across at each transition
    let blend_range = 8usize;
    
    // Find helix regions for dumbbell tapering
    let helix_regions = if settings.fancy_helices {
        find_helix_regions(frames)
    } else {
        Vec::new()
    };
    // Number of frames to taper at helix ends (PyMOL calls this "sampling")
    let helix_taper_frames = 8usize;

    for (i, frame_meta) in frames.iter().enumerate() {
        // Get the base profile for this frame's SS type
        let base_profile = select_profile(frame_meta, settings, i, frames.len(), sheet_termini, &helix_regions, helix_taper_frames);
        
        // Check if we need to blend with a different profile
        let profile = apply_transition_blending(
            i,
            frames,
            &base_profile,
            settings,
            &transition_frames,
            blend_range,
            sheet_termini,
            &helix_regions,
            helix_taper_frames,
        );
        
        // Check if profile type or vertex count changed - need to flush previous segment
        let type_changed = current_profile_type.map_or(false, |t| t != profile.profile_type);
        let len_changed = current_profile_len != 0 && current_profile_len != profile.len();
        
        if type_changed || len_changed {
            // Check for helix-to-loop transition (fancy helices)
            // Extend loop tube backward to helix boundary position
            let is_helix_to_loop = settings.fancy_helices
                && current_profile_type == Some(ProfileType::Flat2Face)
                && profile.profile_type == ProfileType::Round;

            // Check for loop-to-helix transition (fancy helices)
            // Extend loop tube forward to helix start position
            let is_loop_to_helix = settings.fancy_helices
                && current_profile_type == Some(ProfileType::Round)
                && profile.profile_type == ProfileType::Flat2Face;

            if is_loop_to_helix && !segment_ring_starts.is_empty() {
                // Extend the loop forward slightly into the helix taper region
                // Only extend 2 frames - just enough to cover where edge tubes converge to center
                let loop_profile = Profile::circle(settings.loop_radius, settings.quality);
                let loop_extend_frames = 2usize;

                // Find the helix region that starts at frame i
                let helix_end = helix_regions.iter()
                    .find(|&&(start, _)| start == i)
                    .map(|&(_, end)| end)
                    .unwrap_or(i);

                // Extend into helix by a small amount (or to helix end, whichever is smaller)
                let extend_end = (i + loop_extend_frames).min(helix_end + 1);

                for frame_idx in i..extend_end {
                    let extend_frame = &frames[frame_idx];
                    let ring_start = vertices.len();

                    for j in 0..loop_profile.len() {
                        let world_pos = extend_frame.frame.transform_local(loop_profile.points[j]);
                        let world_normal = extend_frame.frame.local_normal(loop_profile.normals[j]);
                        vertices.push(MeshVertex {
                            position: [world_pos.x, world_pos.y, world_pos.z],
                            normal: [world_normal.x, world_normal.y, world_normal.z],
                            color: extend_frame.color,
                        });
                    }

                    segment_ring_starts.push(ring_start);
                    segment_frame_indices.push(frame_idx);
                }
            }

            // Check if this transition involves a helix (Flat2Face profile)
            let involves_helix = current_profile_type == Some(ProfileType::Flat2Face)
                || profile.profile_type == ProfileType::Flat2Face;

            // Only generate overlap geometry for non-helix transitions (e.g., loop-to-sheet)
            // For helix transitions, the loop extension above handles the connection
            if !involves_helix && !segment_ring_starts.is_empty() && current_profile_len > 0 {
                // Get the previous segment's profile type based on the last frame's SS type
                let prev_ss_type = if let Some(&last_idx) = segment_frame_indices.last() {
                    frames[last_idx].ss_type
                } else {
                    frames[i].ss_type
                };

                // Generate the OLD profile for the current frame position
                let old_profile = create_profile_for_ss(prev_ss_type, settings, i, frames.len(), sheet_termini, &helix_regions, helix_taper_frames);

                // Only extend if profile types match (same vertex count)
                if old_profile.len() == current_profile_len {
                    let ring_start = vertices.len();

                    // Generate vertices at frame i with the OLD profile
                    for j in 0..old_profile.len() {
                        let local_pos = old_profile.points[j];
                        let local_normal = old_profile.normals[j];
                        let world_pos = frame_meta.frame.transform_local(local_pos);
                        let world_normal = frame_meta.frame.local_normal(local_normal);

                        vertices.push(MeshVertex {
                            position: [world_pos.x, world_pos.y, world_pos.z],
                            normal: [world_normal.x, world_normal.y, world_normal.z],
                            color: frame_meta.color,
                        });
                    }

                    segment_ring_starts.push(ring_start);
                    segment_frame_indices.push(i);
                }
            }

            // Flush the previous segment with appropriate rendering
            flush_segment(
                &mut indices,
                &segment_ring_starts,
                current_profile_type.unwrap_or(ProfileType::Round),
                current_profile_len,
            );

            segment_ring_starts.clear();
            segment_frame_indices.clear();

            // For helix-to-loop transition, extend the loop backward slightly into the helix taper region
            // Only extend 2 frames - just enough to cover where edge tubes converge to center
            if is_helix_to_loop && i > 0 {
                let loop_profile = Profile::circle(settings.loop_radius, settings.quality);
                let loop_extend_frames = 2usize;

                // Find the helix region that ends at frame i-1
                let helix_start = helix_regions.iter()
                    .find(|&&(_, end)| end == i - 1)
                    .map(|&(start, _)| start)
                    .unwrap_or(i - 1);

                // Extend back into helix by a small amount (or to helix start, whichever is larger)
                let extend_start = (i - 1).saturating_sub(loop_extend_frames - 1).max(helix_start);

                // Add rings from extend_start to i-1 (in forward order for proper mesh connectivity)
                for frame_idx in extend_start..i {
                    let extend_frame = &frames[frame_idx];
                    let ring_start = vertices.len();

                    for j in 0..loop_profile.len() {
                        let world_pos = extend_frame.frame.transform_local(loop_profile.points[j]);
                        let world_normal = extend_frame.frame.local_normal(loop_profile.normals[j]);
                        vertices.push(MeshVertex {
                            position: [world_pos.x, world_pos.y, world_pos.z],
                            normal: [world_normal.x, world_normal.y, world_normal.z],
                            color: extend_frame.color,
                        });
                    }

                    segment_ring_starts.push(ring_start);
                    segment_frame_indices.push(frame_idx);
                }
                // Note: current_profile_type and current_profile_len will be set
                // in the normal loop iteration below when processing frame i
            }
        }
        
        // Calculate dumbbell taper factor for fancy helices
        // PyMOL tapers the lateral extent (Z component) at helix ends using smooth(f, 2)
        // This creates a smooth visual transition where the ribbon narrows to a line
        let dumbbell_taper = if profile.profile_type == ProfileType::Flat2Face && settings.fancy_helices {
            calculate_dumbbell_taper(i, &helix_regions, helix_taper_frames)
        } else {
            1.0
        };
        
        // Record this ring's start and frame index
        let ring_start = vertices.len();
        segment_ring_starts.push(ring_start);
        segment_frame_indices.push(i);
        current_profile_type = Some(profile.profile_type);
        current_profile_len = profile.len();
        
        // Check if this frame is in an arrow taper region (for normal adjustment)
        let is_arrow = settings.fancy_sheets 
            && frame_meta.ss_type == SecondaryStructure::Sheet
            && is_in_arrow_region(i, sheet_termini, settings.arrow_length);

        // Generate vertices for this profile ring
        for j in 0..profile.len() {
            let local_pos = profile.points[j];
            let local_normal = profile.normals[j];
            
            // Apply dumbbell taper to the X (normal) component
            // This matches PyMOL's ExtrudeCGOSurfacePolygonTaper which scales sv[2] (Z/normal)
            let tapered_pos = if dumbbell_taper < 1.0 {
                (local_pos.0 * dumbbell_taper, local_pos.1)
            } else {
                local_pos
            };

            let world_pos = frame_meta.frame.transform_local(tapered_pos);
            
            // Transform the 2D profile normal to world space using the frame
            let world_normal = if is_arrow {
                // Get adjusted normal with tangent tilt for arrow shading
                let (nx, ny, tilt) = adjust_arrow_normal(local_normal, true, 0.4);
                let base_normal = frame_meta.frame.local_normal((nx, ny));
                let tangent_contrib = frame_meta.frame.tangent * tilt;
                let adjusted = base_normal + tangent_contrib;
                normalize_safe(adjusted)
            } else {
                frame_meta.frame.local_normal(local_normal)
            };

            vertices.push(MeshVertex {
                position: [world_pos.x, world_pos.y, world_pos.z],
                normal: [world_normal.x, world_normal.y, world_normal.z],
                color: frame_meta.color,
            });
        }
    }

    // Flush the final segment
    if let Some(profile_type) = current_profile_type {
        flush_segment(
            &mut indices,
            &segment_ring_starts,
            profile_type,
            current_profile_len,
        );
    }

    // Cap the ends
    if frames.len() >= 2 {
        // Compute taper for start cap
        let start_taper = if settings.fancy_helices {
            calculate_dumbbell_taper(0, &helix_regions, helix_taper_frames)
        } else {
            1.0
        };
        // Cap the start (never an arrow tip)
        cap_tube_end(&mut vertices, &mut indices, &frames[0], settings, true, false, start_taper);
        
        // Compute taper for end cap
        let last_idx = frames.len() - 1;
        let end_taper = if settings.fancy_helices {
            calculate_dumbbell_taper(last_idx, &helix_regions, helix_taper_frames)
        } else {
            1.0
        };
        let is_arrow_tip = is_at_arrow_tip(last_idx, sheet_termini);
        cap_tube_end(
            &mut vertices,
            &mut indices,
            &frames[last_idx],
            settings,
            false,
            is_arrow_tip,
            end_taper,
        );
    }

    // Generate arrow back faces at body/head transition points
    if settings.fancy_sheets && settings.arrow_length > 0 {
        for (i, frame_meta) in frames.iter().enumerate() {
            if frame_meta.ss_type == SecondaryStructure::Sheet {
                if is_at_arrow_transition(i, sheet_termini, settings.arrow_length).is_some() {
                    generate_arrow_back_face(&mut vertices, &mut indices, frame_meta, settings);
                }
            }
        }
    }

    // Generate edge tubes for dumbbell helices (if enabled)
    if settings.fancy_helices {
        let helix_regions = find_helix_regions(frames);
        if !helix_regions.is_empty() {
            let (mut edge_vertices, edge_indices) =
                generate_dumbbell_edge_tubes(frames, &helix_regions, settings);

            // Offset indices for edge tube geometry
            let base_index = vertices.len() as u32;
            let offset_indices: Vec<u32> = edge_indices.iter().map(|i| i + base_index).collect();

            // Append edge tube geometry
            vertices.append(&mut edge_vertices);
            indices.extend(offset_indices);
        }
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
    _helix_regions: &[(usize, usize)],
    _helix_taper_frames: usize,
) -> Profile {
    // Ribbon mode: use uniform circular tube for ALL secondary structures
    if settings.uniform_tube {
        return Profile::circle(settings.loop_radius, settings.quality);
    }
    
    match frame_meta.ss_type {
        SecondaryStructure::Helix
        | SecondaryStructure::Helix310
        | SecondaryStructure::HelixPi => {
            // Profile selection for helices:
            // - fancy_helices: flat dumbbell profile (2 parallel faces) with edge tubes
            // - round_helices: elliptical profile (smooth oval cross-section)
            // - default: rectangular profile
            if settings.fancy_helices {
                // Dumbbell profile: flat ribbon with 2 parallel faces
                // Edge tubes are added separately in generate_dumbbell_edge_tubes()
                // Taper is clamped to minimum 0.15 in the main loop to avoid triangular artifacts
                Profile::dumbbell(settings.dumbbell_width, settings.dumbbell_length)
            } else if settings.round_helices {
                // Ellipse profile: smooth oval cross-section
                // helix_height = oval_length (wide), helix_width = oval_width (thin)
                Profile::ellipse(settings.helix_height, settings.helix_width, settings.quality)
            } else {
                Profile::rectangle(settings.helix_height, settings.helix_width, settings.quality)
            }
        }
        SecondaryStructure::Sheet => {
            // Check if this is near a sheet terminus for arrow tapering
            let scale = if settings.fancy_sheets {
                calculate_arrow_scale(
                    frame_idx,
                    sheet_termini,
                    settings.arrow_length,
                    settings.arrow_tip_scale,
                )
            } else {
                1.0 // No arrow tapering when fancy_sheets is disabled
            };
            
            // Use flat_rectangle for PyMOL-style flat sheet appearance
            // This creates distinct flat top/bottom faces with sharp edges
            // 
            // PyMOL coordinate mapping:
            // - First param (thickness): perpendicular to surface (small ~0.2)
            // - Second param (width): lateral ribbon extent (large ~1.4)
            let ribbon_thickness = settings.sheet_width * 0.5;  // ~0.2, perpendicular to surface
            let ribbon_width = settings.sheet_height * scale;   // ~1.4 * scale, lateral extent
            
            if ribbon_width < 0.01 {
                // At or near the arrow tip - create minimal profile
                Profile::flat_rectangle(ribbon_thickness, 0.01)
            } else {
                Profile::flat_rectangle(ribbon_thickness, ribbon_width)
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

/// Check if a frame index is at the arrow body/head transition point
/// Returns the sheet end index if this frame is at a transition, None otherwise
fn is_at_arrow_transition(
    frame_idx: usize,
    sheet_termini: &[(usize, usize)],
    arrow_length: usize,
) -> Option<usize> {
    if arrow_length == 0 {
        return None;
    }
    
    for &(_start, end) in sheet_termini {
        let transition = end.saturating_sub(arrow_length);
        if frame_idx == transition && transition < end {
            return Some(end);
        }
    }
    None
}

/// Check if a frame index is in the arrow taper region of a sheet
fn is_in_arrow_region(
    frame_idx: usize,
    sheet_termini: &[(usize, usize)],
    arrow_length: usize,
) -> bool {
    if arrow_length == 0 {
        return false;
    }
    
    for &(_start, end) in sheet_termini {
        if frame_idx >= end.saturating_sub(arrow_length) && frame_idx <= end {
            return true;
        }
    }
    false
}

/// Calculate the taper factor for dumbbell profiles at helix ends
///
/// PyMOL tapers the Z (normal) component of the dumbbell profile at helix ends
/// using the smooth() easing function. This creates smooth end caps.
///
/// # Arguments
/// * `frame_idx` - Current frame index
/// * `helix_regions` - (start, end) indices for helix regions  
/// * `taper_frames` - Number of frames to taper at each helix end
///
/// # Returns
/// Taper factor in range [0, 1] where 1.0 = no taper (full size), 0.0 = fully tapered
fn calculate_dumbbell_taper(
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
                // Near start of helix - taper in
                let f = (frame_idx - start) as f32 / taper_frames as f32;
                return smooth(f, 2.0);
            } else if frame_idx > start + sub_n {
                // Near end of helix - taper out  
                let f = (end - frame_idx) as f32 / taper_frames as f32;
                return smooth(f, 2.0);
            }
            // Middle of helix - no taper
            return 1.0;
        }
    }
    1.0 // Not in a helix region
}

/// Adjust normals on arrow taper edges for proper shading (PyMOL-style)
///
/// PyMOL tilts normals outward on the slanted arrow edges by adding 0.4 to the
/// tangent component. This makes the tapered faces catch light properly.
///
/// - `local_normal`: The 2D profile normal in local coordinates (x=normal dir, y=binormal dir)
/// - `is_arrow_region`: Whether this vertex is in the arrow taper region
/// - `tilt_factor`: How much to tilt (PyMOL uses 0.4)
///
/// Returns the adjusted normal (still in local 2D coordinates, but will be transformed
/// to include tangent tilt when converted to 3D)
fn adjust_arrow_normal(
    local_normal: (f32, f32),
    is_arrow_region: bool,
    tilt_factor: f32,
) -> (f32, f32, f32) {
    // The third component represents the tangent direction tilt
    // PyMOL checks if the normal has a significant z component (binormal direction)
    // and tilts the x component (tangent direction) outward
    if is_arrow_region && local_normal.1.abs() > 0.001 {
        // Tilt outward along tangent direction
        // The sign of the tilt should match the sign of the binormal component
        let tilt = tilt_factor * local_normal.1.signum();
        (local_normal.0, local_normal.1, tilt)
    } else {
        (local_normal.0, local_normal.1, 0.0)
    }
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
    helix_regions: &[(usize, usize)],
    helix_taper_frames: usize,
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
    
    // Skip blending for helix transitions (creates arrow-like flaring)
    // PyMOL uses sharp transitions at helix boundaries
    let involves_helix = matches!(
        before_ss,
        SecondaryStructure::Helix | SecondaryStructure::Helix310 | SecondaryStructure::HelixPi
    ) || matches!(
        after_ss,
        SecondaryStructure::Helix | SecondaryStructure::Helix310 | SecondaryStructure::HelixPi
    );
    
    if involves_helix {
        return base_profile.clone();
    }
    
    // Skip blending for sheet transitions (flat_rectangle has different vertex count)
    // This matches PyMOL's sharp-edge behavior for sheets
    let involves_sheet = before_ss == SecondaryStructure::Sheet 
        || after_ss == SecondaryStructure::Sheet;
    
    if involves_sheet {
        return base_profile.clone();
    }
    
    // Get the "other" profile (the one we're transitioning to/from)
    let other_ss = if frame_idx < trans_idx { after_ss } else { before_ss };
    
    // Create profile for the other SS type
    let other_profile = create_profile_for_ss(other_ss, settings, frame_idx, frames.len(), sheet_termini, helix_regions, helix_taper_frames);
    
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
    _helix_regions: &[(usize, usize)],
    _helix_taper_frames: usize,
) -> Profile {
    // Ribbon mode: use uniform circular tube for ALL secondary structures
    if settings.uniform_tube {
        return Profile::circle(settings.loop_radius, settings.quality);
    }
    
    match ss_type {
        SecondaryStructure::Helix
        | SecondaryStructure::Helix310
        | SecondaryStructure::HelixPi => {
            // Profile selection for helices:
            // - fancy_helices: flat dumbbell profile (2 parallel faces) with edge tubes
            // - round_helices: elliptical profile (smooth oval cross-section)
            // - default: rectangular profile
            if settings.fancy_helices {
                // Dumbbell profile: flat ribbon with 2 parallel faces
                // Taper is clamped to minimum 0.15 in the main loop to avoid triangular artifacts
                Profile::dumbbell(settings.dumbbell_width, settings.dumbbell_length)
            } else if settings.round_helices {
                Profile::ellipse(settings.helix_height, settings.helix_width, settings.quality)
            } else {
                Profile::rectangle(settings.helix_height, settings.helix_width, settings.quality)
            }
        }
        SecondaryStructure::Sheet => {
            let scale = if settings.fancy_sheets {
                calculate_arrow_scale(
                    frame_idx,
                    sheet_termini,
                    settings.arrow_length,
                    settings.arrow_tip_scale,
                )
            } else {
                1.0
            };
            // First param = thickness (normal dir), second param = width (binormal dir)
            let ribbon_thickness = settings.sheet_width * 0.5;
            let ribbon_width = settings.sheet_height * scale;
            if ribbon_width < 0.01 {
                Profile::flat_rectangle(ribbon_thickness, 0.01)
            } else {
                Profile::flat_rectangle(ribbon_thickness, ribbon_width)
            }
        }
        _ => Profile::circle(settings.loop_radius, settings.quality),
    }
}

/// Flush accumulated rings to indices using the appropriate rendering method
///
/// Dispatches to either connect_rings() for round profiles or generate_face_strips()
/// for flat profiles, depending on the profile type.
fn flush_segment(
    indices: &mut Vec<u32>,
    ring_starts: &[usize],
    profile_type: ProfileType,
    vertices_per_ring: usize,
) {
    if ring_starts.len() < 2 || vertices_per_ring == 0 {
        return;
    }
    
    match profile_type {
        ProfileType::Round => {
            // Use connect_rings() which wraps around all vertices
            for i in 0..(ring_starts.len() - 1) {
                connect_rings(
                    indices,
                    ring_starts[i] as u32,
                    ring_starts[i + 1] as u32,
                    vertices_per_ring as u32,
                );
            }
        }
        ProfileType::Flat4Face => {
            // 8 vertices, 4 face pairs - use face-based rendering
            generate_face_strips(indices, ring_starts, 4, vertices_per_ring);
        }
        ProfileType::Flat2Face => {
            // 4 vertices, 2 face pairs (top/bottom only) - use face-based rendering
            generate_face_strips(indices, ring_starts, 2, vertices_per_ring);
        }
    }
}

/// Connect two profile rings with triangles (for round profiles)
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

/// Generate mesh using face-based strips (PyMOL-style)
///
/// For flat profiles organized as face pairs, this generates separate triangle strips
/// for each face without wrapping around. This is the correct approach for sheets
/// and dumbbells where we want distinct top/bottom/side faces.
///
/// - `indices`: Output index buffer
/// - `ring_starts`: Start vertex index of each ring along the path
/// - `num_faces`: Number of face pairs (4 for sheets with sides, 2 for dumbbells)
/// - `vertices_per_ring`: Total vertices per ring (num_faces * 2)
fn generate_face_strips(
    indices: &mut Vec<u32>,
    ring_starts: &[usize],
    num_faces: usize,
    vertices_per_ring: usize,
) {
    let num_rings = ring_starts.len();
    if num_rings < 2 {
        return;
    }
    
    // Each face pair consists of 2 consecutive profile vertices
    // Face 0: vertices 0,1
    // Face 1: vertices 2,3
    // Face 2: vertices 4,5
    // etc.
    for face in 0..num_faces {
        let v0_offset = (face * 2) % vertices_per_ring;      // First vertex of face pair
        let v1_offset = (face * 2 + 1) % vertices_per_ring;  // Second vertex of face pair
        
        // Generate triangles along the path for this face
        for ring_idx in 0..(num_rings - 1) {
            let ring1 = ring_starts[ring_idx] as u32;
            let ring2 = ring_starts[ring_idx + 1] as u32;
            
            // Quad vertices: (ring1,v0), (ring1,v1), (ring2,v0), (ring2,v1)
            let r1v0 = ring1 + v0_offset as u32;
            let r1v1 = ring1 + v1_offset as u32;
            let r2v0 = ring2 + v0_offset as u32;
            let r2v1 = ring2 + v1_offset as u32;
            
            // Two triangles per quad - CCW winding when viewed from outside
            // For flat profiles: face normal points in +/- normal direction
            // Vertices are at (normal_coord, binormal_coord) where:
            // - v0 is at binormal = -w (left)
            // - v1 is at binormal = +w (right)
            // When viewed from the face normal direction (+normal for top face):
            // - Tangent points "forward" (into screen)
            // - Binormal points "right"
            // So we need: r1v0 (back-left) -> r2v0 (front-left) -> r1v1 (back-right) for CCW
            // Triangle 1: r1v0 -> r2v0 -> r1v1
            // Triangle 2: r1v1 -> r2v0 -> r2v1
            indices.push(r1v0);
            indices.push(r2v0);
            indices.push(r1v1);
            
            indices.push(r1v1);
            indices.push(r2v0);
            indices.push(r2v1);
        }
    }
}

/// Generate the back face of an arrow head (flat quad closing the arrow body)
///
/// For 8-vertex Flat4Face profiles with consistent vertex ordering, corners are:
/// - Index 0 (top-left), 1 (top-right)
/// - Index 4 (bottom-left), 5 (bottom-right)
///
/// The back face closes the gap between the body and the wider arrow head.
fn generate_arrow_back_face(
    vertices: &mut Vec<MeshVertex>,
    indices: &mut Vec<u32>,
    frame_meta: &FrameWithMetadata,
    settings: &CartoonGeometrySettings,
) {
    // The back face normal points backward (negative tangent direction)
    let back_normal = frame_meta.frame.tangent * -1.0;
    
    // Get the arrow profile at the transition point (at 1.5x scale)
    let ribbon_thickness = settings.sheet_width * 0.5;
    // Arrow head starts at 1.5x width
    let ribbon_width = settings.sheet_height * settings.arrow_tip_scale;
    let profile = Profile::flat_rectangle(ribbon_thickness, ribbon_width);
    
    // For Flat4Face profiles (8 vertices), the rectangle corners are:
    // Top-left: vertex 0, Top-right: vertex 1
    // Bottom-left: vertex 4, Bottom-right: vertex 5
    // (vertices 2-3 and 6-7 are the right and left side edges)
    
    // We'll create new vertices for the back face with the back-pointing normal
    let start_idx = vertices.len() as u32;
    
    // Get the four corners - order: top-left, top-right, bottom-right, bottom-left
    let corner_indices = [0, 1, 5, 4];
    for &idx in &corner_indices {
        let local_pos = profile.points[idx];
        let world_pos = frame_meta.frame.transform_local(local_pos);
        
        vertices.push(MeshVertex {
            position: [world_pos.x, world_pos.y, world_pos.z],
            normal: [back_normal.x, back_normal.y, back_normal.z],
            color: frame_meta.color,
        });
    }
    
    // Create two triangles for the quad (CCW winding when looking at the back)
    // Vertices in our array: 0=top-left, 1=top-right, 2=bottom-right, 3=bottom-left
    // Triangle 1: top-left -> bottom-left -> top-right
    indices.push(start_idx + 0);
    indices.push(start_idx + 3);
    indices.push(start_idx + 1);
    
    // Triangle 2: bottom-left -> bottom-right -> top-right
    indices.push(start_idx + 3);
    indices.push(start_idx + 2);
    indices.push(start_idx + 1);
}

/// Add end cap to the tube
///
/// - `is_start`: true for the start of a segment, false for the end
/// - `is_arrow_tip`: true if this is at an arrow tip (skip cap since geometry tapers to a point)
/// - `dumbbell_taper`: taper factor for Flat2Face caps (0.0 = fully tapered, 1.0 = full size)
fn cap_tube_end(
    vertices: &mut Vec<MeshVertex>,
    indices: &mut Vec<u32>,
    frame_meta: &FrameWithMetadata,
    settings: &CartoonGeometrySettings,
    is_start: bool,
    is_arrow_tip: bool,
    _dumbbell_taper: f32,
) {
    // Skip cap for arrow tips - the geometry tapers to a point
    if is_arrow_tip && !is_start {
        return;
    }

    // Ribbon mode: use uniform circular tube for ALL secondary structures
    let profile = if settings.uniform_tube {
        Profile::circle(settings.loop_radius, settings.quality)
    } else {
        match frame_meta.ss_type {
            SecondaryStructure::Helix
            | SecondaryStructure::Helix310
            | SecondaryStructure::HelixPi => {
                // Profile selection for helix caps - must match select_profile()
                if settings.fancy_helices {
                    Profile::dumbbell(settings.dumbbell_width, settings.dumbbell_length)
                } else if settings.round_helices {
                    Profile::ellipse(settings.helix_height, settings.helix_width, settings.quality)
                } else {
                    Profile::rectangle(settings.helix_height, settings.helix_width, settings.quality)
                }
            }
            SecondaryStructure::Sheet => {
                // Use flat_rectangle for sheets (consistent with select_profile)
                // First param = thickness, second param = width (lateral)
                Profile::flat_rectangle(settings.sheet_width * 0.5, settings.sheet_height)
            }
            _ => Profile::circle(settings.loop_radius, settings.quality),
        }
    };

    let normal_dir = if is_start { -1.0 } else { 1.0 };
    let cap_normal = frame_meta.frame.tangent * normal_dir;

    // Different cap generation based on profile type
    match profile.profile_type {
        ProfileType::Round => {
            // Use triangle fan from center for round profiles
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
        ProfileType::Flat4Face => {
            // 8-vertex face-pair profile - use quad cap with corners
            // With new vertex ordering:
            // 0: (t, -w) top-left, 1: (t, w) top-right
            // 4: (-t, -w) bottom-left, 5: (-t, w) bottom-right
            let start_idx = vertices.len() as u32;
            let corner_indices = [0, 1, 5, 4]; // top-left, top-right, bottom-right, bottom-left
            
            for &idx in &corner_indices {
                let world_pos = frame_meta.frame.transform_local(profile.points[idx]);
                vertices.push(MeshVertex {
                    position: [world_pos.x, world_pos.y, world_pos.z],
                    normal: [cap_normal.x, cap_normal.y, cap_normal.z],
                    color: frame_meta.color,
                });
            }
            
            // Two triangles for the quad
            if is_start {
                // Reverse winding for start cap
                indices.push(start_idx + 0);
                indices.push(start_idx + 3);
                indices.push(start_idx + 1);
                
                indices.push(start_idx + 3);
                indices.push(start_idx + 2);
                indices.push(start_idx + 1);
            } else {
                indices.push(start_idx + 0);
                indices.push(start_idx + 1);
                indices.push(start_idx + 3);
                
                indices.push(start_idx + 1);
                indices.push(start_idx + 2);
                indices.push(start_idx + 3);
            }
        }
        ProfileType::Flat2Face => {
            // Skip cap generation for dumbbell helices - tapering provides visual termination
            // The ribbon tapers to a line at helix ends, and edge tubes taper to center
            // This avoids the triangular artifacts that occurred with wedge caps
            // PyMOL also doesn't render explicit caps here - the visual termination comes
            // from the tapered geometry converging at the transition point
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

/// Find helix region indices in a sequence of frames
///
/// Returns a vector of (start_idx, end_idx) tuples for each contiguous helix region.
pub fn find_helix_regions(frames: &[FrameWithMetadata]) -> Vec<(usize, usize)> {
    let mut regions = Vec::new();
    let mut in_helix = false;
    let mut helix_start = 0;

    for (i, frame) in frames.iter().enumerate() {
        let is_helix = matches!(
            frame.ss_type,
            SecondaryStructure::Helix | SecondaryStructure::Helix310 | SecondaryStructure::HelixPi
        );

        if is_helix && !in_helix {
            // Helix starts
            helix_start = i;
            in_helix = true;
        } else if !is_helix && in_helix {
            // Helix ends
            regions.push((helix_start, i - 1));
            in_helix = false;
        }
    }

    // Handle helix at end
    if in_helix {
        regions.push((helix_start, frames.len() - 1));
    }

    regions
}

/// Generate edge tube geometry for dumbbell helices (PyMOL-style)
///
/// Creates two cylindrical tubes along the edges of dumbbell helices.
/// The tubes are offset from the backbone by `±sin(π/4) * dumbbell_length`
/// in the normal direction, with smooth tapering at helix ends.
///
/// Matches PyMOL's ExtrudeDumbbellEdge function behavior.
///
/// # Arguments
/// * `frames` - Reference frames along the backbone
/// * `helix_regions` - (start, end) indices for helix regions
/// * `settings` - Geometry settings containing dumbbell parameters
///
/// # Returns
/// (vertices, indices) for the edge tubes
fn generate_dumbbell_edge_tubes(
    frames: &[FrameWithMetadata],
    helix_regions: &[(usize, usize)],
    settings: &CartoonGeometrySettings,
) -> (Vec<MeshVertex>, Vec<u32>) {
    if !settings.fancy_helices || helix_regions.is_empty() {
        return (Vec::new(), Vec::new());
    }

    let mut vertices = Vec::new();
    let mut indices = Vec::new();
    
    // PyMOL edge tube offset: sin(π/4) * dumbbell_length
    // This positions the tubes at the lateral edges of the dumbbell ribbon
    let sin45 = std::f32::consts::FRAC_1_SQRT_2;
    let base_offset = sin45 * settings.dumbbell_length;
    
    // Tube radius from settings
    let tube_radius = settings.dumbbell_radius;
    
    // Circular profile for edge tubes
    let tube_quality = 16u32;
    let tube_profile = Profile::circle(tube_radius, tube_quality);
    
    // Number of frames to taper at helix ends (PyMOL's "sampling")
    let taper_frames = 8usize;

    for &(start, end) in helix_regions {
        let region_len = end - start + 1;
        let sub_n = region_len.saturating_sub(taper_frames);

        // Generate two edge tubes: one at +offset, one at -offset
        for sign in [-1.0_f32, 1.0_f32] {
            let mut prev_ring_start: Option<usize> = None;

            // Generate tube along helix region
            for i in start..=end {
                let frame_meta = &frames[i];
                let ring_start = vertices.len();

                // Calculate tapered offset using PyMOL's smooth function
                // Edge tubes taper to 0 at helix ends, converging to center point
                // This provides the visual transition while the main ribbon stays at full width
                let frame_in_region = i - start;
                let taper_factor = if frame_in_region < taper_frames {
                    // Near start - taper in from 0 to 1
                    smooth(frame_in_region as f32 / taper_frames as f32, 2.0)
                } else if frame_in_region > sub_n {
                    // Near end - taper out from 1 to 0
                    smooth((end - i) as f32 / taper_frames as f32, 2.0)
                } else {
                    // Middle - full offset
                    1.0
                };

                let offset = base_offset * sign * taper_factor;

                // Offset position along normal direction (lateral edge of ribbon)
                let offset_pos = frame_meta.frame.position + frame_meta.frame.normal * offset;

                // Add flat cap at the center where edge tube starts (first frame of helix)
                if i == start {
                    cap_edge_tube_end(
                        &mut vertices,
                        &mut indices,
                        frame_meta,
                        0.0, // offset = 0 since taper_factor is 0 at start
                        tube_radius,
                        tube_quality,
                        true, // is_start
                    );
                }
                
                // Generate vertices for circular cross-section
                for j in 0..tube_profile.len() {
                    let local_pos = tube_profile.points[j];
                    let local_normal = tube_profile.normals[j];
                    
                    // Transform local position to world space
                    // The tube is centered at offset_pos, oriented along the tangent
                    let world_pos = offset_pos 
                        + frame_meta.frame.normal * local_pos.0 
                        + frame_meta.frame.binormal * local_pos.1;
                    
                    // Transform normal to world space
                    let world_normal = normalize_safe(
                        frame_meta.frame.normal * local_normal.0 
                        + frame_meta.frame.binormal * local_normal.1
                    );
                    
                    vertices.push(MeshVertex {
                        position: [world_pos.x, world_pos.y, world_pos.z],
                        normal: [world_normal.x, world_normal.y, world_normal.z],
                        color: frame_meta.color,
                    });
                }
                
                // Connect to previous ring
                if let Some(prev_start) = prev_ring_start {
                    connect_rings(
                        &mut indices,
                        prev_start as u32,
                        ring_start as u32,
                        tube_profile.len() as u32,
                    );
                }
                
                prev_ring_start = Some(ring_start);

                // Add flat cap at the center where edge tube ends (last frame of helix)
                if i == end {
                    cap_edge_tube_end(
                        &mut vertices,
                        &mut indices,
                        frame_meta,
                        0.0, // offset = 0 since taper_factor is 0 at end
                        tube_radius,
                        tube_quality,
                        false, // is_start = false (end cap)
                    );
                }
            }
        }
    }

    (vertices, indices)
}

/// Add flat disc end cap to an edge tube (PyMOL cCylCap::Flat style)
///
/// Creates a simple flat circular cap using a triangle fan from center to edge vertices.
/// This matches PyMOL's flat cap style for dumbbell edge tubes.
fn cap_edge_tube_end(
    vertices: &mut Vec<MeshVertex>,
    indices: &mut Vec<u32>,
    frame_meta: &FrameWithMetadata,
    normal_offset: f32,
    radius: f32,
    quality: u32,
    is_start: bool,
) {
    // Offset position along normal direction (tube center at this frame)
    let center_pos = frame_meta.frame.position + frame_meta.frame.normal * normal_offset;
    
    // Cap normal direction (pointing away from helix body)
    let cap_dir = if is_start { -1.0 } else { 1.0 };
    let cap_normal = frame_meta.frame.tangent * cap_dir;
    
    // Add center vertex for the flat disc
    let center_idx = vertices.len() as u32;
    vertices.push(MeshVertex {
        position: [center_pos.x, center_pos.y, center_pos.z],
        normal: [cap_normal.x, cap_normal.y, cap_normal.z],
        color: frame_meta.color,
    });
    
    // Add edge vertices around the disc circumference
    let edge_start = vertices.len() as u32;
    let n = quality as usize;
    
    for i in 0..n {
        let theta = 2.0 * std::f32::consts::PI * (i as f32 / n as f32);
        let cos_theta = theta.cos();
        let sin_theta = theta.sin();
        
        // Position on disc edge
        let local_x = radius * cos_theta; // normal direction
        let local_y = radius * sin_theta; // binormal direction
        
        let world_pos = center_pos
            + frame_meta.frame.normal * local_x
            + frame_meta.frame.binormal * local_y;
        
        vertices.push(MeshVertex {
            position: [world_pos.x, world_pos.y, world_pos.z],
            normal: [cap_normal.x, cap_normal.y, cap_normal.z],
            color: frame_meta.color,
        });
    }
    
    // Generate triangle fan indices
    for i in 0..n as u32 {
        let i_next = (i + 1) % n as u32;

        if is_start {
            // Reverse winding for start cap
            indices.push(center_idx);
            indices.push(edge_start + i_next);
            indices.push(edge_start + i);
        } else {
            indices.push(center_idx);
            indices.push(edge_start + i);
            indices.push(edge_start + i_next);
        }
    }
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
    
    #[test]
    fn test_flat_rectangle_profile() {
        // flat_rectangle(thickness, width) - 8-vertex face-pair layout
        let profile = Profile::flat_rectangle(0.2, 1.4);
        
        // Should have 8 vertices (2 per face, 4 faces)
        assert_eq!(profile.len(), 8);
        assert_eq!(profile.profile_type, ProfileType::Flat4Face);
        
        // Check normals are unit vectors
        for (nx, ny) in &profile.normals {
            let len = (nx * nx + ny * ny).sqrt();
            assert!((len - 1.0).abs() < 0.01, "Normal should be unit length, got {}", len);
        }
        
        // Expected scaled values (cos45 ≈ 0.707)
        let cos45 = std::f32::consts::FRAC_1_SQRT_2;
        let t = cos45 * 0.2;   // ~0.14
        let w = cos45 * 1.4;   // ~0.99
        
        // Face 0-1: Top face at (t, -w), (t, w) with normal (1, 0)
        assert!((profile.points[0].0 - t).abs() < 0.01 && (profile.points[0].1 - (-w)).abs() < 0.01);
        assert!((profile.points[1].0 - t).abs() < 0.01 && (profile.points[1].1 - w).abs() < 0.01);
        assert!((profile.normals[0].0 - 1.0).abs() < 0.01);
        assert!((profile.normals[1].0 - 1.0).abs() < 0.01);
        
        // Face 2-3: Right side at (t, w), (-t, w) with normal (0, 1)
        assert!((profile.points[2].0 - t).abs() < 0.01 && (profile.points[2].1 - w).abs() < 0.01);
        assert!((profile.points[3].0 - (-t)).abs() < 0.01 && (profile.points[3].1 - w).abs() < 0.01);
        assert!((profile.normals[2].1 - 1.0).abs() < 0.01);
        assert!((profile.normals[3].1 - 1.0).abs() < 0.01);
        
        // Face 4-5: Bottom face at (-t, -w), (-t, w) with normal (-1, 0) - left to right
        assert!((profile.points[4].0 - (-t)).abs() < 0.01 && (profile.points[4].1 - (-w)).abs() < 0.01);
        assert!((profile.points[5].0 - (-t)).abs() < 0.01 && (profile.points[5].1 - w).abs() < 0.01);
        assert!((profile.normals[4].0 - (-1.0)).abs() < 0.01);
        assert!((profile.normals[5].0 - (-1.0)).abs() < 0.01);
        
        // Face 6-7: Left side at (t, -w), (-t, -w) with normal (0, -1) - top to bottom
        assert!((profile.points[6].0 - t).abs() < 0.01 && (profile.points[6].1 - (-w)).abs() < 0.01);
        assert!((profile.points[7].0 - (-t)).abs() < 0.01 && (profile.points[7].1 - (-w)).abs() < 0.01);
        assert!((profile.normals[6].1 - (-1.0)).abs() < 0.01);
        assert!((profile.normals[7].1 - (-1.0)).abs() < 0.01);
    }
    
    #[test]
    fn test_dumbbell_profile() {
        let profile = Profile::dumbbell(0.17, 1.6);
        
        // Should have 4 vertices (2 face pairs: top and bottom)
        assert_eq!(profile.len(), 4);
        assert_eq!(profile.profile_type, ProfileType::Flat2Face);
        
        // Check normals are unit vectors
        for (nx, ny) in &profile.normals {
            let len = (nx * nx + ny * ny).sqrt();
            assert!((len - 1.0).abs() < 0.01, "Normal should be unit length, got {}", len);
        }
        
        // Expected vertex positions with cos(π/4) ≈ 0.707 scaling
        let w = std::f32::consts::FRAC_1_SQRT_2 * 0.17;  // ~0.12
        let l = std::f32::consts::FRAC_1_SQRT_2 * 1.6;   // ~1.13
        
        // Coordinate mapping to match ellipse orientation (wide in normal direction):
        // - First coord (normal direction): ribbon width, spanning ±l
        // - Second coord (binormal direction): face position at ±w
        // Face 0-1: Top face at (-l, w), (l, w) with normal (0, 1) - left to right in normal
        // Face 2-3: Bottom face at (-l, -w), (l, -w) with normal (0, -1) - left to right in normal
        assert!((profile.points[0].0 - (-l)).abs() < 0.01 && (profile.points[0].1 - w).abs() < 0.01);
        assert!((profile.points[1].0 - l).abs() < 0.01 && (profile.points[1].1 - w).abs() < 0.01);
        assert!((profile.points[2].0 - (-l)).abs() < 0.01 && (profile.points[2].1 - (-w)).abs() < 0.01);
        assert!((profile.points[3].0 - l).abs() < 0.01 && (profile.points[3].1 - (-w)).abs() < 0.01);
        
        // Check that top face (v0, v1) has +Y normals, bottom face (v2, v3) has -Y normals
        assert!((profile.normals[0].0 - 0.0).abs() < 0.01 && (profile.normals[0].1 - 1.0).abs() < 0.01);
        assert!((profile.normals[1].0 - 0.0).abs() < 0.01 && (profile.normals[1].1 - 1.0).abs() < 0.01);
        assert!((profile.normals[2].0 - 0.0).abs() < 0.01 && (profile.normals[2].1 - (-1.0)).abs() < 0.01);
        assert!((profile.normals[3].0 - 0.0).abs() < 0.01 && (profile.normals[3].1 - (-1.0)).abs() < 0.01);
    }
    
    #[test]
    fn test_is_in_arrow_region() {
        let sheet_termini = vec![(1, 10)];
        let arrow_length = 3;
        
        // Before arrow region
        assert!(!is_in_arrow_region(0, &sheet_termini, arrow_length));
        assert!(!is_in_arrow_region(5, &sheet_termini, arrow_length));
        assert!(!is_in_arrow_region(6, &sheet_termini, arrow_length));
        
        // In arrow region (frames 7-10 when arrow_length=3, end=10)
        assert!(is_in_arrow_region(7, &sheet_termini, arrow_length));
        assert!(is_in_arrow_region(8, &sheet_termini, arrow_length));
        assert!(is_in_arrow_region(9, &sheet_termini, arrow_length));
        assert!(is_in_arrow_region(10, &sheet_termini, arrow_length));
        
        // After arrow region
        assert!(!is_in_arrow_region(11, &sheet_termini, arrow_length));
    }
    
    #[test]
    fn test_adjust_arrow_normal() {
        // Non-arrow region should not modify
        let (_nx, _ny, tilt) = adjust_arrow_normal((0.0, 1.0), false, 0.4);
        assert!((tilt - 0.0).abs() < 0.01);
        
        // Arrow region with binormal component should add tilt
        let (_nx, _ny, tilt) = adjust_arrow_normal((0.0, 1.0), true, 0.4);
        assert!((tilt - 0.4).abs() < 0.01);
        
        // Arrow region with negative binormal component should add negative tilt
        let (_nx, _ny, tilt) = adjust_arrow_normal((0.0, -1.0), true, 0.4);
        assert!((tilt - (-0.4)).abs() < 0.01);
        
        // Arrow region with zero binormal component should not add tilt
        let (_nx, _ny, tilt) = adjust_arrow_normal((1.0, 0.0), true, 0.4);
        assert!((tilt - 0.0).abs() < 0.01);
    }
    
    #[test]
    fn test_find_helix_regions() {
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
            make_frame(SecondaryStructure::Helix),
            make_frame(SecondaryStructure::Helix),
            make_frame(SecondaryStructure::Helix),
            make_frame(SecondaryStructure::Loop),
            make_frame(SecondaryStructure::Helix310),
            make_frame(SecondaryStructure::Helix310),
            make_frame(SecondaryStructure::Sheet),
            make_frame(SecondaryStructure::HelixPi),
            make_frame(SecondaryStructure::Loop),
        ];

        let regions = find_helix_regions(&frames);
        assert_eq!(regions.len(), 3);
        assert_eq!(regions[0], (1, 3)); // Helix
        assert_eq!(regions[1], (5, 6)); // Helix310
        assert_eq!(regions[2], (8, 8)); // HelixPi (single frame)
    }
    
    #[test]
    fn test_dumbbell_edge_tubes_disabled() {
        let make_frame = |ss, pos: f32| FrameWithMetadata {
            frame: ReferenceFrame::new(
                Vec3::new(pos, 0.0, 0.0),
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
            make_frame(SecondaryStructure::Helix, 0.0),
            make_frame(SecondaryStructure::Helix, 1.0),
            make_frame(SecondaryStructure::Helix, 2.0),
        ];
        
        let helix_regions = find_helix_regions(&frames);
        assert_eq!(helix_regions.len(), 1);
        
        // With fancy_helices disabled, should return empty
        let settings = CartoonGeometrySettings {
            fancy_helices: false,
            ..Default::default()
        };
        
        let (vertices, indices) = generate_dumbbell_edge_tubes(&frames, &helix_regions, &settings);
        assert!(vertices.is_empty());
        assert!(indices.is_empty());
    }
    
    #[test]
    fn test_dumbbell_edge_tubes_enabled() {
        let make_frame = |ss, pos: f32| FrameWithMetadata {
            frame: ReferenceFrame::new(
                Vec3::new(pos, 0.0, 0.0),
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
            make_frame(SecondaryStructure::Helix, 0.0),
            make_frame(SecondaryStructure::Helix, 1.0),
            make_frame(SecondaryStructure::Helix, 2.0),
            make_frame(SecondaryStructure::Helix, 3.0),
        ];
        
        let helix_regions = find_helix_regions(&frames);
        assert_eq!(helix_regions.len(), 1);
        
        // With fancy_helices enabled, should generate edge tubes
        let settings = CartoonGeometrySettings {
            fancy_helices: true,
            dumbbell_length: 1.6,
            dumbbell_radius: 0.16,
            ..Default::default()
        };
        
        let (vertices, indices) = generate_dumbbell_edge_tubes(&frames, &helix_regions, &settings);
        
        // Should have generated geometry for 2 edge tubes
        assert!(!vertices.is_empty());
        assert!(!indices.is_empty());
        
        // Check that indices are valid
        let max_vertex_idx = vertices.len() as u32;
        for &idx in &indices {
            assert!(idx < max_vertex_idx, "Index {} out of bounds (max: {})", idx, max_vertex_idx);
        }
    }
}
