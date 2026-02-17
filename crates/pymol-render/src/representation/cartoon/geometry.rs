//! Geometry generation for cartoon representation
//!
//! Provides cross-section profiles and mesh extrusion for different
//! secondary structure types (helices, sheets, loops).

use lin_alg::f32::Vec3;
use pymol_mol::SecondaryStructure;
use pymol_settings::SettingResolver;

use super::frame::FrameWithMetadata;
use super::utils::{find_ss_regions, is_helix, normalize_safe, smooth};
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
    /// This creates a flat ribbon with 4 faces organized as face pairs for PyMOL-style
    /// rendering. Each face has 2 vertices with the same normal.
    ///
    /// Called as `flat_rectangle(sheet_height, sheet_width)`:
    /// - `thickness`: Half-extent in the normal direction (ribbon width, LARGE ~1.4)
    /// - `width`: Half-extent in the binormal direction (ribbon thickness, SMALL ~0.4)
    ///
    /// This matches PyMOL's ExtrudeRectangle(width, length) where:
    /// - length (cartoon_rect_length=1.4) → Z/normal direction → our first coord (thickness)
    /// - width (cartoon_rect_width=0.4) → Y/binormal direction → our second coord (width)
    ///
    /// 8 vertices organized as 4 face pairs:
    /// - Face 0-1: Side edge at +t normal (narrow edge of ribbon)
    /// - Face 2-3: Flat top surface at +w binormal (wide visible face)
    /// - Face 4-5: Side edge at -t normal (narrow edge of ribbon)
    /// - Face 6-7: Flat bottom surface at -w binormal (wide visible face)
    pub fn flat_rectangle(thickness: f32, width: f32) -> Self {
        // PyMOL uses cos(π/4) scaling for the rectangular shape
        let cos45 = std::f32::consts::FRAC_1_SQRT_2;
        let t = cos45 * thickness;  // half-extent in normal dir (LARGE = ribbon width)
        let w = cos45 * width;      // half-extent in binormal dir (SMALL = ribbon thickness)

        // 8 vertices: 2 per face, organized as face pairs
        // t = normal direction (first coord, ribbon width), w = binormal direction (second coord, ribbon thickness)
        let points = vec![
            // Face 0-1: Side edge at +t normal, spanning ±w binormal
            (t, -w), (t, w),
            // Face 2-3: Flat top surface at +w binormal, spanning +t to -t normal
            (t, w), (-t, w),
            // Face 4-5: Side edge at -t normal, spanning ±w binormal
            (-t, -w), (-t, w),
            // Face 6-7: Flat bottom surface at -w binormal, spanning +t to -t normal
            (t, -w), (-t, -w),
        ];

        // Flat normals for each face - all vertices in a face pair have same normal
        let normals = vec![
            // Face 0-1: +normal direction (side edge outward)
            (1.0, 0.0), (1.0, 0.0),
            // Face 2-3: +binormal direction (flat top surface outward)
            (0.0, 1.0), (0.0, 1.0),
            // Face 4-5: -normal direction (side edge outward)
            (-1.0, 0.0), (-1.0, 0.0),
            // Face 6-7: -binormal direction (flat bottom surface outward)
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

    /// Get the number of points in the profile
    pub fn len(&self) -> usize {
        self.points.len()
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
            fancy_sheets: settings.get_bool_if_defined(CARTOON_FANCY_SHEETS).unwrap_or(true),
            fancy_helices: settings.get_bool_if_defined(CARTOON_FANCY_HELICES).unwrap_or(false),
            dumbbell_length: settings.get_float_if_defined(CARTOON_DUMBBELL_LENGTH).unwrap_or(1.6),
            dumbbell_width: settings.get_float_if_defined(CARTOON_DUMBBELL_WIDTH).unwrap_or(0.17),
            dumbbell_radius: settings.get_float_if_defined(CARTOON_DUMBBELL_RADIUS).unwrap_or(0.16),
            arrow_tip_scale: 1.5,
            arrow_length: 0, // Will be calculated from subdivisions
            arrow_residues: 1, // Match PyMOL's sampling (~7-8 frames)
            uniform_tube: false, // Default: different profiles for different SS types
        }
    }

    /// Set the arrow length based on spline subdivisions
    ///
    /// The arrow head length matches PyMOL's sampling parameter.
    /// PyMOL uses sampling (~subdivisions) as arrow head length.
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
            fancy_sheets: true,
            fancy_helices: false,
            dumbbell_length: 1.6,
            dumbbell_width: 0.17,
            dumbbell_radius: 0.16,
            arrow_tip_scale: 1.5,
            arrow_length: 8, // Default for subdivisions=7: 1 * (7+1) = 8
            arrow_residues: 1,
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
        // Skip ALL sheet frames — they are rendered by generate_explicit_sheet().
        // No overlap: loop ends exactly at boundary, sheet starts exactly at boundary.
        if frame_meta.ss_type == SecondaryStructure::Sheet && !settings.uniform_tube {
            // Flush any pending segment before skipping
            if let Some(profile_type) = current_profile_type {
                if !segment_ring_starts.is_empty() {
                    flush_segment(
                        &mut indices,
                        &segment_ring_starts,
                        profile_type,
                        current_profile_len,
                    );
                    segment_ring_starts.clear();
                    segment_frame_indices.clear();
                }
            }
            current_profile_type = None;
            current_profile_len = 0;
            continue;
        }

        // Get the base profile for this frame's SS type
        let base_profile = profile_for_ss(frame_meta.ss_type, settings);

        // Check if we need to blend with a different profile
        let profile = apply_transition_blending(
            i,
            frames,
            &base_profile,
            settings,
            &transition_frames,
            blend_range,
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
                let old_profile = profile_for_ss(prev_ss_type, settings);

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

        // Generate vertices for this profile ring
        for j in 0..profile.len() {
            let local_pos = profile.points[j];
            let local_normal = profile.normals[j];

            // Apply dumbbell taper for fancy helices at helix ends
            let tapered_pos = if dumbbell_taper < 1.0 {
                (local_pos.0 * dumbbell_taper, local_pos.1)
            } else {
                local_pos
            };

            let world_pos = frame_meta.frame.transform_local(tapered_pos);
            let world_normal = frame_meta.frame.local_normal(local_normal);

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

    // Generate explicit sheet geometry for each sheet run
    if !settings.uniform_tube {
        for &(start, end) in sheet_termini {
            generate_explicit_sheet(&mut vertices, &mut indices, frames, start, end, settings);
        }
    }

    // Cap the ends (skip if the end frame is a sheet — explicit sheets have their own caps)
    if frames.len() >= 2 {
        let first_is_sheet = frames[0].ss_type == SecondaryStructure::Sheet && !settings.uniform_tube;
        let last_idx = frames.len() - 1;
        let last_is_sheet = frames[last_idx].ss_type == SecondaryStructure::Sheet && !settings.uniform_tube;

        if !first_is_sheet {
            cap_tube_end(&mut vertices, &mut indices, &frames[0], settings, true, false);
        }

        if !last_is_sheet {
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

/// Select the appropriate cross-section profile for a secondary structure type
fn profile_for_ss(ss_type: SecondaryStructure, settings: &CartoonGeometrySettings) -> Profile {
    // Ribbon mode: use uniform circular tube for ALL secondary structures
    if settings.uniform_tube {
        return Profile::circle(settings.loop_radius, settings.quality);
    }

    if is_helix(ss_type) {
        // Profile selection for helices:
        // - fancy_helices: flat dumbbell profile (2 parallel faces) with edge tubes
        // - round_helices: elliptical profile (smooth oval cross-section)
        // - default: rectangular profile
        if settings.fancy_helices {
            Profile::dumbbell(settings.dumbbell_width, settings.dumbbell_length)
        } else if settings.round_helices {
            Profile::ellipse(settings.helix_height, settings.helix_width, settings.quality)
        } else {
            Profile::rectangle(settings.helix_height, settings.helix_width, settings.quality)
        }
    } else if ss_type == SecondaryStructure::Sheet {
        Profile::flat_rectangle(settings.sheet_height, settings.sheet_width)
    } else {
        Profile::circle(settings.loop_radius, settings.quality)
    }
}

/// Check if a frame index is at an arrow tip (end of a sheet segment)
fn is_at_arrow_tip(frame_idx: usize, sheet_termini: &[(usize, usize)]) -> bool {
    sheet_termini.iter().any(|&(_start, end)| frame_idx == end)
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
    let involves_helix = is_helix(before_ss) || is_helix(after_ss);

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
    let other_profile = profile_for_ss(other_ss, settings);

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
pub fn connect_rings(indices: &mut Vec<u32>, ring1_start: u32, ring2_start: u32, profile_len: u32) {
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
pub fn generate_face_strips(
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

/// Generate a tube connector bridging a sheet boundary to a loop region.
///
/// Creates a smooth tube that transitions between full loop radius and the
/// thin sheet thickness. Used at both N-terminal and C-terminal sheet boundaries.
///
/// * `frame_start`..`frame_end` — range of frame indices to generate rings for (exclusive end)
/// * `sheet_normal` — the normal from the sheet boundary frame, for consistent orientation
/// * `expanding` — if false, scale goes 1.0→min (shrinking toward sheet);
///                 if true, scale goes min→1.0 (expanding away from sheet)
fn generate_sheet_connector(
    vertices: &mut Vec<MeshVertex>,
    indices: &mut Vec<u32>,
    frames: &[FrameWithMetadata],
    frame_start: usize,
    frame_end: usize,
    sheet_normal: Vec3,
    settings: &CartoonGeometrySettings,
    half_thick: f32,
    expanding: bool,
) {
    let loop_profile = Profile::circle(settings.loop_radius, settings.quality);
    let profile_len = loop_profile.len();
    let frame_count = frame_end - frame_start;
    if frame_count < 2 {
        return;
    }
    let mut ring_starts: Vec<usize> = Vec::new();
    let min_scale = (half_thick / settings.loop_radius).min(1.0);

    for fi in frame_start..frame_end {
        let f = &frames[fi];
        let ring_start = vertices.len();

        let idx = fi - frame_start;
        let t_frac = idx as f32 / (frame_count - 1).max(1) as f32;
        let scale = if expanding {
            min_scale + (1.0 - min_scale) * t_frac
        } else {
            1.0 - (1.0 - min_scale) * t_frac
        };

        let ring_binormal = normalize_safe(f.frame.tangent.cross(sheet_normal));
        let ring_normal = normalize_safe(ring_binormal.cross(f.frame.tangent));
        for j in 0..profile_len {
            let (lx, ly) = loop_profile.points[j];
            let (nx_l, ny_l) = loop_profile.normals[j];
            let world_pos = f.frame.position + ring_normal * (lx * scale) + ring_binormal * (ly * scale);
            let world_normal = normalize_safe(ring_normal * nx_l + ring_binormal * ny_l);
            vertices.push(MeshVertex {
                position: [world_pos.x, world_pos.y, world_pos.z],
                normal: [world_normal.x, world_normal.y, world_normal.z],
                color: f.color,
            });
        }
        ring_starts.push(ring_start);
    }

    for r in 0..ring_starts.len() - 1 {
        let r0 = ring_starts[r] as u32;
        let r1 = ring_starts[r + 1] as u32;
        let n_verts = profile_len as u32;
        for j in 0..n_verts {
            let j_next = (j + 1) % n_verts;
            indices.extend([
                r0 + j, r0 + j_next, r1 + j_next,
                r0 + j, r1 + j_next, r1 + j,
            ]);
        }
    }
}

/// Generate explicit sheet geometry for a contiguous sheet run.
///
/// Builds 4 quad strips (top, bottom, left side, right side) plus caps.
/// In fancy mode, the C-terminal end has an arrow head instead of a flat cap.
pub fn generate_explicit_sheet(
    vertices: &mut Vec<MeshVertex>,
    indices: &mut Vec<u32>,
    frames: &[FrameWithMetadata],
    sheet_start: usize,
    sheet_end: usize,
    settings: &CartoonGeometrySettings,
) {
    if sheet_start >= sheet_end || sheet_end >= frames.len() {
        return;
    }

    let half_width = settings.sheet_height; // Wide dimension (~1.4), in normal direction
    let half_thick = settings.sheet_width * 0.5; // Thin dimension (~0.2), in binormal direction

    // Determine arrow region
    let sheet_len = sheet_end - sheet_start + 1;
    let (body_end, has_arrow) = if settings.fancy_sheets && settings.arrow_length > 0 {
        let max_arrow = sheet_len / 2;
        let effective_len = settings.arrow_length.min(max_arrow);
        if effective_len >= 2 {
            (sheet_end + 1 - effective_len, true)
        } else {
            (sheet_end, false)
        }
    } else {
        (sheet_end, false)
    };

    // Helper: compute 4 corners for a frame at a given width
    let corners = |frame: &FrameWithMetadata, hw: f32| -> (Vec3, Vec3, Vec3, Vec3) {
        let n = frame.frame.normal;
        let b = frame.frame.binormal;
        let p = frame.frame.position;
        // TL = +normal, +binormal; TR = -normal, +binormal
        // BL = +normal, -binormal; BR = -normal, -binormal
        (
            p + n * hw + b * half_thick, // TL
            p - n * hw + b * half_thick, // TR
            p + n * hw - b * half_thick, // BL
            p - n * hw - b * half_thick, // BR
        )
    };

    // Helper: push a vertex
    let push_vert = |verts: &mut Vec<MeshVertex>, pos: Vec3, normal: Vec3, color: [f32; 4]| -> u32 {
        let idx = verts.len() as u32;
        verts.push(MeshVertex {
            position: [pos.x, pos.y, pos.z],
            normal: [normal.x, normal.y, normal.z],
            color,
        });
        idx
    };

    // Helper: push a quad (two CCW triangles) given 4 vertex indices
    // Vertices should be in order: v0-v1-v2-v3 going around the quad
    let push_quad = |inds: &mut Vec<u32>, v0: u32, v1: u32, v2: u32, v3: u32| {
        inds.extend([v0, v1, v2, v0, v2, v3]);
    };

    // ========== BODY (rectangular tube) ==========
    // Render frames from sheet_start to body_end as quad strips

    // Collect corner positions for body frames
    let body_frame_count = body_end - sheet_start + 1;
    if body_frame_count < 2 {
        return;
    }

    // Generate 4 quad strips for the body
    for fi in sheet_start..body_end {
        let f0 = &frames[fi];
        let f1 = &frames[fi + 1];
        let (tl0, tr0, bl0, br0) = corners(f0, half_width);
        let (tl1, tr1, bl1, br1) = corners(f1, half_width);
        let color0 = f0.color;
        let color1 = f1.color;

        // TOP surface (facing +binormal)
        let top_n0 = f0.frame.binormal;
        let top_n1 = f1.frame.binormal;
        let v0 = push_vert(vertices, tl0, top_n0, color0);
        let v1 = push_vert(vertices, tr0, top_n0, color0);
        let v2 = push_vert(vertices, tr1, top_n1, color1);
        let v3 = push_vert(vertices, tl1, top_n1, color1);
        push_quad(indices, v0, v1, v2, v3);

        // BOTTOM surface (facing -binormal)
        let bot_n0 = f0.frame.binormal * -1.0;
        let bot_n1 = f1.frame.binormal * -1.0;
        let v0 = push_vert(vertices, br0, bot_n0, color0);
        let v1 = push_vert(vertices, bl0, bot_n0, color0);
        let v2 = push_vert(vertices, bl1, bot_n1, color1);
        let v3 = push_vert(vertices, br1, bot_n1, color1);
        push_quad(indices, v0, v1, v2, v3);

        // LEFT surface (facing +normal)
        let left_n0 = f0.frame.normal;
        let left_n1 = f1.frame.normal;
        let v0 = push_vert(vertices, tl0, left_n0, color0);
        let v1 = push_vert(vertices, tl1, left_n1, color1);
        let v2 = push_vert(vertices, bl1, left_n1, color1);
        let v3 = push_vert(vertices, bl0, left_n0, color0);
        push_quad(indices, v0, v1, v2, v3);

        // RIGHT surface (facing -normal)
        let right_n0 = f0.frame.normal * -1.0;
        let right_n1 = f1.frame.normal * -1.0;
        let v0 = push_vert(vertices, tr1, right_n1, color1);
        let v1 = push_vert(vertices, tr0, right_n0, color0);
        let v2 = push_vert(vertices, br0, right_n0, color0);
        let v3 = push_vert(vertices, br1, right_n1, color1);
        push_quad(indices, v0, v1, v2, v3);
    }

    // ========== N-TERMINAL CONNECTOR ==========
    if sheet_start > 0 {
        let connector_start = sheet_start.saturating_sub(2);
        generate_sheet_connector(
            vertices, indices, frames,
            connector_start, sheet_start + 1,
            frames[sheet_start].frame.normal,
            settings, half_thick, false,
        );
    }

    // N-terminal cap (back face)
    {
        let f = &frames[sheet_start];
        let (tl, tr, bl, br) = corners(f, half_width);
        let cap_n = f.frame.tangent * -1.0;
        let color = f.color;
        let v0 = push_vert(vertices, tl, cap_n, color);
        let v1 = push_vert(vertices, bl, cap_n, color);
        let v2 = push_vert(vertices, br, cap_n, color);
        let v3 = push_vert(vertices, tr, cap_n, color);
        push_quad(indices, v0, v1, v2, v3);
    }

    // ========== C-TERMINAL CONNECTOR ==========
    if sheet_end + 1 < frames.len() {
        let connector_end = (sheet_end + 3).min(frames.len());
        generate_sheet_connector(
            vertices, indices, frames,
            sheet_end, connector_end,
            frames[sheet_end].frame.normal,
            settings, half_thick, true,
        );
    }

    if !has_arrow {
        // C-terminal cap (front face) — no arrow
        let f = &frames[sheet_end];
        let (tl, tr, bl, br) = corners(f, half_width);
        let cap_n = f.frame.tangent;
        let color = f.color;
        let v0 = push_vert(vertices, tl, cap_n, color);
        let v1 = push_vert(vertices, tr, cap_n, color);
        let v2 = push_vert(vertices, br, cap_n, color);
        let v3 = push_vert(vertices, bl, cap_n, color);
        push_quad(indices, v0, v1, v2, v3);

        return;
    }

    // ========== ARROW HEAD ==========
    // Arrow REPLACES the last few frames of the sheet — tip at frames[sheet_end].
    // Uses the base frame's normal/binormal (same as body) for orientation, ensuring
    // C0 continuity at the body-arrow junction. Positions are linearly interpolated
    // from base to tip so the straight spine matches the fixed orientation.

    let base = &frames[body_end];
    let arrow_width = half_width * settings.arrow_tip_scale;
    let arrow_len = sheet_end - body_end;

    // Use the base frame's orientation directly — same vectors the body uses via corners().
    // This guarantees the arrow starts with exactly the same orientation as the body end.
    let n = base.frame.normal;
    let b = base.frame.binormal;

    // Linearly interpolated spine from base to tip
    let start_pos = base.frame.position;
    let end_pos = frames[sheet_end].frame.position;

    for ai in 0..arrow_len {
        // Interpolation parameters along the arrow
        let t0 = ai as f32 / arrow_len as f32;
        let t1 = (ai + 1) as f32 / arrow_len as f32;

        // Linear taper: arrow_width at base -> 0 at tip
        let w0 = arrow_width * (1.0 - t0);
        let w1 = arrow_width * (1.0 - t1);

        // Linearly interpolated positions (straight spine)
        let pos0 = start_pos + (end_pos - start_pos) * t0;
        let pos1 = start_pos + (end_pos - start_pos) * t1;

        // Colors from actual frames
        let color0 = frames[body_end + ai].color;
        let color1 = frames[body_end + ai + 1].color;

        let tl0 = pos0 + n * w0 + b * half_thick;
        let tr0 = pos0 - n * w0 + b * half_thick;
        let bl0 = pos0 + n * w0 - b * half_thick;
        let br0 = pos0 - n * w0 - b * half_thick;
        let tl1 = pos1 + n * w1 + b * half_thick;
        let tr1 = pos1 - n * w1 + b * half_thick;
        let bl1 = pos1 + n * w1 - b * half_thick;
        let br1 = pos1 - n * w1 - b * half_thick;

        // TOP surface (facing +binormal)
        let v0 = push_vert(vertices, tl0, b, color0);
        let v1 = push_vert(vertices, tr0, b, color0);
        let v2 = push_vert(vertices, tr1, b, color1);
        let v3 = push_vert(vertices, tl1, b, color1);
        push_quad(indices, v0, v1, v2, v3);

        // BOTTOM surface (facing -binormal)
        let neg_b = b * -1.0;
        let v0 = push_vert(vertices, br0, neg_b, color0);
        let v1 = push_vert(vertices, bl0, neg_b, color0);
        let v2 = push_vert(vertices, bl1, neg_b, color1);
        let v3 = push_vert(vertices, br1, neg_b, color1);
        push_quad(indices, v0, v1, v2, v3);

        // LEFT taper edge (facing +normal)
        let v0 = push_vert(vertices, tl0, n, color0);
        let v1 = push_vert(vertices, tl1, n, color1);
        let v2 = push_vert(vertices, bl1, n, color1);
        let v3 = push_vert(vertices, bl0, n, color0);
        push_quad(indices, v0, v1, v2, v3);

        // RIGHT taper edge (facing -normal)
        let neg_n = n * -1.0;
        let v0 = push_vert(vertices, tr1, neg_n, color1);
        let v1 = push_vert(vertices, tr0, neg_n, color0);
        let v2 = push_vert(vertices, br0, neg_n, color0);
        let v3 = push_vert(vertices, br1, neg_n, color1);
        push_quad(indices, v0, v1, v2, v3);
    }

    // Arrow shoulder caps (back-face where arrow widens beyond body)
    {
        let cap_n = base.frame.tangent * -1.0;
        let color = base.color;
        let (body_tl, body_tr, body_bl, body_br) = corners(base, half_width);
        let base_pos = base.frame.position;
        let arrow_tl = base_pos + n * arrow_width + b * half_thick;
        let arrow_tr = base_pos - n * arrow_width + b * half_thick;
        let arrow_bl = base_pos + n * arrow_width - b * half_thick;
        let arrow_br = base_pos - n * arrow_width - b * half_thick;

        // Left shoulder
        let v0 = push_vert(vertices, body_tl, cap_n, color);
        let v1 = push_vert(vertices, arrow_tl, cap_n, color);
        let v2 = push_vert(vertices, arrow_bl, cap_n, color);
        let v3 = push_vert(vertices, body_bl, cap_n, color);
        push_quad(indices, v0, v1, v2, v3);

        // Right shoulder
        let v0 = push_vert(vertices, arrow_tr, cap_n, color);
        let v1 = push_vert(vertices, body_tr, cap_n, color);
        let v2 = push_vert(vertices, body_br, cap_n, color);
        let v3 = push_vert(vertices, arrow_br, cap_n, color);
        push_quad(indices, v0, v1, v2, v3);
    }

}

/// Add end cap to the tube
///
/// - `is_start`: true for the start of a segment, false for the end
/// - `is_arrow_tip`: true if this is at an arrow tip (skip cap since geometry tapers to a point)
/// - `dumbbell_taper`: taper factor for Flat2Face caps (0.0 = fully tapered, 1.0 = full size)
pub fn cap_tube_end(
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

    let profile = profile_for_ss(frame_meta.ss_type, settings);

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
    find_ss_regions(frames, |ss| ss == SecondaryStructure::Sheet)
}

/// Find helix region indices in a sequence of frames
pub fn find_helix_regions(frames: &[FrameWithMetadata]) -> Vec<(usize, usize)> {
    find_ss_regions(frames, is_helix)
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
pub fn generate_dumbbell_edge_tubes(
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
    fn test_find_sheet_termini() {
        let make_frame = |ss| FrameWithMetadata {
            frame: ReferenceFrame::new(
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
            ),
            color: [1.0, 1.0, 1.0, 1.0],
            ss_type: ss,

            segment_idx: 0,

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

            segment_idx: 0,

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
    fn test_find_helix_regions() {
        let make_frame = |ss| FrameWithMetadata {
            frame: ReferenceFrame::new(
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
            ),
            color: [1.0, 1.0, 1.0, 1.0],
            ss_type: ss,

            segment_idx: 0,

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

            segment_idx: 0,

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

            segment_idx: 0,

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
