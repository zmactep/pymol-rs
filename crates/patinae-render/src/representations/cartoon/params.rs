//! `CartoonParams` — uniform shared by the three cartoon compute stages.
//!
//! Holds sampling, smoothing, profile geometry, and segment flags in a stable
//! binary layout so cached pipelines remain valid across settings changes.

use bytemuck::{Pod, Zeroable};
use patinae_settings::ResolvedSettings;

use crate::render_input::SceneLod;

/// Bit flags packed into [`CartoonParams::flags`].
pub mod flags {
    /// 0 → uniform tube of radius `tube_radius`. 1 → use SS-aware profiles
    /// where helix/sheet/loop classes choose separate cross-sections.
    pub const RESPECT_SS: u32 = 1 << 0;
    pub const FANCY_HELICES: u32 = 1 << 1;
    pub const FANCY_SHEETS: u32 = 1 << 2;
    pub const ROUND_HELICES: u32 = 1 << 3;
}

/// Default sample count when `cartoon_sampling = -1`.
pub const DEFAULT_SAMPLES_PER_RESIDUE: u32 = 7;

/// Number of vertices around the tube profile. Default for `SceneLod::High`
/// (and any explicit `cartoon_quality = high`). Larger structures auto-
/// downgrade through `lod_overrides` below.
pub const PROFILE_SEGMENTS: u32 = 12;

/// Sampling / profile downgrade applied when the host classifies the
/// scene size. Returns `(samples_cap, profile_segments)` — the sampling
/// half is a *cap* over the user's `cartoon_sampling`, not a replacement,
/// so users who explicitly request lower sampling still get it.
fn lod_overrides(lod: SceneLod) -> (u32, u32) {
    match lod {
        // `Auto` and `High` both go through the high-quality path; `Auto`
        // means "host did not classify the scene", in which case we must
        // not silently degrade — leave per-rep heuristics alone.
        SceneLod::Auto | SceneLod::High => (u32::MAX, PROFILE_SEGMENTS),
        SceneLod::Medium => (8, 10),
        SceneLod::Low => (6, 8),
        SceneLod::Minimum => (4, 6),
    }
}

/// Loop / coil tube radius in Angstroms.
pub const DEFAULT_LOOP_RADIUS: f32 = 0.2;

/// Default tube radius for the `cartoon_uniform_tube` mode.
pub const DEFAULT_TUBE_RADIUS: f32 = 0.5;

/// Helix oval **wide** half-axis — placed along the **binormal** in
/// `cartoon_extrude.wgsl`.
/// With round_helix orientation = axis × tangent, the binormal lands along
/// the helix axis, so the wide flat face stays glued to the axis and the
/// ribbon reads as continuous (not as the per-residue spiral we get when
/// the wide face is radial).
pub const DEFAULT_HELIX_RADIUS: f32 = 1.35;

/// Helix oval **narrow** half-axis — placed along the **normal** in the
/// extrusion shader.
pub const DEFAULT_HELIX_WIDTH: f32 = 0.25;

/// Sheet rectangle **wide** half-extent — placed along the **binormal**
/// (= in-plane direction perpendicular to the strand tangent). Mirrors
/// `cartoon_rect_length` (1.4).
pub const DEFAULT_SHEET_RADIUS: f32 = 1.4;

/// Sheet rectangle **narrow** half-extent — placed along the **normal**
/// (= perpendicular to the strand plane). Mirrors `cartoon_rect_width` (0.4).
pub const DEFAULT_SHEET_WIDTH: f32 = 0.4;

/// Multiplier on `sheet_radius` for the leading edge of an arrow head.
/// The arrow's wide barb sits at
/// `sheet_radius * arrow_tip_scale` along the binormal; the trailing tip
/// tapers linearly to 0.
pub const DEFAULT_ARROW_TIP_SCALE: f32 = 1.5;

/// Mirrors `CartoonParams` in `cartoon_*.wgsl`. 64 B (4×16 B chunks). Field
/// order and padding must stay in sync — tests assert the layout.
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Pod, Zeroable)]
pub struct CartoonParams {
    /// Number of `BackboneAtom` records driving the chain (= visible Cαs).
    pub n_atoms: u32,
    /// Samples per residue. Always ≥ 1; comes from `cartoon_sampling`
    /// (with `-1 → DEFAULT_SAMPLES_PER_RESIDUE`).
    pub n_samples_per: u32,
    /// Number of Laplacian smoothing iterations (`cartoon_smooth_cycles`).
    pub smooth_cycles: u32,
    /// Vertices around each tube cross-section.
    pub profile_segments: u32,

    pub loop_radius: f32,
    pub tube_radius: f32,
    /// Helix oval **wide** half-axis (placed along the binormal in the extrude
    /// shader, which lands along the helix axis after round-helix orientation
    /// fixup).
    pub helix_radius: f32,
    /// Sheet rectangle **wide** half-extent (along binormal = in-plane).
    pub sheet_radius: f32,

    /// See [`flags`].
    pub flags: u32,
    /// Global cartoon alpha multiplier (1 - `cartoon_transparency`). Sent
    /// to the cartoon FS via `CartoonParams.alpha_mul`; this slot in the
    /// extrude params is informational only and unused by the compute chain.
    pub alpha_mul: f32,
    /// Helix oval **narrow** half-axis (along the normal).
    pub helix_width: f32,
    /// Sheet rectangle **narrow** half-extent (along normal = out-of-plane).
    pub sheet_width: f32,

    /// Number of axial-Laplacian smoothing cycles run after sampling.
    /// Current GPU code reads it host-side and dispatches that many times.
    pub refine_cycles: u32,
    /// Multiplier on `sheet_radius` for the leading edge of a sheet-arrow
    /// terminus.
    pub arrow_tip_scale: f32,
    // 8 B reservation for future geometry knobs.
    pub _reserved2: u32,
    pub _reserved3: u32,
}

impl CartoonParams {
    pub const SIZE: u64 = std::mem::size_of::<Self>() as u64;
}

const _: () = assert!(std::mem::size_of::<CartoonParams>() == 64);

/// Build `CartoonParams` from resolved settings. The caller supplies `n_atoms`
/// (the count of `BackboneAtom` records — already filtered to visible
/// residues). `respect_ss` lets the host force the legacy "uniform tube"
/// behaviour when the user has it disabled.
pub fn pack_from_settings(
    settings: &ResolvedSettings,
    n_atoms: u32,
    respect_ss: bool,
    lod: SceneLod,
) -> CartoonParams {
    let (sample_cap, profile_segments) = lod_overrides(lod);
    let sampling = settings.cartoon.sampling;
    let raw_samples = if sampling <= 0 {
        DEFAULT_SAMPLES_PER_RESIDUE
    } else {
        sampling as u32
    };
    // LOD only acts as a *cap*: if the user explicitly asked for fewer
    // samples we honour that; if the resolved default is higher than the
    // bucket allows we clamp.
    let n_samples_per = raw_samples.min(sample_cap).max(1);
    let smooth_cycles = settings.cartoon.smooth_cycles.max(0) as u32;
    // Round refine_cycles up to even so the GPU ping-pong lands the final
    // result in the canonical `samples_buf` (= the host's `cartoon_finalize_frames`
    // input). Default 5 → 6; cost is negligible (~0.5 ms on 4hhb).
    let raw_refine = settings.cartoon.refine.max(0) as u32;
    let refine_cycles = raw_refine + (raw_refine & 1u32);
    let alpha_mul = (1.0 - settings.cartoon.transparency).clamp(0.0, 1.0);
    let mut flags_word: u32 = 0;
    if respect_ss {
        flags_word |= flags::RESPECT_SS;
    }
    if settings.cartoon.round_helices {
        flags_word |= flags::ROUND_HELICES;
    }
    CartoonParams {
        n_atoms,
        n_samples_per,
        smooth_cycles,
        profile_segments,
        loop_radius: DEFAULT_LOOP_RADIUS,
        tube_radius: DEFAULT_TUBE_RADIUS,
        helix_radius: DEFAULT_HELIX_RADIUS,
        sheet_radius: DEFAULT_SHEET_RADIUS,
        flags: flags_word,
        alpha_mul,
        helix_width: DEFAULT_HELIX_WIDTH,
        sheet_width: DEFAULT_SHEET_WIDTH,
        refine_cycles,
        arrow_tip_scale: DEFAULT_ARROW_TIP_SCALE,
        _reserved2: 0,
        _reserved3: 0,
    }
}

/// Build `CartoonParams` from resolved settings for the **ribbon**
/// representation. Ribbon is a uniform-tube cartoon: all SS classes render as
/// the same circular cross-section of `ribbon_radius`, no helix oval, no sheet
/// rectangle, no arrow taper. The shader path falls through `RESPECT_SS = 0`.
///
/// The sampling / smoothing pipeline is unchanged — ribbon reuses cartoon's
/// smooth → sample → refine → finalize → extrude compute chain. Ribbon-specific
/// `ribbon_sampling` / `ribbon_radius` settings replace the cartoon defaults;
/// `cartoon_smooth_cycles` and `cartoon_refine` still drive the smoothing.
pub fn pack_from_ribbon_settings(
    settings: &ResolvedSettings,
    n_atoms: u32,
    lod: SceneLod,
) -> CartoonParams {
    let (sample_cap, profile_segments) = lod_overrides(lod);
    // `ribbon_sampling <= 0` uses the default 10 samples per residue;
    // otherwise values below 7 are raised because coarser tubes read as
    // polylines. LOD caps both values and may deliberately go lower for
    // assembly-scale scenes.
    let raw_sampling = settings.ribbon.sampling;
    let raw_samples = if raw_sampling <= 0 {
        10u32
    } else {
        (raw_sampling as u32).max(7)
    };
    let n_samples_per = raw_samples.min(sample_cap).max(1);
    // `ribbon_radius = 0.0` means "auto" -> 0.3 Angstroms.
    let radius = if settings.ribbon.radius > 0.0 {
        settings.ribbon.radius
    } else {
        0.3
    };

    // Ribbon shares cartoon's smoothing config; the chain is the same Laplacian
    // ping-pong + Jacobi refine the cartoon path uses.
    let smooth_cycles = settings.cartoon.smooth_cycles.max(0) as u32;
    let raw_refine = settings.cartoon.refine.max(0) as u32;
    let refine_cycles = raw_refine + (raw_refine & 1u32);

    CartoonParams {
        n_atoms,
        n_samples_per,
        smooth_cycles,
        profile_segments,
        loop_radius: radius,
        tube_radius: radius,
        helix_radius: radius,
        sheet_radius: radius,
        // RESPECT_SS off → every sample renders as a uniform tube of
        // `tube_radius`. ROUND_HELICES is irrelevant when SS is ignored.
        flags: 0,
        // Ribbon transparency is not in the typed settings yet. When it
        // lands, multiply it into alpha just like cartoon does.
        alpha_mul: 1.0,
        helix_width: radius,
        sheet_width: radius,
        refine_cycles,
        arrow_tip_scale: 1.0,
        _reserved2: 0,
        _reserved3: 0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use patinae_settings::{groups::Settings, ResolvedSettings};

    #[test]
    fn layout_is_64_bytes() {
        assert_eq!(std::mem::size_of::<CartoonParams>(), 64);
        // 4-byte alignment is fine for std140 uniform upload — the WGSL side
        // gets 16-byte alignment automatically because all fields are u32/f32.
        assert_eq!(std::mem::align_of::<CartoonParams>(), 4);
    }

    #[test]
    fn default_sampling_resolves_to_seven() {
        let global = Settings::default();
        let settings = ResolvedSettings::resolve(&global, None);
        let p = pack_from_settings(&settings, 100, false, SceneLod::Auto);
        assert_eq!(p.n_samples_per, DEFAULT_SAMPLES_PER_RESIDUE);
        assert_eq!(p.profile_segments, PROFILE_SEGMENTS);
        assert_eq!(p.n_atoms, 100);
        assert!((p.alpha_mul - 1.0).abs() < 1e-6);
        assert_eq!(p.flags & flags::RESPECT_SS, 0);
        // Default cartoon.round_helices = true → bit set.
        assert_ne!(p.flags & flags::ROUND_HELICES, 0);
    }

    #[test]
    fn explicit_sampling_overrides_default() {
        let mut global = Settings::default();
        global.cartoon.sampling = 4;
        let settings = ResolvedSettings::resolve(&global, None);
        let p = pack_from_settings(&settings, 1, false, SceneLod::Auto);
        assert_eq!(p.n_samples_per, 4);
    }

    #[test]
    fn transparency_propagates_to_alpha_mul() {
        let mut global = Settings::default();
        global.cartoon.transparency = 0.25;
        let settings = ResolvedSettings::resolve(&global, None);
        let p = pack_from_settings(&settings, 1, false, SceneLod::Auto);
        assert!((p.alpha_mul - 0.75).abs() < 1e-6);
    }
}
