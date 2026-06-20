//! Renderer memory profile selection.
//!
//! The types in this module describe construction-time memory policy. They
//! intentionally stop before representation budgeting: profile-aware rep
//! downgrade, chunking, and skip decisions belong to a later planning phase.

use std::fmt;
use std::str::FromStr;

use crate::byte_units::{gib_to_bytes, mib_to_bytes, BYTES_PER_MIB};

const LOW_MEMORY_BUFFER_THRESHOLD: u64 = gib_to_bytes(1);
const LOW_MEMORY_STORAGE_BINDING_THRESHOLD: u32 = mib_to_bytes(512) as u32;
/// Desired single-buffer target for capable native adapters.
pub const PERFORMANCE_MAX_BUFFER_SIZE: u64 = gib_to_bytes(4);
/// Desired single storage-buffer binding target for capable adapters.
pub const PERFORMANCE_MAX_STORAGE_BUFFER_BINDING_SIZE: u32 = (gib_to_bytes(2) - 1) as u32;

const DEFAULT_PICKING_SCALE: f32 = 0.5;
const LOW_MEMORY_PICKING_SCALE: f32 = 0.25;
const PERFORMANCE_SHADOW_MAP_SIZE: u32 = 4096;
const BALANCED_SHADOW_MAP_SIZE: u32 = 2048;
const LOW_MEMORY_SHADOW_MAP_SIZE: u32 = 1024;
const PERFORMANCE_ATLAS_DIRECTIONS: u32 = 256;
const BALANCED_ATLAS_DIRECTIONS: u32 = 64;
const LOW_MEMORY_ATLAS_DIRECTIONS: u32 = 16;

/// User-visible renderer memory profile.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RenderMemoryProfile {
    /// Preserve current high-quality, eager-allocation behavior.
    Performance,
    /// Prefer lower baseline memory without disabling normal interaction.
    Balanced,
    /// Disable or deny optional full-resolution screen resources.
    LowMemory,
    /// Apply profile policy from an explicit byte budget.
    Budgeted { bytes: u64 },
}

impl RenderMemoryProfile {
    /// Returns the configured budget, if this profile is budgeted.
    pub const fn budget_bytes(self) -> Option<u64> {
        match self {
            Self::Budgeted { bytes } => Some(bytes),
            Self::Performance | Self::Balanced | Self::LowMemory => None,
        }
    }
}

impl fmt::Display for RenderMemoryProfile {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Performance => f.write_str("performance"),
            Self::Balanced => f.write_str("balanced"),
            Self::LowMemory => f.write_str("low"),
            Self::Budgeted { bytes } => write!(f, "low:{}MiB", bytes / BYTES_PER_MIB),
        }
    }
}

impl FromStr for RenderMemoryProfile {
    type Err = RenderMemoryProfileParseError;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        let value = value.trim().to_ascii_lowercase();
        match value.as_str() {
            "performance" | "perf" => return Ok(Self::Performance),
            "balanced" | "balance" => return Ok(Self::Balanced),
            "low" | "low-memory" | "low_memory" => return Ok(Self::LowMemory),
            _ => {}
        }

        let Some((kind, mib)) = value.split_once(':') else {
            return Err(RenderMemoryProfileParseError);
        };
        if kind != "low" && kind != "budget" && kind != "budgeted" {
            return Err(RenderMemoryProfileParseError);
        }
        let mib = mib
            .parse::<u64>()
            .map_err(|_| RenderMemoryProfileParseError)?;
        if mib == 0 {
            return Err(RenderMemoryProfileParseError);
        }
        Ok(Self::Budgeted {
            bytes: mib_to_bytes(mib),
        })
    }
}

/// Error returned when a memory profile override is invalid.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RenderMemoryProfileParseError;

impl fmt::Display for RenderMemoryProfileParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("expected performance, balanced, low, or low:<MiB>")
    }
}

impl std::error::Error for RenderMemoryProfileParseError {}

/// Hit-test picking allocation policy.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PickingPolicy {
    /// Whether hit-test picking resources may be allocated.
    pub hit_test_enabled: bool,
    /// Fraction of the viewport used for hit-test picking textures.
    pub scale: f32,
    /// Whether the reprojection compute path may be allocated.
    pub reprojection_enabled: bool,
}

/// Per-frame target allocation policy.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FrameTargetPolicy {
    /// Allocate SSAO textures during renderer construction.
    pub ssao_eager: bool,
    /// Allocate color scratch during renderer construction.
    pub color_scratch_eager: bool,
    /// Number of native host textures kept in the viewport handoff ring.
    pub viewport_handoff_ring: usize,
}

/// Visual overlay allocation policy.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct OverlayPolicy {
    /// Whether selection and hover marking overlays may allocate resources.
    pub selection_enabled: bool,
    /// Whether silhouette overlays may allocate resources.
    pub silhouette_enabled: bool,
}

/// Full-screen postprocess policy.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PostprocessPolicy {
    /// Whether SSAO may allocate and run.
    pub ssao_enabled: bool,
    /// Whether FXAA may allocate and run.
    pub fxaa_enabled: bool,
}

/// Shadow and atlas resource caps.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ShadowPolicy {
    /// Maximum directional shadow-map size.
    pub max_shadow_map_size: u32,
    /// Maximum atlas-AO tile size.
    pub max_atlas_tile_size: u32,
    /// Maximum atlas-AO direction count.
    pub max_atlas_directions: u32,
}

/// Representation allocation policy.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RepresentationBudgetPolicy {
    /// Whether representation downgrade/skip decisions are active.
    pub enabled: bool,
}

/// Construction-time renderer memory policy.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RenderMemoryPolicy {
    /// Selected profile.
    pub profile: RenderMemoryProfile,
    /// Explicit byte budget for budgeted profiles.
    pub budget_bytes: Option<u64>,
    /// Hit-test picking resource policy.
    pub picking: PickingPolicy,
    /// Frame target resource policy.
    pub frame_targets: FrameTargetPolicy,
    /// Visual overlay policy.
    pub overlays: OverlayPolicy,
    /// Postprocess policy.
    pub postprocess: PostprocessPolicy,
    /// Shadow and atlas caps.
    pub shadows: ShadowPolicy,
    /// Future representation-budgeting policy.
    pub reps: RepresentationBudgetPolicy,
}

impl RenderMemoryPolicy {
    /// Returns the current high-quality policy.
    pub const fn performance() -> Self {
        Self {
            profile: RenderMemoryProfile::Performance,
            budget_bytes: None,
            picking: PickingPolicy {
                hit_test_enabled: true,
                scale: DEFAULT_PICKING_SCALE,
                reprojection_enabled: true,
            },
            frame_targets: FrameTargetPolicy {
                ssao_eager: true,
                color_scratch_eager: true,
                viewport_handoff_ring: 3,
            },
            overlays: OverlayPolicy {
                selection_enabled: true,
                silhouette_enabled: true,
            },
            postprocess: PostprocessPolicy {
                ssao_enabled: true,
                fxaa_enabled: true,
            },
            shadows: ShadowPolicy {
                max_shadow_map_size: PERFORMANCE_SHADOW_MAP_SIZE,
                max_atlas_tile_size: PERFORMANCE_SHADOW_MAP_SIZE,
                max_atlas_directions: PERFORMANCE_ATLAS_DIRECTIONS,
            },
            reps: RepresentationBudgetPolicy { enabled: false },
        }
    }

    /// Returns a lower-baseline-memory policy.
    pub const fn balanced() -> Self {
        Self {
            profile: RenderMemoryProfile::Balanced,
            budget_bytes: None,
            picking: PickingPolicy {
                hit_test_enabled: true,
                scale: DEFAULT_PICKING_SCALE,
                reprojection_enabled: true,
            },
            frame_targets: FrameTargetPolicy {
                ssao_eager: false,
                color_scratch_eager: false,
                viewport_handoff_ring: 2,
            },
            overlays: OverlayPolicy {
                selection_enabled: true,
                silhouette_enabled: true,
            },
            postprocess: PostprocessPolicy {
                ssao_enabled: true,
                fxaa_enabled: true,
            },
            shadows: ShadowPolicy {
                max_shadow_map_size: BALANCED_SHADOW_MAP_SIZE,
                max_atlas_tile_size: 1024,
                max_atlas_directions: BALANCED_ATLAS_DIRECTIONS,
            },
            reps: RepresentationBudgetPolicy { enabled: true },
        }
    }

    /// Returns the constrained low-memory policy.
    pub const fn low_memory() -> Self {
        Self {
            profile: RenderMemoryProfile::LowMemory,
            budget_bytes: None,
            picking: PickingPolicy {
                hit_test_enabled: true,
                scale: LOW_MEMORY_PICKING_SCALE,
                reprojection_enabled: false,
            },
            frame_targets: FrameTargetPolicy {
                ssao_eager: false,
                color_scratch_eager: false,
                viewport_handoff_ring: 1,
            },
            overlays: OverlayPolicy {
                selection_enabled: false,
                silhouette_enabled: false,
            },
            postprocess: PostprocessPolicy {
                ssao_enabled: false,
                fxaa_enabled: false,
            },
            shadows: ShadowPolicy {
                max_shadow_map_size: LOW_MEMORY_SHADOW_MAP_SIZE,
                max_atlas_tile_size: 512,
                max_atlas_directions: LOW_MEMORY_ATLAS_DIRECTIONS,
            },
            reps: RepresentationBudgetPolicy { enabled: true },
        }
    }

    /// Returns the constrained policy with an explicit representation budget.
    pub const fn budgeted(bytes: u64) -> Self {
        let mut policy = Self::low_memory();
        policy.profile = RenderMemoryProfile::Budgeted { bytes };
        policy.budget_bytes = Some(bytes);
        policy.reps.enabled = true;
        policy
    }

    /// Returns policy for a profile.
    pub const fn from_profile(profile: RenderMemoryProfile) -> Self {
        match profile {
            RenderMemoryProfile::Performance => Self::performance(),
            RenderMemoryProfile::Balanced => Self::balanced(),
            RenderMemoryProfile::LowMemory => Self::low_memory(),
            RenderMemoryProfile::Budgeted { bytes } => Self::budgeted(bytes),
        }
    }

    /// Returns the preferred `wgpu` memory hint for this policy.
    pub const fn wgpu_memory_hints(self) -> wgpu::MemoryHints {
        match self.profile {
            RenderMemoryProfile::Performance => wgpu::MemoryHints::Performance,
            RenderMemoryProfile::Balanced
            | RenderMemoryProfile::LowMemory
            | RenderMemoryProfile::Budgeted { .. } => wgpu::MemoryHints::MemoryUsage,
        }
    }
}

impl Default for RenderMemoryPolicy {
    fn default() -> Self {
        Self::performance()
    }
}

/// Adapter class used by policy selection.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RenderAdapterType {
    /// Unknown or uncategorized adapter.
    Other,
    /// Integrated GPU.
    IntegratedGpu,
    /// Discrete GPU.
    DiscreteGpu,
    /// Virtual GPU.
    VirtualGpu,
    /// CPU fallback adapter.
    Cpu,
}

impl From<wgpu::DeviceType> for RenderAdapterType {
    fn from(value: wgpu::DeviceType) -> Self {
        match value {
            wgpu::DeviceType::IntegratedGpu => Self::IntegratedGpu,
            wgpu::DeviceType::DiscreteGpu => Self::DiscreteGpu,
            wgpu::DeviceType::VirtualGpu => Self::VirtualGpu,
            wgpu::DeviceType::Cpu => Self::Cpu,
            wgpu::DeviceType::Other => Self::Other,
        }
    }
}

/// Backend class used by policy selection.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RenderBackend {
    /// Native Vulkan backend.
    Vulkan,
    /// Native Metal backend.
    Metal,
    /// Native DirectX 12 backend.
    Dx12,
    /// Browser WebGPU backend.
    BrowserWebGpu,
    /// Browser WebGL backend.
    BrowserWebGl,
    /// OpenGL backend.
    Gl,
    /// Unknown backend.
    Other,
}

impl From<wgpu::Backend> for RenderBackend {
    fn from(value: wgpu::Backend) -> Self {
        match value {
            wgpu::Backend::Vulkan => Self::Vulkan,
            wgpu::Backend::Metal => Self::Metal,
            wgpu::Backend::Dx12 => Self::Dx12,
            wgpu::Backend::BrowserWebGpu => Self::BrowserWebGpu,
            wgpu::Backend::Gl => Self::Gl,
            _ => Self::Other,
        }
    }
}

/// Synthetic adapter inputs for profile selection.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RenderMemorySelectionInput {
    /// Adapter type.
    pub adapter_type: RenderAdapterType,
    /// Adapter backend.
    pub backend: RenderBackend,
    /// Adapter maximum buffer size.
    pub max_buffer_size: u64,
    /// Adapter maximum storage-buffer binding size.
    pub max_storage_buffer_binding_size: u32,
    /// Adapter maximum two-dimensional texture dimension.
    pub max_texture_dimension_2d: u32,
    /// Whether this is a browser-hosted renderer.
    pub is_web: bool,
    /// Whether prior runtime state requested a downgrade.
    pub observed_downgrade: bool,
}

impl RenderMemorySelectionInput {
    /// Builds selector input from `wgpu` adapter data.
    pub fn from_wgpu(info: &wgpu::AdapterInfo, limits: &wgpu::Limits, is_web: bool) -> Self {
        Self {
            adapter_type: info.device_type.into(),
            backend: info.backend.into(),
            max_buffer_size: limits.max_buffer_size,
            max_storage_buffer_binding_size: limits.max_storage_buffer_binding_size,
            max_texture_dimension_2d: limits.max_texture_dimension_2d,
            is_web,
            observed_downgrade: false,
        }
    }
}

/// Selects a memory policy from adapter facts and an optional override.
pub fn select_render_memory_policy(
    input: RenderMemorySelectionInput,
    override_profile: Option<RenderMemoryProfile>,
) -> RenderMemoryPolicy {
    if let Some(profile) = override_profile {
        return RenderMemoryPolicy::from_profile(profile);
    }
    if input.observed_downgrade
        || input.adapter_type == RenderAdapterType::Cpu
        || input.adapter_type == RenderAdapterType::VirtualGpu
        || input.max_buffer_size < LOW_MEMORY_BUFFER_THRESHOLD
        || input.max_storage_buffer_binding_size < LOW_MEMORY_STORAGE_BINDING_THRESHOLD
        || input.max_texture_dimension_2d < 8192
    {
        return RenderMemoryPolicy::low_memory();
    }
    if input.is_web
        || matches!(
            input.backend,
            RenderBackend::BrowserWebGpu | RenderBackend::BrowserWebGl
        )
        || input.max_buffer_size < PERFORMANCE_MAX_BUFFER_SIZE
        || input.max_storage_buffer_binding_size < PERFORMANCE_MAX_STORAGE_BUFFER_BINDING_SIZE
    {
        return RenderMemoryPolicy::balanced();
    }
    RenderMemoryPolicy::performance()
}

/// Builds `wgpu` required limits for a memory policy.
pub fn required_limits_for_memory_policy(
    adapter_limits: &wgpu::Limits,
    _policy: RenderMemoryPolicy,
) -> wgpu::Limits {
    wgpu::Limits {
        max_buffer_size: PERFORMANCE_MAX_BUFFER_SIZE.min(adapter_limits.max_buffer_size),
        max_storage_buffer_binding_size: PERFORMANCE_MAX_STORAGE_BUFFER_BINDING_SIZE
            .min(adapter_limits.max_storage_buffer_binding_size),
        ..wgpu::Limits::default()
    }
    .using_resolution(adapter_limits.clone())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::byte_units::{gib_to_bytes, mib_to_bytes};

    fn limits(max_buffer_size: u64, max_storage_buffer_binding_size: u32) -> wgpu::Limits {
        wgpu::Limits {
            max_buffer_size,
            max_storage_buffer_binding_size,
            max_texture_dimension_2d: 16_384,
            ..wgpu::Limits::default()
        }
    }

    fn input(
        adapter_type: RenderAdapterType,
        max_buffer_size: u64,
        max_storage_buffer_binding_size: u32,
    ) -> RenderMemorySelectionInput {
        RenderMemorySelectionInput {
            adapter_type,
            backend: RenderBackend::Metal,
            max_buffer_size,
            max_storage_buffer_binding_size,
            max_texture_dimension_2d: 16_384,
            is_web: false,
            observed_downgrade: false,
        }
    }

    #[test]
    fn profile_override_parses_expected_values() {
        assert_eq!(
            "performance".parse::<RenderMemoryProfile>(),
            Ok(RenderMemoryProfile::Performance)
        );
        assert_eq!(
            "balanced".parse::<RenderMemoryProfile>(),
            Ok(RenderMemoryProfile::Balanced)
        );
        assert_eq!(
            "low".parse::<RenderMemoryProfile>(),
            Ok(RenderMemoryProfile::LowMemory)
        );
        assert_eq!(
            "low:1024".parse::<RenderMemoryProfile>(),
            Ok(RenderMemoryProfile::Budgeted {
                bytes: gib_to_bytes(1)
            })
        );
    }

    #[test]
    fn high_end_adapter_selects_performance() {
        let policy = select_render_memory_policy(
            input(RenderAdapterType::IntegratedGpu, gib_to_bytes(8), u32::MAX),
            None,
        );

        assert_eq!(policy.profile, RenderMemoryProfile::Performance);
        assert_eq!(policy.frame_targets.viewport_handoff_ring, 3);
        assert!(policy.postprocess.fxaa_enabled);
        assert!(policy.postprocess.ssao_enabled);
    }

    #[test]
    fn moderate_adapter_selects_balanced() {
        let policy = select_render_memory_policy(
            input(
                RenderAdapterType::IntegratedGpu,
                gib_to_bytes(2),
                gib_to_bytes(1) as u32,
            ),
            None,
        );

        assert_eq!(policy.profile, RenderMemoryProfile::Balanced);
        assert!(!policy.frame_targets.ssao_eager);
        assert!(!policy.frame_targets.color_scratch_eager);
        assert_eq!(policy.frame_targets.viewport_handoff_ring, 2);
    }

    #[test]
    fn constrained_adapter_selects_low_memory() {
        let policy = select_render_memory_policy(
            input(
                RenderAdapterType::IntegratedGpu,
                mib_to_bytes(512),
                mib_to_bytes(256) as u32,
            ),
            None,
        );

        assert_eq!(policy.profile, RenderMemoryProfile::LowMemory);
        assert!(policy.picking.hit_test_enabled);
        assert!(!policy.picking.reprojection_enabled);
        assert!(!policy.postprocess.fxaa_enabled);
        assert!(!policy.overlays.selection_enabled);
        assert_eq!(policy.frame_targets.viewport_handoff_ring, 1);
    }

    #[test]
    fn forced_low_memory_changes_policy_without_weak_hardware() {
        let policy = select_render_memory_policy(
            input(RenderAdapterType::DiscreteGpu, gib_to_bytes(8), u32::MAX),
            Some(RenderMemoryProfile::LowMemory),
        );

        assert_eq!(policy.profile, RenderMemoryProfile::LowMemory);
        assert!(!policy.postprocess.ssao_enabled);
        assert!(policy.reps.enabled);
    }

    #[test]
    fn performance_keeps_rep_budgeting_pass_through() {
        let policy = RenderMemoryPolicy::performance();

        assert!(!policy.reps.enabled);
    }

    #[test]
    fn budgeted_profile_enforces_rep_budgeting_even_on_large_budget() {
        let policy = RenderMemoryPolicy::budgeted(gib_to_bytes(8));

        assert_eq!(
            policy.profile,
            RenderMemoryProfile::Budgeted {
                bytes: gib_to_bytes(8)
            }
        );
        assert!(policy.reps.enabled);
        assert!(!policy.postprocess.fxaa_enabled);
        assert!(!policy.postprocess.ssao_enabled);
        assert!(policy.picking.hit_test_enabled);
        assert!(!policy.picking.reprojection_enabled);
        assert!(!policy.overlays.selection_enabled);
    }

    #[test]
    fn required_limits_never_exceed_adapter_limits() {
        let adapter = limits(mib_to_bytes(512), mib_to_bytes(384) as u32);
        let required =
            required_limits_for_memory_policy(&adapter, RenderMemoryPolicy::low_memory());

        assert!(required.max_buffer_size <= adapter.max_buffer_size);
        assert!(
            required.max_storage_buffer_binding_size <= adapter.max_storage_buffer_binding_size
        );
        assert_eq!(required.max_texture_dimension_2d, 16_384);
    }

    #[test]
    fn required_limits_preserve_performance_targets_on_capable_adapter() {
        let adapter = limits(gib_to_bytes(8), u32::MAX);
        let required =
            required_limits_for_memory_policy(&adapter, RenderMemoryPolicy::performance());

        assert_eq!(required.max_buffer_size, PERFORMANCE_MAX_BUFFER_SIZE);
        assert_eq!(
            required.max_storage_buffer_binding_size,
            PERFORMANCE_MAX_STORAGE_BUFFER_BINDING_SIZE
        );
    }
}
