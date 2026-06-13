//! patinae-render — GPU-first molecular renderer.
//!
//! Core renderer contracts include fixed representation bind groups,
//! indirect draws, settings-as-uniforms, and GPU-first compute criteria.
//!
//! Public API:
//! - `RenderState::new` / `resize` / `sync` / `render` / `pick`
//! - `RenderInput` + `RenderObjectInput` for pull-based host integration
//! - `RepKind`, `ObjectId`, `PickHit` for picking round-trips

pub mod capture;
#[allow(dead_code)]
pub(crate) mod compute;
mod context;
pub mod frame;
mod geometry_export;
pub(crate) mod lut_buffer;
mod map_contour;
pub(crate) mod passes;
pub mod picking;
pub(crate) mod pipelines;
pub(crate) mod postprocess;
mod render_input;
mod render_state;
#[allow(dead_code)]
pub(crate) mod representations;
pub mod scene_store;
mod shader_source;
#[cfg(feature = "stats")]
mod stats;
/// `FrameStatsHistory` is compiled unconditionally so hosts can surface
/// fps + percentile tail-latency in a HUD even when the `stats` cargo
/// feature (per-pass GPU timestamps) is off.
mod stats_history;
mod uniforms;

pub use context::RenderContext;
pub use frame::FrameTargets;
pub use geometry_export::{
    DisplayedGeometry, DisplayedMaterial, DisplayedMesh, DisplayedMeshVertex,
    DisplayedObjectGeometry, DisplayedPrimitive, GeometryExportError, GeometryExportOptions,
};
pub use picking::{ObjectId, PackedId, PickHit, PickingMode, RenderConfig, RepKind};
pub use render_input::{
    pack_rep_rgb8, ColorLutEntry, MarkerUpdate, RenderInput, RenderMapInput, RenderMapMode,
    RenderObjectInput, RepColorLutEntry, SceneLod, REP_COLOR_INHERIT,
};
pub use render_state::{RenderState, RenderSyncTimings};
pub use representations::sphere::SphereLodDiagnostics;
pub use representations::stick::StickLodDiagnostics;
pub use scene_store::SceneStoreFragmentationStats;
#[cfg(feature = "stats")]
pub use stats::{FrameStats, FrameStatsCollector};
pub use stats_history::FrameStatsHistory;
pub use uniforms::FrameUniforms;
