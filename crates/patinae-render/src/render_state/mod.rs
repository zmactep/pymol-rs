//! `RenderState` — public entry point for the crate.
//!
//! It owns the long-lived GPU resources and coordinates scene sync,
//! per-frame compute, draw passes, postprocess, and picking.

mod artifact_flow;
mod compute_build_flow;
mod construction;
mod culling_flow;
mod frame_flow;
mod geometry_export_flow;
mod math;
mod picking_flow;
mod resize;
mod screen_flow;
mod settings;
mod shadow_flow;
pub(crate) mod state;
mod sync_flow;
mod visible_flow;

pub use state::{RenderState, RenderSyncTimings};
