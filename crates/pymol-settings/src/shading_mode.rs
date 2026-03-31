//! Shading mode abstraction.
//!
//! Each shading mode defines its own lighting pipeline and parameters,
//! ensuring that switching between modes has no side effects.

use serde::{Deserialize, Serialize};

use crate::impl_setting_enum;

/// Available shading modes.
///
/// Each mode reads only its own set of parameters:
/// - **Classic**: multi-light PyMOL model (ambient, direct, reflect, specular, lights).
/// - **Skripkin**: pure ambient + multi-directional shadow AO (own ambient/specular/shininess).
/// - **Full**: Classic lighting + per-light directional shadow maps.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[repr(i32)]
pub enum ShadingMode {
    /// Classic PyMOL lighting: multi-light, ambient + direct + reflect + specular.
    Classic = 0,
    /// Skripkin mode: pure ambient + multi-directional shadow AO, no directional lights.
    Skripkin = 1,
    /// Full mode: Classic multi-light lighting + per-light directional shadow maps.
    Full = 2,
}

impl Default for ShadingMode {
    fn default() -> Self {
        Self::Classic
    }
}

impl_setting_enum! {
    ShadingMode {
        Classic = 0 => "classic",
        Skripkin = 1 => "skripkin",
        Full = 2 => "full",
    }
    default: Classic
}

