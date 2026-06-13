//! Shading mode abstraction.
//!
//! Each shading mode defines its own lighting pipeline and parameters,
//! ensuring that switching between modes has no side effects.

use serde::{Deserialize, Serialize};

use crate::impl_setting_enum;

/// Available shading modes.
///
/// Each mode reads only its own set of parameters:
/// - **Classic**: multi-light model with ambient, direct, reflect, and specular terms.
/// - **Skripkin**: pure ambient lighting with AO-controlled contrast.
/// - **Full**: Classic lighting + a directional shadow map.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default, Serialize, Deserialize)]
#[repr(i32)]
pub enum ShadingMode {
    /// Classic multi-light shading: ambient + direct + reflect + specular.
    #[default]
    Classic = 0,
    /// Skripkin mode: pure ambient lighting with AO-controlled contrast.
    Skripkin = 1,
    /// Full mode: Classic multi-light lighting + a directional shadow map.
    Full = 2,
}

impl_setting_enum! {
    ShadingMode {
        Classic = 0 => "classic",
        Skripkin = 1 => "skripkin",
        Full = 2 => "full",
    }
    default: Classic
}
