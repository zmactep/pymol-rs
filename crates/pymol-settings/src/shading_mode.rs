//! Shading mode abstraction.
//!
//! Each shading mode defines its own lighting pipeline and parameters,
//! ensuring that switching between modes has no side effects.

use serde::{Deserialize, Serialize};

use crate::store::GlobalSettings;

/// Available shading modes.
///
/// Each mode reads only its own set of parameters — Classic ignores
/// skripkin_*, Skripkin ignores light_count/ambient/direct/reflect/specular.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[repr(i32)]
pub enum ShadingMode {
    /// Classic PyMOL lighting: multi-light, ambient + direct + reflect + specular.
    Classic = 0,
    /// Skripkin mode: pure ambient + multi-directional shadow AO, no directional lights.
    Skripkin = 1,
}

impl ShadingMode {
    /// Parse from a string alias (case-insensitive).
    pub fn from_str_alias(s: &str) -> Option<Self> {
        match s.to_lowercase().as_str() {
            "classic" | "0" => Some(ShadingMode::Classic),
            "skripkin" | "1" => Some(ShadingMode::Skripkin),
            _ => None,
        }
    }

    /// Get the current shading mode from global settings.
    pub fn from_settings(settings: &GlobalSettings) -> Self {
        Self::from(settings.get_int(crate::definitions::id::shading_mode))
    }

    /// Display name for UI/feedback.
    pub fn name(&self) -> &'static str {
        match self {
            ShadingMode::Classic => "classic",
            ShadingMode::Skripkin => "skripkin",
        }
    }
}

impl From<i32> for ShadingMode {
    fn from(v: i32) -> Self {
        match v {
            1 => ShadingMode::Skripkin,
            _ => ShadingMode::Classic,
        }
    }
}

impl From<ShadingMode> for i32 {
    fn from(m: ShadingMode) -> i32 {
        m as i32
    }
}
