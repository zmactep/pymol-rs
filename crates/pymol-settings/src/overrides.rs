//! Per-object override containers and resolved settings.
//!
//! `ObjectOverrides` holds `Option<T>` wrappers for every object-overridable
//! group. `ResolvedSettings` is the fully-merged result of global + overrides,
//! built once per object per rebuild.

use crate::groups::{
    CartoonOverrides, CartoonSettings, DotOverrides, DotSettings, LineOverrides, LineSettings,
    MeshOverrides, MeshSettings, RibbonOverrides, RibbonSettings, SphereOverrides, SphereSettings,
    StickOverrides, StickSettings, SurfaceOverrides, SurfaceSettings,
};
use crate::macros::Merge;

/// Per-object overrides — only object-overridable groups appear here.
///
/// Global-only groups (shading, ray, ui, movie, behavior) are absent,
/// enforcing at the type level that they cannot be overridden per-object.
#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct ObjectOverrides {
    pub cartoon: CartoonOverrides,
    pub stick: StickOverrides,
    pub sphere: SphereOverrides,
    pub surface: SurfaceOverrides,
    pub ribbon: RibbonOverrides,
    pub line: LineOverrides,
    pub dot: DotOverrides,
    pub mesh: MeshOverrides,
}

/// Fully resolved settings for a single object (global + overrides merged).
///
/// Built once in `MoleculeObject::ensure_representations()` via `Merge::with_overrides`.
#[derive(Debug, Clone)]
pub struct ResolvedSettings {
    pub cartoon: CartoonSettings,
    pub stick: StickSettings,
    pub sphere: SphereSettings,
    pub surface: SurfaceSettings,
    pub ribbon: RibbonSettings,
    pub line: LineSettings,
    pub dot: DotSettings,
    pub mesh: MeshSettings,
}

impl ResolvedSettings {
    /// Build resolved settings from global defaults + per-object overrides.
    pub fn resolve(global: &crate::groups::Settings, overrides: Option<&ObjectOverrides>) -> Self {
        match overrides {
            Some(o) => Self {
                cartoon: global.cartoon.with_overrides(&o.cartoon),
                stick: global.stick.with_overrides(&o.stick),
                sphere: global.sphere.with_overrides(&o.sphere),
                surface: global.surface.with_overrides(&o.surface),
                ribbon: global.ribbon.with_overrides(&o.ribbon),
                line: global.line.with_overrides(&o.line),
                dot: global.dot.with_overrides(&o.dot),
                mesh: global.mesh.with_overrides(&o.mesh),
            },
            None => Self {
                cartoon: global.cartoon.clone(),
                stick: global.stick.clone(),
                sphere: global.sphere.clone(),
                surface: global.surface.clone(),
                ribbon: global.ribbon.clone(),
                line: global.line.clone(),
                dot: global.dot.clone(),
                mesh: global.mesh.clone(),
            },
        }
    }
}
