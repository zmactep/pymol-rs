//! Typed settings groups.
//!
//! Each group is defined via [`define_settings_group!`] and contains a plain
//! data struct with concrete Rust types. Groups marked `group` are
//! object-overridable; `group_global` groups are global-only.

pub mod behavior;
pub mod cartoon;
pub mod dot;
pub mod line;
pub mod mesh;
pub mod movie;
pub mod ribbon;
pub mod shading;
pub mod sphere;
pub mod stick;
pub mod surface;
pub mod ui;

pub use behavior::BehaviorSettings;
pub use cartoon::{CartoonOverrides, CartoonSettings};
pub use dot::{DotOverrides, DotSettings};
pub use line::{LineOverrides, LineSettings};
pub use mesh::{MeshOverrides, MeshSettings};
pub use movie::MovieSettings;
pub use ribbon::{RibbonOverrides, RibbonSettings};
pub use shading::{
    ClassicShadingSettings, CommonShadingSettings, FullShadingSettings, ShadingSettings,
    SkripkinShadingSettings,
};
pub use sphere::{SphereOverrides, SphereSettings};
pub use stick::{StickOverrides, StickSettings};
pub use surface::{SurfaceOverrides, SurfaceSettings};
pub use ui::UiSettings;

/// Top-level settings container — holds all typed groups.
///
/// Global-only groups live directly on `Settings`.
/// Object-overridable groups also live here (as the global defaults).
#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct Settings {
    // Global-only
    pub shading: ShadingSettings,
    pub ui: UiSettings,
    pub movie: MovieSettings,
    pub behavior: BehaviorSettings,

    // Object-overridable (global defaults)
    pub cartoon: CartoonSettings,
    pub stick: StickSettings,
    pub sphere: SphereSettings,
    pub surface: SurfaceSettings,
    pub ribbon: RibbonSettings,
    pub line: LineSettings,
    pub dot: DotSettings,
    pub mesh: MeshSettings,
}

