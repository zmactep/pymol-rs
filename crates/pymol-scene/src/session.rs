//! Session — pure scene state
//!
//! [`Session`] holds all domain/scene data that can exist independently of
//! rendering: object registry, camera, selections, settings, colors, etc.
//!
//! This is the single source of truth that both the GUI (`pymol-gui`) and
//! any headless adapter can own. GPU resources live elsewhere (e.g., in
//! `pymol-render::RenderContext`).

use pymol_color::{ChainColors, ElementColors, NamedColors};
use pymol_settings::GlobalSettings;

use crate::camera::Camera;
use crate::movie::Movie;
use crate::object::ObjectRegistry;
use crate::scene::SceneManager;
use crate::selection::SelectionManager;
use crate::view::ViewManager;
use crate::viewer_trait::RaytracedImage;

/// Pure scene state — no GPU resources, no window, no event loop.
///
/// Owns all molecular objects, camera state, named selections, scenes,
/// views, animation, settings, and color tables.
pub struct Session {
    // =========================================================================
    // Scene
    // =========================================================================
    /// Object registry (molecules, surfaces, maps, CGO, etc.)
    pub registry: ObjectRegistry,
    /// Camera for view control
    pub camera: Camera,
    /// Named selections manager
    pub selections: SelectionManager,
    /// Scene manager for named snapshots (camera + object state)
    pub scenes: SceneManager,
    /// Named views (camera state only — simpler than scenes)
    pub views: ViewManager,
    /// Movie player for frame-based animation
    pub movie: Movie,

    // =========================================================================
    // Settings and Colors
    // =========================================================================
    /// Global rendering settings
    pub settings: GlobalSettings,
    /// Named colors table (e.g., "red", "carbon")
    pub named_colors: NamedColors,
    /// Per-element color defaults
    pub element_colors: ElementColors,
    /// Chain-based coloring (unit struct)
    pub chain_colors: ChainColors,

    // =========================================================================
    // Visual Properties
    // =========================================================================
    /// Background (clear) color as linear RGB floats
    pub clear_color: [f32; 3],

    // =========================================================================
    // Raytraced Image Overlay
    // =========================================================================
    /// Stored raytraced image for display (from `ray` command without filename)
    pub raytraced_image: Option<RaytracedImage>,
}

impl Default for Session {
    fn default() -> Self {
        Self::new()
    }
}

impl Session {
    /// Create a new session with default values.
    pub fn new() -> Self {
        Self {
            registry: ObjectRegistry::new(),
            camera: Camera::new(),
            selections: SelectionManager::new(),
            scenes: SceneManager::new(),
            views: ViewManager::new(),
            movie: Movie::new(),
            settings: GlobalSettings::new(),
            named_colors: NamedColors::default(),
            element_colors: ElementColors::default(),
            chain_colors: ChainColors,
            clear_color: [0.0, 0.0, 0.0],
            raytraced_image: None,
        }
    }
}
