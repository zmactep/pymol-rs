//! Application State
//!
//! Contains all domain/scene data that can exist independently of rendering:
//! object registry, camera, selections, settings, colors, and command executor.

use pymol_cmd::{ArgHint, CommandExecutor};
use pymol_color::{ChainColors, ElementColors, NamedColors};
use pymol_scene::{Camera, Movie, ObjectRegistry, RaytracedImage, SceneManager, SelectionManager, ViewManager};
use pymol_settings::GlobalSettings;

/// Application state containing all scene/domain data
pub struct AppState {
    // =========================================================================
    // Scene
    // =========================================================================
    /// Object registry (molecules, etc.)
    pub registry: ObjectRegistry,
    /// Camera for view control
    pub camera: Camera,
    /// Named selections manager
    pub selections: SelectionManager,
    /// Scene manager for named view snapshots
    pub scenes: SceneManager,
    /// Named views (simpler than scenes - just camera state)
    pub views: ViewManager,
    /// Movie player for frame-based animation
    pub movie: Movie,

    // =========================================================================
    // Settings and Colors
    // =========================================================================
    /// Global settings
    pub settings: GlobalSettings,
    /// Named colors
    pub named_colors: NamedColors,
    /// Element colors
    pub element_colors: ElementColors,
    /// Chain colors
    pub chain_colors: ChainColors,

    // =========================================================================
    // Command System
    // =========================================================================
    /// Command executor
    pub executor: CommandExecutor,

    // =========================================================================
    // Visual Properties
    // =========================================================================
    /// Clear/background color
    pub clear_color: [f32; 3],

    // =========================================================================
    // Raytraced Image Overlay
    // =========================================================================
    /// Stored raytraced image for display (from `ray` command without filename)
    pub raytraced_image: Option<RaytracedImage>,
}

impl Default for AppState {
    fn default() -> Self {
        Self::new()
    }
}

impl AppState {
    /// Create a new application state with default values
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
            executor: CommandExecutor::new(),
            clear_color: [0.0, 0.0, 0.0],
            raytraced_image: None,
        }
    }

    /// Get all command names for autocomplete (queries executor registry)
    pub fn command_names(&self) -> impl Iterator<Item = &str> {
        self.executor.registry().all_names()
    }

    /// Get commands that take file paths as first argument
    pub fn path_commands(&self) -> Vec<&str> {
        self.executor.registry().commands_with_hint(ArgHint::Path, 0)
    }
}
