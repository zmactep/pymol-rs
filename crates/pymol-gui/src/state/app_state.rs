//! Application State
//!
//! Contains all domain/scene data that can exist independently of rendering:
//! object registry, camera, selections, settings, colors, and command executor.

use std::collections::HashMap;

use pymol_cmd::{ArgHint, CommandExecutor};
use pymol_color::{ChainColors, ElementColors, NamedColors};
use pymol_scene::{Camera, ObjectRegistry, SelectionEntry};
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
    /// Named selections (name -> entry with expression and visibility)
    pub selections: HashMap<String, SelectionEntry>,

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
            selections: HashMap::new(),
            settings: GlobalSettings::new(),
            named_colors: NamedColors::default(),
            element_colors: ElementColors::default(),
            chain_colors: ChainColors,
            executor: CommandExecutor::new(),
            clear_color: [0.0, 0.0, 0.0],
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
