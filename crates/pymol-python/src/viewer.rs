//! Viewer wrapper for Python
//!
//! Provides access to the Viewer for headless and GUI modes.

// Allow dead code - these methods are part of the public API for future GUI mode
#![allow(dead_code)]

use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;

use parking_lot::RwLock;
use pymol_color::NamedColors;
use pymol_mol::ObjectMolecule;
use pymol_scene::{ObjectRegistry, Viewer};
use pymol_settings::GlobalSettings;

/// Headless viewer state for scripting without GUI
///
/// This provides a lightweight interface for molecular operations
/// without requiring window creation or GPU initialization.
pub struct HeadlessViewer {
    /// Object registry for named objects
    pub objects: ObjectRegistry,
    /// Global settings
    pub settings: GlobalSettings,
    /// Named colors
    pub named_colors: NamedColors,
    /// Named selections (name -> selection expression)
    pub selections: std::collections::HashMap<String, String>,
    /// Current state index
    pub current_state: usize,
}

impl Default for HeadlessViewer {
    fn default() -> Self {
        Self::new()
    }
}

impl HeadlessViewer {
    /// Create a new headless viewer
    pub fn new() -> Self {
        HeadlessViewer {
            objects: ObjectRegistry::new(),
            settings: GlobalSettings::new(),
            named_colors: NamedColors::default(),
            selections: std::collections::HashMap::new(),
            current_state: 0,
        }
    }

    /// Add a molecule to the viewer
    pub fn add_molecule(&mut self, mol: ObjectMolecule) -> &str {
        use pymol_scene::MoleculeObject;
        let obj = MoleculeObject::new(mol);
        self.objects.add(obj)
    }

    /// Get a molecule by name
    pub fn get_molecule(&self, name: &str) -> Option<&ObjectMolecule> {
        self.objects.get_molecule(name).map(|m| m.molecule())
    }

    /// Get a mutable molecule by name
    pub fn get_molecule_mut(&mut self, name: &str) -> Option<&mut ObjectMolecule> {
        self.objects.get_molecule_mut(name).map(|m| m.molecule_mut())
    }

    /// List all object names
    pub fn get_names(&self) -> Vec<String> {
        self.objects.names().map(|s| s.to_string()).collect()
    }

    /// Remove an object
    pub fn delete(&mut self, name: &str) -> bool {
        self.objects.remove(name).is_some()
    }

    /// Delete all objects
    pub fn delete_all(&mut self) {
        self.objects.clear();
        self.selections.clear();
    }

    /// Register a named selection
    pub fn define_selection(&mut self, name: &str, selection: &str) {
        self.selections.insert(name.to_string(), selection.to_string());
    }

    /// Get a named selection expression
    pub fn get_selection(&self, name: &str) -> Option<&str> {
        self.selections.get(name).map(|s| s.as_str())
    }

    /// Compute bounding box of all objects
    pub fn extent(&self) -> Option<(lin_alg::f32::Vec3, lin_alg::f32::Vec3)> {
        self.objects.extent()
    }

    /// Count atoms in all molecules or matching a selection
    pub fn count_atoms(&self, selection: Option<&str>) -> usize {
        match selection {
            None => {
                // Count all atoms
                self.objects
                    .names()
                    .filter_map(|name| self.objects.get_molecule(name))
                    .map(|m| m.molecule().atom_count())
                    .sum()
            }
            Some(sel) => {
                // Count atoms matching selection
                let mut count = 0;
                for name in self.objects.names() {
                    if let Some(mol_obj) = self.objects.get_molecule(name) {
                        if let Ok(result) = pymol_select::select(mol_obj.molecule(), sel) {
                            count += result.count();
                        }
                    }
                }
                count
            }
        }
    }
}

/// Viewer mode enum
pub enum ViewerMode {
    /// Headless mode - no window, no GPU
    Headless(HeadlessViewer),
    /// GUI mode - with window and GPU rendering
    Gui(Viewer),
}

impl ViewerMode {
    /// Check if this is headless mode
    pub fn is_headless(&self) -> bool {
        matches!(self, ViewerMode::Headless(_))
    }

    /// Get object registry (works for both modes)
    pub fn objects(&self) -> &ObjectRegistry {
        match self {
            ViewerMode::Headless(h) => &h.objects,
            ViewerMode::Gui(v) => v.objects(),
        }
    }

    /// Get mutable object registry
    pub fn objects_mut(&mut self) -> &mut ObjectRegistry {
        match self {
            ViewerMode::Headless(h) => &mut h.objects,
            ViewerMode::Gui(v) => v.objects_mut(),
        }
    }

    /// Get settings
    pub fn settings(&self) -> &GlobalSettings {
        match self {
            ViewerMode::Headless(h) => &h.settings,
            ViewerMode::Gui(v) => v.settings(),
        }
    }

    /// Get mutable settings
    pub fn settings_mut(&mut self) -> &mut GlobalSettings {
        match self {
            ViewerMode::Headless(h) => &mut h.settings,
            ViewerMode::Gui(v) => v.settings_mut(),
        }
    }

    /// Get named colors
    pub fn named_colors(&self) -> &NamedColors {
        match self {
            ViewerMode::Headless(h) => &h.named_colors,
            ViewerMode::Gui(v) => v.named_colors(),
        }
    }

    /// Add a molecule
    pub fn add_molecule(&mut self, mol: ObjectMolecule) -> &str {
        match self {
            ViewerMode::Headless(h) => h.add_molecule(mol),
            ViewerMode::Gui(v) => v.add_molecule(mol),
        }
    }

    /// Zoom to fit all objects (GUI mode only, no-op in headless)
    pub fn zoom_all(&mut self) {
        if let ViewerMode::Gui(v) = self {
            v.zoom_all();
        }
    }

    /// Center on all objects (GUI mode only, no-op in headless)
    pub fn center_all(&mut self) {
        if let ViewerMode::Gui(v) = self {
            v.center_all();
        }
    }

    /// Reset view (GUI mode only, no-op in headless)
    pub fn reset_view(&mut self) {
        if let ViewerMode::Gui(v) = self {
            v.reset_view();
        }
    }

    /// Request redraw (GUI mode only, no-op in headless)
    pub fn request_redraw(&mut self) {
        if let ViewerMode::Gui(v) = self {
            v.request_redraw();
        }
    }
}

// =============================================================================
// Shared Viewer State (for non-blocking GUI mode)
// =============================================================================

/// Shared state that can be accessed from both Python and the viewer thread
pub struct SharedViewerState {
    /// Object registry (molecules, surfaces, etc.)
    pub objects: ObjectRegistry,
    /// Global settings
    pub settings: GlobalSettings,
    /// Named colors
    pub named_colors: NamedColors,
    /// Named selections (name -> selection expression)
    pub selections: std::collections::HashMap<String, String>,
    /// Flag indicating state has changed and needs redraw
    pub dirty: bool,
}

impl Default for SharedViewerState {
    fn default() -> Self {
        Self::new()
    }
}

impl SharedViewerState {
    /// Create a new shared viewer state
    pub fn new() -> Self {
        SharedViewerState {
            objects: ObjectRegistry::new(),
            settings: GlobalSettings::new(),
            named_colors: NamedColors::default(),
            selections: std::collections::HashMap::new(),
            dirty: false,
        }
    }

    /// Add a molecule to the registry
    pub fn add_molecule(&mut self, mol: ObjectMolecule) -> String {
        use pymol_scene::MoleculeObject;
        let obj = MoleculeObject::new(mol);
        let name = self.objects.add(obj).to_string();
        self.dirty = true;
        name
    }

    /// Mark as needing redraw
    pub fn mark_dirty(&mut self) {
        self.dirty = true;
    }

    /// Check and clear dirty flag
    pub fn take_dirty(&mut self) -> bool {
        let was_dirty = self.dirty;
        self.dirty = false;
        was_dirty
    }
}

/// Handle to a running viewer
///
/// This allows checking if the viewer is still running and provides
/// a way to request it to close.
pub struct ViewerHandle {
    /// Whether the viewer is still running
    running: Arc<AtomicBool>,
    /// Shared state with the viewer
    state: Arc<RwLock<SharedViewerState>>,
    /// Thread handle (if running in background)
    #[allow(dead_code)]
    thread: Option<std::thread::JoinHandle<()>>,
}

impl ViewerHandle {
    /// Create a new viewer handle
    pub fn new(
        running: Arc<AtomicBool>,
        state: Arc<RwLock<SharedViewerState>>,
        thread: Option<std::thread::JoinHandle<()>>,
    ) -> Self {
        ViewerHandle {
            running,
            state,
            thread,
        }
    }

    /// Check if the viewer is still running
    pub fn is_running(&self) -> bool {
        self.running.load(Ordering::SeqCst)
    }

    /// Request the viewer to close
    pub fn request_close(&self) {
        self.running.store(false, Ordering::SeqCst);
    }

    /// Get access to the shared state
    pub fn state(&self) -> &Arc<RwLock<SharedViewerState>> {
        &self.state
    }

    /// Wait for the viewer to close (blocks)
    pub fn wait(self) {
        if let Some(thread) = self.thread {
            let _ = thread.join();
        }
    }
}
