//! Viewer abstraction trait
//!
//! This module defines the `ViewerLike` trait that abstracts the viewer interface,
//! allowing commands to work with different viewer implementations (e.g., `Viewer`,
//! GUI adapters).

use std::path::Path;

use crate::camera::Camera;
use crate::object::ObjectRegistry;
use pymol_settings::GlobalSettings;

/// Stored raytraced image for display in the viewport
///
/// When `ray` is called without a filename, the raytraced image is stored
/// here for display. It persists until camera or scene changes occur.
#[derive(Clone)]
pub struct RaytracedImage {
    /// RGBA image data, row-major, top-to-bottom
    pub data: Vec<u8>,
    /// Image width in pixels
    pub width: u32,
    /// Image height in pixels
    pub height: u32,
}

/// Trait for types that can serve as a viewer backend for command execution
///
/// This trait abstracts the viewer interface, allowing commands to work with
/// different viewer implementations (e.g., `Viewer`, GUI adapters).
pub trait ViewerLike {
    // =========================================================================
    // Core - Objects, Camera, and Settings
    // =========================================================================

    /// Get a reference to the object registry
    fn objects(&self) -> &ObjectRegistry;

    /// Get a mutable reference to the object registry
    fn objects_mut(&mut self) -> &mut ObjectRegistry;

    /// Get a reference to the camera
    fn camera(&self) -> &Camera;

    /// Get a mutable reference to the camera
    fn camera_mut(&mut self) -> &mut Camera;

    /// Get a reference to the global settings
    fn settings(&self) -> &GlobalSettings;

    /// Get a mutable reference to the global settings
    fn settings_mut(&mut self) -> &mut GlobalSettings;

    /// Zoom to fit all objects while preserving rotation
    fn zoom_all(&mut self);

    /// Zoom to fit a specific object while preserving rotation
    fn zoom_on(&mut self, name: &str);

    /// Center on all objects without changing zoom or rotation
    fn center_all(&mut self);

    /// Center on a specific object without changing zoom or rotation
    fn center_on(&mut self, name: &str);

    /// Reset the camera to default view
    fn reset_view(&mut self);

    /// Request a redraw
    fn request_redraw(&mut self);

    /// Get a color index by name
    fn color_index(&self, name: &str) -> Option<u32>;

    /// Set the background color
    fn set_background_color(&mut self, r: f32, g: f32, b: f32);

    // =========================================================================
    // Named Selections
    // =========================================================================

    /// Get a named selection expression by name
    fn get_selection(&self, name: &str) -> Option<&str>;

    /// Define (store) a named selection expression
    fn define_selection(&mut self, name: &str, selection: &str);

    /// Remove a named selection, returns true if it existed
    fn remove_selection(&mut self, name: &str) -> bool;

    /// Get all named selection names
    fn selection_names(&self) -> Vec<String>;

    /// Set the visibility of a named selection's indicators
    fn set_selection_visible(&mut self, name: &str, visible: bool);

    /// Check if a named selection's indicators are visible
    fn is_selection_visible(&self, name: &str) -> bool;

    /// Set the indicated selection (shows pink indicators in the 3D view)
    /// For backwards compatibility - creates/updates an "indicate" selection
    fn indicate_selection(&mut self, selection: &str);

    /// Clear the indicated selection (hides all selection indicators)
    fn clear_indication(&mut self);

    /// Get the currently indicated selection expression
    /// For backwards compatibility - returns first visible selection
    fn indicated_selection(&self) -> Option<&str>;

    // =========================================================================
    // Screenshot Capture
    // =========================================================================

    /// Capture a PNG screenshot (optional - not all viewers support this)
    ///
    /// Returns an error by default. Override for viewers that support screenshots.
    fn capture_png(
        &mut self,
        _path: &Path,
        _width: Option<u32>,
        _height: Option<u32>,
    ) -> Result<(), String> {
        Err("Screenshot capture not supported by this viewer".to_string())
    }

    // =========================================================================
    // Raytracing
    // =========================================================================

    /// Perform raytracing and return the image data as RGBA bytes
    ///
    /// # Arguments
    /// * `width` - Image width (None = use viewport width)
    /// * `height` - Image height (None = use viewport height)
    /// * `antialias` - Antialiasing level (1 = no AA, 2-4 = supersampling)
    ///
    /// # Returns
    /// * `Ok(Vec<u8>)` - RGBA image data, row-major
    /// * `Err(String)` - If raytracing fails or is not supported
    fn raytrace(
        &mut self,
        _width: Option<u32>,
        _height: Option<u32>,
        _antialias: u32,
    ) -> Result<Vec<u8>, String> {
        Err("Raytracing not supported by this viewer".to_string())
    }

    /// Perform raytracing and save to a PNG file
    ///
    /// # Arguments
    /// * `path` - Output file path
    /// * `width` - Image width (None = use viewport width)
    /// * `height` - Image height (None = use viewport height)
    /// * `antialias` - Antialiasing level (1 = no AA, 2-4 = supersampling)
    ///
    /// # Returns
    /// * `Ok((u32, u32))` - Actual image dimensions (width, height)
    /// * `Err(String)` - If raytracing or file save fails
    fn raytrace_to_file(
        &mut self,
        _path: &Path,
        _width: Option<u32>,
        _height: Option<u32>,
        _antialias: u32,
    ) -> Result<(u32, u32), String> {
        Err("Raytracing not supported by this viewer".to_string())
    }

    /// Store a raytraced image for display in the viewport
    ///
    /// Called by the `ray` command when no filename is provided.
    /// The image will be displayed as an overlay until cleared.
    fn set_raytraced_image(&mut self, _image: Option<RaytracedImage>) {
        // Default: do nothing (viewer doesn't support raytraced overlay)
    }

    /// Get the stored raytraced image, if any
    ///
    /// Returns the raytraced image that should be displayed as an overlay.
    fn get_raytraced_image(&self) -> Option<&RaytracedImage> {
        None
    }

    /// Clear the stored raytraced image
    ///
    /// Should be called when camera or scene changes invalidate the image.
    fn clear_raytraced_image(&mut self) {
        self.set_raytraced_image(None);
    }

    // =========================================================================
    // Async Operations
    // =========================================================================

    /// Request an async fetch operation from RCSB PDB (non-blocking)
    ///
    /// This method allows commands to request background fetch operations without
    /// blocking the main thread. The viewer implementation decides how to handle
    /// the async operation.
    ///
    /// # Arguments
    ///
    /// * `code` - PDB ID to fetch (e.g., "1ubq")
    /// * `name` - Object name to use when adding to registry
    /// * `format` - Format code: 0 = CIF (default), 1 = PDB
    ///
    /// # Returns
    ///
    /// * `true` - Request was accepted; fetch is running in background
    /// * `false` - Async fetch not supported; caller should use sync fallback
    ///
    /// The default implementation returns `false`, indicating async fetch is not
    /// supported. GUI viewers should override this to spawn background tasks.
    fn request_async_fetch(&mut self, _code: &str, _name: &str, _format: u8) -> bool {
        false
    }

    // =========================================================================
    // Scene Management
    // =========================================================================

    /// Store the current view as a named scene
    ///
    /// # Arguments
    /// * `key` - Scene name/key
    /// * `storemask` - Bitmask of what to store (use SceneStoreMask::ALL.bits() for everything)
    fn scene_store(&mut self, _key: &str, _storemask: u32) {
        // Default: no-op (scene management not supported)
    }

    /// Recall a named scene
    ///
    /// # Arguments
    /// * `key` - Scene name to recall
    /// * `animate` - Whether to animate the transition
    /// * `duration` - Animation duration in seconds
    fn scene_recall(&mut self, _key: &str, _animate: bool, _duration: f32) -> Result<(), String> {
        Err("Scene management not supported by this viewer".to_string())
    }

    /// Delete a named scene
    ///
    /// Returns true if the scene existed and was deleted
    fn scene_delete(&mut self, _key: &str) -> bool {
        false
    }

    /// Get list of all scene names
    fn scene_list(&self) -> Vec<String> {
        Vec::new()
    }

    /// Rename a scene
    fn scene_rename(&mut self, _old_key: &str, _new_key: &str) -> Result<(), String> {
        Err("Scene management not supported by this viewer".to_string())
    }

    /// Clear all scenes
    fn scene_clear(&mut self) {
        // Default: no-op
    }

    // =========================================================================
    // Movie Control
    // =========================================================================

    /// Start movie playback
    fn movie_play(&mut self) {
        // Default: no-op
    }

    /// Stop movie and reset to first frame
    fn movie_stop(&mut self) {
        // Default: no-op
    }

    /// Pause movie playback
    fn movie_pause(&mut self) {
        // Default: no-op
    }

    /// Toggle play/pause
    fn movie_toggle(&mut self) {
        // Default: no-op
    }

    /// Go to a specific frame (0-indexed)
    fn movie_goto(&mut self, _frame: usize) {
        // Default: no-op
    }

    /// Advance to next frame
    fn movie_next(&mut self) {
        // Default: no-op
    }

    /// Go to previous frame
    fn movie_prev(&mut self) {
        // Default: no-op
    }

    /// Set frames per second
    fn movie_set_fps(&mut self, _fps: f32) {
        // Default: no-op
    }

    /// Get total number of frames
    fn movie_frame_count(&self) -> usize {
        0
    }

    /// Get current frame index (0-indexed)
    fn movie_current_frame(&self) -> usize {
        0
    }

    /// Check if movie is currently playing
    fn movie_is_playing(&self) -> bool {
        false
    }

    /// Set loop mode: 0 = once, 1 = loop, 2 = swing
    fn movie_set_loop_mode(&mut self, _mode: u8) {
        // Default: no-op
    }

    /// Set number of movie frames
    fn movie_set_frame_count(&mut self, _count: usize) {
        // Default: no-op
    }

    // =========================================================================
    // Rock Animation
    // =========================================================================

    /// Toggle rock mode (Y-axis oscillation)
    fn rock_toggle(&mut self) {
        // Default: no-op
    }

    /// Set rock mode explicitly
    fn rock_set(&mut self, _enabled: bool) {
        // Default: no-op
    }

    /// Check if rock is enabled
    fn rock_is_enabled(&self) -> bool {
        false
    }

    // =========================================================================
    // Named View Storage (simpler than scenes - just camera state)
    // =========================================================================

    /// Store the current camera view under a name
    ///
    /// Unlike scenes, views only store the camera state (18 values),
    /// not colors, representations, or frame state.
    fn view_store(&mut self, _key: &str) {
        // Default: no-op
    }

    /// Recall a named view
    ///
    /// # Arguments
    /// * `key` - View name to recall
    /// * `animate` - Animation duration in seconds (0 = instant)
    fn view_recall(&mut self, _key: &str, _animate: f32) -> Result<(), String> {
        Err("View storage not supported by this viewer".to_string())
    }

    /// Delete a named view
    ///
    /// Returns true if the view existed and was deleted
    fn view_delete(&mut self, _key: &str) -> bool {
        false
    }

    /// Get list of all stored view names
    fn view_list(&self) -> Vec<String> {
        Vec::new()
    }

    /// Clear all stored views
    fn view_clear(&mut self) {
        // Default: no-op
    }

    // =========================================================================
    // Viewport / Window Size
    // =========================================================================

    /// Get the current viewport size (width, height)
    fn viewport_size(&self) -> (u32, u32) {
        (0, 0) // Default: unknown
    }

    /// Request a viewport resize
    ///
    /// Note: This may not work in all contexts (e.g., embedded viewers)
    fn viewport_set_size(&mut self, _width: u32, _height: u32) {
        // Default: no-op
    }

    // =========================================================================
    // Fullscreen Mode
    // =========================================================================

    /// Check if fullscreen mode is active
    fn is_fullscreen(&self) -> bool {
        false
    }

    /// Set fullscreen mode
    ///
    /// # Arguments
    /// * `enabled` - true to enable fullscreen, false to disable
    fn set_fullscreen(&mut self, _enabled: bool) {
        // Default: no-op
    }

    /// Toggle fullscreen mode
    fn toggle_fullscreen(&mut self) {
        // Default: no-op
    }
}
