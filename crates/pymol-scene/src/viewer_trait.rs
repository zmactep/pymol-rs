//! Viewer abstraction trait
//!
//! This module defines the `ViewerLike` trait that abstracts the viewer interface,
//! allowing commands to work with different viewer implementations (e.g., `Viewer`,
//! GUI adapters).
//!
//! The trait uses accessor methods (`movie()`, `scenes()`, `selections()`, etc.)
//! to provide default implementations for ~50 methods, so implementors only need
//! to provide the accessors plus a few GPU/window-specific overrides.

use std::path::Path;

use pymol_color::NamedColors;
use serde::{Deserialize, Serialize};
use pymol_settings::GlobalSettings;

use crate::camera::Camera;
use crate::movie::{LoopMode, Movie};
use crate::object::{Object, ObjectRegistry};
use crate::scene::SceneManager;
use crate::selection::SelectionManager;
use crate::session::Session;
use crate::view::ViewManager;

/// Image for display as a viewport overlay
///
/// Used by the `ray` command to show raytraced output, but any source
/// (plugins, computed images, etc.) can provide an image for display.
/// The overlay persists until camera or scene changes occur.
#[derive(Clone, Serialize, Deserialize)]
pub struct ViewportImage {
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
///
/// # Required methods
///
/// Implementors must provide accessor methods for the core state managers
/// and a few GPU/window-specific operations. Most command-facing methods
/// have default implementations that use these accessors.
pub trait ViewerLike {
    // =========================================================================
    // Required Accessors — Core State
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

    /// Request a redraw
    fn request_redraw(&mut self);

    /// Get a reference to the full session
    fn session(&self) -> &Session;

    /// Get a mutable reference to the full session
    fn session_mut(&mut self) -> &mut Session;

    /// Replace the entire session (for loading .pse/.prs files)
    ///
    /// Implementors should mark all objects as dirty so representations rebuild.
    fn replace_session(&mut self, session: Session);

    // =========================================================================
    // Required Accessors — Sub-managers
    // =========================================================================

    /// Get a reference to the movie player
    fn movie(&self) -> &Movie;

    /// Get a mutable reference to the movie player
    fn movie_mut(&mut self) -> &mut Movie;

    /// Get a reference to the scene manager
    fn scenes(&self) -> &SceneManager;

    /// Get a mutable reference to the scene manager
    fn scenes_mut(&mut self) -> &mut SceneManager;

    /// Get a reference to the view manager
    fn views(&self) -> &ViewManager;

    /// Get a mutable reference to the view manager
    fn views_mut(&mut self) -> &mut ViewManager;

    /// Get a reference to the selection manager
    fn selections(&self) -> &SelectionManager;

    /// Get a mutable reference to the selection manager
    fn selections_mut(&mut self) -> &mut SelectionManager;

    /// Get a reference to the named colors table
    fn named_colors(&self) -> &NamedColors;
    /// Get a mutable reference to the named colors table
    fn named_colors_mut(&mut self) -> &mut NamedColors;

    // =========================================================================
    // Required Accessors — Simple State
    // =========================================================================

    /// Get the clear (background) color
    fn clear_color(&self) -> [f32; 3];

    /// Set the clear (background) color
    fn set_clear_color(&mut self, color: [f32; 3]);

    /// Get a reference to the stored viewport image overlay, if any
    fn viewport_image_ref(&self) -> Option<&ViewportImage>;

    /// Set or clear the stored viewport image overlay
    fn set_viewport_image_internal(&mut self, image: Option<ViewportImage>);

    // =========================================================================
    // Camera / View — Default Implementations
    // =========================================================================

    /// Zoom to fit all objects while preserving rotation
    fn zoom_all(&mut self) {
        if let Some((min, max)) = self.objects().extent() {
            self.camera_mut().zoom_to(min, max);
            self.set_viewport_image_internal(None);
            self.request_redraw();
        }
    }

    /// Zoom to fit a specific object while preserving rotation
    ///
    /// For groups, uses the combined extent of all children.
    fn zoom_on(&mut self, name: &str) {
        let extent = if self.objects().get_group(name).is_some() {
            self.objects().group_extent(name)
        } else {
            self.objects().get(name).and_then(|o| o.extent())
        };
        if let Some((min, max)) = extent {
            self.camera_mut().zoom_to(min, max);
            self.set_viewport_image_internal(None);
            self.request_redraw();
        }
    }

    /// Center on all objects without changing zoom or rotation
    fn center_all(&mut self) {
        if let Some((min, max)) = self.objects().extent() {
            self.camera_mut().center_to(min, max);
            self.set_viewport_image_internal(None);
            self.request_redraw();
        }
    }

    /// Center on a specific object without changing zoom or rotation
    ///
    /// For groups, uses the combined extent of all children.
    fn center_on(&mut self, name: &str) {
        let extent = if self.objects().get_group(name).is_some() {
            self.objects().group_extent(name)
        } else {
            self.objects().get(name).and_then(|o| o.extent())
        };
        if let Some((min, max)) = extent {
            self.camera_mut().center_to(min, max);
            self.set_viewport_image_internal(None);
            self.request_redraw();
        }
    }

    /// Reset the camera to default view
    fn reset_view(&mut self) {
        *self.camera_mut() = Camera::new();
        self.set_viewport_image_internal(None);
        if let Some((min, max)) = self.objects().extent() {
            self.camera_mut().reset_view(min, max);
            self.request_redraw();
        }
    }

    // =========================================================================
    // Color — Default Implementations
    // =========================================================================

    /// Get a color index by name
    fn color_index(&self, name: &str) -> Option<u32> {
        self.named_colors().get_by_name(name).map(|(idx, _)| idx)
    }

    /// Set the background color
    fn set_background_color(&mut self, r: f32, g: f32, b: f32) {
        self.set_clear_color([r, g, b]);
        self.request_redraw();
    }

    // =========================================================================
    // Named Selections — Default Implementations
    // =========================================================================

    /// Get a named selection expression by name
    fn get_selection(&self, name: &str) -> Option<&str> {
        self.selections().get_expression(name)
    }

    /// Define (store) a named selection expression
    fn define_selection(&mut self, name: &str, selection: &str) {
        self.selections_mut().define(name, selection);
        self.request_redraw();
    }

    /// Remove a named selection, returns true if it existed
    fn remove_selection(&mut self, name: &str) -> bool {
        let removed = self.selections_mut().remove(name);
        if removed {
            self.request_redraw();
        }
        removed
    }

    /// Rename a named selection, returns true if it existed and was renamed
    fn rename_selection(&mut self, old_name: &str, new_name: &str) -> bool {
        let renamed = self.selections_mut().rename(old_name, new_name);
        if renamed {
            self.request_redraw();
        }
        renamed
    }

    /// Get all named selection names
    fn selection_names(&self) -> Vec<String> {
        self.selections().names()
    }

    /// Set the visibility of a named selection's indicators
    fn set_selection_visible(&mut self, name: &str, visible: bool) {
        self.selections_mut().set_visible(name, visible);
        self.request_redraw();
    }

    /// Check if a named selection's indicators are visible
    fn is_selection_visible(&self, name: &str) -> bool {
        self.selections().is_visible(name)
    }

    /// Set the indicated selection (shows pink indicators in the 3D view)
    fn indicate_selection(&mut self, selection: &str) {
        self.selections_mut().indicate(selection);
        self.request_redraw();
    }

    /// Clear the indicated selection (hides all selection indicators)
    fn clear_indication(&mut self) {
        self.selections_mut().clear_indication();
        self.request_redraw();
    }

    /// Get the currently indicated selection expression
    fn indicated_selection(&self) -> Option<&str> {
        self.selections().indicated_selection()
    }

    // =========================================================================
    // Screenshot Capture — Optional Overrides
    // =========================================================================

    /// Capture a PNG screenshot (optional - not all viewers support this)
    fn capture_png(
        &mut self,
        _path: &Path,
        _width: Option<u32>,
        _height: Option<u32>,
    ) -> Result<(), String> {
        Err("Screenshot capture not supported by this viewer".to_string())
    }

    // =========================================================================
    // Raytracing — Optional Overrides
    // =========================================================================

    /// Perform raytracing and return the image data as RGBA bytes
    fn raytrace(
        &mut self,
        _width: Option<u32>,
        _height: Option<u32>,
        _antialias: u32,
    ) -> Result<Vec<u8>, String> {
        Err("Raytracing not supported by this viewer".to_string())
    }

    /// Perform raytracing and save to a PNG file
    fn raytrace_to_file(
        &mut self,
        _path: &Path,
        _width: Option<u32>,
        _height: Option<u32>,
        _antialias: u32,
    ) -> Result<(u32, u32), String> {
        Err("Raytracing not supported by this viewer".to_string())
    }

    // =========================================================================
    // Viewport Image Overlay — Default Implementations
    // =========================================================================

    /// Store an image for display as a viewport overlay
    fn set_viewport_image(&mut self, image: Option<ViewportImage>) {
        self.set_viewport_image_internal(image);
        self.request_redraw();
    }

    /// Get the stored viewport image overlay, if any
    fn get_viewport_image(&self) -> Option<&ViewportImage> {
        self.viewport_image_ref()
    }

    /// Clear the viewport image overlay
    fn clear_viewport_image(&mut self) {
        self.set_viewport_image_internal(None);
    }

    // =========================================================================
    // Async Operations — Optional Override
    // =========================================================================

    /// Request an async fetch operation from RCSB PDB (non-blocking)
    ///
    /// Returns `true` if the request was accepted, `false` if async fetch
    /// is not supported (caller should use sync fallback).
    fn request_async_fetch(&mut self, _code: &str, _name: &str, _format: u8) -> bool {
        false
    }

    // =========================================================================
    // Scene Management
    //
    // scene_store and scene_recall MUST be per-impl because they need
    // simultaneous borrows on different parts of self (e.g. &mut scenes + &camera).
    // The simpler methods have default implementations.
    // =========================================================================

    /// Store the current view as a named scene
    fn scene_store(&mut self, _key: &str, _storemask: u32);

    /// Recall a named scene
    fn scene_recall(&mut self, _key: &str, _animate: bool, _duration: f32) -> Result<(), String>;

    /// Delete a named scene
    fn scene_delete(&mut self, key: &str) -> bool {
        self.scenes_mut().delete(key).is_some()
    }

    /// Get list of all scene names
    fn scene_list(&self) -> Vec<String> {
        self.scenes().list().to_vec()
    }

    /// Rename a scene
    fn scene_rename(&mut self, old_key: &str, new_key: &str) -> Result<(), String> {
        self.scenes_mut().rename(old_key, new_key).map_err(|e| e.to_string())
    }

    /// Clear all scenes
    fn scene_clear(&mut self) {
        self.scenes_mut().clear();
    }

    // =========================================================================
    // Movie Control — Default Implementations
    // =========================================================================

    /// Start movie playback
    fn movie_play(&mut self) {
        let loop_setting = self.settings().get_bool(pymol_settings::id::movie_loop);
        if loop_setting {
            self.movie_mut().set_loop_mode(LoopMode::Loop);
        }
        self.movie_mut().play();
        self.request_redraw();
    }

    /// Stop movie and reset to first frame
    fn movie_stop(&mut self) {
        self.movie_mut().stop();
        self.sync_movie_to_objects();
        self.request_redraw();
    }

    /// Pause movie playback
    fn movie_pause(&mut self) {
        self.movie_mut().pause();
    }

    /// Toggle play/pause
    fn movie_toggle(&mut self) {
        self.movie_mut().toggle();
        self.request_redraw();
    }

    /// Go to a specific frame (0-indexed)
    fn movie_goto(&mut self, frame: usize) {
        self.movie_mut().goto_frame(frame);
        self.sync_movie_to_objects();
        self.request_redraw();
    }

    /// Advance to next frame
    fn movie_next(&mut self) {
        self.movie_mut().next_frame();
        self.sync_movie_to_objects();
        self.request_redraw();
    }

    /// Go to previous frame
    fn movie_prev(&mut self) {
        self.movie_mut().prev_frame();
        self.sync_movie_to_objects();
        self.request_redraw();
    }

    /// Update `Movie::n_object_states` from the current max across all objects.
    /// Call after loading/deleting objects.
    fn update_movie_state_count(&mut self) {
        let max_states = self.objects().iter()
            .map(|obj| obj.n_states())
            .max()
            .unwrap_or(1);
        self.movie_mut().set_n_object_states(max_states);
    }

    /// Sync the movie's current frame to all objects' display states and camera.
    fn sync_movie_to_objects(&mut self) {
        let state_index = self.movie().frame_to_state(self.movie().current_frame());
        let view = self.movie().interpolated_view();

        let names: Vec<String> = self.objects().names().map(|s| s.to_string()).collect();
        for name in &names {
            if let Some(obj) = self.objects_mut().get_molecule_mut(name) {
                obj.set_display_state(state_index);
            }
        }

        self.session_mut().apply_movie_object_transforms();

        if let Some(view) = view {
            self.camera_mut().set_view(view);
        }
    }

    /// Set frames per second
    fn movie_set_fps(&mut self, fps: f32) {
        self.movie_mut().set_fps(fps);
    }

    /// Get total number of frames (effective: movie frames or max object states)
    fn movie_frame_count(&self) -> usize {
        self.movie().effective_frame_count()
    }

    /// Get current frame index (0-indexed)
    fn movie_current_frame(&self) -> usize {
        self.movie().current_frame()
    }

    /// Check if movie is currently playing
    fn movie_is_playing(&self) -> bool {
        self.movie().is_playing()
    }

    /// Set loop mode: 0 = once, 1 = loop, 2 = swing
    fn movie_set_loop_mode(&mut self, mode: u8) {
        let loop_mode = match mode {
            1 => LoopMode::Loop,
            2 => LoopMode::Swing,
            _ => LoopMode::Once,
        };
        self.movie_mut().set_loop_mode(loop_mode);
    }

    /// Set number of movie frames
    fn movie_set_frame_count(&mut self, count: usize) {
        self.movie_mut().set_frame_count(count);
    }

    /// Set movie frames from a specification (state index per frame)
    fn movie_set_from_spec(&mut self, states: Vec<usize>) {
        self.movie_mut().set_from_spec(states);
    }

    /// Append movie frames from a specification
    fn movie_append_from_spec(&mut self, states: Vec<usize>) {
        self.movie_mut().append_from_spec(states);
    }

    /// Store a camera keyframe at the given frame (0-indexed)
    fn movie_store_view(&mut self, frame: usize) {
        let view = self.camera().current_view();
        self.movie_mut().store_camera_keyframe(frame, view);
        let loop_movie = self.settings().get_bool(pymol_settings::id::movie_loop);
        if self.settings().get_bool(pymol_settings::id::movie_auto_interpolate) {
            self.movie_mut().interpolate_keyframes(loop_movie);
        }
        self.request_redraw();
    }

    /// Store a scene keyframe at the given frame (0-indexed)
    fn movie_store_scene(&mut self, frame: usize, scene_name: &str) {
        let view = self.scenes().get(scene_name)
            .map(|s| s.view.clone())
            .unwrap_or_else(|| self.camera().current_view());
        self.movie_mut().store_scene_keyframe(frame, scene_name, view);
        let loop_movie = self.settings().get_bool(pymol_settings::id::movie_loop);
        if self.settings().get_bool(pymol_settings::id::movie_auto_interpolate) {
            self.movie_mut().interpolate_keyframes(loop_movie);
        }
        self.request_redraw();
    }

    /// Store an object keyframe at the given frame
    fn movie_store_object(&mut self, frame: usize, object: &str, state: Option<usize>) {
        // Capture the object's current TTT transform matrix
        let transform = self
            .objects()
            .get_molecule(object)
            .map(|mol| mol.state().transform.clone());
        self.movie_mut()
            .store_object_keyframe(frame, object, state, transform);
        // Auto-interpolate object transforms between keyframes
        let loop_movie = self.settings().get_bool(pymol_settings::id::movie_loop);
        if self.settings().get_bool(pymol_settings::id::movie_auto_interpolate) {
            self.movie_mut().interpolate_object_keyframes(loop_movie);
        }
        self.request_redraw();
    }

    /// Clear camera keyframe(s). None = clear all.
    fn movie_clear_view(&mut self, frame: Option<usize>) {
        self.movie_mut().clear_camera_keyframes(frame);
    }

    /// Trigger keyframe interpolation
    fn movie_interpolate(&mut self) {
        let loop_movie = self.settings().get_bool(pymol_settings::id::movie_loop);
        self.movie_mut().interpolate_keyframes(loop_movie);
        self.request_redraw();
    }

    /// Render a specific frame to PNG (for mpng)
    fn capture_frame_png(
        &mut self,
        _frame: usize,
        _path: &Path,
        _width: Option<u32>,
        _height: Option<u32>,
    ) -> Result<(), String> {
        Err("Not supported".into())
    }

    // =========================================================================
    // Rock Animation — Default Implementations
    // =========================================================================

    /// Toggle rock mode (Y-axis oscillation)
    fn rock_toggle(&mut self) {
        let current = self.movie().is_rock_enabled();
        self.movie_mut().set_rock(!current);
        self.request_redraw();
    }

    /// Set rock mode explicitly
    fn rock_set(&mut self, enabled: bool) {
        self.movie_mut().set_rock(enabled);
        self.request_redraw();
    }

    /// Check if rock is enabled
    fn rock_is_enabled(&self) -> bool {
        self.movie().is_rock_enabled()
    }

    // =========================================================================
    // Named View Storage — Default Implementations
    //
    // view_store uses camera().current_view() then views_mut() — no conflict.
    // view_recall MUST be per-impl (needs &views + &mut camera simultaneously).
    // =========================================================================

    /// Store the current camera view under a name
    fn view_store(&mut self, key: &str) {
        let view = self.camera().current_view();
        self.views_mut().store_view(key, view);
    }

    /// Recall a named view (MUST be per-impl due to borrow conflict)
    fn view_recall(&mut self, _key: &str, _animate: f32) -> Result<(), String>;

    /// Delete a named view
    fn view_delete(&mut self, key: &str) -> bool {
        self.views_mut().delete(key).is_some()
    }

    /// Get list of all stored view names
    fn view_list(&self) -> Vec<String> {
        self.views().list().to_vec()
    }

    /// Clear all stored views
    fn view_clear(&mut self) {
        self.views_mut().clear();
    }

    // =========================================================================
    // Viewport / Window Size — Optional Overrides
    // =========================================================================

    /// Get the current viewport size (width, height)
    fn viewport_size(&self) -> (u32, u32) {
        (0, 0)
    }

    /// Request a viewport resize
    fn viewport_set_size(&mut self, _width: u32, _height: u32) {
        // Default: no-op
    }

    // =========================================================================
    // Fullscreen Mode — Optional Overrides
    // =========================================================================

    /// Check if fullscreen mode is active
    fn is_fullscreen(&self) -> bool {
        false
    }

    /// Set fullscreen mode
    fn set_fullscreen(&mut self, _enabled: bool) {
        // Default: no-op
    }

    /// Toggle fullscreen mode
    fn toggle_fullscreen(&mut self) {
        // Default: no-op
    }
}
