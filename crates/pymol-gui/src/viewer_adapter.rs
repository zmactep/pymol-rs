//! ViewerAdapter - bridges App state to the ViewerLike trait for commands
//!
//! This module provides the adapter that allows the command system to interact
//! with the GUI's state through the `ViewerLike` trait.

use std::path::Path;

use pymol_cmd::ViewerLike;
use pymol_render::RenderContext;
use pymol_scene::{capture_png_to_file, Camera, LoopMode, ObjectRegistry, SceneStoreMask, SelectionEntry};

use crate::async_tasks::TaskRunner;
use crate::fetch::FetchTask;
use crate::state::AppState;

/// Adapter that wraps AppState to implement ViewerLike for command execution
///
/// This struct borrows a reference to the application state and
/// implements the `ViewerLike` trait, allowing commands to interact with
/// the application state in a uniform way.
pub struct ViewerAdapter<'a> {
    /// Reference to the application state
    pub state: &'a mut AppState,
    /// Reference to task runner for background operations
    pub task_runner: &'a TaskRunner,
    /// Reference to render context for GPU operations (None if not yet initialized)
    pub render_context: Option<&'a RenderContext>,
    /// Default size for PNG capture when width/height not specified
    pub default_size: (u32, u32),
    /// Flag to indicate if a redraw is needed (external to state for borrow checker)
    pub needs_redraw: &'a mut bool,
}

impl<'a> ViewerLike for ViewerAdapter<'a> {
    fn objects(&self) -> &ObjectRegistry {
        &self.state.registry
    }

    fn objects_mut(&mut self) -> &mut ObjectRegistry {
        &mut self.state.registry
    }

    fn camera(&self) -> &Camera {
        &self.state.camera
    }

    fn camera_mut(&mut self) -> &mut Camera {
        &mut self.state.camera
    }

    fn zoom_all(&mut self) {
        if let Some((min, max)) = self.state.registry.extent() {
            self.state.camera.zoom_to(min, max);
            *self.needs_redraw = true;
        }
    }

    fn zoom_on(&mut self, name: &str) {
        if let Some(obj) = self.state.registry.get(name) {
            if let Some((min, max)) = obj.extent() {
                self.state.camera.zoom_to(min, max);
                *self.needs_redraw = true;
            }
        }
    }

    fn center_all(&mut self) {
        if let Some((min, max)) = self.state.registry.extent() {
            self.state.camera.center_to(min, max);
            *self.needs_redraw = true;
        }
    }

    fn center_on(&mut self, name: &str) {
        if let Some(obj) = self.state.registry.get(name) {
            if let Some((min, max)) = obj.extent() {
                self.state.camera.center_to(min, max);
                *self.needs_redraw = true;
            }
        }
    }

    fn reset_view(&mut self) {
        self.state.camera = Camera::new();
        if let Some((min, max)) = self.state.registry.extent() {
            self.state.camera.reset_view(min, max);
            *self.needs_redraw = true;
        }
    }

    fn request_redraw(&mut self) {
        *self.needs_redraw = true;
    }

    fn color_index(&self, name: &str) -> Option<u32> {
        self.state.named_colors.get_by_name(name).map(|(idx, _)| idx)
    }

    fn set_background_color(&mut self, r: f32, g: f32, b: f32) {
        self.state.clear_color[0] = r;
        self.state.clear_color[1] = g;
        self.state.clear_color[2] = b;
        *self.needs_redraw = true;
    }

    fn capture_png(
        &mut self,
        path: &Path,
        width: Option<u32>,
        height: Option<u32>,
    ) -> Result<(), String> {
        let context = self.render_context.ok_or(
            "No render context available. GUI may not be fully initialized yet."
        )?;

        capture_png_to_file(
            path,
            width,
            height,
            context,
            &mut self.state.camera,
            &mut self.state.registry,
            &self.state.settings,
            &self.state.named_colors,
            &self.state.element_colors,
            &self.state.chain_colors,
            self.state.clear_color,
            self.default_size,
        ).map_err(|e| e.to_string())
    }

    fn get_selection(&self, name: &str) -> Option<&str> {
        self.state.selections.get(name).map(|e| e.expression.as_str())
    }

    fn define_selection(&mut self, name: &str, selection: &str) {
        self.state.selections.insert(name.to_string(), SelectionEntry::new(selection.to_string()));
        *self.needs_redraw = true;
    }

    fn remove_selection(&mut self, name: &str) -> bool {
        let removed = self.state.selections.remove(name).is_some();
        if removed {
            *self.needs_redraw = true;
        }
        removed
    }

    fn selection_names(&self) -> Vec<String> {
        self.state.selections.keys().cloned().collect()
    }

    fn set_selection_visible(&mut self, name: &str, visible: bool) {
        if let Some(entry) = self.state.selections.get_mut(name) {
            entry.visible = visible;
            *self.needs_redraw = true;
        }
    }

    fn is_selection_visible(&self, name: &str) -> bool {
        self.state.selections.get(name).map(|e| e.visible).unwrap_or(false)
    }

    fn indicate_selection(&mut self, selection: &str) {
        // For backwards compatibility: create or update "indicate" selection and make it visible
        self.state.selections.insert("indicate".to_string(), SelectionEntry::new(selection.to_string()));
        *self.needs_redraw = true;
    }

    fn clear_indication(&mut self) {
        // Hide all selection indicators
        for entry in self.state.selections.values_mut() {
            entry.visible = false;
        }
        *self.needs_redraw = true;
    }

    fn indicated_selection(&self) -> Option<&str> {
        // Return the first visible selection expression (for backwards compat)
        self.state.selections.values()
            .find(|e| e.visible)
            .map(|e| e.expression.as_str())
    }

    fn request_async_fetch(&mut self, code: &str, name: &str, format: u8) -> bool {
        let fmt = match format {
            1 => pymol_io::FetchFormat::Pdb,
            _ => pymol_io::FetchFormat::Cif,
        };
        self.task_runner.spawn(FetchTask::new(code.to_string(), name.to_string(), fmt));
        true
    }

    // =========================================================================
    // Scene Management
    // =========================================================================

    fn scene_store(&mut self, key: &str, storemask: u32) {
        let mask = SceneStoreMask::from_bits_truncate(storemask);
        self.state.scenes.store(key, mask, &self.state.camera, &self.state.registry);
        *self.needs_redraw = true;
    }

    fn scene_recall(&mut self, key: &str, animate: bool, duration: f32) -> Result<(), String> {
        self.state
            .scenes
            .recall(key, &mut self.state.camera, &mut self.state.registry, animate, duration)
            .map_err(|e| e.to_string())?;
        *self.needs_redraw = true;
        Ok(())
    }

    fn scene_delete(&mut self, key: &str) -> bool {
        self.state.scenes.delete(key).is_some()
    }

    fn scene_list(&self) -> Vec<String> {
        self.state.scenes.list().to_vec()
    }

    fn scene_rename(&mut self, old_key: &str, new_key: &str) -> Result<(), String> {
        self.state.scenes.rename(old_key, new_key).map_err(|e| e.to_string())
    }

    fn scene_clear(&mut self) {
        self.state.scenes.clear();
    }

    // =========================================================================
    // Movie Control
    // =========================================================================

    fn movie_play(&mut self) {
        self.state.movie.play();
        *self.needs_redraw = true;
    }

    fn movie_stop(&mut self) {
        self.state.movie.stop();
        *self.needs_redraw = true;
    }

    fn movie_pause(&mut self) {
        self.state.movie.pause();
    }

    fn movie_toggle(&mut self) {
        self.state.movie.toggle();
        *self.needs_redraw = true;
    }

    fn movie_goto(&mut self, frame: usize) {
        self.state.movie.goto_frame(frame);
        *self.needs_redraw = true;
    }

    fn movie_next(&mut self) {
        self.state.movie.next_frame();
        *self.needs_redraw = true;
    }

    fn movie_prev(&mut self) {
        self.state.movie.prev_frame();
        *self.needs_redraw = true;
    }

    fn movie_set_fps(&mut self, fps: f32) {
        self.state.movie.set_fps(fps);
    }

    fn movie_frame_count(&self) -> usize {
        self.state.movie.frame_count()
    }

    fn movie_current_frame(&self) -> usize {
        self.state.movie.current_frame()
    }

    fn movie_is_playing(&self) -> bool {
        self.state.movie.is_playing()
    }

    fn movie_set_loop_mode(&mut self, mode: u8) {
        let loop_mode = match mode {
            1 => LoopMode::Loop,
            2 => LoopMode::Swing,
            _ => LoopMode::Once,
        };
        self.state.movie.set_loop_mode(loop_mode);
    }

    fn movie_set_frame_count(&mut self, count: usize) {
        self.state.movie.set_frame_count(count);
    }

    // =========================================================================
    // Rock Animation
    // =========================================================================

    fn rock_toggle(&mut self) {
        let current = self.state.movie.is_rock_enabled();
        self.state.movie.set_rock(!current);
        *self.needs_redraw = true;
    }

    fn rock_set(&mut self, enabled: bool) {
        self.state.movie.set_rock(enabled);
        *self.needs_redraw = true;
    }

    fn rock_is_enabled(&self) -> bool {
        self.state.movie.is_rock_enabled()
    }
}
