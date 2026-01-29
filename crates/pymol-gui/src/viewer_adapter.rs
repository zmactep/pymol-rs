//! ViewerAdapter - bridges App state to the ViewerLike trait for commands
//!
//! This module provides the adapter that allows the command system to interact
//! with the GUI's state through the `ViewerLike` trait.

use std::collections::HashMap;
use std::path::Path;

use pymol_cmd::ViewerLike;
use pymol_color::{ChainColors, ElementColors, NamedColors};
use pymol_render::RenderContext;
use pymol_scene::{capture_png_to_file, Camera, ObjectRegistry, SelectionEntry};
use pymol_settings::GlobalSettings;

use crate::async_tasks::TaskRunner;
use crate::fetch::FetchTask;

/// Adapter that wraps App's fields to implement ViewerLike for command execution
///
/// This struct borrows references to the various components of the App and
/// implements the `ViewerLike` trait, allowing commands to interact with
/// the application state in a uniform way.
pub struct ViewerAdapter<'a> {
    /// Reference to the object registry
    pub registry: &'a mut ObjectRegistry,
    /// Reference to the camera
    pub camera: &'a mut Camera,
    /// Reference to named colors
    pub named_colors: &'a NamedColors,
    /// Reference to the background/clear color
    pub clear_color: &'a mut [f32; 3],
    /// Flag to indicate if a redraw is needed
    pub needs_redraw: &'a mut bool,
    /// Reference to named selections (with visibility state)
    pub selections: &'a mut HashMap<String, SelectionEntry>,
    /// Reference to task runner for background operations
    pub task_runner: &'a TaskRunner,
    // --- Fields for PNG capture ---
    /// Reference to render context for GPU operations (None if not yet initialized)
    pub render_context: Option<&'a RenderContext>,
    /// Reference to global settings
    pub settings: &'a GlobalSettings,
    /// Reference to element colors
    pub element_colors: &'a ElementColors,
    /// Reference to chain colors
    pub chain_colors: &'a ChainColors,
    /// Default size for PNG capture when width/height not specified
    pub default_size: (u32, u32),
}

impl<'a> ViewerLike for ViewerAdapter<'a> {
    fn objects(&self) -> &ObjectRegistry {
        self.registry
    }

    fn objects_mut(&mut self) -> &mut ObjectRegistry {
        self.registry
    }

    fn camera(&self) -> &Camera {
        self.camera
    }

    fn camera_mut(&mut self) -> &mut Camera {
        self.camera
    }

    fn zoom_all(&mut self) {
        if let Some((min, max)) = self.registry.extent() {
            self.camera.zoom_to(min, max);
            *self.needs_redraw = true;
        }
    }

    fn zoom_on(&mut self, name: &str) {
        if let Some(obj) = self.registry.get(name) {
            if let Some((min, max)) = obj.extent() {
                self.camera.zoom_to(min, max);
                *self.needs_redraw = true;
            }
        }
    }

    fn center_all(&mut self) {
        if let Some((min, max)) = self.registry.extent() {
            self.camera.center_to(min, max);
            *self.needs_redraw = true;
        }
    }

    fn center_on(&mut self, name: &str) {
        if let Some(obj) = self.registry.get(name) {
            if let Some((min, max)) = obj.extent() {
                self.camera.center_to(min, max);
                *self.needs_redraw = true;
            }
        }
    }

    fn reset_view(&mut self) {
        *self.camera = Camera::new();
        if let Some((min, max)) = self.registry.extent() {
            self.camera.reset_view(min, max);
            *self.needs_redraw = true;
        }
    }

    fn request_redraw(&mut self) {
        *self.needs_redraw = true;
    }

    fn color_index(&self, name: &str) -> Option<u32> {
        self.named_colors.get_by_name(name).map(|(idx, _)| idx)
    }

    fn set_background_color(&mut self, r: f32, g: f32, b: f32) {
        self.clear_color[0] = r;
        self.clear_color[1] = g;
        self.clear_color[2] = b;
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
            self.camera,
            self.registry,
            self.settings,
            self.named_colors,
            self.element_colors,
            self.chain_colors,
            *self.clear_color,
            self.default_size,
        ).map_err(|e| e.to_string())
    }

    fn get_selection(&self, name: &str) -> Option<&str> {
        self.selections.get(name).map(|e| e.expression.as_str())
    }

    fn define_selection(&mut self, name: &str, selection: &str) {
        self.selections.insert(name.to_string(), SelectionEntry::new(selection.to_string()));
        *self.needs_redraw = true;
    }

    fn remove_selection(&mut self, name: &str) -> bool {
        let removed = self.selections.remove(name).is_some();
        if removed {
            *self.needs_redraw = true;
        }
        removed
    }

    fn selection_names(&self) -> Vec<String> {
        self.selections.keys().cloned().collect()
    }

    fn set_selection_visible(&mut self, name: &str, visible: bool) {
        if let Some(entry) = self.selections.get_mut(name) {
            entry.visible = visible;
            *self.needs_redraw = true;
        }
    }

    fn is_selection_visible(&self, name: &str) -> bool {
        self.selections.get(name).map(|e| e.visible).unwrap_or(false)
    }

    fn indicate_selection(&mut self, selection: &str) {
        // For backwards compatibility: create or update "indicate" selection and make it visible
        self.selections.insert("indicate".to_string(), SelectionEntry::new(selection.to_string()));
        *self.needs_redraw = true;
    }

    fn clear_indication(&mut self) {
        // Hide all selection indicators
        for entry in self.selections.values_mut() {
            entry.visible = false;
        }
        *self.needs_redraw = true;
    }

    fn indicated_selection(&self) -> Option<&str> {
        // Return the first visible selection expression (for backwards compat)
        self.selections.values()
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
}
