//! ViewerAdapter - bridges App state to the ViewerLike trait for commands
//!
//! This module provides the adapter that allows the command system to interact
//! with the GUI's state through the `ViewerLike` trait.

use std::path::Path;

use pymol_cmd::ViewerLike;
use pymol_render::RenderContext;
use pymol_render::ColorResolver;
use pymol_scene::{
    capture_png_to_file, raytrace_scene, raytrace_to_file, Camera, Movie, NamedColors,
    ObjectRegistry, RaytracedImage, RaytraceInput, SceneManager, SceneStoreMask,
    SelectionManager, ViewManager,
};

use pymol_scene::Session;

use crate::async_tasks::TaskRunner;
use crate::fetch::FetchTask;

/// Adapter that wraps Session to implement ViewerLike for command execution
///
/// Borrows scene state and provides GUI-specific overrides (async fetch,
/// render context access) on top of the standard `ViewerLike` interface.
pub struct ViewerAdapter<'a> {
    /// Reference to the scene state
    pub state: &'a mut Session,
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
    // =========================================================================
    // Required Accessors
    // =========================================================================

    fn objects(&self) -> &ObjectRegistry { &self.state.registry }
    fn objects_mut(&mut self) -> &mut ObjectRegistry { &mut self.state.registry }
    fn camera(&self) -> &Camera { &self.state.camera }
    fn camera_mut(&mut self) -> &mut Camera { &mut self.state.camera }
    fn settings(&self) -> &pymol_settings::GlobalSettings { &self.state.settings }
    fn settings_mut(&mut self) -> &mut pymol_settings::GlobalSettings { &mut self.state.settings }
    fn movie(&self) -> &Movie { &self.state.movie }
    fn movie_mut(&mut self) -> &mut Movie { &mut self.state.movie }
    fn scenes(&self) -> &SceneManager { &self.state.scenes }
    fn scenes_mut(&mut self) -> &mut SceneManager { &mut self.state.scenes }
    fn views(&self) -> &ViewManager { &self.state.views }
    fn views_mut(&mut self) -> &mut ViewManager { &mut self.state.views }
    fn selections(&self) -> &SelectionManager { &self.state.selections }
    fn selections_mut(&mut self) -> &mut SelectionManager { &mut self.state.selections }
    fn named_colors(&self) -> &NamedColors { &self.state.named_colors }
    fn clear_color(&self) -> [f32; 3] { self.state.clear_color }
    fn set_clear_color(&mut self, color: [f32; 3]) { self.state.clear_color = color; }
    fn raytraced_image_ref(&self) -> Option<&RaytracedImage> { self.state.raytraced_image.as_ref() }
    fn set_raytraced_image_internal(&mut self, image: Option<RaytracedImage>) { self.state.raytraced_image = image; }

    fn request_redraw(&mut self) {
        *self.needs_redraw = true;
    }

    // =========================================================================
    // GPU-specific overrides
    // =========================================================================

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

    fn raytrace(
        &mut self,
        width: Option<u32>,
        height: Option<u32>,
        antialias: u32,
    ) -> Result<Vec<u8>, String> {
        self.prepare_representations_for_raytrace()?;

        let context = self.render_context.ok_or(
            "No render context available. GUI may not be fully initialized yet."
        )?;

        let mut input = RaytraceInput {
            device: context.device(),
            queue: context.queue(),
            camera: &mut self.state.camera,
            registry: &self.state.registry,
            settings: &self.state.settings,
            named_colors: &self.state.named_colors,
            element_colors: &self.state.element_colors,
            chain_colors: &self.state.chain_colors,
            clear_color: self.state.clear_color,
            default_size: self.default_size,
        };

        raytrace_scene(&mut input, width, height, antialias)
            .map(|(data, _, _)| data)
            .map_err(|e| e.to_string())
    }

    fn raytrace_to_file(
        &mut self,
        path: &Path,
        width: Option<u32>,
        height: Option<u32>,
        antialias: u32,
    ) -> Result<(u32, u32), String> {
        self.prepare_representations_for_raytrace()?;

        let context = self.render_context.ok_or(
            "No render context available. GUI may not be fully initialized yet."
        )?;

        let mut input = RaytraceInput {
            device: context.device(),
            queue: context.queue(),
            camera: &mut self.state.camera,
            registry: &self.state.registry,
            settings: &self.state.settings,
            named_colors: &self.state.named_colors,
            element_colors: &self.state.element_colors,
            chain_colors: &self.state.chain_colors,
            clear_color: self.state.clear_color,
            default_size: self.default_size,
        };

        raytrace_to_file(&mut input, path, width, height, antialias)
            .map_err(|e| e.to_string())
    }

    // =========================================================================
    // GUI-specific override
    // =========================================================================

    fn request_async_fetch(&mut self, code: &str, name: &str, format: u8) -> bool {
        let fmt = match format {
            1 => pymol_io::FetchFormat::Pdb,
            _ => pymol_io::FetchFormat::Cif,
        };
        self.task_runner.spawn(FetchTask::new(code.to_string(), name.to_string(), fmt));
        true
    }

    // =========================================================================
    // Per-impl methods (borrow checker constraints)
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
            .map_err(|e: pymol_scene::SceneError| e.to_string())?;
        *self.needs_redraw = true;
        Ok(())
    }

    fn view_recall(&mut self, key: &str, animate: f32) -> Result<(), String> {
        self.state
            .views
            .recall(key, &mut self.state.camera, animate)
            .map_err(|e| e.to_string())?;
        *self.needs_redraw = true;
        Ok(())
    }

    fn capture_frame_png(
        &mut self,
        frame: usize,
        path: &Path,
        width: Option<u32>,
        height: Option<u32>,
    ) -> Result<(), String> {
        self.state.movie.goto_frame(frame);
        if let Some(view) = self.state.movie.interpolated_view() {
            self.state.camera.set_view(view);
        }
        if let Some(scene_name) = self.state.movie.current_scene_name().map(|s| s.to_string()) {
            if let Some(scene) = self.state.scenes.get(&scene_name) {
                scene.apply(&mut self.state.camera, &mut self.state.registry, false, 0.0);
                if let Some(view) = self.state.movie.interpolated_view() {
                    self.state.camera.set_view(view);
                }
            }
        }
        self.capture_png(path, width, height)
    }

    // =========================================================================
    // Viewport override
    // =========================================================================

    fn viewport_size(&self) -> (u32, u32) {
        self.default_size
    }
}

impl<'a> ViewerAdapter<'a> {
    /// Ensure all representations are built before raytracing.
    ///
    /// This is needed when raytracing from scripts where the render loop hasn't run.
    fn prepare_representations_for_raytrace(&mut self) -> Result<(), String> {
        let context = self.render_context.ok_or(
            "No render context available. GUI may not be fully initialized yet."
        )?;

        let names: Vec<String> = self.state.registry.names().map(|s| s.to_string()).collect();
        for name in &names {
            let color_resolver = ColorResolver::new(
                &self.state.named_colors,
                &self.state.element_colors,
                &self.state.chain_colors,
            );
            if let Some(mol_obj) = self.state.registry.get_molecule_mut(name) {
                mol_obj.prepare_render(context, &color_resolver, &self.state.settings);
            }
        }

        Ok(())
    }
}
