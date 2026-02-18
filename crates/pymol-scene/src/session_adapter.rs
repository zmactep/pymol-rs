//! SessionAdapter — bridges Session to the ViewerLike trait
//!
//! [`SessionAdapter`] borrows a [`Session`] and an optional [`RenderContext`]
//! and implements [`ViewerLike`], allowing the command system to operate on
//! pure scene state with optional GPU support.

use std::path::Path;

use pymol_color::NamedColors;
use pymol_render::{ColorResolver, RenderContext};
use pymol_settings::GlobalSettings;

use crate::camera::Camera;
use crate::capture::capture_png_to_file;
use crate::movie::Movie;
use crate::object::ObjectRegistry;
use crate::raytrace::{raytrace_scene, raytrace_to_file, RaytraceInput};
use crate::scene::{SceneManager, SceneStoreMask};
use crate::selection::SelectionManager;
use crate::session::Session;
use crate::view::ViewManager;
use crate::viewer_trait::{RaytracedImage, ViewerLike};

/// Adapter that wraps a [`Session`] to implement [`ViewerLike`] for command execution.
///
/// Borrows a mutable reference to a `Session` (scene state) plus an optional
/// `RenderContext` (GPU resources) and a `needs_redraw` flag.
pub struct SessionAdapter<'a> {
    /// The scene state
    pub session: &'a mut Session,
    /// GPU render context (None if headless / not yet initialized)
    pub render_context: Option<&'a RenderContext>,
    /// Default viewport size for capture when width/height not specified
    pub default_size: (u32, u32),
    /// Redraw flag — set to true when a re-render is needed
    pub needs_redraw: &'a mut bool,
}

impl<'a> ViewerLike for SessionAdapter<'a> {
    // =========================================================================
    // Required Accessors
    // =========================================================================

    fn objects(&self) -> &ObjectRegistry { &self.session.registry }
    fn objects_mut(&mut self) -> &mut ObjectRegistry { &mut self.session.registry }
    fn camera(&self) -> &Camera { &self.session.camera }
    fn camera_mut(&mut self) -> &mut Camera { &mut self.session.camera }
    fn settings(&self) -> &GlobalSettings { &self.session.settings }
    fn settings_mut(&mut self) -> &mut GlobalSettings { &mut self.session.settings }
    fn movie(&self) -> &Movie { &self.session.movie }
    fn movie_mut(&mut self) -> &mut Movie { &mut self.session.movie }
    fn scenes(&self) -> &SceneManager { &self.session.scenes }
    fn scenes_mut(&mut self) -> &mut SceneManager { &mut self.session.scenes }
    fn views(&self) -> &ViewManager { &self.session.views }
    fn views_mut(&mut self) -> &mut ViewManager { &mut self.session.views }
    fn selections(&self) -> &SelectionManager { &self.session.selections }
    fn selections_mut(&mut self) -> &mut SelectionManager { &mut self.session.selections }
    fn named_colors(&self) -> &NamedColors { &self.session.named_colors }
    fn clear_color(&self) -> [f32; 3] { self.session.clear_color }
    fn set_clear_color(&mut self, color: [f32; 3]) { self.session.clear_color = color; }
    fn raytraced_image_ref(&self) -> Option<&RaytracedImage> { self.session.raytraced_image.as_ref() }
    fn set_raytraced_image_internal(&mut self, image: Option<RaytracedImage>) { self.session.raytraced_image = image; }

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
        let context = self.render_context.ok_or_else(|| {
            "No render context available".to_string()
        })?;

        capture_png_to_file(
            path,
            width,
            height,
            context,
            &mut self.session.camera,
            &mut self.session.registry,
            &self.session.settings,
            &self.session.named_colors,
            &self.session.element_colors,
            &self.session.chain_colors,
            self.session.clear_color,
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

        let context = self.render_context.ok_or_else(|| {
            "No render context available".to_string()
        })?;

        let mut input = RaytraceInput {
            device: context.device(),
            queue: context.queue(),
            camera: &mut self.session.camera,
            registry: &self.session.registry,
            settings: &self.session.settings,
            named_colors: &self.session.named_colors,
            element_colors: &self.session.element_colors,
            chain_colors: &self.session.chain_colors,
            clear_color: self.session.clear_color,
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

        let context = self.render_context.ok_or_else(|| {
            "No render context available".to_string()
        })?;

        let mut input = RaytraceInput {
            device: context.device(),
            queue: context.queue(),
            camera: &mut self.session.camera,
            registry: &self.session.registry,
            settings: &self.session.settings,
            named_colors: &self.session.named_colors,
            element_colors: &self.session.element_colors,
            chain_colors: &self.session.chain_colors,
            clear_color: self.session.clear_color,
            default_size: self.default_size,
        };

        raytrace_to_file(&mut input, path, width, height, antialias)
            .map_err(|e| e.to_string())
    }

    fn capture_frame_png(
        &mut self,
        frame: usize,
        path: &Path,
        width: Option<u32>,
        height: Option<u32>,
    ) -> Result<(), String> {
        self.session.movie.goto_frame(frame);
        if let Some(view) = self.session.movie.interpolated_view() {
            self.session.camera.set_view(view);
        }
        if let Some(scene_name) = self.session.movie.current_scene_name().map(|s| s.to_string()) {
            if let Some(scene) = self.session.scenes.get(&scene_name) {
                scene.apply(&mut self.session.camera, &mut self.session.registry, false, 0.0);
                if let Some(view) = self.session.movie.interpolated_view() {
                    self.session.camera.set_view(view);
                }
            }
        }
        self.capture_png(path, width, height)
    }

    // =========================================================================
    // Per-impl methods (borrow checker constraints)
    // =========================================================================

    fn scene_store(&mut self, key: &str, storemask: u32) {
        let mask = SceneStoreMask::from_bits_truncate(storemask);
        self.session.scenes.store(key, mask, &self.session.camera, &self.session.registry);
        *self.needs_redraw = true;
    }

    fn scene_recall(&mut self, key: &str, animate: bool, duration: f32) -> Result<(), String> {
        self.session
            .scenes
            .recall(key, &mut self.session.camera, &mut self.session.registry, animate, duration)
            .map_err(|e| e.to_string())?;
        *self.needs_redraw = true;
        Ok(())
    }

    fn view_recall(&mut self, key: &str, animate: f32) -> Result<(), String> {
        self.session
            .views
            .recall(key, &mut self.session.camera, animate)
            .map_err(|e| e.to_string())?;
        *self.needs_redraw = true;
        Ok(())
    }

    fn viewport_size(&self) -> (u32, u32) {
        self.default_size
    }
}

impl<'a> SessionAdapter<'a> {
    /// Ensure all representations are built before raytracing.
    fn prepare_representations_for_raytrace(&mut self) -> Result<(), String> {
        let context = self.render_context.ok_or_else(|| {
            "No render context available".to_string()
        })?;

        let names: Vec<String> = self.session.registry.names().map(|s| s.to_string()).collect();
        for name in &names {
            let color_resolver = ColorResolver::new(
                &self.session.named_colors,
                &self.session.element_colors,
                &self.session.chain_colors,
            );
            if let Some(mol_obj) = self.session.registry.get_molecule_mut(name) {
                mol_obj.prepare_render(context, &color_resolver, &self.session.settings);
            }
        }

        Ok(())
    }
}
