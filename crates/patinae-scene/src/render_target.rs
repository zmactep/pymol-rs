//! `CaptureRenderer` — abstract render-target trait used by `SessionAdapter`
//! and the command kernel.
//!
//! Hosts plug in their own implementation: the interactive viewport uses
//! `patinae_render::RenderState`, while headless callers can provide any
//! renderer that can produce a PNG from the current scene state.

use std::path::Path;
use std::sync::Arc;

use patinae_color::{NamedPalette, ThemedPalette};
use patinae_settings::Settings;

use crate::camera::Camera;
use crate::error::ViewerError;
use crate::object::ObjectRegistry;

#[cfg(feature = "render-bridge")]
use patinae_render::{DisplayedGeometry, GeometryExportOptions};

/// Headless-renderable target. The host owns the GPU resources and
/// implements `capture_png` however it wishes.
pub trait CaptureRenderer {
    /// Shared device, exposed for commands that need raw GPU handles.
    fn gpu_device(&self) -> &Arc<wgpu::Device>;
    /// Shared queue, exposed for commands that need raw GPU handles.
    fn gpu_queue(&self) -> &Arc<wgpu::Queue>;

    /// Render the current scene state into a PNG file. Width/height come
    /// from the caller (defaults from the viewport when unset upstream).
    /// Implementations are free to mutate `camera` (animation interpolation
    /// already happens upstream in `capture_frame_png`) and `registry`
    /// (rep rebuilds on dirty flags).
    ///
    /// # Errors
    ///
    /// Returns an error if the renderer cannot capture or write the image.
    #[allow(clippy::too_many_arguments)]
    fn capture_png(
        &mut self,
        path: &Path,
        width: u32,
        height: u32,
        camera: &mut Camera,
        registry: &mut ObjectRegistry,
        settings: &Settings,
        named: &NamedPalette,
        themed: &ThemedPalette,
        clear_color: [f32; 3],
    ) -> Result<(), ViewerError>;

    /// Export the current displayed scene as renderer-neutral geometry.
    ///
    /// # Errors
    ///
    /// Returns an error if the renderer does not support geometry export.
    #[cfg(feature = "render-bridge")]
    #[allow(clippy::too_many_arguments)]
    fn export_displayed_geometry(
        &mut self,
        camera: &mut Camera,
        registry: &mut ObjectRegistry,
        settings: &Settings,
        named: &NamedPalette,
        themed: &ThemedPalette,
        clear_color: [f32; 3],
        options: &GeometryExportOptions,
    ) -> Result<DisplayedGeometry, ViewerError> {
        let _ = (
            camera,
            registry,
            settings,
            named,
            themed,
            clear_color,
            options,
        );
        Err(ViewerError::capture_error(
            "Displayed geometry export not supported by this renderer",
        ))
    }
}
