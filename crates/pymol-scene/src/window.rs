//! Window management for molecular visualization
//!
//! Wraps winit window with wgpu surface management.

use std::sync::Arc;
use winit::dpi::PhysicalSize;
use winit::event_loop::ActiveEventLoop;
use winit::window::WindowAttributes;

use crate::error::WindowError;

/// A window for molecular visualization
///
/// Combines a winit window with a wgpu surface for GPU rendering.
pub struct Window {
    /// The underlying winit window
    window: Arc<winit::window::Window>,
    /// The wgpu surface for rendering
    surface: wgpu::Surface<'static>,
    /// Surface configuration
    surface_config: wgpu::SurfaceConfiguration,
    /// Depth texture view
    depth_view: wgpu::TextureView,
    /// Current window size
    size: (u32, u32),
}

impl Window {
    /// Create a new window
    ///
    /// This creates a winit window and a wgpu surface.
    pub fn new(
        event_loop: &ActiveEventLoop,
        title: &str,
        size: (u32, u32),
        instance: &wgpu::Instance,
        adapter: &wgpu::Adapter,
        device: &wgpu::Device,
    ) -> Result<Self, WindowError> {
        // Create window attributes
        let window_attrs = WindowAttributes::default()
            .with_title(title)
            .with_inner_size(PhysicalSize::new(size.0, size.1));

        // Create the window
        let window = event_loop
            .create_window(window_attrs)
            .map_err(WindowError::OsError)?;
        let window = Arc::new(window);

        // Get actual size (may differ from requested)
        let inner_size = window.inner_size();
        let (width, height) = (inner_size.width.max(1), inner_size.height.max(1));

        // Create surface
        let surface = instance
            .create_surface(window.clone())
            .map_err(|e| WindowError::SurfaceFailed(e.to_string()))?;

        // Get surface capabilities
        let caps = surface.get_capabilities(adapter);
        let format = caps
            .formats
            .iter()
            .find(|f| **f == wgpu::TextureFormat::Bgra8Unorm)
            .or_else(|| caps.formats.iter().find(|f| **f == wgpu::TextureFormat::Rgba8Unorm))
            .copied()
            .unwrap_or(caps.formats[0]);

        // Configure surface
        let surface_config = wgpu::SurfaceConfiguration {
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT,
            format,
            width,
            height,
            present_mode: wgpu::PresentMode::AutoVsync,
            alpha_mode: caps.alpha_modes[0],
            view_formats: vec![],
            desired_maximum_frame_latency: 2,
        };
        surface.configure(device, &surface_config);

        // Create depth texture
        let depth_view = create_depth_texture(device, width, height);

        Ok(Self {
            window,
            surface,
            surface_config,
            depth_view,
            size: (width, height),
        })
    }

    /// Get a reference to the underlying winit window
    pub fn window(&self) -> &winit::window::Window {
        &self.window
    }

    /// Get the window ID
    pub fn id(&self) -> winit::window::WindowId {
        self.window.id()
    }

    /// Get the current window size
    pub fn size(&self) -> (u32, u32) {
        self.size
    }

    /// Get the width
    pub fn width(&self) -> u32 {
        self.size.0
    }

    /// Get the height
    pub fn height(&self) -> u32 {
        self.size.1
    }

    /// Get the aspect ratio
    pub fn aspect_ratio(&self) -> f32 {
        self.size.0 as f32 / self.size.1 as f32
    }

    /// Get the surface format
    pub fn surface_format(&self) -> wgpu::TextureFormat {
        self.surface_config.format
    }

    /// Get the depth texture view
    pub fn depth_view(&self) -> &wgpu::TextureView {
        &self.depth_view
    }

    /// Handle window resize
    pub fn resize(&mut self, new_size: (u32, u32), device: &wgpu::Device) {
        let (width, height) = (new_size.0.max(1), new_size.1.max(1));

        if width == self.size.0 && height == self.size.1 {
            return;
        }

        self.size = (width, height);
        self.surface_config.width = width;
        self.surface_config.height = height;
        self.surface.configure(device, &self.surface_config);
        self.depth_view = create_depth_texture(device, width, height);
    }

    /// Get the current surface texture for rendering
    pub fn get_current_texture(&self) -> Result<wgpu::SurfaceTexture, wgpu::SurfaceError> {
        self.surface.get_current_texture()
    }

    /// Request a redraw
    pub fn request_redraw(&self) {
        self.window.request_redraw();
    }

    /// Set the window title
    pub fn set_title(&self, title: &str) {
        self.window.set_title(title);
    }

    /// Check if the window is minimized
    pub fn is_minimized(&self) -> bool {
        self.window
            .is_minimized()
            .unwrap_or(self.size.0 == 0 || self.size.1 == 0)
    }

    /// Set cursor visibility
    pub fn set_cursor_visible(&self, visible: bool) {
        self.window.set_cursor_visible(visible);
    }

    /// Set cursor grab mode
    pub fn set_cursor_grab(
        &self,
        mode: winit::window::CursorGrabMode,
    ) -> Result<(), winit::error::ExternalError> {
        self.window.set_cursor_grab(mode)
    }

    /// Get the scale factor (for HiDPI displays)
    pub fn scale_factor(&self) -> f64 {
        self.window.scale_factor()
    }

    /// Convert physical position to logical position
    pub fn physical_to_logical(&self, pos: (f32, f32)) -> (f32, f32) {
        let scale = self.scale_factor() as f32;
        (pos.0 / scale, pos.1 / scale)
    }

    /// Get normalized device coordinates from physical position
    ///
    /// Returns (x, y) in range [-1, 1] with y-up.
    pub fn to_ndc(&self, pos: (f32, f32)) -> (f32, f32) {
        let x = (pos.0 / self.size.0 as f32) * 2.0 - 1.0;
        let y = 1.0 - (pos.1 / self.size.1 as f32) * 2.0;
        (x, y)
    }
}

/// Create a depth texture for the given dimensions
fn create_depth_texture(device: &wgpu::Device, width: u32, height: u32) -> wgpu::TextureView {
    let texture = device.create_texture(&wgpu::TextureDescriptor {
        label: Some("Depth Texture"),
        size: wgpu::Extent3d {
            width,
            height,
            depth_or_array_layers: 1,
        },
        mip_level_count: 1,
        sample_count: 1,
        dimension: wgpu::TextureDimension::D2,
        format: wgpu::TextureFormat::Depth32Float,
        usage: wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::TEXTURE_BINDING,
        view_formats: &[],
    });

    texture.create_view(&wgpu::TextureViewDescriptor::default())
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_to_ndc() {
        // Mock window with size 800x600
        // Since we can't create a real Window in tests, test the calculation directly
        let width = 800u32;
        let height = 600u32;

        // Helper function that matches Window::to_ndc
        let to_ndc = |pos: (f32, f32)| -> (f32, f32) {
            let x = (pos.0 / width as f32) * 2.0 - 1.0;
            let y = 1.0 - (pos.1 / height as f32) * 2.0;
            (x, y)
        };

        // Center should be (0, 0)
        let (x, y) = to_ndc((400.0, 300.0));
        assert!((x - 0.0).abs() < 0.001);
        assert!((y - 0.0).abs() < 0.001);

        // Top-left should be (-1, 1)
        let (x, y) = to_ndc((0.0, 0.0));
        assert!((x - (-1.0)).abs() < 0.001);
        assert!((y - 1.0).abs() < 0.001);

        // Bottom-right should be (1, -1)
        let (x, y) = to_ndc((800.0, 600.0));
        assert!((x - 1.0).abs() < 0.001);
        assert!((y - (-1.0)).abs() < 0.001);
    }
}
