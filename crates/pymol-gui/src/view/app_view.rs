//! Application View
//!
//! Contains all GPU and windowing resources for rendering:
//! wgpu instance, window, surface, egui integration, and UI configuration.

use std::sync::Arc;

use egui::ViewportId;
use pymol_render::RenderContext;
use winit::dpi::PhysicalSize;
use winit::window::Window;

use super::UiConfig;

/// Application view containing all rendering resources
pub struct AppView {
    // =========================================================================
    // GPU Resources
    // =========================================================================
    /// wgpu instance
    pub instance: wgpu::Instance,
    /// Render context for molecular rendering (owns device and queue)
    pub render_context: Option<RenderContext>,

    // =========================================================================
    // Window State
    // =========================================================================
    /// The main window
    pub window: Option<Arc<Window>>,
    /// Window surface
    pub surface: Option<wgpu::Surface<'static>>,
    /// Surface configuration
    pub surface_config: Option<wgpu::SurfaceConfiguration>,
    /// Depth texture view
    pub depth_view: Option<wgpu::TextureView>,

    // =========================================================================
    // egui Integration
    // =========================================================================
    /// egui context
    pub egui_ctx: egui::Context,
    /// egui-winit state
    pub egui_state: Option<egui_winit::State>,
    /// egui-wgpu renderer
    pub egui_renderer: Option<egui_wgpu::Renderer>,
    /// Viewport rect for 3D rendering (in physical pixels)
    pub viewport_rect: Option<egui::Rect>,

    // =========================================================================
    // UI Configuration
    // =========================================================================
    /// UI layout configuration
    pub ui_config: UiConfig,
}

impl Default for AppView {
    fn default() -> Self {
        Self::new()
    }
}

impl AppView {
    /// Create a new application view with default values
    pub fn new() -> Self {
        let instance = wgpu::Instance::new(&wgpu::InstanceDescriptor {
            backends: wgpu::Backends::all(),
            ..Default::default()
        });

        Self {
            instance,
            render_context: None,
            window: None,
            surface: None,
            surface_config: None,
            depth_view: None,
            egui_ctx: egui::Context::default(),
            egui_state: None,
            egui_renderer: None,
            viewport_rect: None,
            ui_config: UiConfig::new(),
        }
    }

    /// Initialize GPU resources
    pub async fn init_gpu(&mut self, window: Arc<Window>) -> Result<(), String> {
        // Create surface
        let surface = self
            .instance
            .create_surface(window.clone())
            .map_err(|e| format!("Failed to create surface: {}", e))?;

        // Request adapter
        let adapter = self
            .instance
            .request_adapter(&wgpu::RequestAdapterOptions {
                power_preference: wgpu::PowerPreference::HighPerformance,
                compatible_surface: Some(&surface),
                force_fallback_adapter: false,
            })
            .await
            .map_err(|e| format!("No suitable GPU adapter found: {}", e))?;

        // Request device
        let (device, queue) = adapter
            .request_device(&wgpu::DeviceDescriptor {
                label: Some("PyMOL-RS Device"),
                required_features: wgpu::Features::empty(),
                required_limits: wgpu::Limits::default(),
                memory_hints: wgpu::MemoryHints::Performance,
                experimental_features: wgpu::ExperimentalFeatures::default(),
                trace: wgpu::Trace::Off,
            })
            .await
            .map_err(|e| format!("Failed to create device: {}", e))?;

        // Get surface capabilities
        let caps = surface.get_capabilities(&adapter);
        let format = caps
            .formats
            .iter()
            .find(|f| **f == wgpu::TextureFormat::Bgra8Unorm)
            .or_else(|| caps.formats.iter().find(|f| **f == wgpu::TextureFormat::Rgba8Unorm))
            .copied()
            .unwrap_or(caps.formats[0]);

        let size = window.inner_size();
        let config = wgpu::SurfaceConfiguration {
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT,
            format,
            width: size.width.max(1),
            height: size.height.max(1),
            present_mode: wgpu::PresentMode::AutoVsync,
            alpha_mode: caps.alpha_modes[0],
            view_formats: vec![],
            desired_maximum_frame_latency: 2,
        };
        surface.configure(&device, &config);

        // Create depth texture
        let depth_view = Self::create_depth_texture(&device, size.width, size.height);

        // Initialize egui
        let egui_state = egui_winit::State::new(
            self.egui_ctx.clone(),
            ViewportId::ROOT,
            &window,
            Some(window.scale_factor() as f32),
            None,
            None,
        );

        let egui_renderer = egui_wgpu::Renderer::new(&device, format, egui_wgpu::RendererOptions {
            depth_stencil_format: None,
            msaa_samples: 1,
            dithering: false,
            ..Default::default()
        });

        // Create render context (takes ownership of device and queue)
        let render_context = RenderContext::new(device, queue, format);

        // Store everything
        self.window = Some(window);
        self.surface = Some(surface);
        self.surface_config = Some(config);
        self.depth_view = Some(depth_view);
        self.render_context = Some(render_context);
        self.egui_state = Some(egui_state);
        self.egui_renderer = Some(egui_renderer);

        Ok(())
    }

    /// Get adapter info (returns device name and backend)
    pub fn adapter_info(&self) -> Option<(String, String)> {
        // This is called after init_gpu, so we can get the info from the instance
        // However, the adapter is consumed. For now, we'll return None
        // The info logging should be done during init_gpu
        None
    }

    /// Create depth texture
    pub fn create_depth_texture(device: &wgpu::Device, width: u32, height: u32) -> wgpu::TextureView {
        let texture = device.create_texture(&wgpu::TextureDescriptor {
            label: Some("Depth Texture"),
            size: wgpu::Extent3d {
                width: width.max(1),
                height: height.max(1),
                depth_or_array_layers: 1,
            },
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: wgpu::TextureFormat::Depth32Float,
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT,
            view_formats: &[],
        });
        texture.create_view(&wgpu::TextureViewDescriptor::default())
    }

    /// Handle window resize
    pub fn resize(&mut self, new_size: PhysicalSize<u32>) {
        if new_size.width == 0 || new_size.height == 0 {
            return;
        }

        if let (Some(surface), Some(config), Some(context)) =
            (&self.surface, &mut self.surface_config, &self.render_context)
        {
            let device = context.device();
            config.width = new_size.width;
            config.height = new_size.height;
            surface.configure(device, config);

            // Recreate depth texture
            self.depth_view = Some(Self::create_depth_texture(device, new_size.width, new_size.height));
        }
    }

    /// Get the current window size
    pub fn window_size(&self) -> Option<PhysicalSize<u32>> {
        self.window.as_ref().map(|w| w.inner_size())
    }

    /// Get the scale factor
    pub fn scale_factor(&self) -> f32 {
        self.window.as_ref().map(|w| w.scale_factor()).unwrap_or(1.0) as f32
    }

    /// Show the window
    pub fn show_window(&self) {
        if let Some(window) = &self.window {
            window.set_visible(true);
        }
    }

    /// Hide the window
    pub fn hide_window(&self) {
        if let Some(window) = &self.window {
            window.set_visible(false);
        }
    }

    /// Check if window is visible
    pub fn is_window_visible(&self) -> Option<bool> {
        self.window.as_ref().and_then(|w| w.is_visible())
    }

    /// Request a window redraw
    pub fn request_redraw(&self) {
        if let Some(window) = &self.window {
            window.request_redraw();
        }
    }

    /// Check if mouse position is over the viewport
    pub fn is_over_viewport(&self, mouse_pos: (f32, f32)) -> bool {
        self.viewport_rect
            .map(|rect| rect.contains(egui::pos2(mouse_pos.0, mouse_pos.1)))
            .unwrap_or(true) // Default to true if no viewport rect yet
    }
}
