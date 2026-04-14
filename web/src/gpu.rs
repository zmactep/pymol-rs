//! Async WebGPU initialization from an HTML canvas element.

use wasm_bindgen::JsCast;
use web_sys::HtmlCanvasElement;

use pymol_render::RenderContext;
use pymol_render::silhouette::SilhouettePipeline;

/// GPU resources initialized from a canvas element.
pub struct GpuState {
    pub render_context: RenderContext,
    pub surface: wgpu::Surface<'static>,
    pub surface_config: wgpu::SurfaceConfiguration,
    pub depth_view: wgpu::TextureView,
    pub silhouette: SilhouettePipeline,
}

impl GpuState {
    /// Create GPU resources from a `<canvas>` element ID.
    #[allow(unreachable_code, unused_variables)]
    pub async fn from_canvas(canvas_id: &str) -> Result<Self, String> {
        // Get the canvas element
        let window = web_sys::window().ok_or("No global window")?;
        let document = window.document().ok_or("No document")?;
        let element = document
            .get_element_by_id(canvas_id)
            .ok_or_else(|| format!("Canvas '{}' not found", canvas_id))?;
        let canvas: HtmlCanvasElement = element
            .dyn_into()
            .map_err(|_| "Element is not a canvas")?;

        let width = canvas.client_width().max(1) as u32;
        let height = canvas.client_height().max(1) as u32;

        #[cfg(not(target_arch = "wasm32"))]
        {
            let _ = (canvas, width, height);
            return Err("WebViewer only runs on wasm32".into());
        }

        // Create wgpu instance with WebGPU backend
        let instance = wgpu::Instance::new(&wgpu::InstanceDescriptor {
            backends: wgpu::Backends::BROWSER_WEBGPU,
            ..Default::default()
        });

        // Create surface from canvas (wasm-only API)
        #[cfg(target_arch = "wasm32")]
        let surface = instance
            .create_surface(wgpu::SurfaceTarget::Canvas(canvas))
            .map_err(|e| format!("Failed to create surface: {}", e))?;

        // Request adapter
        let adapter = instance
            .request_adapter(&wgpu::RequestAdapterOptions {
                power_preference: wgpu::PowerPreference::HighPerformance,
                compatible_surface: Some(&surface),
                force_fallback_adapter: false,
            })
            .await
            .map_err(|e| format!("No suitable GPU adapter: {}", e))?;

        log::info!("GPU adapter: {}", adapter.get_info().name);

        // Request device
        let (device, queue) = adapter
            .request_device(
                &wgpu::DeviceDescriptor {
                    label: Some("pymol-web"),
                    required_features: wgpu::Features::empty(),
                    required_limits: adapter.limits(),
                    memory_hints: wgpu::MemoryHints::default(),
                    experimental_features: wgpu::ExperimentalFeatures::default(),
                    trace: wgpu::Trace::Off,
                },
            )
            .await
            .map_err(|e| format!("Failed to create device: {}", e))?;

        // Configure surface
        let caps = surface.get_capabilities(&adapter);
        let format = caps
            .formats
            .iter()
            .find(|f| **f == wgpu::TextureFormat::Bgra8Unorm)
            .or_else(|| {
                caps.formats
                    .iter()
                    .find(|f| **f == wgpu::TextureFormat::Rgba8Unorm)
            })
            .copied()
            .unwrap_or(caps.formats[0]);

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
        surface.configure(&device, &surface_config);

        // Create depth texture
        let depth_view = create_depth_texture(&device, width, height);

        // Create silhouette pipeline
        let silhouette = SilhouettePipeline::new(&device, format);

        // Create render context
        let render_context = RenderContext::new(device, queue, format);

        Ok(Self {
            render_context,
            surface,
            surface_config,
            depth_view,
            silhouette,
        })
    }

    /// Handle canvas resize.
    pub fn resize(&mut self, width: u32, height: u32) {
        let width = width.max(1);
        let height = height.max(1);

        if width == self.surface_config.width && height == self.surface_config.height {
            return;
        }

        self.surface_config.width = width;
        self.surface_config.height = height;
        self.surface
            .configure(self.render_context.device(), &self.surface_config);
        self.depth_view = create_depth_texture(self.render_context.device(), width, height);
    }
}

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
