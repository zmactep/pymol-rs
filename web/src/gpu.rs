//! Async WebGPU initialization from an HTML canvas element.
//!
//! Hosts `patinae_render::RenderState` for the frame pipeline. The surface
//! owns the canvas-backed swapchain.

#[cfg(target_arch = "wasm32")]
use std::sync::Arc;

use wasm_bindgen::JsCast;
use web_sys::HtmlCanvasElement;

#[cfg(target_arch = "wasm32")]
use patinae_render::{
    required_limits_for_memory_policy, select_render_memory_policy, RenderMemorySelectionInput,
};
use patinae_render::{RenderConfig, RenderState};

/// GPU resources initialized from a canvas element. The renderer is the
/// shared `patinae_render::RenderState`; the surface + config drive the
/// canvas swapchain.
pub struct GpuState {
    pub state: RenderState,
    pub surface: wgpu::Surface<'static>,
    pub surface_config: wgpu::SurfaceConfiguration,
}

impl GpuState {
    /// Create GPU resources from a `<canvas>` element ID.
    #[allow(unreachable_code, unused_variables)]
    pub async fn from_canvas(canvas_id: &str, config: RenderConfig) -> Result<Self, String> {
        // Get the canvas element
        let window = web_sys::window().ok_or("No global window")?;
        let document = window.document().ok_or("No document")?;
        let element = document
            .get_element_by_id(canvas_id)
            .ok_or_else(|| format!("Canvas '{}' not found", canvas_id))?;
        let canvas: HtmlCanvasElement = element
            .dyn_into()
            .map_err(|_| "Element is not a canvas")?;

        let width = canvas.width().max(1);
        let height = canvas.height().max(1);

        #[cfg(not(target_arch = "wasm32"))]
        {
            let _ = (canvas, width, height);
            return Err("WebViewer only runs on wasm32".into());
        }

        #[cfg(target_arch = "wasm32")]
        {
        let mut config = config;
        // Create wgpu instance with WebGPU backend
        let instance = wgpu::Instance::new(&wgpu::InstanceDescriptor {
            backends: wgpu::Backends::BROWSER_WEBGPU,
            ..Default::default()
        });

        // Create surface from canvas (wasm-only API)
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

        let adapter_info = adapter.get_info();
        log::info!("GPU adapter: {}", adapter_info.name);
        let adapter_limits = adapter.limits();
        let memory_policy = select_render_memory_policy(
            RenderMemorySelectionInput::from_wgpu(&adapter_info, &adapter_limits, true),
            None,
        );
        config.memory = memory_policy;
        log::info!(
            "render memory profile: profile={} budget={:?}",
            memory_policy.profile,
            memory_policy.budget_bytes
        );

        // Request device
        let (device, queue) = adapter
            .request_device(&wgpu::DeviceDescriptor {
                label: Some("patinae-web"),
                required_features: wgpu::Features::empty(),
                required_limits: required_limits_for_memory_policy(
                    &adapter_limits,
                    memory_policy,
                ),
                memory_hints: memory_policy.wgpu_memory_hints(),
                experimental_features: wgpu::ExperimentalFeatures::default(),
                trace: wgpu::Trace::Off,
            })
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

        // patinae-render takes Arc<Device>/Arc<Queue>; wrap the wgpu handles
        // (which are themselves Arc-backed internally — clone is cheap).
        let device_arc = Arc::new(device);
        let queue_arc = Arc::new(queue);
        let state = RenderState::with_config(
            device_arc,
            queue_arc,
            format,
            (width, height),
            config,
        );

        Ok(Self {
            state,
            surface,
            surface_config,
        })
        }
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
            .configure(&self.state.ctx.device, &self.surface_config);
        self.state.resize((width, height));
    }
}
