//! Main viewer and render loop
//!
//! The [`Viewer`] struct is the main application entry point, coordinating
//! window management, input handling, and rendering.

use std::path::Path;
use std::sync::Arc;
use std::time::Instant;

use pymol_color::{ChainColors, ElementColors, NamedColors};
use pymol_render::{ColorResolver, GlobalUniforms, RenderContext};
use pymol_select::{EvalContext, SelectionResult};
use pymol_settings::GlobalSettings;
use winit::application::ApplicationHandler;
use winit::event::{KeyEvent, WindowEvent};
use winit::event_loop::{ActiveEventLoop, ControlFlow, EventLoop};
use winit::keyboard::PhysicalKey;
use winit::window::WindowId;

use pymol_mol::RepMask;

use crate::camera::Camera;
use crate::error::{ViewerError, WindowError};
use crate::keybindings::{KeyBinding, KeyBindings};
use crate::object::ObjectRegistry;

// ============================================================================
// PNG Capture - Single source of truth for screenshot functionality
// ============================================================================

/// Capture the current scene to a PNG file
///
/// This is the single source of truth for PNG capture logic, used by both
/// the `Viewer` and GUI's `ViewerAdapter`.
///
/// # Arguments
///
/// * `path` - Output file path
/// * `width` - Optional width in pixels
/// * `height` - Optional height in pixels  
/// * `context` - Render context with GPU device/queue
/// * `camera` - Camera for view/projection matrices (aspect ratio temporarily modified)
/// * `registry` - Object registry containing molecules to render
/// * `settings` - Global settings for lighting, fog, etc.
/// * `named_colors` - Named color definitions
/// * `element_colors` - Element color definitions
/// * `chain_colors` - Chain color definitions
/// * `clear_color` - Background color RGB
/// * `default_size` - Fallback size when width/height not specified
///
/// # Example
///
/// ```ignore
/// capture_png_to_file(
///     Path::new("screenshot.png"),
///     Some(1920), Some(1080),
///     &context, &mut camera, &mut registry,
///     &settings, &named_colors, &element_colors, &chain_colors,
///     [0.0, 0.0, 0.0], (1024, 768),
/// )?;
/// ```
pub fn capture_png_to_file(
    path: &Path,
    width: Option<u32>,
    height: Option<u32>,
    context: &RenderContext,
    camera: &mut Camera,
    registry: &mut ObjectRegistry,
    settings: &GlobalSettings,
    named_colors: &NamedColors,
    element_colors: &ElementColors,
    chain_colors: &ChainColors,
    clear_color: [f32; 3],
    default_size: (u32, u32),
) -> Result<(), ViewerError> {
    // Determine output size
    let (output_width, output_height) = match (width, height) {
        (Some(w), Some(h)) => (w, h),
        (Some(w), None) => {
            // Maintain aspect ratio based on width
            let aspect = camera.aspect();
            (w, (w as f32 / aspect) as u32)
        }
        (None, Some(h)) => {
            // Maintain aspect ratio based on height
            let aspect = camera.aspect();
            ((h as f32 * aspect) as u32, h)
        }
        (None, None) => default_size,
    };

    let device = context.device();
    let queue = context.queue();

    // Create offscreen color texture using the same format as render pipelines
    let texture_format = context.surface_format();
    let color_texture = device.create_texture(&wgpu::TextureDescriptor {
        label: Some("Screenshot Color Texture"),
        size: wgpu::Extent3d {
            width: output_width,
            height: output_height,
            depth_or_array_layers: 1,
        },
        mip_level_count: 1,
        sample_count: 1,
        dimension: wgpu::TextureDimension::D2,
        format: texture_format,
        usage: wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::COPY_SRC,
        view_formats: &[],
    });
    let color_view = color_texture.create_view(&wgpu::TextureViewDescriptor::default());

    // Create offscreen depth texture
    let depth_texture = device.create_texture(&wgpu::TextureDescriptor {
        label: Some("Screenshot Depth Texture"),
        size: wgpu::Extent3d {
            width: output_width,
            height: output_height,
            depth_or_array_layers: 1,
        },
        mip_level_count: 1,
        sample_count: 1,
        dimension: wgpu::TextureDimension::D2,
        format: wgpu::TextureFormat::Depth32Float,
        usage: wgpu::TextureUsages::RENDER_ATTACHMENT,
        view_formats: &[],
    });
    let depth_view = depth_texture.create_view(&wgpu::TextureViewDescriptor::default());

    // Update uniforms for the capture resolution
    let mut uniforms = GlobalUniforms::new();

    // Temporarily update camera aspect ratio for correct projection
    let original_aspect = camera.aspect();
    camera.set_aspect(output_width as f32 / output_height as f32);

    uniforms.set_camera(camera.view_matrix(), camera.projection_matrix());
    uniforms.set_background(clear_color);
    uniforms.set_viewport(output_width as f32, output_height as f32);

    let camera_pos = camera.world_position();
    uniforms.camera_pos = [camera_pos.x, camera_pos.y, camera_pos.z, 1.0];

    // Lighting settings
    let ambient = settings.get_float(pymol_settings::id::ambient);
    let direct = settings.get_float(pymol_settings::id::direct);
    let reflect = settings.get_float(pymol_settings::id::reflect);
    let specular = settings.get_float(pymol_settings::id::specular);
    let shininess = settings.get_float(pymol_settings::id::shininess);
    let spec_direct = settings.get_float(pymol_settings::id::spec_direct);
    let spec_direct_power = settings.get_float(pymol_settings::id::spec_direct_power);
    uniforms.set_lighting(ambient, direct, reflect, specular, shininess, spec_direct, spec_direct_power);

    // Clip planes
    let current_view = camera.current_view();
    uniforms.set_clip_planes(current_view.clip_front, current_view.clip_back);

    // Fog parameters
    let depth_cue_enabled = settings.get_bool(pymol_settings::id::depth_cue);
    let fog_density = settings.get_float(pymol_settings::id::fog);
    if depth_cue_enabled && fog_density > 0.0 {
        let fog_start_ratio = settings.get_float(pymol_settings::id::fog_start);
        let fog_start_actual = (current_view.clip_back - current_view.clip_front) * fog_start_ratio + current_view.clip_front;
        let fog_end_actual = if (fog_density - 1.0).abs() < 0.001 {
            current_view.clip_back
        } else {
            fog_start_actual + (current_view.clip_back - fog_start_actual) / fog_density
        };
        uniforms.set_fog(fog_start_actual, fog_end_actual, fog_density, clear_color);
        uniforms.set_depth_cue(1.0);
    }

    context.update_uniforms(&uniforms);

    // Restore original aspect ratio
    camera.set_aspect(original_aspect);

    // Prepare molecules for rendering
    let names: Vec<_> = registry.names().map(|s| s.to_string()).collect();
    for name in &names {
        let color_resolver = ColorResolver::new(named_colors, element_colors, chain_colors);
        if let Some(mol_obj) = registry.get_molecule_mut(name) {
            mol_obj.prepare_render(context, &color_resolver, settings);
        }
    }

    // Create command encoder
    let mut encoder = device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
        label: Some("Screenshot Encoder"),
    });

    // Render pass
    {
        let mut render_pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
            label: Some("Screenshot Render Pass"),
            color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                view: &color_view,
                resolve_target: None,
                ops: wgpu::Operations {
                    load: wgpu::LoadOp::Clear(wgpu::Color {
                        r: clear_color[0] as f64,
                        g: clear_color[1] as f64,
                        b: clear_color[2] as f64,
                        a: 1.0,
                    }),
                    store: wgpu::StoreOp::Store,
                },
                depth_slice: None,
            })],
            depth_stencil_attachment: Some(wgpu::RenderPassDepthStencilAttachment {
                view: &depth_view,
                depth_ops: Some(wgpu::Operations {
                    load: wgpu::LoadOp::Clear(1.0),
                    store: wgpu::StoreOp::Store,
                }),
                stencil_ops: None,
            }),
            timestamp_writes: None,
            occlusion_query_set: None,
        });

        // Render all enabled objects
        for name in &names {
            if let Some(mol_obj) = registry.get_molecule(name) {
                if mol_obj.is_enabled() {
                    mol_obj.render(&mut render_pass, context);
                }
            }
        }
    }

    // Calculate buffer size with proper alignment
    // wgpu requires rows to be aligned to 256 bytes
    let bytes_per_pixel = 4u32; // RGBA8
    let unpadded_bytes_per_row = output_width * bytes_per_pixel;
    let align = wgpu::COPY_BYTES_PER_ROW_ALIGNMENT;
    let padded_bytes_per_row = (unpadded_bytes_per_row + align - 1) / align * align;
    let buffer_size = (padded_bytes_per_row * output_height) as u64;

    // Create output buffer
    let output_buffer = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("Screenshot Buffer"),
        size: buffer_size,
        usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
        mapped_at_creation: false,
    });

    // Copy texture to buffer
    encoder.copy_texture_to_buffer(
        wgpu::TexelCopyTextureInfo {
            texture: &color_texture,
            mip_level: 0,
            origin: wgpu::Origin3d::ZERO,
            aspect: wgpu::TextureAspect::All,
        },
        wgpu::TexelCopyBufferInfo {
            buffer: &output_buffer,
            layout: wgpu::TexelCopyBufferLayout {
                offset: 0,
                bytes_per_row: Some(padded_bytes_per_row),
                rows_per_image: Some(output_height),
            },
        },
        wgpu::Extent3d {
            width: output_width,
            height: output_height,
            depth_or_array_layers: 1,
        },
    );

    // Submit and wait
    queue.submit(std::iter::once(encoder.finish()));

    // Map buffer and read data
    let buffer_slice = output_buffer.slice(..);
    let (tx, rx) = std::sync::mpsc::channel();
    buffer_slice.map_async(wgpu::MapMode::Read, move |result| {
        tx.send(result).unwrap();
    });
    device.poll(wgpu::PollType::Wait { submission_index: None, timeout: None }).ok();
    rx.recv()
        .map_err(|e| ViewerError::GpuInitFailed(format!("Failed to receive map result: {}", e)))?
        .map_err(|e| ViewerError::GpuInitFailed(format!("Failed to map buffer: {:?}", e)))?;

    // Read pixel data (removing padding)
    let data = buffer_slice.get_mapped_range();
    let mut pixels = Vec::with_capacity((output_width * output_height * bytes_per_pixel) as usize);
    for row in 0..output_height {
        let start = (row * padded_bytes_per_row) as usize;
        let end = start + (output_width * bytes_per_pixel) as usize;
        pixels.extend_from_slice(&data[start..end]);
    }
    drop(data);
    output_buffer.unmap();

    // Convert BGRA to RGBA if needed (common on many platforms)
    let is_bgra = matches!(
        texture_format,
        wgpu::TextureFormat::Bgra8Unorm | wgpu::TextureFormat::Bgra8UnormSrgb
    );
    if is_bgra {
        // Swap B and R channels: BGRA -> RGBA
        for chunk in pixels.chunks_exact_mut(4) {
            chunk.swap(0, 2);
        }
    }

    // Save as PNG using image crate
    let img: image::RgbaImage = image::ImageBuffer::from_raw(output_width, output_height, pixels)
        .ok_or_else(|| ViewerError::GpuInitFailed("Failed to create image buffer".to_string()))?;

    img.save(path)
        .map_err(|e| ViewerError::IoError(format!("Failed to save PNG: {}", e)))?;

    Ok(())
}

// ============================================================================
// Selection Entry - stores selection expression and visibility state
// ============================================================================

/// Entry for a named selection
#[derive(Debug, Clone)]
pub struct SelectionEntry {
    /// The selection expression
    pub expression: String,
    /// Whether the selection indicators are visible
    pub visible: bool,
}

impl SelectionEntry {
    /// Create a new selection entry with visibility enabled by default
    pub fn new(expression: String) -> Self {
        Self {
            expression,
            visible: true,
        }
    }
}

/// Trait for types that can serve as a viewer backend for command execution
///
/// This trait abstracts the viewer interface, allowing commands to work with
/// different viewer implementations (e.g., `Viewer`, GUI adapters).
pub trait ViewerLike {
    /// Get a reference to the object registry
    fn objects(&self) -> &ObjectRegistry;

    /// Get a mutable reference to the object registry
    fn objects_mut(&mut self) -> &mut ObjectRegistry;

    /// Get a reference to the camera
    fn camera(&self) -> &Camera;

    /// Get a mutable reference to the camera
    fn camera_mut(&mut self) -> &mut Camera;

    /// Zoom to fit all objects while preserving rotation
    fn zoom_all(&mut self);

    /// Zoom to fit a specific object while preserving rotation
    fn zoom_on(&mut self, name: &str);

    /// Center on all objects without changing zoom or rotation
    fn center_all(&mut self);

    /// Center on a specific object without changing zoom or rotation
    fn center_on(&mut self, name: &str);

    /// Reset the camera to default view
    fn reset_view(&mut self);

    /// Request a redraw
    fn request_redraw(&mut self);

    /// Get a color index by name
    fn color_index(&self, name: &str) -> Option<u32>;

    /// Set the background color
    fn set_background_color(&mut self, r: f32, g: f32, b: f32);

    // =========================================================================
    // Named Selections
    // =========================================================================

    /// Get a named selection expression by name
    fn get_selection(&self, name: &str) -> Option<&str>;

    /// Define (store) a named selection expression
    fn define_selection(&mut self, name: &str, selection: &str);

    /// Remove a named selection, returns true if it existed
    fn remove_selection(&mut self, name: &str) -> bool;

    /// Get all named selection names
    fn selection_names(&self) -> Vec<String>;

    /// Set the visibility of a named selection's indicators
    fn set_selection_visible(&mut self, name: &str, visible: bool);

    /// Check if a named selection's indicators are visible
    fn is_selection_visible(&self, name: &str) -> bool;

    /// Set the indicated selection (shows pink indicators in the 3D view)
    /// For backwards compatibility - creates/updates an "indicate" selection
    fn indicate_selection(&mut self, selection: &str);

    /// Clear the indicated selection (hides all selection indicators)
    fn clear_indication(&mut self);

    /// Get the currently indicated selection expression
    /// For backwards compatibility - returns first visible selection
    fn indicated_selection(&self) -> Option<&str>;

    /// Capture a PNG screenshot (optional - not all viewers support this)
    ///
    /// Returns an error by default. Override for viewers that support screenshots.
    fn capture_png(
        &mut self,
        _path: &Path,
        _width: Option<u32>,
        _height: Option<u32>,
    ) -> Result<(), String> {
        Err("Screenshot capture not supported by this viewer".to_string())
    }

    // =========================================================================
    // Async Operations
    // =========================================================================

    /// Request an async fetch operation from RCSB PDB (non-blocking)
    ///
    /// This method allows commands to request background fetch operations without
    /// blocking the main thread. The viewer implementation decides how to handle
    /// the async operation.
    ///
    /// # Arguments
    ///
    /// * `code` - PDB ID to fetch (e.g., "1ubq")
    /// * `name` - Object name to use when adding to registry
    /// * `format` - Format code: 0 = CIF (default), 1 = PDB
    ///
    /// # Returns
    ///
    /// * `true` - Request was accepted; fetch is running in background
    /// * `false` - Async fetch not supported; caller should use sync fallback
    ///
    /// The default implementation returns `false`, indicating async fetch is not
    /// supported. GUI viewers should override this to spawn background tasks.
    fn request_async_fetch(&mut self, _code: &str, _name: &str, _format: u8) -> bool {
        false
    }
}

/// Type alias for key action callbacks
///
/// Key actions are closures that receive a mutable reference to the Viewer
/// and can perform any operation on it.
pub type KeyAction = Arc<dyn Fn(&mut Viewer) + Send + Sync>;
use crate::input::{CameraDelta, InputState};
use crate::object::{MoleculeObject, Object};
use crate::scene::SceneManager;
use crate::window::Window;

/// Main molecular visualization viewer
///
/// Coordinates all aspects of the visualization application:
/// - GPU resource management (wgpu device, queue, context)
/// - Window and surface management
/// - Camera and input handling
/// - Object and scene management
/// - Rendering
pub struct Viewer {
    // =========================================================================
    // GPU Resources
    // =========================================================================
    /// wgpu instance
    instance: wgpu::Instance,
    /// Render context (owns device, queue, pipelines, shared buffers)
    render_context: Option<RenderContext>,

    // =========================================================================
    // Window State
    // =========================================================================
    /// The main window (created on resume)
    window: Option<Window>,
    /// Pending window creation request
    pending_window: Option<(String, u32, u32)>,

    // =========================================================================
    // Scene State
    // =========================================================================
    /// Camera for view control
    camera: Camera,
    /// Object registry
    registry: ObjectRegistry,
    /// Scene manager
    scenes: SceneManager,
    /// Named selections (name -> entry with expression and visibility)
    selections: ahash::AHashMap<String, SelectionEntry>,

    // =========================================================================
    // Settings and Colors
    // =========================================================================
    /// Global settings
    settings: GlobalSettings,
    /// Named colors table
    named_colors: NamedColors,
    /// Element colors table
    element_colors: ElementColors,
    /// Chain colors (unit struct)
    chain_colors: ChainColors,

    // =========================================================================
    // Input State
    // =========================================================================
    /// Input handler
    input: InputState,

    // =========================================================================
    // Frame Timing
    // =========================================================================
    /// Last frame timestamp
    last_frame: Instant,
    /// Whether a redraw is needed
    needs_redraw: bool,

    // =========================================================================
    // Background Color
    // =========================================================================
    /// Clear color
    clear_color: [f32; 3],

    // =========================================================================
    // Key Bindings
    // =========================================================================
    /// Keyboard shortcuts mapped to action callbacks
    key_bindings: KeyBindings<KeyAction>,
}

impl Default for Viewer {
    fn default() -> Self {
        Self::new()
    }
}

impl Viewer {
    /// Create a new viewer
    ///
    /// The viewer is created without a window. Use [`create_window`] to
    /// create a window after the event loop is running.
    pub fn new() -> Self {
        let instance = wgpu::Instance::new(&wgpu::InstanceDescriptor {
            backends: wgpu::Backends::all(),
            ..Default::default()
        });

        Self {
            instance,
            render_context: None,
            window: None,
            pending_window: Some(("PyMOL-RS".to_string(), 1024, 768)),
            camera: Camera::new(),
            registry: ObjectRegistry::new(),
            scenes: SceneManager::new(),
            selections: ahash::AHashMap::new(),
            settings: GlobalSettings::new(),
            named_colors: NamedColors::default(),
            element_colors: ElementColors::default(),
            chain_colors: ChainColors,
            input: InputState::new(),
            last_frame: Instant::now(),
            needs_redraw: true,
            clear_color: [0.0, 0.0, 0.0],
            key_bindings: KeyBindings::new(),
        }
    }

    /// Initialize GPU and create window
    async fn init_gpu_and_window(
        &mut self,
        event_loop: &ActiveEventLoop,
        title: &str,
        width: u32,
        height: u32,
    ) -> Result<(), ViewerError> {
        // Request adapter
        let adapter = self
            .instance
            .request_adapter(&wgpu::RequestAdapterOptions {
                power_preference: wgpu::PowerPreference::HighPerformance,
                compatible_surface: None,
                force_fallback_adapter: false,
            })
            .await
            .map_err(|e| ViewerError::GpuInitFailed(format!("No suitable adapter found: {}", e)))?;

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
            .map_err(|e| ViewerError::GpuInitFailed(e.to_string()))?;

        // Create window
        let window = Window::new(
            event_loop,
            title,
            (width, height),
            &self.instance,
            &adapter,
            &device,
        ).map_err(ViewerError::Window)?;

        let surface_format = window.surface_format();

        // Create render context - it takes ownership of device and queue
        let render_context = RenderContext::new(device, queue, surface_format);

        // Set camera aspect ratio
        self.camera.set_aspect(window.aspect_ratio());

        self.window = Some(window);
        self.render_context = Some(render_context);
        self.needs_redraw = true;

        Ok(())
    }

    // =========================================================================
    // Public API
    // =========================================================================

    /// Get a reference to the camera
    pub fn camera(&self) -> &Camera {
        &self.camera
    }

    /// Get a mutable reference to the camera
    pub fn camera_mut(&mut self) -> &mut Camera {
        &mut self.camera
    }

    /// Get a reference to the object registry
    pub fn objects(&self) -> &ObjectRegistry {
        &self.registry
    }

    /// Get a mutable reference to the object registry
    pub fn objects_mut(&mut self) -> &mut ObjectRegistry {
        &mut self.registry
    }

    /// Get a reference to the scene manager
    pub fn scenes(&self) -> &SceneManager {
        &self.scenes
    }

    /// Get a mutable reference to the scene manager
    pub fn scenes_mut(&mut self) -> &mut SceneManager {
        &mut self.scenes
    }

    /// Get a reference to the global settings
    pub fn settings(&self) -> &GlobalSettings {
        &self.settings
    }

    /// Get a mutable reference to the global settings
    pub fn settings_mut(&mut self) -> &mut GlobalSettings {
        &mut self.settings
    }

    /// Add a molecule object to the scene
    pub fn add_molecule(&mut self, mol: pymol_mol::ObjectMolecule) -> &str {
        let obj = MoleculeObject::new(mol);
        self.registry.add(obj)
    }

    /// Zoom the camera to fit all objects while preserving rotation
    pub fn zoom_all(&mut self) {
        if let Some((min, max)) = self.registry.extent() {
            self.camera.zoom_to(min, max);
            self.needs_redraw = true;
        }
    }

    /// Zoom the camera to fit a specific object while preserving rotation
    pub fn zoom_on(&mut self, name: &str) {
        if let Some(obj) = self.registry.get(name) {
            if let Some((min, max)) = obj.extent() {
                self.camera.zoom_to(min, max);
                self.needs_redraw = true;
            }
        }
    }

    /// Center the camera on all objects without changing zoom or rotation
    pub fn center_all(&mut self) {
        if let Some((min, max)) = self.registry.extent() {
            self.camera.center_to(min, max);
            self.needs_redraw = true;
        }
    }

    /// Center the camera on a specific object without changing zoom or rotation
    pub fn center_on(&mut self, name: &str) {
        if let Some(obj) = self.registry.get(name) {
            if let Some((min, max)) = obj.extent() {
                self.camera.center_to(min, max);
                self.needs_redraw = true;
            }
        }
    }

    /// Reset the camera to default view (resets everything including rotation)
    pub fn reset_view(&mut self) {
        self.camera = Camera::new();
        if let Some((min, max)) = self.registry.extent() {
            self.camera.reset_view(min, max);
            self.needs_redraw = true;
        }
    }

    /// Set the background color
    pub fn set_background_color(&mut self, r: f32, g: f32, b: f32) {
        self.clear_color = [r, g, b];
        self.needs_redraw = true;
    }

    /// Request a redraw
    pub fn request_redraw(&mut self) {
        self.needs_redraw = true;
        if let Some(window) = &self.window {
            window.request_redraw();
        }
    }

    // =========================================================================
    // Color Management
    // =========================================================================

    /// Get a color index by name
    ///
    /// Returns the color index for a named color (e.g., "red", "green", "blue").
    /// This index can be assigned to atoms via `atom.color`.
    ///
    /// Returns `None` if the color name is not found.
    pub fn color_index(&self, name: &str) -> Option<u32> {
        self.named_colors.get_by_name(name).map(|(idx, _)| idx)
    }

    /// Get a reference to the named colors registry
    pub fn named_colors(&self) -> &NamedColors {
        &self.named_colors
    }

    // =========================================================================
    // Selection Indication
    // =========================================================================

    /// Set the indicated selection (shows pink indicators in the 3D view)
    ///
    /// This evaluates the selection and displays pink indicator dots at
    /// the positions of all selected atoms.
    ///
    /// # Arguments
    /// * `selection` - A selection expression (e.g., "chain A", "organic", "sele")
    pub fn indicate_selection(&mut self, selection: &str) {
        // Create or update "indicate" selection and make it visible
        self.selections.insert("indicate".to_string(), SelectionEntry::new(selection.to_string()));
        self.needs_redraw = true;
    }

    /// Clear the indicated selection
    ///
    /// Hides all selection indicators in the 3D view.
    pub fn clear_indication(&mut self) {
        // Hide all selection indicators
        for entry in self.selections.values_mut() {
            entry.visible = false;
        }
        // Clear indicators from all molecules
        let names: Vec<_> = self.registry.names().map(|s| s.to_string()).collect();
        for name in names {
            if let Some(mol) = self.registry.get_molecule_mut(&name) {
                mol.clear_selection_indicator();
            }
        }
        self.needs_redraw = true;
    }

    /// Get the currently indicated selection expression
    /// Returns the first visible selection expression for backwards compatibility
    pub fn indicated_selection(&self) -> Option<&str> {
        self.selections.values()
            .find(|e| e.visible)
            .map(|e| e.expression.as_str())
    }

    /// Evaluate all visible selections for all molecules
    ///
    /// Returns a vector of (object_name, SelectionResult) pairs where the
    /// SelectionResult is the union of all visible selections.
    fn evaluate_visible_selections(
        &self,
        names: &[String],
    ) -> Vec<(String, SelectionResult)> {
        // Collect visible selection expressions
        let visible_selections: Vec<(&str, &str)> = self.selections
            .iter()
            .filter(|(_, entry)| entry.visible)
            .map(|(name, entry)| (name.as_str(), entry.expression.as_str()))
            .collect();

        if visible_selections.is_empty() {
            return Vec::new();
        }

        let mut results: Vec<(String, SelectionResult)> = Vec::new();

        // Evaluate for each molecule
        for mol_name in names {
            if let Some(mol_obj) = self.registry.get_molecule(mol_name) {
                let mol = mol_obj.molecule();
                
                // Build context for this molecule with named selections
                let mut ctx = EvalContext::single(mol);
                
                // Add all named selections to context (for reference resolution)
                for (sel_name, entry) in &self.selections {
                    if let Ok(sel_ast) = pymol_select::parse(&entry.expression) {
                        if let Ok(sel_result) = pymol_select::evaluate(&sel_ast, &ctx) {
                            ctx.add_selection(sel_name.clone(), sel_result);
                        }
                    }
                }

                // Evaluate and combine all visible selections
                let mut combined_result: Option<SelectionResult> = None;

                for (sel_name, sel_expr) in &visible_selections {
                    let expr = match pymol_select::parse(sel_expr) {
                        Ok(e) => e,
                        Err(e) => {
                            log::debug!("Failed to parse selection '{}': {:?}", sel_name, e);
                            continue;
                        }
                    };

                    match pymol_select::evaluate(&expr, &ctx) {
                        Ok(result) => {
                            combined_result = match combined_result {
                                Some(existing) => Some(existing.union(&result)),
                                None => Some(result),
                            };
                        }
                        Err(e) => {
                            log::debug!("Failed to evaluate selection '{}' for '{}': {:?}", sel_name, mol_name, e);
                        }
                    }
                }

                // Add combined result if any atoms matched
                if let Some(result) = combined_result {
                    if result.any() {
                        results.push((mol_name.clone(), result));
                    }
                }
            }
        }

        results
    }

    // =========================================================================
    // Representation Control
    // =========================================================================

    /// Toggle a representation for all molecules
    ///
    /// If the representation is currently visible on any molecule, it will be hidden.
    /// Otherwise, it will be shown on all molecules.
    pub fn toggle_representation(&mut self, rep: u32) {
        let names: Vec<_> = self.registry.names().map(|s| s.to_string()).collect();
        for name in names {
            if let Some(mol) = self.registry.get_molecule_mut(&name) {
                mol.toggle(rep);
            }
        }
        self.needs_redraw = true;
    }

    /// Show a representation for all molecules
    pub fn show_representation(&mut self, rep: u32) {
        let names: Vec<_> = self.registry.names().map(|s| s.to_string()).collect();
        for name in names {
            if let Some(mol) = self.registry.get_molecule_mut(&name) {
                mol.show(rep);
            }
        }
        self.needs_redraw = true;
    }

    /// Hide a representation for all molecules
    pub fn hide_representation(&mut self, rep: u32) {
        let names: Vec<_> = self.registry.names().map(|s| s.to_string()).collect();
        for name in names {
            if let Some(mol) = self.registry.get_molecule_mut(&name) {
                mol.hide(rep);
            }
        }
        self.needs_redraw = true;
    }

    /// Hide all representations for all molecules
    pub fn hide_all_representations(&mut self) {
        let names: Vec<_> = self.registry.names().map(|s| s.to_string()).collect();
        for name in names {
            if let Some(mol) = self.registry.get_molecule_mut(&name) {
                mol.hide_all();
            }
        }
        self.needs_redraw = true;
    }

    /// Show default representations (lines + sticks) for all molecules
    pub fn show_default_representations(&mut self) {
        let names: Vec<_> = self.registry.names().map(|s| s.to_string()).collect();
        for name in names {
            if let Some(mol) = self.registry.get_molecule_mut(&name) {
                mol.hide_all();
                mol.show(RepMask::LINES);
                mol.show(RepMask::STICKS);
            }
        }
        self.needs_redraw = true;
    }

    // =========================================================================
    // Surface Quality Control
    // =========================================================================

    /// Set surface quality for all molecules (-4 to 4)
    ///
    /// Lower values are faster but coarser, higher values are smoother but slower.
    /// - -4: Very coarse (fastest)
    /// - 0: Default quality
    /// - 4: Very fine (slowest)
    pub fn set_surface_quality(&mut self, quality: i32) {
        let names: Vec<_> = self.registry.names().map(|s| s.to_string()).collect();
        for name in names {
            if let Some(mol) = self.registry.get_molecule_mut(&name) {
                mol.set_surface_quality(quality);
            }
        }
        self.needs_redraw = true;
    }

    /// Get current surface quality (from first molecule, or default 0)
    pub fn surface_quality(&self) -> i32 {
        self.registry
            .names()
            .next()
            .and_then(|n| self.registry.get_molecule(n))
            .map(|m| m.surface_quality())
            .unwrap_or(0)
    }

    /// Increase surface quality (finer mesh, slower)
    pub fn increase_surface_quality(&mut self) {
        let current = self.surface_quality();
        let new_quality = (current + 1).min(4);
        if new_quality != current {
            self.set_surface_quality(new_quality);
            log::info!("Surface quality: {}", new_quality);
        }
    }

    /// Decrease surface quality (coarser mesh, faster)
    pub fn decrease_surface_quality(&mut self) {
        let current = self.surface_quality();
        let new_quality = (current - 1).max(-4);
        if new_quality != current {
            self.set_surface_quality(new_quality);
            log::info!("Surface quality: {}", new_quality);
        }
    }

    // =========================================================================
    // Key Bindings
    // =========================================================================

    /// Bind a key or key combination to an action callback
    ///
    /// When the specified key (with optional modifiers) is pressed, the provided
    /// closure will be called with a mutable reference to the Viewer. If the key
    /// was already bound, the previous binding is replaced.
    ///
    /// The key parameter accepts either a `KeyCode` for simple bindings or a
    /// `KeyBinding` for combinations with modifiers.
    ///
    /// # Example
    ///
    /// ```ignore
    /// use pymol_scene::{Viewer, KeyCode, KeyBinding};
    /// use pymol_mol::RepMask;
    ///
    /// let mut viewer = Viewer::new();
    ///
    /// // Bind simple keys (no modifiers)
    /// viewer.bind_key(KeyCode::KeyR, |v| v.reset_view());
    /// viewer.bind_key(KeyCode::Digit1, |v| v.toggle_representation(RepMask::LINES));
    ///
    /// // Bind key combinations with modifiers
    /// viewer.bind_key(KeyBinding::new(KeyCode::KeyR).ctrl(), |v| {
    ///     log::info!("Ctrl+R pressed!");
    ///     v.reset_view();
    /// });
    /// viewer.bind_key(KeyBinding::new(KeyCode::KeyS).ctrl().shift(), |v| {
    ///     log::info!("Ctrl+Shift+S pressed!");
    /// });
    /// ```
    pub fn bind_key<K, F>(&mut self, key: K, action: F)
    where
        K: Into<KeyBinding>,
        F: Fn(&mut Viewer) + Send + Sync + 'static,
    {
        self.key_bindings.bind(key, Arc::new(action));
    }

    /// Remove a key binding
    ///
    /// Returns true if a binding was removed.
    pub fn unbind_key<K: Into<KeyBinding>>(&mut self, key: K) -> bool {
        self.key_bindings.unbind(key)
    }

    /// Check if a key is bound to an action
    pub fn is_key_bound<K: Into<KeyBinding>>(&self, key: K) -> bool {
        self.key_bindings.is_bound(key)
    }

    // =========================================================================
    // Rendering
    // =========================================================================

    /// Render a single frame
    fn render(&mut self) -> Result<(), ViewerError> {
        let window = self.window.as_ref().ok_or(ViewerError::Window(
            WindowError::NotAvailable,
        ))?;

        if window.is_minimized() {
            return Ok(());
        }

        let context = self.render_context.as_ref().ok_or(ViewerError::Window(
            WindowError::NotAvailable,
        ))?;

        // Get surface texture
        let output = match window.get_current_texture() {
            Ok(output) => output,
            Err(wgpu::SurfaceError::Lost | wgpu::SurfaceError::Outdated) => {
                // Reconfigure surface - need mutable access
                if let Some(window) = &mut self.window {
                    if let Some(context) = &self.render_context {
                        let size = window.size();
                        window.resize(size, context.device());
                    }
                }
                return Ok(());
            }
            Err(wgpu::SurfaceError::OutOfMemory) => {
                return Err(ViewerError::GpuInitFailed("Out of memory".to_string()));
            }
            Err(wgpu::SurfaceError::Timeout) => {
                return Ok(());
            }
            Err(wgpu::SurfaceError::Other) => {
                return Err(ViewerError::GpuInitFailed("Surface error".to_string()));
            }
        };

        let view = output
            .texture
            .create_view(&wgpu::TextureViewDescriptor::default());

        // Update uniforms
        let mut uniforms = GlobalUniforms::new();
        uniforms.set_camera(self.camera.view_matrix(), self.camera.projection_matrix());
        uniforms.set_background(self.clear_color);
        uniforms.set_viewport(window.width() as f32, window.height() as f32);

        let camera_pos = self.camera.world_position();
        uniforms.camera_pos = [camera_pos.x, camera_pos.y, camera_pos.z, 1.0];

        // Get lighting from settings (PyMOL dual-light model)
        let ambient = self.settings.get_float(pymol_settings::id::ambient);
        let direct = self.settings.get_float(pymol_settings::id::direct);
        let reflect = self.settings.get_float(pymol_settings::id::reflect);
        let specular = self.settings.get_float(pymol_settings::id::specular);
        let shininess = self.settings.get_float(pymol_settings::id::shininess);
        let spec_direct = self.settings.get_float(pymol_settings::id::spec_direct);
        let spec_direct_power = self.settings.get_float(pymol_settings::id::spec_direct_power);
        uniforms.set_lighting(ambient, direct, reflect, specular, shininess, spec_direct, spec_direct_power);

        // Set clip planes from camera
        let current_view = self.camera.current_view();
        let clip_front = current_view.clip_front;
        let clip_back = current_view.clip_back;
        uniforms.set_clip_planes(clip_front, clip_back);

        // Compute fog parameters (PyMOL-compatible)
        // depth_cue controls whether fog is enabled
        let depth_cue_enabled = self.settings.get_bool(pymol_settings::id::depth_cue);
        let fog_density = self.settings.get_float(pymol_settings::id::fog);

        if depth_cue_enabled && fog_density > 0.0 {
            // fog_start setting is a ratio (0.0-1.0) between clip planes
            let fog_start_ratio = self.settings.get_float(pymol_settings::id::fog_start);

            // Compute actual fog distances based on clip planes
            // FogStart = (back - front) * fog_start_ratio + front
            let fog_start_actual = (clip_back - clip_front) * fog_start_ratio + clip_front;

            // FogEnd depends on fog density
            let fog_end_actual = if (fog_density - 1.0).abs() < 0.001 {
                // If density is ~1.0, fog ends at back clip plane
                clip_back
            } else {
                // Otherwise scale the fog range by density
                fog_start_actual + (clip_back - fog_start_actual) / fog_density
            };

            // Set fog with background color (so objects fade to background)
            uniforms.set_fog(fog_start_actual, fog_end_actual, fog_density, self.clear_color);

            // Enable depth cue (darkening based on depth)
            uniforms.set_depth_cue(1.0);
        }
        // else: fog_density and depth_cue remain at 0.0 (disabled) from defaults

        context.update_uniforms(&uniforms);

        // Prepare molecules for rendering
        // Create color resolver with borrowed data - need to avoid borrowing self during registry iteration
        let names: Vec<_> = self.registry.names().map(|s| s.to_string()).collect();
        
        // Evaluate all visible selections
        let selection_results = self.evaluate_visible_selections(&names);
        
        // Get selection indicator size from settings (selection_width, ID 80)
        // Use a minimum size to ensure visibility
        let selection_width = self.settings.get_float(pymol_settings::id::selection_width).max(6.0);
        
        for name in &names {
            // Create color resolver fresh each iteration to avoid borrow conflicts
            let color_resolver = ColorResolver::new(&self.named_colors, &self.element_colors, &self.chain_colors);
            if let Some(mol_obj) = self.registry.get_molecule_mut(name) {
                mol_obj.prepare_render(context, &color_resolver, &self.settings);
                
                // Update selection indicator if we have results for this molecule
                if let Some((_, selection_result)) = selection_results.iter().find(|(n, _)| n == name) {
                    mol_obj.set_selection_indicator_with_size(selection_result, context, Some(selection_width));
                } else {
                    // Clear indicator - no visible selection matches this molecule
                    mol_obj.clear_selection_indicator();
                }
            }
        }

        // Create command encoder
        let mut encoder = context.device().create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("Render Encoder"),
        });

        // Begin render pass
        {
            let mut render_pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("Main Render Pass"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &view,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Clear(wgpu::Color {
                            r: self.clear_color[0] as f64,
                            g: self.clear_color[1] as f64,
                            b: self.clear_color[2] as f64,
                            a: 1.0,
                        }),
                        store: wgpu::StoreOp::Store,
                    },
                    depth_slice: None,
                })],
                depth_stencil_attachment: Some(wgpu::RenderPassDepthStencilAttachment {
                    view: window.depth_view(),
                    depth_ops: Some(wgpu::Operations {
                        load: wgpu::LoadOp::Clear(1.0),
                        store: wgpu::StoreOp::Store,
                    }),
                    stencil_ops: None,
                }),
                timestamp_writes: None,
                occlusion_query_set: None,
            });

            // Render all enabled objects
            for name in self.registry.names().collect::<Vec<_>>() {
                if let Some(mol_obj) = self.registry.get_molecule(&name) {
                    if mol_obj.is_enabled() {
                        mol_obj.render(&mut render_pass, context);
                    }
                }
            }
        }

        // Submit commands
        context.queue().submit(std::iter::once(encoder.finish()));
        output.present();

        Ok(())
    }

    // =========================================================================
    // Screenshot Capture
    // =========================================================================

    /// Capture the current view and save as a PNG file
    ///
    /// # Arguments
    ///
    /// * `path` - Output file path
    /// * `width` - Optional width in pixels (uses current window width if None)
    /// * `height` - Optional height in pixels (uses current window height if None)
    ///
    /// # Example
    ///
    /// ```ignore
    /// viewer.capture_png("screenshot.png", None, None)?;
    /// viewer.capture_png("hires.png", Some(1920), Some(1080))?;
    /// ```
    pub fn capture_png(
        &mut self,
        path: impl AsRef<std::path::Path>,
        width: Option<u32>,
        height: Option<u32>,
    ) -> Result<(), ViewerError> {
        let context = self.render_context.as_ref().ok_or(ViewerError::Window(
            WindowError::NotAvailable,
        ))?;

        // Get default size from window or use fallback
        let default_size = self.window.as_ref()
            .map(|w| w.size())
            .unwrap_or((1024, 768));

        // Delegate to the shared capture function
        capture_png_to_file(
            path.as_ref(),
            width,
            height,
            context,
            &mut self.camera,
            &mut self.registry,
            &self.settings,
            &self.named_colors,
            &self.element_colors,
            &self.chain_colors,
            self.clear_color,
            default_size,
        )
    }

    /// Update the viewer state
    fn update(&mut self, dt: f32) {
        // Process input deltas
        let deltas = self.input.take_camera_deltas();
        for delta in deltas {
            match delta {
                CameraDelta::Rotate { x, y } => {
                    self.camera.rotate_x(x);
                    self.camera.rotate_y(y);
                    self.needs_redraw = true;
                }
                CameraDelta::Translate(v) => {
                    self.camera.translate(v);
                    self.needs_redraw = true;
                }
                CameraDelta::Zoom(factor) => {
                    self.camera.zoom(factor);
                    self.needs_redraw = true;
                }
                CameraDelta::Clip { front, back } => {
                    let view = self.camera.view_mut();
                    view.clip_front = (view.clip_front + front).max(0.01);
                    view.clip_back = (view.clip_back + back).max(view.clip_front + 0.01);
                    self.needs_redraw = true;
                }
            }
        }

        // Update camera animation
        if self.camera.update(dt) {
            self.needs_redraw = true;
        }
    }

    /// Handle a keyboard event
    ///
    /// Looks up the pressed key (with current modifiers) in the key bindings map
    /// and executes the associated action callback if found.
    fn handle_key(&mut self, event: KeyEvent) {
        if event.state != winit::event::ElementState::Pressed {
            return;
        }

        if let PhysicalKey::Code(key) = event.physical_key {
            // Create a KeyBinding with current modifier state
            let binding = KeyBinding {
                key,
                ctrl: self.input.ctrl_held(),
                shift: self.input.shift_held(),
                alt: self.input.alt_held(),
            };

            // Clone the Arc to avoid borrow conflicts when executing
            if let Some(action) = self.key_bindings.get_cloned(&binding) {
                action(self);
            }
        }
    }
}

impl ApplicationHandler for Viewer {
    fn resumed(&mut self, event_loop: &ActiveEventLoop) {
        // Create window if pending (combines GPU init and window creation)
        if let Some((title, width, height)) = self.pending_window.take() {
            if let Err(e) = pollster::block_on(self.init_gpu_and_window(event_loop, &title, width, height)) {
                log::error!("Failed to initialize: {}", e);
                event_loop.exit();
                return;
            }
        }

        // Request initial redraw
        if let Some(window) = &self.window {
            window.request_redraw();
        }
    }

    fn window_event(&mut self, event_loop: &ActiveEventLoop, window_id: WindowId, event: WindowEvent) {
        // Ignore events for other windows
        if self.window.as_ref().map(|w| w.id()) != Some(window_id) {
            return;
        }

        match event {
            WindowEvent::CloseRequested => {
                event_loop.exit();
            }

            WindowEvent::Resized(size) => {
                if let (Some(context), Some(window)) = (&self.render_context, &mut self.window) {
                    window.resize((size.width, size.height), context.device());
                    self.camera.set_aspect(window.aspect_ratio());
                    self.needs_redraw = true;
                }
            }

            WindowEvent::RedrawRequested => {
                // Calculate delta time
                let now = Instant::now();
                let dt = (now - self.last_frame).as_secs_f32();
                self.last_frame = now;

                // Update state
                self.update(dt);

                // Render if needed
                if self.needs_redraw {
                    if let Err(e) = self.render() {
                        log::error!("Render error: {}", e);
                    }
                    self.needs_redraw = false;
                }

                // Request continuous redraw if animating or input is active
                if self.camera.is_animating() || self.input.any_button_pressed() {
                    if let Some(window) = &self.window {
                        window.request_redraw();
                    }
                }
            }

            WindowEvent::MouseInput { state, button, .. } => {
                self.input.handle_mouse_button(state, button);
                if let Some(window) = &self.window {
                    window.request_redraw();
                }
            }

            WindowEvent::CursorMoved { position, .. } => {
                self.input.handle_mouse_motion((position.x, position.y));
                if self.input.any_button_pressed() {
                    self.needs_redraw = true;
                    if let Some(window) = &self.window {
                        window.request_redraw();
                    }
                }
            }

            WindowEvent::MouseWheel { delta, .. } => {
                self.input.handle_scroll(delta);
                self.needs_redraw = true;
                if let Some(window) = &self.window {
                    window.request_redraw();
                }
            }

            WindowEvent::ModifiersChanged(modifiers) => {
                self.input.handle_modifiers(modifiers.state());
            }

            WindowEvent::KeyboardInput { event, .. } => {
                self.handle_key(event);
                if let Some(window) = &self.window {
                    window.request_redraw();
                }
            }

            _ => {}
        }
    }
}

// ============================================================================
// ViewerLike implementation
// ============================================================================

impl ViewerLike for Viewer {
    fn objects(&self) -> &ObjectRegistry {
        &self.registry
    }

    fn objects_mut(&mut self) -> &mut ObjectRegistry {
        &mut self.registry
    }

    fn camera(&self) -> &Camera {
        &self.camera
    }

    fn camera_mut(&mut self) -> &mut Camera {
        &mut self.camera
    }

    fn zoom_all(&mut self) {
        Viewer::zoom_all(self)
    }

    fn zoom_on(&mut self, name: &str) {
        Viewer::zoom_on(self, name)
    }

    fn center_all(&mut self) {
        Viewer::center_all(self)
    }

    fn center_on(&mut self, name: &str) {
        Viewer::center_on(self, name)
    }

    fn reset_view(&mut self) {
        Viewer::reset_view(self)
    }

    fn request_redraw(&mut self) {
        Viewer::request_redraw(self)
    }

    fn color_index(&self, name: &str) -> Option<u32> {
        Viewer::color_index(self, name)
    }

    fn set_background_color(&mut self, r: f32, g: f32, b: f32) {
        Viewer::set_background_color(self, r, g, b)
    }

    fn get_selection(&self, name: &str) -> Option<&str> {
        self.selections.get(name).map(|e| e.expression.as_str())
    }

    fn define_selection(&mut self, name: &str, selection: &str) {
        self.selections.insert(name.to_string(), SelectionEntry::new(selection.to_string()));
        self.needs_redraw = true;
    }

    fn remove_selection(&mut self, name: &str) -> bool {
        let removed = self.selections.remove(name).is_some();
        if removed {
            self.needs_redraw = true;
        }
        removed
    }

    fn selection_names(&self) -> Vec<String> {
        self.selections.keys().cloned().collect()
    }

    fn set_selection_visible(&mut self, name: &str, visible: bool) {
        if let Some(entry) = self.selections.get_mut(name) {
            entry.visible = visible;
            self.needs_redraw = true;
        }
    }

    fn is_selection_visible(&self, name: &str) -> bool {
        self.selections.get(name).map(|e| e.visible).unwrap_or(false)
    }

    fn indicate_selection(&mut self, selection: &str) {
        Viewer::indicate_selection(self, selection)
    }

    fn clear_indication(&mut self) {
        Viewer::clear_indication(self)
    }

    fn indicated_selection(&self) -> Option<&str> {
        Viewer::indicated_selection(self)
    }

    fn capture_png(
        &mut self,
        path: &Path,
        width: Option<u32>,
        height: Option<u32>,
    ) -> Result<(), String> {
        Viewer::capture_png(self, path, width, height).map_err(|e| e.to_string())
    }
}

/// Run the viewer with an event loop
pub fn run(mut viewer: Viewer) -> Result<(), ViewerError> {
    let event_loop = EventLoop::new().map_err(|e| ViewerError::Window(WindowError::CreationFailed(e.to_string())))?;
    event_loop.set_control_flow(ControlFlow::Wait);
    event_loop.run_app(&mut viewer).map_err(|e| ViewerError::Window(WindowError::CreationFailed(e.to_string())))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_viewer_creation() {
        // Just test that we can create a viewer without panicking
        let viewer = Viewer::new();
        assert!(viewer.window.is_none());
        assert!(viewer.registry.is_empty());
    }
}
