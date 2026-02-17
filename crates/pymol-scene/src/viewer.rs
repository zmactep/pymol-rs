//! Main viewer and render loop
//!
//! The [`Viewer`] struct is the main application entry point, coordinating
//! window management, input handling, and rendering.

use std::path::Path;
use std::sync::Arc;
use std::time::Instant;

use pymol_color::{ChainColors, ElementColors, NamedColors};
use pymol_render::silhouette::SilhouettePipeline;
use pymol_render::{ColorResolver, RenderContext};
use pymol_settings::GlobalSettings;
use winit::application::ApplicationHandler;
use winit::event::{KeyEvent, WindowEvent};
use winit::event_loop::{ActiveEventLoop, ControlFlow, EventLoop};
use winit::keyboard::PhysicalKey;
use winit::window::WindowId;

use pymol_mol::RepMask;

use crate::camera::Camera;
use crate::capture::capture_png_to_file;
use crate::error::{ViewerError, WindowError};
use crate::input::{CameraDelta, InputState};
use crate::keybindings::{KeyBinding, KeyBindings};
use crate::movie::Movie;
use crate::object::{MoleculeObject, Object, ObjectRegistry};
use crate::raytrace::{raytrace_scene, RaytraceInput};
use crate::scene::{SceneManager, SceneStoreMask};
use crate::selection::SelectionManager;
use crate::uniform::setup_uniforms;
use crate::view::ViewManager;
use crate::viewer_trait::{RaytracedImage, ViewerLike};
use crate::window::Window;

/// Type alias for key action callbacks
///
/// Key actions are closures that receive a mutable reference to the Viewer
/// and can perform any operation on it.
pub type KeyAction = Arc<dyn Fn(&mut Viewer) + Send + Sync>;

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
    /// Silhouette edge rendering pipeline
    silhouette_pipeline: Option<SilhouettePipeline>,

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
    /// Movie player for frame-based animation
    movie: Movie,
    /// Named selections manager
    selections: SelectionManager,
    /// Named views (simpler than scenes - just camera state)
    views: ViewManager,

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
    // Raytraced Image Overlay
    // =========================================================================
    /// Stored raytraced image for display (from `ray` command without filename)
    raytraced_image: Option<RaytracedImage>,

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
            silhouette_pipeline: None,
            window: None,
            pending_window: Some(("PyMOL-RS".to_string(), 1024, 768)),
            camera: Camera::new(),
            registry: ObjectRegistry::new(),
            scenes: SceneManager::new(),
            movie: Movie::new(),
            selections: SelectionManager::new(),
            views: ViewManager::new(),
            settings: GlobalSettings::new(),
            named_colors: NamedColors::default(),
            element_colors: ElementColors::default(),
            chain_colors: ChainColors,
            input: InputState::new(),
            last_frame: Instant::now(),
            needs_redraw: true,
            clear_color: [0.0, 0.0, 0.0],
            raytraced_image: None,
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

        // Request device with adapter's actual limits (not conservative defaults)
        let (device, queue) = adapter
            .request_device(&wgpu::DeviceDescriptor {
                label: Some("PyMOL-RS Device"),
                required_features: wgpu::Features::empty(),
                required_limits: adapter.limits(), // Use hardware limits, not WebGL2 defaults
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

        // Create silhouette pipeline
        let silhouette_pipeline = SilhouettePipeline::new(render_context.device(), surface_format);

        // Set camera aspect ratio
        self.camera.set_aspect(window.aspect_ratio());

        self.window = Some(window);
        self.render_context = Some(render_context);
        self.silhouette_pipeline = Some(silhouette_pipeline);
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
            self.raytraced_image = None; // Clear raytraced overlay on camera change
            self.needs_redraw = true;
        }
    }

    /// Zoom the camera to fit a specific object while preserving rotation
    pub fn zoom_on(&mut self, name: &str) {
        if let Some(obj) = self.registry.get(name) {
            if let Some((min, max)) = obj.extent() {
                self.camera.zoom_to(min, max);
                self.raytraced_image = None; // Clear raytraced overlay on camera change
                self.needs_redraw = true;
            }
        }
    }

    /// Center the camera on all objects without changing zoom or rotation
    pub fn center_all(&mut self) {
        if let Some((min, max)) = self.registry.extent() {
            self.camera.center_to(min, max);
            self.raytraced_image = None; // Clear raytraced overlay on camera change
            self.needs_redraw = true;
        }
    }

    /// Center the camera on a specific object without changing zoom or rotation
    pub fn center_on(&mut self, name: &str) {
        if let Some(obj) = self.registry.get(name) {
            if let Some((min, max)) = obj.extent() {
                self.camera.center_to(min, max);
                self.raytraced_image = None; // Clear raytraced overlay on camera change
                self.needs_redraw = true;
            }
        }
    }

    /// Reset the camera to default view (resets everything including rotation)
    pub fn reset_view(&mut self) {
        self.camera = Camera::new();
        self.raytraced_image = None; // Clear raytraced overlay on camera change
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
        self.selections.indicate(selection);
        self.needs_redraw = true;
    }

    /// Clear the indicated selection
    ///
    /// Hides all selection indicators in the 3D view.
    pub fn clear_indication(&mut self) {
        self.selections.clear_indication();
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
        self.selections.indicated_selection()
    }

    // =========================================================================
    // Representation Control
    // =========================================================================

    /// Toggle a representation for all molecules
    ///
    /// If the representation is currently visible on any molecule, it will be hidden.
    /// Otherwise, it will be shown on all molecules.
    pub fn toggle_representation(&mut self, rep: pymol_mol::RepMask) {
        let names: Vec<_> = self.registry.names().map(|s| s.to_string()).collect();
        for name in names {
            if let Some(mol) = self.registry.get_molecule_mut(&name) {
                mol.toggle(rep);
            }
        }
        self.needs_redraw = true;
    }

    /// Show a representation for all molecules
    pub fn show_representation(&mut self, rep: pymol_mol::RepMask) {
        let names: Vec<_> = self.registry.names().map(|s| s.to_string()).collect();
        for name in names {
            if let Some(mol) = self.registry.get_molecule_mut(&name) {
                mol.show(rep);
            }
        }
        self.needs_redraw = true;
    }

    /// Hide a representation for all molecules
    pub fn hide_representation(&mut self, rep: pymol_mol::RepMask) {
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

        // Update uniforms using shared setup function
        let uniforms = setup_uniforms(
            &self.camera,
            &self.settings,
            self.clear_color,
            (window.width() as f32, window.height() as f32),
        );
        context.update_uniforms(&uniforms);

        // Prepare molecules for rendering
        // Create color resolver with borrowed data - need to avoid borrowing self during registry iteration
        let names: Vec<_> = self.registry.names().map(|s| s.to_string()).collect();

        // Evaluate all visible selections using the SelectionManager
        let selection_results = self.selections.evaluate_visible(&self.registry);

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

        // Silhouette pass (after main scene, before submit)
        if self.settings.get_bool(pymol_settings::id::silhouettes) {
            if let Some(silhouette) = &self.silhouette_pipeline {
                let width = window.width();
                let height = window.height();
                let thickness = self.settings.get_float(pymol_settings::id::silhouette_width);
                let depth_jump = self.settings.get_float(pymol_settings::id::silhouette_depth_jump);
                let color_rgb = self.settings.get_float3(pymol_settings::id::silhouette_color);
                let color = [color_rgb[0], color_rgb[1], color_rgb[2], 1.0];
                silhouette.render(
                    &mut encoder,
                    context.queue(),
                    context.device(),
                    &view,
                    window.depth_view(),
                    width,
                    height,
                    thickness,
                    depth_jump,
                    color,
                    None,
                );
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

    /// Perform raytracing and return the image data as RGBA bytes
    ///
    /// # Arguments
    /// * `width` - Image width (None = use viewport width)
    /// * `height` - Image height (None = use viewport height)
    /// * `antialias` - Antialiasing level (1 = no AA, 2-4 = supersampling)
    ///
    /// # Returns
    /// * RGBA image data, row-major
    pub fn raytrace(
        &mut self,
        width: Option<u32>,
        height: Option<u32>,
        antialias: u32,
    ) -> Result<(Vec<u8>, u32, u32), ViewerError> {
        self.prepare_representations_for_raytrace()?;

        let context = self.render_context.as_ref().ok_or(ViewerError::Window(
            WindowError::NotAvailable,
        ))?;
        let default_size = self.window.as_ref()
            .map(|w| w.size())
            .unwrap_or((1024, 768));

        let mut input = RaytraceInput {
            device: context.device(),
            queue: context.queue(),
            camera: &mut self.camera,
            registry: &self.registry,
            settings: &self.settings,
            named_colors: &self.named_colors,
            element_colors: &self.element_colors,
            chain_colors: &self.chain_colors,
            clear_color: self.clear_color,
            default_size,
        };

        raytrace_scene(&mut input, width, height, antialias)
            .map_err(|e| ViewerError::GpuInitFailed(e.to_string()))
    }

    /// Perform raytracing and save to a PNG file
    pub fn raytrace_to_file(
        &mut self,
        path: impl AsRef<std::path::Path>,
        width: Option<u32>,
        height: Option<u32>,
        antialias: u32,
    ) -> Result<(u32, u32), ViewerError> {
        self.prepare_representations_for_raytrace()?;

        let context = self.render_context.as_ref().ok_or(ViewerError::Window(
            WindowError::NotAvailable,
        ))?;
        let default_size = self.window.as_ref()
            .map(|w| w.size())
            .unwrap_or((1024, 768));

        let mut input = RaytraceInput {
            device: context.device(),
            queue: context.queue(),
            camera: &mut self.camera,
            registry: &self.registry,
            settings: &self.settings,
            named_colors: &self.named_colors,
            element_colors: &self.element_colors,
            chain_colors: &self.chain_colors,
            clear_color: self.clear_color,
            default_size,
        };

        crate::raytrace::raytrace_to_file(&mut input, path, width, height, antialias)
            .map_err(|e| ViewerError::GpuInitFailed(e.to_string()))
    }

    /// Ensure all representations are built before raytracing.
    ///
    /// This is needed when raytracing from scripts where the render loop hasn't run.
    fn prepare_representations_for_raytrace(&mut self) -> Result<(), ViewerError> {
        let context = self.render_context.as_ref().ok_or(ViewerError::Window(
            WindowError::NotAvailable,
        ))?;

        let names: Vec<String> = self.registry.names().map(|s| s.to_string()).collect();
        for name in &names {
            let color_resolver = ColorResolver::new(&self.named_colors, &self.element_colors, &self.chain_colors);
            if let Some(mol_obj) = self.registry.get_molecule_mut(name) {
                mol_obj.prepare_render(context, &color_resolver, &self.settings);
            }
        }

        Ok(())
    }

    /// Update the viewer state
    fn update(&mut self, dt: f32) {
        // Process input deltas
        let deltas = self.input.take_camera_deltas();
        for delta in deltas {
            // Clear raytraced overlay on any camera change
            if self.raytraced_image.is_some() {
                self.raytraced_image = None;
            }
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

        // Update movie playback
        if self.movie.is_playing() {
            if self.movie.update() {
                // Frame changed: apply interpolated view
                if let Some(view) = self.movie.interpolated_view() {
                    self.camera.set_view(view);
                }
                self.needs_redraw = true;
            }
        }

        // Update rock animation
        if self.movie.is_rock_enabled() {
            let delta = self.movie.update_rock(dt, 3.0_f32.to_radians(), 1.0);
            if delta.abs() > 1e-6 {
                self.camera.rotate_y(delta);
                self.needs_redraw = true;
            }
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

                // Request continuous redraw if animating, playing, rocking, or input is active
                if self.camera.is_animating()
                    || self.input.any_button_pressed()
                    || self.movie.is_playing()
                    || self.movie.is_rock_enabled()
                {
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
    // =========================================================================
    // Required Accessors
    // =========================================================================

    fn objects(&self) -> &ObjectRegistry { &self.registry }
    fn objects_mut(&mut self) -> &mut ObjectRegistry { &mut self.registry }
    fn camera(&self) -> &Camera { &self.camera }
    fn camera_mut(&mut self) -> &mut Camera { &mut self.camera }
    fn settings(&self) -> &GlobalSettings { &self.settings }
    fn settings_mut(&mut self) -> &mut GlobalSettings { &mut self.settings }
    fn movie(&self) -> &Movie { &self.movie }
    fn movie_mut(&mut self) -> &mut Movie { &mut self.movie }
    fn scenes(&self) -> &SceneManager { &self.scenes }
    fn scenes_mut(&mut self) -> &mut SceneManager { &mut self.scenes }
    fn views(&self) -> &ViewManager { &self.views }
    fn views_mut(&mut self) -> &mut ViewManager { &mut self.views }
    fn selections(&self) -> &SelectionManager { &self.selections }
    fn selections_mut(&mut self) -> &mut SelectionManager { &mut self.selections }
    fn named_colors(&self) -> &pymol_color::NamedColors { &self.named_colors }
    fn clear_color(&self) -> [f32; 3] { self.clear_color }
    fn set_clear_color(&mut self, color: [f32; 3]) { self.clear_color = color; }
    fn raytraced_image_ref(&self) -> Option<&RaytracedImage> { self.raytraced_image.as_ref() }
    fn set_raytraced_image_internal(&mut self, image: Option<RaytracedImage>) { self.raytraced_image = image; }

    fn request_redraw(&mut self) {
        self.needs_redraw = true;
        if let Some(window) = &self.window {
            window.request_redraw();
        }
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
        Viewer::capture_png(self, path, width, height).map_err(|e| e.to_string())
    }

    fn raytrace(
        &mut self,
        width: Option<u32>,
        height: Option<u32>,
        antialias: u32,
    ) -> Result<Vec<u8>, String> {
        Viewer::raytrace(self, width, height, antialias)
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
        Viewer::raytrace_to_file(self, path, width, height, antialias).map_err(|e| e.to_string())
    }

    // =========================================================================
    // Per-impl methods (borrow checker constraints)
    // =========================================================================

    fn scene_store(&mut self, key: &str, storemask: u32) {
        let mask = SceneStoreMask::from_bits_truncate(storemask);
        self.scenes.store(key, mask, &self.camera, &self.registry);
        self.needs_redraw = true;
    }

    fn scene_recall(&mut self, key: &str, animate: bool, duration: f32) -> Result<(), String> {
        self.scenes
            .recall(key, &mut self.camera, &mut self.registry, animate, duration)
            .map_err(|e| e.to_string())?;
        self.needs_redraw = true;
        Ok(())
    }

    fn view_recall(&mut self, key: &str, animate: f32) -> Result<(), String> {
        self.views
            .recall(key, &mut self.camera, animate)
            .map_err(|e| e.to_string())?;
        self.needs_redraw = true;
        Ok(())
    }

    fn capture_frame_png(
        &mut self,
        frame: usize,
        path: &Path,
        width: Option<u32>,
        height: Option<u32>,
    ) -> Result<(), String> {
        self.movie.goto_frame(frame);
        if let Some(view) = self.movie.interpolated_view() {
            self.camera.set_view(view);
        }
        if let Some(scene_name) = self.movie.current_scene_name().map(|s| s.to_string()) {
            if let Some(scene) = self.scenes.get(&scene_name) {
                scene.apply(&mut self.camera, &mut self.registry, false, 0.0);
                if let Some(view) = self.movie.interpolated_view() {
                    self.camera.set_view(view);
                }
            }
        }
        Viewer::capture_png(self, path, width, height).map_err(|e| e.to_string())
    }

    // =========================================================================
    // Window-specific overrides
    // =========================================================================

    fn viewport_size(&self) -> (u32, u32) {
        self.window.as_ref().map(|w| w.size()).unwrap_or((0, 0))
    }

    fn viewport_set_size(&mut self, width: u32, height: u32) {
        if let Some(window) = &self.window {
            let _ = window.window().request_inner_size(winit::dpi::LogicalSize::new(width, height));
            self.needs_redraw = true;
        }
    }

    fn is_fullscreen(&self) -> bool {
        self.window.as_ref().map_or(false, |w| w.window().fullscreen().is_some())
    }

    fn set_fullscreen(&mut self, enabled: bool) {
        if let Some(window) = &self.window {
            if enabled {
                window.window().set_fullscreen(Some(winit::window::Fullscreen::Borderless(None)));
            } else {
                window.window().set_fullscreen(None);
            }
            self.needs_redraw = true;
        }
    }

    fn toggle_fullscreen(&mut self) {
        let is_full = self.is_fullscreen();
        self.set_fullscreen(!is_full);
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
