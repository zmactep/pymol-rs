//! Main Application
//!
//! The App struct is the main entry point that combines the Viewer, CommandExecutor,
//! and egui UI into a single application.

use std::sync::Arc;
use std::time::Instant;

use egui::ViewportId;
use pymol_cmd::{ArgHint, CmdError, CommandExecutor, MessageKind, ParsedCommand, parse_command};
use pymol_color::{ChainColors, ElementColors, NamedColors};
use pymol_mol::RepMask;
use pymol_render::{ColorResolver, GlobalUniforms, RenderContext};
use pymol_scene::{MoleculeObject, Object, ObjectRegistry};
use pymol_select::{EvalContext, SelectionResult};
use pymol_scene::{Camera, KeyBinding, KeyBindings};
use pymol_scene::{CameraDelta, InputState};
// Re-export SelectionEntry for use in UI
pub use pymol_scene::SelectionEntry;
use pymol_settings::GlobalSettings;
use winit::application::ApplicationHandler;
use winit::dpi::PhysicalSize;
use winit::event::{ElementState, KeyEvent, WindowEvent};
use winit::event_loop::ActiveEventLoop;
use winit::keyboard::PhysicalKey;
use winit::window::{Window, WindowId};

use crate::async_tasks::{TaskContext, TaskRunner};
use crate::ipc::{
    ExternalCommandRegistry, IpcCallbackTask, IpcRequest, IpcResponse, IpcServer, OutputKind,
};
use crate::state::GuiState;
use crate::ui::command::CommandAction;
use crate::ui::objects::ObjectAction;
use crate::ui::{CommandLinePanel, NotificationOverlay, ObjectListPanel, OutputPanel};
use crate::viewer_adapter::ViewerAdapter;

/// Type alias for key action callbacks
pub type KeyAction = Arc<dyn Fn(&mut App) + Send + Sync>;

/// Main application state
pub struct App {
    // =========================================================================
    // GPU Resources
    // =========================================================================
    /// wgpu instance
    instance: wgpu::Instance,
    /// Render context for molecular rendering (owns device and queue)
    render_context: Option<RenderContext>,

    // =========================================================================
    // Window State
    // =========================================================================
    /// The main window
    window: Option<Arc<Window>>,
    /// Window surface
    surface: Option<wgpu::Surface<'static>>,
    /// Surface configuration
    surface_config: Option<wgpu::SurfaceConfiguration>,
    /// Depth texture view
    depth_view: Option<wgpu::TextureView>,

    // =========================================================================
    // egui Integration
    // =========================================================================
    /// egui context
    egui_ctx: egui::Context,
    /// egui-winit state
    egui_state: Option<egui_winit::State>,
    /// egui-wgpu renderer
    egui_renderer: Option<egui_wgpu::Renderer>,
    /// Viewport rect for 3D rendering (in physical pixels)
    viewport_rect: Option<egui::Rect>,

    // =========================================================================
    // Scene State
    // =========================================================================
    /// Camera for view control
    camera: Camera,
    /// Object registry
    registry: ObjectRegistry,
    /// Named selections (name -> entry with expression and visibility)
    selections: std::collections::HashMap<String, SelectionEntry>,
    /// Global settings
    settings: GlobalSettings,
    /// Named colors
    named_colors: NamedColors,
    /// Element colors
    element_colors: ElementColors,
    /// Chain colors
    chain_colors: ChainColors,

    // =========================================================================
    // Command System
    // =========================================================================
    /// Command executor (for future use with full command parsing)
    #[allow(dead_code)]
    executor: CommandExecutor,

    // =========================================================================
    // GUI State
    // =========================================================================
    /// GUI state (output, command line, etc.)
    gui_state: GuiState,
    /// Object list panel (stateful for menu handling)
    object_list_panel: ObjectListPanel,

    // =========================================================================
    // Input State
    // =========================================================================
    /// Input handler (from pymol-scene, handles mouse with proper sensitivity)
    input: InputState,

    // =========================================================================
    // Frame Timing
    // =========================================================================
    /// Last frame timestamp
    last_frame: Instant,
    /// Whether a redraw is needed
    needs_redraw: bool,
    /// Frame counter for initial setup (egui needs a few frames to layout properly)
    frame_count: u32,

    // =========================================================================
    // Appearance
    // =========================================================================
    /// Clear/background color
    clear_color: [f32; 3],

    // =========================================================================
    // Key Bindings
    // =========================================================================
    /// Keyboard shortcuts
    key_bindings: KeyBindings<KeyAction>,

    // =========================================================================
    // Async Task System
    // =========================================================================
    /// Task runner for background operations (fetch, etc.)
    task_runner: TaskRunner,

    // =========================================================================
    // IPC (Inter-Process Communication)
    // =========================================================================
    /// Optional IPC server for external control (e.g., from Python)
    ipc_server: Option<IpcServer>,
    /// Registry for external commands registered via IPC
    external_commands: ExternalCommandRegistry,

    // =========================================================================
    // Headless Mode
    // =========================================================================
    /// Whether the application is running in headless mode (no window displayed)
    headless: bool,

    // =========================================================================
    // Pending Actions (initialization)
    // =========================================================================
    /// File path to load after GPU initialization
    pending_load_file: Option<String>,

    // =========================================================================
    // Application Lifecycle
    // =========================================================================
    /// Whether the quit command was issued
    quit_requested: bool,
}

// ============================================================================
// TaskContext implementation for App
// ============================================================================

impl TaskContext for App {
    fn add_molecule(&mut self, name: &str, mol: pymol_mol::ObjectMolecule) {
        self.registry.add(MoleculeObject::with_name(mol, name));
    }

    fn execute_command(&mut self, cmd: &str) {
        self.execute_command(cmd);
    }

    fn print_info(&mut self, msg: String) {
        self.gui_state.print_info(msg);
    }

    fn print_warning(&mut self, msg: String) {
        self.gui_state.print_warning(msg);
    }

    fn print_error(&mut self, msg: String) {
        self.gui_state.print_error(msg);
    }
}

impl Default for App {
    fn default() -> Self {
        Self::new(false)
    }
}

impl App {
    /// Create a new application
    ///
    /// # Arguments
    /// * `headless` - If true, the window will not be displayed initially
    pub fn new(headless: bool) -> Self {
        let instance = wgpu::Instance::new(&wgpu::InstanceDescriptor {
            backends: wgpu::Backends::all(),
            ..Default::default()
        });

        // Create executor and extract command names for autocomplete
        let executor = CommandExecutor::new();
        let registry = executor.registry();
        
        // Get all command names (including aliases)
        let command_names: Vec<String> = registry
            .all_names()
            .map(String::from)
            .collect();
        
        // Get commands that take file paths as first argument (position 0)
        let path_commands: Vec<String> = registry
            .commands_with_hint(ArgHint::Path, 0)
            .into_iter()
            .map(String::from)
            .collect();

        let mut gui_state = GuiState::new();
        gui_state.command_names = command_names;
        gui_state.path_commands = path_commands;

        let mut app = Self {
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
            camera: Camera::new(),
            registry: ObjectRegistry::new(),
            selections: std::collections::HashMap::new(),
            settings: GlobalSettings::new(),
            named_colors: NamedColors::default(),
            element_colors: ElementColors::default(),
            chain_colors: ChainColors,
            executor,
            gui_state,
            object_list_panel: ObjectListPanel::new(),
            input: InputState::new(),
            last_frame: Instant::now(),
            needs_redraw: true,
            frame_count: 0,
            clear_color: [0.0, 0.0, 0.0],
            key_bindings: KeyBindings::new(),
            task_runner: TaskRunner::new(),
            ipc_server: None,
            external_commands: ExternalCommandRegistry::new(),
            headless,
            pending_load_file: None,
            quit_requested: false,
        };

        // Set up default key bindings
        app.setup_default_key_bindings();

        app
    }

    /// Create a new application with IPC server enabled
    ///
    /// # Arguments
    /// * `socket_path` - Path to the Unix domain socket for IPC
    /// * `headless` - If true, the window will not be displayed initially
    pub fn with_ipc(socket_path: &std::path::Path, headless: bool) -> Result<Self, String> {
        let mut app = Self::new(headless);
        
        let server = IpcServer::bind(socket_path)
            .map_err(|e| format!("Failed to create IPC server: {}", e))?;
        
        app.ipc_server = Some(server);
        app.gui_state.print_info("IPC server enabled".to_string());
        
        Ok(app)
    }

    /// Check if IPC is enabled
    pub fn ipc_enabled(&self) -> bool {
        self.ipc_server.is_some()
    }

    /// Check if the application is running in headless mode
    pub fn is_headless(&self) -> bool {
        self.headless
    }

    /// Show the window (makes it visible)
    ///
    /// This has no effect if the window hasn't been created yet.
    ///
    /// # Platform-specific
    /// - Android / Wayland / Web: Unsupported
    pub fn show_window(&self) {
        if let Some(window) = &self.window {
            window.set_visible(true);
        }
    }

    /// Hide the window (makes it invisible)
    ///
    /// This has no effect if the window hasn't been created yet.
    ///
    /// # Platform-specific
    /// - Android / Wayland / Web: Unsupported
    pub fn hide_window(&self) {
        if let Some(window) = &self.window {
            window.set_visible(false);
        }
    }

    /// Returns whether the window is currently visible
    ///
    /// Returns `None` if the window hasn't been created yet or if
    /// visibility cannot be determined on the current platform.
    ///
    /// # Platform-specific
    /// - X11: Not implemented (always returns `None`)
    /// - Wayland / iOS / Android / Web: Unsupported (always returns `None`)
    pub fn is_window_visible(&self) -> Option<bool> {
        self.window.as_ref().and_then(|w| w.is_visible())
    }

    /// Setup default keyboard shortcuts
    fn setup_default_key_bindings(&mut self) {
        // TODO: Implement default key bindings
    }

    /// Bind a key to an action
    pub fn bind_key<K, F>(&mut self, key: K, action: F)
    where
        K: Into<KeyBinding>,
        F: Fn(&mut App) + Send + Sync + 'static,
    {
        self.key_bindings.bind(key, Arc::new(action));
    }

    /// Queue a file to load after initialization
    pub fn queue_load_file(&mut self, path: String) {
        self.pending_load_file = Some(path);
    }

    /// Initialize GPU resources
    async fn init_gpu(&mut self, window: Arc<Window>) -> Result<(), String> {
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

        // Log adapter info
        let info = adapter.get_info();
        self.gui_state.print_info(format!("DEVICE: {}", info.name));
        self.gui_state.print_info(format!("ENGINE:  {}", info.backend));

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

        // Initialize egui first (before render context takes ownership of device/queue)
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

        // Set camera aspect ratio
        self.camera.set_aspect(size.width as f32 / size.height as f32);

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

    /// Create depth texture
    fn create_depth_texture(device: &wgpu::Device, width: u32, height: u32) -> wgpu::TextureView {
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
    fn resize(&mut self, new_size: PhysicalSize<u32>) {
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

            // Update camera aspect ratio
            self.camera.set_aspect(new_size.width as f32 / new_size.height as f32);

            self.needs_redraw = true;
        }
    }

    /// Render a frame
    fn render(&mut self) -> Result<(), wgpu::SurfaceError> {
        // Extract config values before any mutable borrows
        let (width, height) = {
            let config = self.surface_config.as_ref().unwrap();
            (config.width, config.height)
        };

        let scale_factor = self.window.as_ref().map(|w| w.scale_factor()).unwrap_or(1.0) as f32;

        // Run egui UI first to get output (this borrows self mutably for draw_ui)
        let egui_output = self.run_egui_ui();

        // Update camera aspect ratio based on viewport dimensions
        let (viewport_width, viewport_height) = if let Some(viewport) = &self.viewport_rect {
            (viewport.width().max(1.0), viewport.height().max(1.0))
        } else {
            (width as f32, height as f32)
        };
        self.camera.set_aspect(viewport_width / viewport_height);

        // Now get immutable references for rendering
        let surface = self.surface.as_ref().unwrap();
        let depth_view = self.depth_view.as_ref().unwrap();
        let context = self.render_context.as_ref().unwrap();
        let device = context.device();
        let queue = context.queue();

        // Get surface texture
        let output = surface.get_current_texture()?;
        let view = output.texture.create_view(&wgpu::TextureViewDescriptor::default());

        // Update global uniforms
        let mut uniforms = GlobalUniforms::new();
        uniforms.set_camera(self.camera.view_matrix(), self.camera.projection_matrix());
        uniforms.set_background(self.clear_color);
        uniforms.set_viewport(viewport_width, viewport_height);

        let camera_pos = self.camera.world_position();
        uniforms.camera_pos = [camera_pos.x, camera_pos.y, camera_pos.z, 1.0];

        // Lighting settings
        let ambient = self.settings.get_float(pymol_settings::id::ambient);
        let direct = self.settings.get_float(pymol_settings::id::direct);
        let reflect = self.settings.get_float(pymol_settings::id::reflect);
        let specular = self.settings.get_float(pymol_settings::id::specular);
        let shininess = self.settings.get_float(pymol_settings::id::shininess);
        let spec_direct = self.settings.get_float(pymol_settings::id::spec_direct);
        let spec_direct_power = self.settings.get_float(pymol_settings::id::spec_direct_power);
        uniforms.set_lighting(ambient, direct, reflect, specular, shininess, spec_direct, spec_direct_power);

        // Clip planes
        let current_view = self.camera.current_view();
        uniforms.set_clip_planes(current_view.clip_front, current_view.clip_back);

        // Fog
        let depth_cue_enabled = self.settings.get_bool(pymol_settings::id::depth_cue);
        let fog_density = self.settings.get_float(pymol_settings::id::fog);
        if depth_cue_enabled && fog_density > 0.0 {
            let fog_start_ratio = self.settings.get_float(pymol_settings::id::fog_start);
            let fog_start = (current_view.clip_back - current_view.clip_front) * fog_start_ratio
                + current_view.clip_front;
            let fog_end = if (fog_density - 1.0).abs() < 0.001 {
                current_view.clip_back
            } else {
                fog_start + (current_view.clip_back - fog_start) / fog_density
            };
            uniforms.set_fog(fog_start, fog_end, fog_density, self.clear_color);
            uniforms.set_depth_cue(1.0);
        }

        context.update_uniforms(&uniforms);

        // Prepare molecules
        let names: Vec<String> = self.registry.names().map(|s: &str| s.to_string()).collect();
        
        // Evaluate all visible selections
        let selection_results = self.evaluate_visible_selections(&names);
        
        // Get selection indicator size from settings (selection_width, ID 80)
        // Use a minimum size to ensure visibility
        let selection_width = self.settings.get_float(pymol_settings::id::selection_width).max(6.0);
        
        for name in &names {
            let color_resolver = ColorResolver::new(&self.named_colors, &self.element_colors, &self.chain_colors);
            if let Some(mol_obj) = self.registry.get_molecule_mut(name) {
                mol_obj.prepare_render(context, &color_resolver, &self.settings);
                
                // Update selection indicator if we have results for this molecule
                if let Some((_, selection_result)) = selection_results.iter().find(|(n, _)| n == name) {
                    log::debug!("Setting selection indicator for '{}' with {} atoms", name, selection_result.count());
                    mol_obj.set_selection_indicator_with_size(selection_result, context, Some(selection_width));
                } else {
                    // Clear indicator - no visible selection matches this molecule
                    mol_obj.clear_selection_indicator();
                }
            }
        }

        // Get context again after mutable borrow ends
        let context = self.render_context.as_ref().unwrap();

        // Create encoder
        let mut encoder = device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("Render Encoder"),
        });

        // Render 3D scene
        {
            let mut render_pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("3D Render Pass"),
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
                    view: depth_view,
                    depth_ops: Some(wgpu::Operations {
                        load: wgpu::LoadOp::Clear(1.0),
                        store: wgpu::StoreOp::Store,
                    }),
                    stencil_ops: None,
                }),
                timestamp_writes: None,
                occlusion_query_set: None,
            });

            // Set viewport and scissor rect to clip 3D rendering to the central area
            if let Some(viewport) = &self.viewport_rect {
                let vp_x = viewport.min.x.max(0.0);
                let vp_y = viewport.min.y.max(0.0);
                let vp_w = viewport.width().max(1.0);
                let vp_h = viewport.height().max(1.0);
                
                render_pass.set_viewport(vp_x, vp_y, vp_w, vp_h, 0.0, 1.0);
                render_pass.set_scissor_rect(
                    vp_x as u32,
                    vp_y as u32,
                    vp_w as u32,
                    vp_h as u32,
                );
            }

            // Render all enabled objects
            for name in &names {
                if let Some(mol_obj) = self.registry.get_molecule(&name) {
                    if mol_obj.is_enabled() {
                        mol_obj.render(&mut render_pass, context);
                    }
                }
            }
        }

        // Render egui output
        if let (Some((clipped_primitives, textures_delta)), Some(egui_renderer)) =
            (egui_output, &mut self.egui_renderer)
        {
            for (id, image_delta) in &textures_delta.set {
                egui_renderer.update_texture(device, queue, *id, image_delta);
            }

            let screen_descriptor = egui_wgpu::ScreenDescriptor {
                size_in_pixels: [width, height],
                pixels_per_point: scale_factor,
            };

            egui_renderer.update_buffers(device, queue, &mut encoder, &clipped_primitives, &screen_descriptor);

            {
                let render_pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                    label: Some("egui Render Pass"),
                    color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                        view: &view,
                        resolve_target: None,
                        ops: wgpu::Operations {
                            load: wgpu::LoadOp::Load,
                            store: wgpu::StoreOp::Store,
                        },
                        depth_slice: None,
                    })],
                    depth_stencil_attachment: None,
                    timestamp_writes: None,
                    occlusion_query_set: None,
                });

                // Convert to 'static lifetime as required by egui-wgpu
                let mut render_pass = render_pass.forget_lifetime();
                egui_renderer.render(&mut render_pass, &clipped_primitives, &screen_descriptor);
            }

            for id in &textures_delta.free {
                egui_renderer.free_texture(id);
            }
        }

        queue.submit(std::iter::once(encoder.finish()));
        output.present();

        Ok(())
    }

    /// Run the egui UI and return the output for rendering later
    fn run_egui_ui(&mut self) -> Option<(Vec<egui::ClippedPrimitive>, egui::TexturesDelta)> {
        // Take egui input first
        let raw_input = {
            let (egui_state, window) = match (&mut self.egui_state, &self.window) {
                (Some(state), Some(window)) => (state, window),
                _ => return None,
            };
            egui_state.take_egui_input(window)
        };

        // Get scale factor for converting to physical pixels
        let scale_factor = self.window.as_ref().map(|w| w.scale_factor()).unwrap_or(1.0) as f32;

        // Collect data needed for UI without borrowing self in the closure
        let mut commands_to_execute = Vec::new();
        let mut object_actions = Vec::new();
        let mut viewport_rect_logical = egui::Rect::NOTHING;

        // Check for pending async tasks and get their messages before entering the closure
        let pending_messages = self.task_runner.pending_messages();

        // Run egui with a closure that captures only what it needs
        let full_output = {
            let gui_state = &mut self.gui_state;
            let registry = &self.registry;
            let selections = &self.selections;
            let object_list_panel = &mut self.object_list_panel;
            
            self.egui_ctx.run(raw_input, |ctx| {
                // Top panel - output and command line
                egui::TopBottomPanel::top("top_panel").show(ctx, |ui| {
                    if gui_state.show_output_panel {
                        OutputPanel::show(ui, gui_state);
                        ui.separator();
                    }

                    match CommandLinePanel::show(ui, gui_state) {
                        CommandAction::Execute(cmd) => {
                            commands_to_execute.push(cmd);
                        }
                        CommandAction::None => {}
                    }
                });

                // Right panel - object list only
                if gui_state.show_control_panel {
                    egui::SidePanel::right("right_panel")
                        .default_width(gui_state.right_panel_width)
                        .show(ctx, |ui| {
                            egui::ScrollArea::vertical().show(ui, |ui| {
                                // Object list with selections (stateful panel)
                                object_actions = object_list_panel.show(ui, registry, selections);
                            });
                        });
                }

                // Central panel to capture the viewport rect (area for 3D rendering)
                let central_response = egui::CentralPanel::default()
                    .frame(egui::Frame::NONE)
                    .show(ctx, |ui| {
                        // Make the central panel interactive to detect hover
                        ui.allocate_response(ui.available_size(), egui::Sense::hover())
                    });
                viewport_rect_logical = central_response.response.rect;
                
                // Track if mouse is over the viewport
                gui_state.viewport_hovered = central_response.inner.hovered();

                // Show notification overlay when async tasks are in progress
                if !pending_messages.is_empty() {
                    NotificationOverlay::show(ctx, &pending_messages);
                }
            })
        };

        // Convert viewport rect from logical to physical pixels
        let viewport_rect_physical = egui::Rect::from_min_max(
            egui::pos2(
                viewport_rect_logical.min.x * scale_factor,
                viewport_rect_logical.min.y * scale_factor,
            ),
            egui::pos2(
                viewport_rect_logical.max.x * scale_factor,
                viewport_rect_logical.max.y * scale_factor,
            ),
        );
        self.viewport_rect = Some(viewport_rect_physical);

        // Handle platform output
        if let (Some(egui_state), Some(window)) = (&mut self.egui_state, &self.window) {
            egui_state.handle_platform_output(window, full_output.platform_output);
        }

        // Process commands and actions
        let has_actions = !commands_to_execute.is_empty() || !object_actions.is_empty();
        for cmd in commands_to_execute {
            self.execute_command(&cmd);
        }
        for action in object_actions {
            self.handle_object_action(action);
        }

        // Check if egui needs a repaint (e.g., for popup menus, animations, etc.)
        // repaint_delay of Duration::ZERO means immediate repaint needed
        let egui_needs_repaint = full_output
            .viewport_output
            .values()
            .any(|v| v.repaint_delay.is_zero());

        // Request window redraw if we processed any actions OR if egui needs repaint
        // This is needed because the event loop is in Wait mode and won't redraw
        // unless explicitly requested
        if has_actions || egui_needs_repaint {
            if let Some(window) = &self.window {
                window.request_redraw();
            }
        }

        let clipped_primitives = self.egui_ctx.tessellate(full_output.shapes, full_output.pixels_per_point);
        
        Some((clipped_primitives, full_output.textures_delta))
    }

    /// Execute a command using the CommandExecutor with ViewerAdapter
    ///
    /// This method is IPC-aware - if the command is an external command
    /// registered via IPC, it will send a callback request to the client.
    fn execute_command(&mut self, cmd: &str) {
        let cmd = cmd.trim();
        if cmd.is_empty() || cmd.starts_with('#') {
            return;
        }

        // Parse using pymol-cmd parser for proper handling of quotes, commas, etc.
        let parsed = match parse_command(cmd) {
            Ok(p) => p,
            Err(e) => {
                self.gui_state.print_error(format!("Parse error: {}", e));
                return;
            }
        };

        // Handle "help <external_command>" - check if help target is an external command
        if parsed.name == "help" {
            if let Some(target_cmd) = parsed.get_str(0) {
                if let Some(help_text) = self.external_commands.help(target_cmd) {
                    self.gui_state.print_info(help_text);
                    return;
                } else if self.external_commands.contains(target_cmd) {
                    // External command without help text
                    self.gui_state.print_info(format!(" '{}' - external command (no help available)", target_cmd));
                    return;
                }
                // Fall through to built-in help command
            }
        }

        // Check if it's an external command (registered via IPC)
        if self.external_commands.contains(&parsed.name) {
            self.execute_external_command(&parsed);
            return;
        }

        // Execute as built-in command using CommandExecutor
        self.execute_builtin_command_internal(cmd);
    }

    /// Execute an external command via IPC callback
    ///
    /// This spawns an async task that waits for the callback response.
    fn execute_external_command(&mut self, parsed: &ParsedCommand) {
        let server = match self.ipc_server.as_mut() {
            Some(s) => s,
            None => {
                self.gui_state.print_error("No IPC server available");
                return;
            }
        };

        let id = server.next_callback_id();
        // Convert parsed arguments to strings using ArgValue's Display impl
        let args: Vec<String> = parsed.args
            .iter()
            .map(|(_, v)| v.to_string())
            .collect();

        // Create oneshot channel for the callback response
        let (tx, rx) = tokio::sync::oneshot::channel();

        // Register the callback with the server
        server.register_callback(id, tx);

        // Send the CallbackRequest to the client
        let response = IpcResponse::CallbackRequest {
            id,
            name: parsed.name.clone(),
            args,
        };

        if let Err(e) = server.send(response) {
            self.gui_state.print_error(format!("IPC error: {}", e));
            return;
        }

        // Spawn the async task that waits for the response
        let task = IpcCallbackTask::new(id, parsed.name.clone(), rx);
        self.task_runner.spawn(task);
    }

    /// Execute a built-in command (internal helper)
    fn execute_builtin_command_internal(&mut self, cmd: &str) {
        // Calculate default size for PNG capture from viewport or window
        let default_size = self.viewport_rect
            .map(|r| (r.width().max(1.0) as u32, r.height().max(1.0) as u32))
            .unwrap_or((1024, 768));

        // Create a ViewerAdapter that wraps our state
        let mut adapter = ViewerAdapter {
            registry: &mut self.registry,
            camera: &mut self.camera,
            named_colors: &self.named_colors,
            clear_color: &mut self.clear_color,
            needs_redraw: &mut self.needs_redraw,
            selections: &mut self.selections,
            task_runner: &self.task_runner,
            // PNG capture fields
            render_context: self.render_context.as_ref(),
            settings: &self.settings,
            element_colors: &self.element_colors,
            chain_colors: &self.chain_colors,
            default_size,
        };

        // Use a fresh executor to avoid borrowing issues with self
        // The executor is stateless for single command execution
        let mut executor = CommandExecutor::new();
        
        // Execute the command using the proper CommandExecutor
        // Use do_with_options to capture output messages
        match executor.do_with_options(&mut adapter, cmd, true, false) {
            Ok(output) => {
                // Display any output messages from the command with appropriate styling
                for msg in output.messages {
                    match msg.kind {
                        MessageKind::Info => self.gui_state.print_info(msg.text),
                        MessageKind::Warning => self.gui_state.print_warning(msg.text),
                        MessageKind::Error => self.gui_state.print_error(msg.text),
                    }
                }
            }
            Err(CmdError::Aborted) => {
                // Quit/exit command was issued - signal application to close
                self.quit_requested = true;
            }
            Err(e) => {
                // Print the error to the GUI output
                self.gui_state.print_error(format!("{}", e));
            }
        }
    }

    /// Request a redraw of the window
    fn request_redraw(&mut self) {
        self.needs_redraw = true;
        if let Some(window) = &self.window {
            window.request_redraw();
        }
    }

    /// Process completed async tasks
    ///
    /// This should be called each frame to handle results from background operations
    /// like fetch, save, etc.
    fn process_async_tasks(&mut self) {
        let mut had_results = false;
        while let Some(result) = self.task_runner.poll() {
            // Each result knows how to apply itself via the TaskResult trait
            result.apply(self);
            had_results = true;
        }
        // Request window redraw if any tasks completed
        if had_results {
            self.request_redraw();
        }
    }

    /// Process IPC requests
    ///
    /// This should be called each frame to handle incoming IPC requests
    /// from external clients (e.g., pymol-python).
    fn process_ipc(&mut self) {
        // Take the server temporarily to avoid borrow conflicts
        let mut server = match self.ipc_server.take() {
            Some(s) => s,
            None => return,
        };

        // Process all pending requests
        while let Some(request) = server.poll() {
            let response = self.handle_ipc_request(&request, &mut server);
            if let Some(resp) = response {
                if let Err(e) = server.send(resp) {
                    log::error!("Failed to send IPC response: {}", e);
                }
            }
        }

        // Put the server back
        self.ipc_server = Some(server);
    }

    /// Handle a single IPC request
    fn handle_ipc_request(&mut self, request: &IpcRequest, server: &mut IpcServer) -> Option<IpcResponse> {
        match request {
            IpcRequest::Execute { id, command } => {
                log::debug!("IPC Execute: {}", command);
                // Execute the command - errors are printed to GUI output
                self.execute_command(command);
                // Always return Ok since errors are handled internally
                Some(IpcResponse::Ok { id: *id })
            }

            IpcRequest::RegisterCommand { name, help } => {
                log::info!("IPC RegisterCommand: {}", name);
                self.external_commands.register(name.clone(), help.clone());
                // Update command names for autocomplete
                self.gui_state.command_names.push(name.clone());
                self.gui_state.command_names.sort();
                None // No response needed
            }

            IpcRequest::UnregisterCommand { name } => {
                log::info!("IPC UnregisterCommand: {}", name);
                self.external_commands.unregister(name);
                // Remove from autocomplete
                self.gui_state.command_names.retain(|n| n != name);
                None // No response needed
            }

            IpcRequest::CallbackResponse { id, success, error, output } => {
                log::debug!("IPC CallbackResponse: id={}, success={}", id, success);
                // This is handled by the async task system via the server's pending callbacks
                let _ = server.handle_callback_response(request);
                
                // Also display output immediately
                for msg in output {
                    match msg.kind {
                        OutputKind::Info => self.gui_state.print_info(msg.text.clone()),
                        OutputKind::Warning => self.gui_state.print_warning(msg.text.clone()),
                        OutputKind::Error => self.gui_state.print_error(msg.text.clone()),
                    }
                }
                
                if !success {
                    if let Some(err) = error {
                        self.gui_state.print_error(err.clone());
                    }
                }
                
                self.needs_redraw = true;
                None // No response needed
            }

            IpcRequest::GetState { id } => {
                // TODO: Implement state serialization
                Some(IpcResponse::Value { 
                    id: *id, 
                    value: serde_json::json!({}) 
                })
            }

            IpcRequest::GetNames { id } => {
                let names: Vec<String> = self.registry.names()
                    .map(|s| s.to_string())
                    .collect();
                Some(IpcResponse::Value { 
                    id: *id, 
                    value: serde_json::json!(names) 
                })
            }

            IpcRequest::CountAtoms { id, selection } => {
                let mut count = 0;
                for name in self.registry.names() {
                    if let Some(mol_obj) = self.registry.get_molecule(name) {
                        if let Ok(result) = pymol_select::select(mol_obj.molecule(), selection) {
                            count += result.count();
                        }
                    }
                }
                Some(IpcResponse::Value { 
                    id: *id, 
                    value: serde_json::json!(count) 
                })
            }

            IpcRequest::Quit => {
                log::info!("IPC Quit received");
                self.quit_requested = true;
                Some(IpcResponse::Closing)
            }

            IpcRequest::Ping { id } => {
                Some(IpcResponse::Pong { id: *id })
            }

            IpcRequest::ShowWindow { id } => {
                log::info!("IPC ShowWindow received");
                self.show_window();
                self.headless = false;
                self.needs_redraw = true;
                Some(IpcResponse::Ok { id: *id })
            }

            IpcRequest::HideWindow { id } => {
                log::info!("IPC HideWindow received");
                self.hide_window();
                self.headless = true;
                Some(IpcResponse::Ok { id: *id })
            }
        }
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

        log::debug!("Evaluating {} visible selections", visible_selections.len());

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
                            log::warn!("Failed to parse selection '{}': {:?}", sel_name, e);
                            continue;
                        }
                    };

                    match pymol_select::evaluate(&expr, &ctx) {
                        Ok(result) => {
                            log::debug!(
                                "Selection '{}' matched {} atoms in '{}'",
                                sel_name,
                                result.count(),
                                mol_name
                            );
                            combined_result = match combined_result {
                                Some(existing) => Some(existing.union(&result)),
                                None => Some(result),
                            };
                        }
                        Err(e) => {
                            log::warn!("Failed to evaluate selection '{}' for '{}': {:?}", sel_name, mol_name, e);
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

        log::debug!("Visible selections evaluation returned {} molecule results", results.len());
        results
    }

    /// Handle object list actions
    fn handle_object_action(&mut self, action: ObjectAction) {
        match action {
            ObjectAction::None => {}
            ObjectAction::ToggleEnabled(name) => {
                if let Some(obj) = self.registry.get_mut(&name) {
                    let enabled = obj.is_enabled();
                    if enabled {
                        obj.disable();
                    } else {
                        obj.enable();
                    }
                    self.needs_redraw = true;
                }
            }
            ObjectAction::EnableAll => {
                let names: Vec<String> = self.registry.names().map(|s| s.to_string()).collect();
                for name in names {
                    if let Some(obj) = self.registry.get_mut(&name) {
                        obj.enable();
                    }
                }
                self.needs_redraw = true;
            }
            ObjectAction::DisableAll => {
                let names: Vec<String> = self.registry.names().map(|s| s.to_string()).collect();
                for name in names {
                    if let Some(obj) = self.registry.get_mut(&name) {
                        obj.disable();
                    }
                }
                self.needs_redraw = true;
            }
            ObjectAction::ShowAll(name) => {
                if let Some(mol) = self.registry.get_molecule_mut(&name) {
                    mol.show(RepMask::ALL.0);
                    self.needs_redraw = true;
                }
            }
            ObjectAction::HideAll(name) => {
                if let Some(mol) = self.registry.get_molecule_mut(&name) {
                    mol.hide_all();
                    self.needs_redraw = true;
                }
            }
            ObjectAction::ShowRep(name, rep) => {
                if let Some(mol) = self.registry.get_molecule_mut(&name) {
                    mol.show(rep);
                    self.needs_redraw = true;
                }
            }
            ObjectAction::HideRep(name, rep) => {
                if let Some(mol) = self.registry.get_molecule_mut(&name) {
                    mol.hide(rep);
                    self.needs_redraw = true;
                }
            }
            ObjectAction::SetColor(name, color) => {
                if let Some(color_idx) = self.named_colors.get_by_name(&color).map(|(idx, _)| idx) {
                    if let Some(mol_obj) = self.registry.get_molecule_mut(&name) {
                        for atom in mol_obj.molecule_mut().atoms_mut() {
                            atom.color = color_idx as i32;
                        }
                        self.needs_redraw = true;
                    }
                }
            }
            ObjectAction::Delete(name) => {
                self.registry.remove(&name);
                self.needs_redraw = true;
            }
            ObjectAction::ZoomTo(name) => {
                self.execute_command(&format!("zoom {}", name));
            }
            ObjectAction::CenterOn(name) => {
                self.execute_command(&format!("center {}", name));
            }
            ObjectAction::DeleteSelection(name) => {
                self.selections.remove(&name);
                self.gui_state.print_info(format!("Deleted selection \"{}\"", name));
                self.needs_redraw = true;
            }
            ObjectAction::ToggleSelectionEnabled(name) => {
                if let Some(entry) = self.selections.get_mut(&name) {
                    entry.visible = !entry.visible;
                    let state = if entry.visible { "shown" } else { "hidden" };
                    log::debug!("Selection '{}' indicators {}", name, state);
                    self.needs_redraw = true;
                }
            }
        }
    }

    /// Process accumulated input and update camera
    /// 
    /// This processes mouse deltas accumulated in InputState and applies them
    /// to the camera. Called once per frame in update().
    fn process_input(&mut self) {
        // Check if mouse is over viewport using stored rect and current mouse position
        // This avoids timing issues with gui_state.viewport_hovered which is updated after this runs
        let mouse_pos = self.input.mouse_position();
        let over_viewport = self.viewport_rect
            .map(|rect| rect.contains(egui::pos2(mouse_pos.0, mouse_pos.1)))
            .unwrap_or(true); // Default to true if no viewport rect yet (initial frames)

        if !over_viewport {
            // Still need to consume the deltas to avoid accumulation
            self.input.take_camera_deltas();
            return;
        }

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
    }

    /// Handle keyboard input
    fn handle_key(&mut self, event: KeyEvent) {
        if event.state != ElementState::Pressed {
            return;
        }

        // Don't process app key bindings when command input has focus
        // Let egui handle the keyboard events for text input
        if self.gui_state.command_has_focus {
            return;
        }

        if let PhysicalKey::Code(key) = event.physical_key {
            let binding = KeyBinding {
                key,
                ctrl: self.input.ctrl_held(),
                shift: self.input.shift_held(),
                alt: self.input.alt_held(),
            };

            if let Some(action) = self.key_bindings.get_cloned(&binding) {
                // Need to work around the borrow checker here
                // Clone the action Arc and execute
                let action = action.clone();
                action(self);
            }
        }
    }

}

impl ApplicationHandler for App {
    fn resumed(&mut self, event_loop: &ActiveEventLoop) {
        if self.window.is_some() {
            return;
        }

        // Create window (hidden if headless mode is enabled)
        let window_attrs = Window::default_attributes()
            .with_title("PyMOL-RS")
            .with_inner_size(PhysicalSize::new(1024, 768))
            .with_visible(!self.headless);

        let window = Arc::new(
            event_loop
                .create_window(window_attrs)
                .expect("Failed to create window"),
        );

        // Initialize GPU
        if let Err(e) = pollster::block_on(self.init_gpu(window.clone())) {
            log::error!("GPU initialization failed: {}", e);
            event_loop.exit();
            return;
        }

        // Load pending file if any (use command executor for consistency)
        if let Some(path) = self.pending_load_file.take() {
            // Quote the path to handle spaces and special characters
            self.execute_command(&format!("load \"{}\"", path));
        }

        window.request_redraw();
    }

    fn window_event(&mut self, event_loop: &ActiveEventLoop, _window_id: WindowId, event: WindowEvent) {
        // IMPORTANT: Always update InputState for mouse/modifier events BEFORE egui processing
        // This ensures button states are tracked correctly even when egui consumes events
        match &event {
            WindowEvent::MouseInput { state, button, .. } => {
                self.input.handle_mouse_button(*state, *button);
            }
            WindowEvent::CursorMoved { position, .. } => {
                self.input.handle_mouse_motion((position.x, position.y));
            }
            WindowEvent::MouseWheel { delta, .. } => {
                self.input.handle_scroll(*delta);
            }
            WindowEvent::ModifiersChanged(modifiers) => {
                self.input.handle_modifiers(modifiers.state());
            }
            _ => {}
        }

        // Pass events to egui
        let egui_wants_input = if let Some(egui_state) = &mut self.egui_state {
            if let Some(window) = &self.window {
                let response = egui_state.on_window_event(window, &event);
                response.consumed
            } else {
                false
            }
        } else {
            false
        };

        match event {
            WindowEvent::CloseRequested => {
                event_loop.exit();
            }

            WindowEvent::Resized(size) => {
                self.resize(size);
            }

            WindowEvent::RedrawRequested => {
                let now = Instant::now();
                let dt = (now - self.last_frame).as_secs_f32();
                self.last_frame = now;
                self.frame_count = self.frame_count.saturating_add(1);

                // Process any completed async tasks (fetch, etc.)
                self.process_async_tasks();

                // Process IPC requests from external clients
                self.process_ipc();

                // Process input and update camera (only if egui doesn't want input)
                self.process_input();

                // Update camera animation
                if self.camera.update(dt) {
                    self.needs_redraw = true;
                }

                // Always render on RedrawRequested to ensure UI is up to date
                match self.render() {
                    Ok(()) => {}
                    Err(wgpu::SurfaceError::Lost | wgpu::SurfaceError::Outdated) => {
                        if let Some(config) = &self.surface_config {
                            self.resize(PhysicalSize::new(config.width, config.height));
                        }
                    }
                    Err(wgpu::SurfaceError::OutOfMemory) => {
                        log::error!("Out of GPU memory");
                        event_loop.exit();
                    }
                    Err(e) => {
                        log::warn!("Surface error: {:?}", e);
                    }
                }
                self.needs_redraw = false;

                // Check if quit was requested by a command (quit/exit)
                if self.quit_requested {
                    event_loop.exit();
                }

                // Request continuous redraw if needed
                // Also request redraws for the first few frames to let egui layout properly
                if self.frame_count < 5
                    || self.camera.is_animating()
                    || self.input.any_button_pressed()
                {
                    self.request_redraw();
                }
            }

            WindowEvent::MouseInput { .. } => {
                // Input already handled above, just request redraw
                self.request_redraw();
            }

            WindowEvent::CursorMoved { .. } => {
                // Input already handled above
                // Always request redraw on cursor move so egui hover states update
                self.request_redraw();
                // Mark 3D scene as needing redraw only when dragging
                if self.input.any_button_pressed() {
                    self.needs_redraw = true;
                }
            }

            WindowEvent::MouseWheel { .. } => {
                // Input already handled above
                if self.gui_state.viewport_hovered && !egui_wants_input {
                    self.needs_redraw = true;
                }
                self.request_redraw();
            }

            WindowEvent::ModifiersChanged(_) => {
                // Input already handled above
            }

            WindowEvent::KeyboardInput { event, .. } => {
                // Only handle key bindings if egui didn't consume the event
                if !egui_wants_input {
                    self.handle_key(event);
                }
                self.request_redraw();
            }

            _ => {}
        }
    }

    fn about_to_wait(&mut self, event_loop: &ActiveEventLoop) {
        // Process IPC requests - this is essential for headless mode where
        // no user interaction generates RedrawRequested events.
        // We must ALWAYS call process_ipc() to accept new connections,
        // not just when a client is already connected.
        if self.ipc_server.is_some() {
            self.process_ipc();
            
            // In headless mode with IPC, use a timeout-based wake-up
            // to periodically check for connections and requests.
            // This avoids 100% CPU usage while still being responsive.
            if self.headless {
                use std::time::{Duration, Instant};
                use winit::event_loop::ControlFlow;
                
                // Wake up every 50ms to check for IPC requests
                event_loop.set_control_flow(ControlFlow::WaitUntil(
                    Instant::now() + Duration::from_millis(50)
                ));
                
                // Check for quit request (in headless mode, no window events trigger this check)
                if self.quit_requested {
                    event_loop.exit();
                }
            } else {
                // In GUI mode, request a redraw if IPC commands modified state
                // This ensures the window updates even when not focused
                if self.needs_redraw {
                    self.request_redraw();
                }
            }
        }

        // Process async tasks
        self.process_async_tasks();
    }
}
