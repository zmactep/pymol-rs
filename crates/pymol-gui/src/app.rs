//! Main Application
//!
//! The App struct is the main entry point that combines the Viewer, CommandExecutor,
//! and egui UI into a single application.

use std::path::Path;
use std::sync::Arc;
use std::time::Instant;

use egui::ViewportId;
use pymol_cmd::{CommandExecutor, ViewerLike};
use pymol_color::{ChainColors, ElementColors, NamedColors};
use pymol_mol::RepMask;
use pymol_render::{ColorResolver, GlobalUniforms, RenderContext};
use pymol_scene::{MoleculeObject, Object, ObjectRegistry};
use pymol_scene::{Camera, KeyBinding, KeyBindings};
use pymol_scene::{CameraDelta, InputState};
use pymol_settings::GlobalSettings;
use winit::application::ApplicationHandler;
use winit::dpi::PhysicalSize;
use winit::event::{ElementState, KeyEvent, WindowEvent};
use winit::event_loop::ActiveEventLoop;
use winit::keyboard::PhysicalKey;
use winit::window::{Window, WindowId};

use crate::state::GuiState;
use crate::ui::command::CommandAction;
use crate::ui::objects::ObjectAction;
use crate::ui::{CommandLinePanel, ObjectListPanel, OutputPanel};

// Unused imports for future use
#[allow(unused_imports)]
use crate::ui::controls::ControlAction;
#[allow(unused_imports)]
use crate::ui::state_bar::StateBarAction;
#[allow(unused_imports)]
use crate::ui::{ControlsPanel, MouseInfoPanel, StateBar};

/// Type alias for key action callbacks
pub type KeyAction = Arc<dyn Fn(&mut App) + Send + Sync>;

// ============================================================================
// ViewerAdapter - bridges App state to the ViewerLike trait for commands
// ============================================================================

/// Adapter that wraps App's fields to implement ViewerLike for command execution
pub struct ViewerAdapter<'a> {
    /// Reference to the object registry
    pub registry: &'a mut ObjectRegistry,
    /// Reference to the camera
    pub camera: &'a mut Camera,
    /// Reference to named colors
    pub named_colors: &'a NamedColors,
    /// Reference to the background/clear color
    pub clear_color: &'a mut [f32; 3],
    /// Flag to indicate if a redraw is needed
    pub needs_redraw: &'a mut bool,
}

impl<'a> ViewerLike for ViewerAdapter<'a> {
    fn objects(&self) -> &ObjectRegistry {
        self.registry
    }

    fn objects_mut(&mut self) -> &mut ObjectRegistry {
        self.registry
    }

    fn camera(&self) -> &Camera {
        self.camera
    }

    fn camera_mut(&mut self) -> &mut Camera {
        self.camera
    }

    fn zoom_all(&mut self) {
        if let Some((min, max)) = self.registry.extent() {
            self.camera.zoom_to(min, max);
            *self.needs_redraw = true;
        }
    }

    fn zoom_on(&mut self, name: &str) {
        if let Some(obj) = self.registry.get(name) {
            if let Some((min, max)) = obj.extent() {
                self.camera.zoom_to(min, max);
                *self.needs_redraw = true;
            }
        }
    }

    fn center_all(&mut self) {
        if let Some((min, max)) = self.registry.extent() {
            self.camera.center_to(min, max);
            *self.needs_redraw = true;
        }
    }

    fn center_on(&mut self, name: &str) {
        if let Some(obj) = self.registry.get(name) {
            if let Some((min, max)) = obj.extent() {
                self.camera.center_to(min, max);
                *self.needs_redraw = true;
            }
        }
    }

    fn reset_view(&mut self) {
        *self.camera = Camera::new();
        if let Some((min, max)) = self.registry.extent() {
            self.camera.reset_view(min, max);
            *self.needs_redraw = true;
        }
    }

    fn request_redraw(&mut self) {
        *self.needs_redraw = true;
    }

    fn color_index(&self, name: &str) -> Option<u32> {
        self.named_colors.get_by_name(name).map(|(idx, _)| idx)
    }

    fn set_background_color(&mut self, r: f32, g: f32, b: f32) {
        self.clear_color[0] = r;
        self.clear_color[1] = g;
        self.clear_color[2] = b;
        *self.needs_redraw = true;
    }

    fn capture_png(
        &mut self,
        _path: &Path,
        _width: Option<u32>,
        _height: Option<u32>,
    ) -> Result<(), String> {
        // GUI doesn't support screenshot capture through this interface
        Err("Screenshot capture not yet implemented in GUI".to_string())
    }
}

/// Main application state
pub struct App {
    // =========================================================================
    // GPU Resources
    // =========================================================================
    /// wgpu instance
    instance: wgpu::Instance,
    /// GPU adapter (kept for info)
    #[allow(dead_code)]
    adapter: Option<wgpu::Adapter>,
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
}

impl Default for App {
    fn default() -> Self {
        Self::new()
    }
}

impl App {
    /// Create a new application
    pub fn new() -> Self {
        let instance = wgpu::Instance::new(wgpu::InstanceDescriptor {
            backends: wgpu::Backends::all(),
            ..Default::default()
        });

        let mut app = Self {
            instance,
            adapter: None,
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
            settings: GlobalSettings::new(),
            named_colors: NamedColors::default(),
            element_colors: ElementColors::default(),
            chain_colors: ChainColors,
            executor: CommandExecutor::new(),
            gui_state: GuiState::new(),
            input: InputState::new(),
            last_frame: Instant::now(),
            needs_redraw: true,
            frame_count: 0,
            clear_color: [0.0, 0.0, 0.0],
            key_bindings: KeyBindings::new(),
        };

        // Set up default key bindings
        app.setup_default_key_bindings();

        app
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
        self.gui_state.pending_load_file = Some(path);
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
            .ok_or_else(|| "No suitable GPU adapter found".to_string())?;

        // Log adapter info
        let info = adapter.get_info();
        self.gui_state.print_info(format!("DEVICE: {}", info.name));
        self.gui_state.print_info(format!("ENGINE:  {}", info.backend));

        // Request device
        let (device, queue) = adapter
            .request_device(
                &wgpu::DeviceDescriptor {
                    label: Some("PyMOL-RS Device"),
                    required_features: wgpu::Features::empty(),
                    required_limits: wgpu::Limits::default(),
                    memory_hints: wgpu::MemoryHints::Performance,
                },
                None,
            )
            .await
            .map_err(|e| format!("Failed to create device: {}", e))?;

        // Get surface capabilities
        let caps = surface.get_capabilities(&adapter);
        let format = caps
            .formats
            .iter()
            .find(|f| f.is_srgb())
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

        let egui_renderer = egui_wgpu::Renderer::new(&device, format, None, 1, false);

        // Create render context (takes ownership of device and queue)
        let render_context = RenderContext::new(device, queue, format);

        // Set camera aspect ratio
        self.camera.set_aspect(size.width as f32 / size.height as f32);

        // Store everything
        self.window = Some(window);
        self.surface = Some(surface);
        self.surface_config = Some(config);
        self.depth_view = Some(depth_view);
        self.adapter = Some(adapter);
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
        for name in &names {
            let color_resolver = ColorResolver::new(&self.named_colors, &self.element_colors, &self.chain_colors);
            if let Some(mol_obj) = self.registry.get_molecule_mut(name) {
                mol_obj.prepare_render(context, &color_resolver, &self.settings);
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

        // Run egui with a closure that captures only what it needs
        let full_output = {
            let gui_state = &mut self.gui_state;
            let registry = &self.registry;
            
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
                                // Object list
                                object_actions = ObjectListPanel::show(ui, registry);
                            });
                        });
                }

                // Central panel to capture the viewport rect (area for 3D rendering)
                let central_response = egui::CentralPanel::default()
                    .frame(egui::Frame::none())
                    .show(ctx, |ui| {
                        // Make the central panel interactive to detect hover
                        ui.allocate_response(ui.available_size(), egui::Sense::hover())
                    });
                viewport_rect_logical = central_response.response.rect;
                
                // Track if mouse is over the viewport
                gui_state.viewport_hovered = central_response.inner.hovered();
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
        for cmd in commands_to_execute {
            self.execute_command(&cmd);
        }
        for action in object_actions {
            self.handle_object_action(action);
        }

        let clipped_primitives = self.egui_ctx.tessellate(full_output.shapes, full_output.pixels_per_point);
        
        Some((clipped_primitives, full_output.textures_delta))
    }

    /// Draw the egui UI (kept for future use/debugging)
    #[allow(dead_code)]
    fn draw_ui(&mut self, ctx: &egui::Context) {
        let mut commands_to_execute = Vec::new();
        let mut control_actions = Vec::new();
        let mut object_actions = Vec::new();
        let mut state_bar_actions = Vec::new();

        // Top panel - output and command line
        egui::TopBottomPanel::top("top_panel").show(ctx, |ui| {
            if self.gui_state.show_output_panel {
                OutputPanel::show(ui, &mut self.gui_state);
                ui.separator();
            }

            match CommandLinePanel::show(ui, &mut self.gui_state) {
                CommandAction::Execute(cmd) => {
                    commands_to_execute.push(cmd);
                }
                CommandAction::None => {}
            }
        });

        // Bottom panel - state bar
        if self.gui_state.show_state_bar {
            egui::TopBottomPanel::bottom("bottom_panel").show(ctx, |ui| {
                state_bar_actions = StateBar::show(ui, &self.gui_state);
            });
        }

        // Right panel - controls, objects, mouse info
        if self.gui_state.show_control_panel {
            egui::SidePanel::right("right_panel")
                .default_width(self.gui_state.right_panel_width)
                .show(ctx, |ui| {
                    egui::ScrollArea::vertical().show(ui, |ui| {
                        // Controls
                        control_actions = ControlsPanel::show(
                            ui,
                            self.gui_state.is_playing,
                            self.gui_state.is_rocking,
                        );

                        ui.separator();

                        // Object list
                        object_actions = ObjectListPanel::show(ui, &self.registry);

                        ui.separator();

                        // Mouse info
                        MouseInfoPanel::show(ui, &mut self.gui_state);
                    });
                });
        }

        // Process commands
        for cmd in commands_to_execute {
            self.execute_command(&cmd);
        }

        // Process control actions
        for action in control_actions {
            self.handle_control_action(action);
        }

        // Process object actions
        for action in object_actions {
            self.handle_object_action(action);
        }

        // Process state bar actions
        for action in state_bar_actions {
            self.handle_state_bar_action(action);
        }
    }

    /// Execute a command using the CommandExecutor with ViewerAdapter
    fn execute_command(&mut self, cmd: &str) {
        let cmd = cmd.trim();
        if cmd.is_empty() {
            return;
        }

        // Create a ViewerAdapter that wraps our state
        let mut adapter = ViewerAdapter {
            registry: &mut self.registry,
            camera: &mut self.camera,
            named_colors: &self.named_colors,
            clear_color: &mut self.clear_color,
            needs_redraw: &mut self.needs_redraw,
        };

        // Use a fresh executor to avoid borrowing issues with self
        // The executor is stateless for single command execution
        let mut executor = CommandExecutor::new();
        
        // Execute the command using the proper CommandExecutor
        match executor.do_(&mut adapter, cmd) {
            Ok(()) => {
                // Command succeeded
            }
            Err(e) => {
                // Print the error to the GUI output
                self.gui_state.print_error(format!("{}", e));
            }
        }
    }

    /// Load a molecular file
    fn load_file(&mut self, path: &str) {
        match pymol_io::read_file(std::path::Path::new(path)) {
            Ok(mol) => {
                let name = mol.name.clone();
                self.gui_state.print_info(format!("CmdLoad: \"{}\" loaded as \"{}\"", path, name));
                let obj = MoleculeObject::new(mol);
                self.registry.add(obj);
                self.zoom_all();
                self.needs_redraw = true;
            }
            Err(e) => {
                self.gui_state.print_error(format!("Failed to load {}: {}", path, e));
            }
        }
    }

    /// Zoom to fit all objects
    pub fn zoom_all(&mut self) {
        if let Some((min, max)) = self.registry.extent() {
            self.camera.zoom_to(min, max);
            self.needs_redraw = true;
        }
    }

    /// Center on all objects
    pub fn center_all(&mut self) {
        if let Some((min, max)) = self.registry.extent() {
            self.camera.center_to(min, max);
            self.needs_redraw = true;
        }
    }

    /// Reset view
    pub fn reset_view(&mut self) {
        self.camera = Camera::new();
        if let Some((min, max)) = self.registry.extent() {
            self.camera.reset_view(min, max);
            self.needs_redraw = true;
        }
    }

    /// Toggle a representation (used by key bindings)
    #[allow(dead_code)]
    fn toggle_representation(&mut self, rep: u32) {
        let names: Vec<_> = self.registry.names().map(|s| s.to_string()).collect();
        for name in names {
            if let Some(mol) = self.registry.get_molecule_mut(&name) {
                mol.toggle(rep);
            }
        }
        self.needs_redraw = true;
    }

    /// Hide all representations (used by key bindings)
    #[allow(dead_code)]
    fn hide_all_representations(&mut self) {
        let names: Vec<_> = self.registry.names().map(|s| s.to_string()).collect();
        for name in names {
            if let Some(mol) = self.registry.get_molecule_mut(&name) {
                mol.hide_all();
            }
        }
        self.needs_redraw = true;
    }

    /// Show default representations (used by key bindings)
    #[allow(dead_code)]
    fn show_default_representations(&mut self) {
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

    /// Increase surface quality (used by key bindings)
    #[allow(dead_code)]
    fn increase_surface_quality(&mut self) {
        let names: Vec<_> = self.registry.names().map(|s| s.to_string()).collect();
        for name in names {
            if let Some(mol) = self.registry.get_molecule_mut(&name) {
                let q = mol.surface_quality();
                mol.set_surface_quality((q + 1).min(4));
            }
        }
        self.needs_redraw = true;
    }

    /// Decrease surface quality (used by key bindings)
    #[allow(dead_code)]
    fn decrease_surface_quality(&mut self) {
        let names: Vec<_> = self.registry.names().map(|s| s.to_string()).collect();
        for name in names {
            if let Some(mol) = self.registry.get_molecule_mut(&name) {
                let q = mol.surface_quality();
                mol.set_surface_quality((q - 1).max(-4));
            }
        }
        self.needs_redraw = true;
    }

    /// Set background color by name (legacy method, now use bg_color command)
    #[allow(dead_code)]
    fn set_bg_color(&mut self, name: &str) {
        let color = match name.to_lowercase().as_str() {
            "black" => [0.0, 0.0, 0.0],
            "white" => [1.0, 1.0, 1.0],
            "gray" | "grey" => [0.5, 0.5, 0.5],
            "red" => [1.0, 0.0, 0.0],
            "green" => [0.0, 1.0, 0.0],
            "blue" => [0.0, 0.0, 1.0],
            _ => {
                self.gui_state.print_error(format!("Unknown color: {}", name));
                return;
            }
        };
        self.clear_color = color;
        self.gui_state.print_info(format!("Background set to {}", name));
        self.needs_redraw = true;
    }

    /// Show representation by name (legacy method, now use show command)
    #[allow(dead_code)]
    fn show_representation_by_name(&mut self, name: &str) {
        let rep = match name.to_lowercase().as_str() {
            "lines" => RepMask::LINES,
            "sticks" => RepMask::STICKS,
            "spheres" => RepMask::SPHERES,
            "cartoon" => RepMask::CARTOON,
            "surface" => RepMask::SURFACE,
            "mesh" => RepMask::MESH,
            "dots" => RepMask::DOTS,
            "ribbon" => RepMask::RIBBON,
            _ => {
                self.gui_state.print_error(format!("Unknown representation: {}", name));
                return;
            }
        };
        let names: Vec<_> = self.registry.names().map(|s| s.to_string()).collect();
        for name in names {
            if let Some(mol) = self.registry.get_molecule_mut(&name) {
                mol.show(rep);
            }
        }
        self.gui_state.print_info(format!("Showing {}", name));
        self.needs_redraw = true;
    }

    /// Hide representation by name (legacy method, now use hide command)
    #[allow(dead_code)]
    fn hide_representation_by_name(&mut self, name: &str) {
        let rep = match name.to_lowercase().as_str() {
            "lines" => RepMask::LINES,
            "sticks" => RepMask::STICKS,
            "spheres" => RepMask::SPHERES,
            "cartoon" => RepMask::CARTOON,
            "surface" => RepMask::SURFACE,
            "mesh" => RepMask::MESH,
            "dots" => RepMask::DOTS,
            "ribbon" => RepMask::RIBBON,
            "everything" | "all" => {
                self.hide_all_representations();
                self.gui_state.print_info("Hiding all representations");
                return;
            }
            _ => {
                self.gui_state.print_error(format!("Unknown representation: {}", name));
                return;
            }
        };
        let names: Vec<_> = self.registry.names().map(|s| s.to_string()).collect();
        for name in names {
            if let Some(mol) = self.registry.get_molecule_mut(&name) {
                mol.hide(rep);
            }
        }
        self.gui_state.print_info(format!("Hiding {}", name));
        self.needs_redraw = true;
    }

    /// Color all objects (legacy method, now use color command)
    #[allow(dead_code)]
    fn color_all(&mut self, color: &str) {
        if let Some(color_idx) = self.named_colors.get_by_name(color).map(|(idx, _)| idx) {
            let names: Vec<String> = self.registry.names().map(|s: &str| s.to_string()).collect();
            for name in names {
                if let Some(mol_obj) = self.registry.get_molecule_mut(&name) {
                    for atom in mol_obj.molecule_mut().atoms_mut() {
                        atom.color = color_idx as i32;
                    }
                }
            }
            self.gui_state.print_info(format!("Colored all objects {}", color));
            self.needs_redraw = true;
        } else {
            self.gui_state.print_error(format!("Unknown color: {}", color));
        }
    }

    /// Print help (legacy method, now use help command)
    #[allow(dead_code)]
    fn print_help(&mut self) {
        self.gui_state.print_info("Available commands:");
        self.gui_state.print("  load <file>     - Load a molecule file");
        self.gui_state.print("  show <rep>      - Show representation (lines, sticks, cartoon, etc.)");
        self.gui_state.print("  hide <rep>      - Hide representation");
        self.gui_state.print("  color <color>   - Color all objects");
        self.gui_state.print("  bg_color <col>  - Set background color");
        self.gui_state.print("  zoom            - Zoom to fit all objects");
        self.gui_state.print("  center          - Center on all objects");
        self.gui_state.print("  reset           - Reset view");
        self.gui_state.print("  help            - Show this help");
    }

    /// Handle control panel actions (for future use)
    #[allow(dead_code)]
    fn handle_control_action(&mut self, action: ControlAction) {
        match action {
            ControlAction::None => {}
            ControlAction::Reset => self.reset_view(),
            ControlAction::Zoom => self.zoom_all(),
            ControlAction::Orient => {
                // TODO: implement orient
                self.gui_state.print_info("Orient not yet implemented");
            }
            ControlAction::ToggleOrtho => {
                self.camera.toggle_projection();
                self.needs_redraw = true;
            }
            ControlAction::Screenshot => {
                self.gui_state.print_info("Screenshot not yet implemented");
            }
            ControlAction::Ray => {
                self.gui_state.print_info("Ray tracing not yet implemented");
            }
            ControlAction::Unpick => {
                self.gui_state.print_info("Unpick not yet implemented");
            }
            ControlAction::Deselect => {
                self.gui_state.print_info("Deselect not yet implemented");
            }
            ControlAction::Rock => {
                self.gui_state.is_rocking = !self.gui_state.is_rocking;
                self.gui_state.print_info(if self.gui_state.is_rocking {
                    "Rocking enabled"
                } else {
                    "Rocking disabled"
                });
            }
            ControlAction::GetView => {
                let view = self.camera.current_view();
                self.gui_state.print_info(format!("View: {:?}", view));
            }
            ControlAction::FirstState => {
                self.gui_state.current_state = 1;
            }
            ControlAction::PreviousState => {
                if self.gui_state.current_state > 1 {
                    self.gui_state.current_state -= 1;
                }
            }
            ControlAction::Stop => {
                self.gui_state.is_playing = false;
            }
            ControlAction::Play => {
                self.gui_state.is_playing = !self.gui_state.is_playing;
            }
            ControlAction::NextState => {
                if self.gui_state.current_state < self.gui_state.total_states {
                    self.gui_state.current_state += 1;
                }
            }
            ControlAction::LastState => {
                self.gui_state.current_state = self.gui_state.total_states;
            }
            ControlAction::MovieClear => {
                self.gui_state.print_info("Movie clear not yet implemented");
            }
            ControlAction::Builder => {
                self.gui_state.print_info("Builder not yet implemented");
            }
            ControlAction::Properties => {
                self.gui_state.print_info("Properties not yet implemented");
            }
            ControlAction::Rebuild => {
                self.needs_redraw = true;
                self.gui_state.print_info("Rebuilding representations...");
            }
        }
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
                if let Some(obj) = self.registry.get(&name) {
                    if let Some((min, max)) = obj.extent() {
                        self.camera.zoom_to(min, max);
                        self.needs_redraw = true;
                    }
                }
            }
            ObjectAction::CenterOn(name) => {
                if let Some(obj) = self.registry.get(&name) {
                    if let Some((min, max)) = obj.extent() {
                        self.camera.center_to(min, max);
                        self.needs_redraw = true;
                    }
                }
            }
        }
    }

    /// Handle state bar actions (for future use)
    #[allow(dead_code)]
    fn handle_state_bar_action(&mut self, action: StateBarAction) {
        match action {
            StateBarAction::None => {}
            StateBarAction::FirstState => {
                self.gui_state.current_state = 1;
            }
            StateBarAction::PreviousState => {
                if self.gui_state.current_state > 1 {
                    self.gui_state.current_state -= 1;
                }
            }
            StateBarAction::Stop => {
                self.gui_state.is_playing = false;
            }
            StateBarAction::TogglePlay => {
                self.gui_state.is_playing = !self.gui_state.is_playing;
            }
            StateBarAction::NextState => {
                if self.gui_state.current_state < self.gui_state.total_states {
                    self.gui_state.current_state += 1;
                }
            }
            StateBarAction::LastState => {
                self.gui_state.current_state = self.gui_state.total_states;
            }
            StateBarAction::SceneCycle => {
                self.gui_state.print_info("Scene cycle not yet implemented");
            }
            StateBarAction::Fullscreen => {
                if let Some(window) = &self.window {
                    let is_fullscreen = window.fullscreen().is_some();
                    if is_fullscreen {
                        window.set_fullscreen(None);
                    } else {
                        window.set_fullscreen(Some(winit::window::Fullscreen::Borderless(None)));
                    }
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

    /// Request a redraw
    fn request_redraw(&self) {
        if let Some(window) = &self.window {
            window.request_redraw();
        }
    }
}

impl ApplicationHandler for App {
    fn resumed(&mut self, event_loop: &ActiveEventLoop) {
        if self.window.is_some() {
            return;
        }

        // Create window
        let window_attrs = Window::default_attributes()
            .with_title("PyMOL-RS")
            .with_inner_size(PhysicalSize::new(1280, 800));

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

        // Load pending file if any
        if let Some(path) = self.gui_state.pending_load_file.take() {
            self.load_file(&path);
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

                // Request continuous redraw if needed
                // Also request redraws for the first few frames to let egui layout properly
                if self.frame_count < 5
                    || self.camera.is_animating()
                    || self.input.any_button_pressed()
                    || self.gui_state.is_playing
                    || self.gui_state.is_rocking
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
                // Request redraw when dragging - viewport check is done in process_input()
                if self.input.any_button_pressed() {
                    self.needs_redraw = true;
                    self.request_redraw();
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
}
