//! Main Application
//!
//! The App struct is the main entry point that combines the Viewer, CommandExecutor,
//! and egui UI into a single application.

use std::sync::Arc;
use std::time::Instant;

use pymol_cmd::{CmdError, CommandExecutor, MessageKind, ParsedCommand, parse_command};
use pymol_mol::RepMask;
use pymol_render::ColorResolver;
use pymol_scene::{MoleculeObject, Object};
use pymol_scene::{expand_pick_to_selection, normalize_matrix, pick_expression_for_hit, PickHit, Picker};
use pymol_select::{build_sele_command, SelectionResult};
use pymol_scene::{CameraDelta, KeyBinding, KeyBindings, InputState, setup_uniforms};
use pymol_render::picking::Ray;
// Re-export SelectionEntry for use in UI
pub use pymol_scene::SelectionEntry;
use winit::application::ApplicationHandler;
use winit::dpi::PhysicalSize;
use winit::event::{ElementState, KeyEvent, WindowEvent};
use winit::event_loop::ActiveEventLoop;
use winit::keyboard::PhysicalKey;
use winit::window::{Window, WindowId};

use crate::async_tasks::{TaskContext, TaskRunner};
use pymol_render::ShadingManager;
use crate::ipc::{
    ExternalCommandRegistry, IpcCallbackTask, IpcRequest, IpcResponse, IpcServer,
    OutputKind as IpcOutputKind,
};
use pymol_scene::Session;
use crate::state::{CommandLineState, OutputBufferState};
use crate::view::AppView;
use crate::ui::command::CommandAction;
use crate::ui::objects::ObjectAction;
use crate::ui::sequence::{ResidueRef, SequenceAction};
use crate::ui::{
    CommandLinePanel, DragDropOverlay, NotificationOverlay, ObjectListPanel, OutputPanel,
    SequencePanel,
};
use pymol_scene::SessionAdapter;

/// Type alias for key action callbacks
pub type KeyAction = Arc<dyn Fn(&mut App) + Send + Sync>;

/// Main application state
pub struct App {
    // =========================================================================
    // Core Components
    // =========================================================================
    /// Application state (scene, camera, settings, colors)
    pub state: Session,
    /// Application view (GPU, window, egui)
    pub view: AppView,
    /// Command executor (registry of all built-in commands)
    executor: CommandExecutor,

    // =========================================================================
    // UI State (separated)
    // =========================================================================
    /// Output buffer state (log messages)
    pub output: OutputBufferState,
    /// Command line state (input, history, autocomplete)
    pub command_line: CommandLineState,

    // =========================================================================
    // UI Panels (stateful)
    // =========================================================================
    /// Object list panel (stateful for menu handling)
    object_list_panel: ObjectListPanel,
    /// Sequence viewer panel (bottom panel)
    sequence_panel: SequencePanel,

    // =========================================================================
    // Shading
    // =========================================================================
    /// Shading mode manager (Classic / Skripkin pipelines)
    shading: ShadingManager,

    // =========================================================================
    // Input State
    // =========================================================================
    /// Input handler (from pymol-scene, handles mouse with proper sensitivity)
    input: InputState,

    // =========================================================================
    // Picking / Hover
    // =========================================================================
    /// CPU ray-casting picker for hover detection
    picker: Picker,
    /// Current hover hit (atom under cursor)
    hover_hit: Option<PickHit>,
    /// Current sequence viewer hover
    sequence_hover: Option<ResidueRef>,
    /// Mouse position at left button press (for click vs drag detection)
    click_start_pos: Option<(f32, f32)>,

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
    /// File currently being dragged over the window (for hover UI hint).
    /// None when no drag is in progress.
    drag_hover_path: Option<std::path::PathBuf>,

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
    fn add_molecule(&mut self, name: &str, mut mol: pymol_mol::ObjectMolecule) {
        // Apply DSS (Define Secondary Structure) if auto_dss is enabled
        let auto_dss = self.state.settings.get_bool(pymol_settings::id::auto_dss);
        if auto_dss {
            use pymol_mol::dss::{assign_secondary_structure, DssSettings};
            let settings = DssSettings::default();
            assign_secondary_structure(&mut mol, 0, &settings);
        }

        self.state.registry.add(MoleculeObject::with_name(mol, name));
    }

    fn execute_command(&mut self, cmd: &str) {
        // Ignore the result - TaskContext doesn't propagate errors
        let _ = self.execute_command(cmd);
    }

    fn print_info(&mut self, msg: String) {
        self.output.print_info(msg);
    }

    fn print_warning(&mut self, msg: String) {
        self.output.print_warning(msg);
    }

    fn print_error(&mut self, msg: String) {
        self.output.print_error(msg);
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
        let mut app = Self {
            state: Session::new(),
            view: AppView::new(),
            executor: CommandExecutor::new(),
            output: OutputBufferState::new(),
            command_line: CommandLineState::new(),
            object_list_panel: ObjectListPanel::new(),
            sequence_panel: SequencePanel::new(),
            shading: ShadingManager::new(),
            input: InputState::new(),
            picker: Picker::new(),
            hover_hit: None,
            sequence_hover: None,
            click_start_pos: None,
            last_frame: Instant::now(),
            needs_redraw: true,
            frame_count: 0,
            key_bindings: KeyBindings::new(),
            task_runner: TaskRunner::new(),
            ipc_server: None,
            external_commands: ExternalCommandRegistry::new(),
            headless,
            pending_load_file: None,
            drag_hover_path: None,
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
        app.output.print_info("IPC server enabled".to_string());

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
        self.view.show_window();
    }

    /// Hide the window (makes it invisible)
    ///
    /// This has no effect if the window hasn't been created yet.
    ///
    /// # Platform-specific
    /// - Android / Wayland / Web: Unsupported
    pub fn hide_window(&self) {
        self.view.hide_window();
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
        self.view.is_window_visible()
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
    async fn init_gpu(&mut self, window: Arc<Window>) -> Result<(String, String), String> {
        // We need to get adapter info before init_gpu consumes it
        // For now, we'll do a quick adapter request just for info
        let adapter = self.view.instance
            .request_adapter(&wgpu::RequestAdapterOptions {
                power_preference: wgpu::PowerPreference::HighPerformance,
                compatible_surface: None,
                force_fallback_adapter: false,
            })
            .await
            .map_err(|e| format!("No suitable GPU adapter found: {}", e))?;

        let info = adapter.get_info();
        let device_name = info.name.clone();
        let backend = format!("{:?}", info.backend);

        // Drop adapter so init_gpu can create its own
        drop(adapter);

        // Initialize the view's GPU resources
        self.view.init_gpu(window.clone()).await?;

        // Set camera aspect ratio
        let size = window.inner_size();
        self.state.camera.set_aspect(size.width as f32 / size.height as f32);

        Ok((device_name, backend))
    }

    /// Handle window resize
    fn resize(&mut self, new_size: PhysicalSize<u32>) {
        if new_size.width == 0 || new_size.height == 0 {
            return;
        }

        self.view.resize(new_size);

        // Update camera aspect ratio
        self.state.camera.set_aspect(new_size.width as f32 / new_size.height as f32);
        self.needs_redraw = true;
    }

    /// Render a frame — thin orchestrator.
    fn render(&mut self) -> Result<(), wgpu::SurfaceError> {
        // Extract config values before any mutable borrows.
        let (width, height) = {
            let config = self.view.surface_config.as_ref().unwrap();
            (config.width, config.height)
        };
        let scale_factor = self.view.scale_factor();

        // Run egui UI first (mutably borrows self to process commands/actions).
        let egui_output = self.run_egui_ui();

        // Update camera aspect ratio based on viewport dimensions.
        let (viewport_width, viewport_height) = if let Some(vp) = &self.view.viewport_rect {
            (vp.width().max(1.0), vp.height().max(1.0))
        } else {
            (width as f32, height as f32)
        };
        self.state.camera.set_aspect(viewport_width / viewport_height);

        // Get surface texture.
        let surface = self.view.surface.as_ref().unwrap();
        let output = surface.get_current_texture()?;
        let output_view = output.texture.create_view(&wgpu::TextureViewDescriptor::default());

        // Upload scene uniforms and set the active shading mode.
        let shading_mode = pymol_settings::ShadingMode::from_settings(&self.state.settings);
        {
            let context = self.view.render_context.as_mut().unwrap();
            self.shading.set_mode(shading_mode, context);
            let uniforms = setup_uniforms(
                &self.state.camera,
                &self.state.settings,
                self.state.clear_color,
                (viewport_width, viewport_height),
            );
            context.update_uniforms(&uniforms);
        }

        // Prepare molecule GPU geometry; detect geometry changes.
        let names: Vec<String> = self.state.registry.names().map(|s: &str| s.to_string()).collect();
        self.prepare_molecules(&names);

        // Create command encoder.
        let mut encoder = {
            let context = self.view.render_context.as_ref().unwrap();
            context.device().create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("Render Encoder"),
            })
        };

        // Give the active pipeline the current scene bounding sphere so it can
        // compute pre-pass matrices (e.g. shadow projections). Modes that don't
        // need it ignore this call.
        if let Some((min, max)) = self.state.registry.extent() {
            let center = [
                (min.x + max.x) * 0.5,
                (min.y + max.y) * 0.5,
                (min.z + max.z) * 0.5,
            ];
            let dx = max.x - min.x;
            let dy = max.y - min.y;
            let dz = max.z - min.z;
            let radius = (dx * dx + dy * dy + dz * dz).sqrt() * 0.5;
            self.shading.set_scene_bounds(center, radius.max(1.0));
        }

        // Run the active pipeline's prepare step. Returns true when the pipeline
        // requires additional render passes this frame (e.g. shadow depth passes).
        let need_shadow_passes = {
            let context = self.view.render_context.as_mut().unwrap();
            self.shading.prepare(context, &self.state.settings)
        };

        // Execute pre-passes if required, then notify the pipeline they're done.
        if need_shadow_passes {
            self.render_shadow_passes(&mut encoder, &names);
            self.shading.finish_shadow_passes();
        }

        // 3D scene pass (opaque + transparent).
        self.render_molecules(&mut encoder, &output_view, &names);

        // Post-process: silhouette edge detection.
        self.render_silhouettes(&mut encoder, &output_view, width, height);

        // egui UI overlay.
        self.render_egui(&mut encoder, &output_view, egui_output, width, height, scale_factor);

        // Submit and present.
        let context = self.view.render_context.as_ref().unwrap();
        context.queue().submit(std::iter::once(encoder.finish()));
        output.present();

        Ok(())
    }

    /// Render shadow depth passes requested by the active shading pipeline.
    ///
    /// Uses `ShadowPassState` from the pipeline — fully mode-agnostic.
    fn render_shadow_passes(&mut self, encoder: &mut wgpu::CommandEncoder, names: &[String]) {
        let context = self.view.render_context.as_ref().unwrap();

        let state = match self.shading.shadow_pass_state() {
            Some(s) => s,
            None => return,
        };

        let bind_group = state.pipelines.create_bind_group(context.device());

        for (i, matrix) in state.matrices.iter().enumerate() {
            state.pipelines.update_uniform(context.queue(), matrix);
            let (vp_x, vp_y, vp_w, vp_h) = state.atlas.tile_viewport(i as u32);

            let mut shadow_pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("Shadow Depth Pass"),
                color_attachments: &[],
                depth_stencil_attachment: Some(wgpu::RenderPassDepthStencilAttachment {
                    view: &state.atlas.depth_view,
                    depth_ops: Some(wgpu::Operations {
                        load: if i == 0 { wgpu::LoadOp::Clear(1.0) } else { wgpu::LoadOp::Load },
                        store: wgpu::StoreOp::Store,
                    }),
                    stencil_ops: None,
                }),
                timestamp_writes: None,
                occlusion_query_set: None,
            });

            shadow_pass.set_viewport(vp_x as f32, vp_y as f32, vp_w as f32, vp_h as f32, 0.0, 1.0);
            shadow_pass.set_scissor_rect(vp_x, vp_y, vp_w, vp_h);

            for name in names {
                if let Some(mol_obj) = self.state.registry.get_molecule(name) {
                    mol_obj.render_shadow_depth(&mut shadow_pass, state.pipelines, &bind_group, context);
                }
            }
        }
    }

    /// Prepare molecule GPU geometry for the frame.
    ///
    /// Calls `prepare_render` on each molecule, updates selection indicators,
    /// and marks shadows dirty if any geometry changed.
    fn prepare_molecules(&mut self, names: &[String]) {
        let selection_results = self.evaluate_visible_selections();
        let selection_width = self.state.settings.get_float(pymol_settings::id::selection_width).max(6.0);
        let mouse_selection_mode = self.state.settings.get_int(pymol_settings::id::mouse_selection_mode);

        // Snapshot hover state to avoid borrow issues
        let hover_hit = self.hover_hit.clone();
        let sequence_hover = self.sequence_hover.clone();

        // Sequence hover uses Ctrl/Cmd for exclusion mode
        let ctrl_held = self.input.ctrl_or_cmd_held();

        const COLOR_YELLOW: [f32; 4] = [1.0, 1.0, 0.5, 0.7];
        const COLOR_RED: [f32; 4] = [1.0, 0.0, 0.0, 1.0];

        let mut geometry_changed = false;
        let context = self.view.render_context.as_ref().unwrap();
        for name in names {
            let color_resolver = ColorResolver::new(
                &self.state.named_colors,
                &self.state.element_colors,
            );
            if let Some(mol_obj) = self.state.registry.get_molecule_mut(name) {
                if mol_obj.is_dirty() {
                    geometry_changed = true;
                }
                mol_obj.prepare_render(context, color_resolver, &self.state.settings);

                // Selection indicator
                let sele_for_obj = selection_results.iter().find(|(n, _)| n == name).map(|(_, s)| s);
                if let Some(sel) = sele_for_obj {
                    log::debug!("Setting selection indicator for '{}' with {} atoms", name, sel.count());
                    mol_obj.set_selection_indicator_with_size(sel, context, Some(selection_width));
                } else {
                    mol_obj.clear_selection_indicator();
                }

                // Hover indicator: 3D viewport pick takes priority, then sequence hover.
                if let Some(ref hit) = hover_hit {
                    if hit.object_name == *name {
                        let sel = expand_pick_to_selection(hit, mouse_selection_mode, mol_obj.molecule());
                        // 3D viewport: red if hovered atoms overlap sele, yellow otherwise
                        let overlaps_sele = sele_for_obj
                            .is_some_and(|sele| sel.intersection(sele).any());
                        let color = if overlaps_sele { COLOR_RED } else { COLOR_YELLOW };
                        mol_obj.set_hover_indicator(&sel, context, selection_width, color);
                    } else {
                        mol_obj.clear_hover_indicator();
                    }
                } else if let Some(ref hover) = sequence_hover {
                    if hover.object_name == *name {
                        let mol = mol_obj.molecule();
                        let sel = SelectionResult::from_indices(
                            mol.atom_count(),
                            mol.atoms_indexed().filter_map(|(idx, atom)| {
                                if atom.residue.key.chain == hover.chain_id
                                    && atom.residue.key.resv == hover.resv
                                {
                                    Some(idx)
                                } else {
                                    None
                                }
                            }),
                        );
                        // Sequence hover: red with Ctrl (exclusion), yellow without
                        let color = if ctrl_held { COLOR_RED } else { COLOR_YELLOW };
                        mol_obj.set_hover_indicator(&sel, context, selection_width, color);
                    } else {
                        mol_obj.clear_hover_indicator();
                    }
                } else {
                    mol_obj.clear_hover_indicator();
                }
            }
        }

        // Prepare measurement objects
        for name in names {
            if let Some(meas_obj) = self.state.registry.get_measurement_mut(name) {
                meas_obj.prepare_render(context);
            }
        }

        if geometry_changed {
            self.shading.invalidate_shadows();
        }
    }

    /// 3D scene render pass — opaque + transparent molecules.
    fn render_molecules(
        &mut self,
        encoder: &mut wgpu::CommandEncoder,
        output_view: &wgpu::TextureView,
        names: &[String],
    ) {
        let context = self.view.render_context.as_ref().unwrap();
        let depth_view = self.view.depth_view.as_ref().unwrap();

        let mut render_pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
            label: Some("3D Render Pass"),
            color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                view: output_view,
                resolve_target: None,
                ops: wgpu::Operations {
                    load: wgpu::LoadOp::Clear(wgpu::Color {
                        r: self.state.clear_color[0] as f64,
                        g: self.state.clear_color[1] as f64,
                        b: self.state.clear_color[2] as f64,
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

        // Restrict 3D rendering to the central viewport area, clamped to surface bounds.
        if let Some(vp) = &self.view.viewport_rect {
            let (surf_w, surf_h) = self
                .view
                .surface_config
                .as_ref()
                .map(|c| (c.width, c.height))
                .unwrap_or((1, 1));
            let vp_x = (vp.min.x.max(0.0) as u32).min(surf_w.saturating_sub(1));
            let vp_y = (vp.min.y.max(0.0) as u32).min(surf_h.saturating_sub(1));
            let vp_w = (vp.width().max(1.0) as u32).min(surf_w - vp_x);
            let vp_h = (vp.height().max(1.0) as u32).min(surf_h - vp_y);
            render_pass.set_viewport(vp_x as f32, vp_y as f32, vp_w as f32, vp_h as f32, 0.0, 1.0);
            render_pass.set_scissor_rect(vp_x, vp_y, vp_w, vp_h);
        }

        for name in names {
            if let Some(mol_obj) = self.state.registry.get_molecule(name) {
                if mol_obj.is_enabled() {
                    mol_obj.render(&mut render_pass, context);
                }
            }
        }

        // Render measurement objects (dashed lines)
        for name in names {
            if let Some(meas_obj) = self.state.registry.get_measurement(name) {
                if meas_obj.is_enabled() {
                    meas_obj.render(&mut render_pass, context);
                }
            }
        }
    }

    /// Post-process silhouette edge detection pass.
    fn render_silhouettes(
        &mut self,
        encoder: &mut wgpu::CommandEncoder,
        output_view: &wgpu::TextureView,
        width: u32,
        height: u32,
    ) {
        if !self.state.settings.get_bool(pymol_settings::id::silhouettes) {
            return;
        }
        if let (Some(silhouette), Some(depth_view)) =
            (&self.view.silhouette_pipeline, &self.view.depth_view)
        {
            let context = self.view.render_context.as_ref().unwrap();
            let thickness = self.state.settings.get_float(pymol_settings::id::silhouette_width);
            let depth_jump = self.state.settings.get_float(pymol_settings::id::silhouette_depth_jump);
            let color_int = self.state.settings.get_color(pymol_settings::id::silhouette_color);
            let color = pymol_color::Color::from_packed_rgb(color_int).to_rgba(1.0);
            silhouette.render(
                encoder,
                context.queue(),
                context.device(),
                output_view,
                depth_view,
                width,
                height,
                thickness,
                depth_jump,
                color,
                None,
            );
        }
    }

    /// egui UI overlay pass.
    fn render_egui(
        &mut self,
        encoder: &mut wgpu::CommandEncoder,
        output_view: &wgpu::TextureView,
        egui_output: Option<(Vec<egui::ClippedPrimitive>, egui::TexturesDelta)>,
        width: u32,
        height: u32,
        scale_factor: f32,
    ) {
        let (clipped_primitives, textures_delta) = match egui_output {
            Some(o) => o,
            None => return,
        };
        let egui_renderer = match &mut self.view.egui_renderer {
            Some(r) => r,
            None => return,
        };

        let context = self.view.render_context.as_ref().unwrap();
        let device = context.device();
        let queue = context.queue();

        for (id, image_delta) in &textures_delta.set {
            egui_renderer.update_texture(device, queue, *id, image_delta);
        }

        let screen_descriptor = egui_wgpu::ScreenDescriptor {
            size_in_pixels: [width, height],
            pixels_per_point: scale_factor,
        };
        egui_renderer.update_buffers(device, queue, encoder, &clipped_primitives, &screen_descriptor);

        {
            let render_pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("egui Render Pass"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: output_view,
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
            let mut render_pass = render_pass.forget_lifetime();
            egui_renderer.render(&mut render_pass, &clipped_primitives, &screen_descriptor);
        }

        for id in &textures_delta.free {
            egui_renderer.free_texture(id);
        }
    }

    /// Run the egui UI and return the output for rendering later
    fn run_egui_ui(&mut self) -> Option<(Vec<egui::ClippedPrimitive>, egui::TexturesDelta)> {
        // Take egui input first
        let raw_input = {
            let (egui_state, window) = match (&mut self.view.egui_state, &self.view.window) {
                (Some(state), Some(window)) => (state, window),
                _ => return None,
            };
            egui_state.take_egui_input(window)
        };

        // Get scale factor for converting to physical pixels
        let scale_factor = self.view.scale_factor();

        // Collect data needed for UI without borrowing self in the closure
        let mut commands_to_execute = Vec::new();
        let mut object_actions = Vec::new();
        let mut sequence_actions = Vec::new();
        let mut viewport_rect_logical = egui::Rect::NOTHING;

        // Check for pending async tasks and get their messages before entering the closure
        let pending_messages = self.task_runner.pending_messages();

        // Get command names for autocomplete (combine built-in + external)
        let builtin_names: Vec<&str> = self.executor.registry().all_names().collect();
        let external_names: Vec<&str> = self.external_commands.names().collect();
        let all_command_names: Vec<&str> = builtin_names.into_iter().chain(external_names).collect();

        // Build completion context
        let setting_names = pymol_settings::setting_names();
        let color_names: Vec<String> = self.state.named_colors.names().iter().map(|s| s.to_string()).collect();
        let object_names: Vec<String> = self.state.registry.names().map(|s| s.to_string()).collect();
        let selection_names: Vec<String> = self.state.selections.names();

        // Update raytraced overlay texture if there's a new image
        // This must happen outside the egui closure to avoid borrow conflicts
        let ray_overlay_info = if let Some(ray_img) = &self.state.raytraced_image {
            // Check if we need to update the texture (new image or size changed)
            let needs_update = self.view.ray_overlay_size != Some((ray_img.width, ray_img.height));
            if needs_update {
                self.view.update_ray_overlay(&ray_img.data, ray_img.width, ray_img.height);
            }
            self.view.ray_overlay_texture_id.map(|id| (id, ray_img.width, ray_img.height))
        } else {
            // Clear overlay if raytraced image was removed
            if self.view.ray_overlay_texture_id.is_some() {
                self.view.clear_ray_overlay();
            }
            None
        };

        // Run egui with a closure that captures only what it needs
        let full_output = {
            let output = &mut self.output;
            let command_line = &mut self.command_line;
            let ui_config = &mut self.view.ui_config;
            let registry = &self.state.registry;
            let camera = &self.state.camera;
            let selections = &self.state.selections;
            let named_colors = &self.state.named_colors;
            let object_list_panel = &mut self.object_list_panel;
            let sequence_panel = &mut self.sequence_panel;
            let cmd_registry = self.executor.registry();
            let setting_names_refs: Vec<&str> = setting_names.iter().copied().collect();

            let completion_ctx = crate::ui::completion::CompletionContext {
                command_names: &all_command_names,
                registry: &cmd_registry,
                setting_names: &setting_names_refs,
                color_names: &color_names,
                object_names: &object_names,
                selection_names: &selection_names,
            };

            self.view.egui_ctx.run(raw_input, |ctx| {
                // Top panel - output and command line
                egui::TopBottomPanel::top("top_panel").show(ctx, |ui| {
                    if ui_config.show_output_panel {
                        OutputPanel::show(ui, output, ui_config.output_panel_height);

                        // Draggable resize handle between output and command line
                        let resize_id = ui.id().with("output_resize_handle");
                        let sep_response = ui.separator();
                        let handle_rect = sep_response.rect.expand2(egui::vec2(0.0, 2.0));
                        let response =
                            ui.interact(handle_rect, resize_id, egui::Sense::drag());
                        if response.dragged() {
                            ui_config.output_panel_height =
                                (ui_config.output_panel_height + response.drag_delta().y)
                                    .clamp(30.0, 600.0);
                        }
                        if response.hovered() || response.dragged() {
                            ui.ctx().set_cursor_icon(egui::CursorIcon::ResizeVertical);
                        }
                    }

                    match CommandLinePanel::show(ui, command_line, &completion_ctx) {
                        CommandAction::Execute(cmd) => {
                            commands_to_execute.push(cmd);
                        }
                        CommandAction::None => {}
                    }
                });

                // Bottom panel - sequence viewer (docked or floating)
                if ui_config.show_sequence_panel {
                    let has_sele = selections.contains("sele");

                    if ui_config.sequence_panel_floating {
                        // Floating mode: egui::Window
                        let mut open = true;
                        egui::Window::new("Sequence")
                            .id(egui::Id::new("sequence_panel_window"))
                            .open(&mut open)
                            .default_size(ui_config.sequence_window_size)
                            .min_size([200.0, 60.0])
                            .resizable(true)
                            .collapsible(true)
                            .show(ctx, |ui| {
                                ui.horizontal(|ui| {
                                    ui.with_layout(
                                        egui::Layout::right_to_left(egui::Align::Center),
                                        |ui| {
                                            if ui
                                                .small_button("\u{2B07}")
                                                .on_hover_text("Dock panel")
                                                .clicked()
                                            {
                                                ui_config.sequence_panel_floating = false;
                                                sequence_panel.clear_drag();
                                            }
                                        },
                                    );
                                });
                                ui.separator();
                                sequence_actions =
                                    sequence_panel.show(ui, registry, named_colors, has_sele);
                            });
                        if !open {
                            ui_config.show_sequence_panel = false;
                        }
                    } else {
                        // Docked mode: TopBottomPanel at bottom
                        let panel_height = sequence_panel.desired_height(registry);
                        egui::TopBottomPanel::bottom("sequence_panel")
                            .resizable(true)
                            .default_height(panel_height)
                            .height_range(30.0..=200.0)
                            .show(ctx, |ui| {
                                // Float button in top-right corner
                                let rect = ui.max_rect();
                                let btn_size = egui::vec2(20.0, 16.0);
                                let btn_pos = egui::pos2(
                                    rect.right() - btn_size.x - 4.0,
                                    rect.top() + 2.0,
                                );
                                let btn_rect = egui::Rect::from_min_size(btn_pos, btn_size);
                                if ui
                                    .put(
                                        btn_rect,
                                        egui::Button::new(
                                            egui::RichText::new("\u{2197}").size(12.0),
                                        )
                                        .frame(false),
                                    )
                                    .on_hover_text("Float panel")
                                    .clicked()
                                {
                                    ui_config.sequence_panel_floating = true;
                                    sequence_panel.clear_drag();
                                }

                                sequence_actions =
                                    sequence_panel.show(ui, registry, named_colors, has_sele);
                            });
                    }
                }

                // Right panel - object list only
                if ui_config.show_control_panel {
                    egui::SidePanel::right("right_panel")
                        .default_width(ui_config.right_panel_width)
                        .show(ctx, |ui| {
                            let panel_x_range = ui.max_rect().x_range();

                            // Bottom toolbar (pinned below the scrollable list)
                            egui::TopBottomPanel::bottom("right_panel_toolbar")
                                .show_separator_line(false)
                                .show_inside(ui, |ui| {
                                    // Full-width separator using the outer panel's x range
                                    let y = ui.max_rect().top();
                                    let stroke = ui.visuals().widgets.noninteractive.bg_stroke;
                                    ui.painter().hline(panel_x_range, y, stroke);
                                    ui.add_space(4.0);

                                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                                        let size = ui.spacing().interact_size.y;
                                        let btn = egui::Button::new("S")
                                            .min_size(egui::vec2(size, size))
                                            .selected(ui_config.show_sequence_panel);
                                        if ui.add(btn).on_hover_text("Sequence viewer").clicked() {
                                            ui_config.show_sequence_panel = !ui_config.show_sequence_panel;
                                        }
                                    });
                                });

                            // Scrollable object list fills remaining space
                            egui::ScrollArea::vertical().show(ui, |ui| {
                                object_actions = object_list_panel.show(ui, registry, selections);
                            });
                        });
                }

                // Central panel to capture the viewport rect (area for 3D rendering)
                let central_response = egui::CentralPanel::default()
                    .frame(egui::Frame::NONE)
                    .show(ctx, |ui| {
                        let available = ui.available_size();

                        // If we have a raytraced overlay, render it scaled to fit the viewport
                        if let Some((texture_id, img_width, img_height)) = ray_overlay_info {
                            // Calculate scaling to fit the image in the viewport while maintaining aspect ratio
                            let img_aspect = img_width as f32 / img_height as f32;
                            let vp_aspect = available.x / available.y;

                            let (display_width, display_height) = if img_aspect > vp_aspect {
                                // Image is wider than viewport - fit to width
                                (available.x, available.x / img_aspect)
                            } else {
                                // Image is taller than viewport - fit to height
                                (available.y * img_aspect, available.y)
                            };

                            // Center the image in the viewport
                            let offset_x = (available.x - display_width) / 2.0;
                            let offset_y = (available.y - display_height) / 2.0;

                            // Allocate space for interaction
                            let (rect, response) = ui.allocate_exact_size(available, egui::Sense::hover());

                            let image_rect = egui::Rect::from_min_size(
                                egui::pos2(rect.min.x + offset_x, rect.min.y + offset_y),
                                egui::vec2(display_width, display_height),
                            );

                            // Render the raytraced image
                            ui.painter().image(
                                texture_id,
                                image_rect,
                                egui::Rect::from_min_max(egui::pos2(0.0, 0.0), egui::pos2(1.0, 1.0)),
                                egui::Color32::WHITE,
                            );

                            response
                        } else {
                            // No overlay - just allocate the space for 3D rendering
                            ui.allocate_response(available, egui::Sense::hover())
                        }
                    });
                viewport_rect_logical = central_response.response.rect;

                // Paint atom labels as 2D overlay on the 3D viewport
                {
                    let vp = viewport_rect_logical;
                    let vp_tuple = (vp.min.x, vp.min.y, vp.width(), vp.height());
                    let painter = ctx.layer_painter(egui::LayerId::new(
                        egui::Order::Foreground,
                        egui::Id::new("labels_overlay"),
                    ));
                    let font = egui::FontId::new(14.0, egui::FontFamily::Proportional);
                    let label_color = egui::Color32::WHITE;

                    for name in registry.names() {
                        let mol_obj = match registry.get_molecule(name) {
                            Some(m) if m.is_enabled() => m,
                            _ => continue,
                        };

                        for (pos, text) in mol_obj.collect_labels() {
                            if let Some((sx, sy)) = camera.project_to_screen(pos, vp_tuple) {
                                if vp.contains(egui::pos2(sx, sy)) {
                                    painter.text(
                                        egui::pos2(sx, sy),
                                        egui::Align2::LEFT_BOTTOM,
                                        text,
                                        font.clone(),
                                        label_color,
                                    );
                                }
                            }
                        }
                    }

                    // Paint measurement labels
                    let meas_font = egui::FontId::new(13.0, egui::FontFamily::Proportional);
                    let meas_color = egui::Color32::from_rgb(255, 255, 0);
                    for name in registry.names() {
                        if let Some(meas_obj) = registry.get_measurement(name) {
                            if !meas_obj.is_enabled() {
                                continue;
                            }
                            for (pos, text) in meas_obj.collect_labels() {
                                if let Some((sx, sy)) = camera.project_to_screen(pos, vp_tuple) {
                                    if vp.contains(egui::pos2(sx, sy)) {
                                        painter.text(
                                            egui::pos2(sx, sy),
                                            egui::Align2::CENTER_CENTER,
                                            text,
                                            meas_font.clone(),
                                            meas_color,
                                        );
                                    }
                                }
                            }
                        }
                    }
                }

                // Show notification overlay when async tasks are in progress
                if !pending_messages.is_empty() {
                    NotificationOverlay::show(ctx, &pending_messages);
                }

                // Show drag-and-drop hint when a file is being dragged over the window
                if let Some(ref path) = self.drag_hover_path {
                    DragDropOverlay::show(ctx, path);
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
        self.view.viewport_rect = Some(viewport_rect_physical);

        // Handle platform output
        if let (Some(egui_state), Some(window)) = (&mut self.view.egui_state, &self.view.window) {
            egui_state.handle_platform_output(window, full_output.platform_output);
        }

        // Process commands and actions
        let has_actions = !commands_to_execute.is_empty()
            || !object_actions.is_empty()
            || !sequence_actions.is_empty();
        for cmd in commands_to_execute {
            // Errors are displayed in the GUI output, no need to propagate
            let _ = self.execute_command(&cmd);
        }
        for action in object_actions {
            self.handle_object_action(action);
        }
        for action in sequence_actions {
            match action {
                SequenceAction::Execute(cmd) => {
                    let _ = self.execute_command(&cmd);
                }
                SequenceAction::Notify(msg) => {
                    self.output.print_info(msg);
                }
                SequenceAction::Hover(info) => {
                    self.sequence_hover = info;
                    self.needs_redraw = true;
                }
            }
        }

        // Sync 3D selection highlights to sequence panel
        if has_actions {
            self.sync_selection_to_sequence();
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
            self.view.request_redraw();
        }

        let clipped_primitives = self.view.egui_ctx.tessellate(full_output.shapes, full_output.pixels_per_point);

        Some((clipped_primitives, full_output.textures_delta))
    }

    /// Execute a command using the CommandExecutor with SessionAdapter
    ///
    /// This method is IPC-aware - if the command is an external command
    /// registered via IPC, it will send a callback request to the client.
    ///
    /// Returns Ok(()) on success, Err(message) on failure.
    fn execute_command(&mut self, cmd: &str) -> Result<(), String> {
        let cmd = cmd.trim();
        if cmd.is_empty() || cmd.starts_with('#') {
            return Ok(());
        }

        // Echo the command to output
        self.output.print_command(format!("PyMOL> {}", cmd));

        // Parse using pymol-cmd parser for proper handling of quotes, commas, etc.
        let parsed = match parse_command(cmd) {
            Ok(p) => p,
            Err(e) => {
                let msg = format!("Parse error: {}", e);
                self.output.print_error(msg.clone());
                return Err(msg);
            }
        };

        // Handle "help <external_command>" - check if help target is an external command
        if parsed.name == "help" {
            if let Some(target_cmd) = parsed.get_str(0) {
                if let Some(help_text) = self.external_commands.help(target_cmd) {
                    self.output.print_info(help_text);
                    return Ok(());
                } else if self.external_commands.contains(target_cmd) {
                    // External command without help text
                    self.output.print_info(format!(" '{}' - external command (no help available)", target_cmd));
                    return Ok(());
                }
                // Fall through to built-in help command
            }
        }

        // Check if it's an external command (registered via IPC)
        if self.external_commands.contains(&parsed.name) {
            self.execute_external_command(&parsed);
            // External commands are async, errors are handled separately via callback
            return Ok(());
        }

        // Execute as built-in command using CommandExecutor
        self.execute_builtin_command_internal(cmd)
    }

    /// Execute an external command via IPC callback
    ///
    /// This spawns an async task that waits for the callback response.
    fn execute_external_command(&mut self, parsed: &ParsedCommand) {
        let server = match self.ipc_server.as_mut() {
            Some(s) => s,
            None => {
                self.output.print_error("No IPC server available");
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
            self.output.print_error(format!("IPC error: {}", e));
            return;
        }

        // Spawn the async task that waits for the response
        let task = IpcCallbackTask::new(id, parsed.name.clone(), rx);
        self.task_runner.spawn(task);
    }

    /// Execute a built-in command (internal helper)
    ///
    /// Returns Ok(()) on success, Err(message) on failure.
    fn execute_builtin_command_internal(&mut self, cmd: &str) -> Result<(), String> {
        // Calculate default size for PNG capture from viewport or window
        let default_size = self.view.viewport_rect
            .map(|r| (r.width().max(1.0) as u32, r.height().max(1.0) as u32))
            .unwrap_or((1024, 768));

        // Create a SessionAdapter that wraps our state
        let task_runner = &self.task_runner;
        let mut adapter = SessionAdapter {
            session: &mut self.state,
            render_context: self.view.render_context.as_ref(),
            default_size,
            needs_redraw: &mut self.needs_redraw,
            async_fetch_fn: Some(Box::new(|code: &str, name: &str, format: u8| {
                let fmt = match format {
                    1 => pymol_io::FetchFormat::Pdb,
                    0 => pymol_io::FetchFormat::Cif,
                    _ => pymol_io::FetchFormat::Bcif,
                };
                task_runner.spawn(crate::fetch::FetchTask::new(
                    code.to_string(),
                    name.to_string(),
                    fmt,
                ));
                true
            })),
        };

        // Temporarily take the executor out to avoid a borrow conflict
        // (adapter already borrows self.state and self.needs_redraw mutably).
        let mut executor = std::mem::replace(&mut self.executor, CommandExecutor::new());
        let result = executor.do_with_options(&mut adapter, cmd, true, false);
        self.executor = executor;

        match result {
            Ok(output) => {
                // Display any output messages from the command with appropriate styling
                for msg in output.messages {
                    match msg.kind {
                        MessageKind::Info => self.output.print_info(msg.text),
                        MessageKind::Warning => self.output.print_warning(msg.text),
                        MessageKind::Error => self.output.print_error(msg.text),
                    }
                }
                Ok(())
            }
            Err(CmdError::Aborted) => {
                // Quit/exit command was issued - signal application to close
                self.quit_requested = true;
                Ok(())  // Not an error to the caller
            }
            Err(e) => {
                // Print the error to the GUI output
                let msg = format!("{}", e);
                self.output.print_error(msg.clone());
                Err(msg)
            }
        }
    }

    /// Request a redraw of the window
    fn request_redraw(&mut self) {
        self.needs_redraw = true;
        self.view.request_redraw();
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
                // Execute the command and propagate errors to the client
                match self.execute_command(command) {
                    Ok(()) => Some(IpcResponse::Ok { id: *id }),
                    Err(message) => Some(IpcResponse::Error { id: *id, message }),
                }
            }

            IpcRequest::RegisterCommand { name, help } => {
                log::info!("IPC RegisterCommand: {}", name);
                self.external_commands.register(name.clone(), help.clone());
                // Note: We no longer need to update a cached list - command names are queried on demand
                None // No response needed
            }

            IpcRequest::UnregisterCommand { name } => {
                log::info!("IPC UnregisterCommand: {}", name);
                self.external_commands.unregister(name);
                // Note: We no longer need to update a cached list - command names are queried on demand
                None // No response needed
            }

            IpcRequest::CallbackResponse { id, success, error, output } => {
                log::debug!("IPC CallbackResponse: id={}, success={}", id, success);
                // This is handled by the async task system via the server's pending callbacks
                let _ = server.handle_callback_response(request);

                // Also display output immediately
                for msg in output {
                    match msg.kind {
                        IpcOutputKind::Info => self.output.print_info(msg.text.clone()),
                        IpcOutputKind::Warning => self.output.print_warning(msg.text.clone()),
                        IpcOutputKind::Error => self.output.print_error(msg.text.clone()),
                    }
                }

                if !success {
                    if let Some(err) = error {
                        self.output.print_error(err.clone());
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
                let names: Vec<String> = self.state.registry.names()
                    .map(|s| s.to_string())
                    .collect();
                Some(IpcResponse::Value {
                    id: *id,
                    value: serde_json::json!(names)
                })
            }

            IpcRequest::CountAtoms { id, selection } => {
                let mut count = 0;
                for name in self.state.registry.names() {
                    if let Some(mol_obj) = self.state.registry.get_molecule(name) {
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

            IpcRequest::GetView { id } => {
                // Get current view as 18 floats
                let view = self.state.camera.current_view();
                let r = &view.rotation;

                // Build the 18-value array:
                // [0-8]: 3x3 rotation matrix (row-major)
                // [9-11]: Camera position
                // [12-14]: Origin
                // [15]: Front clip, [16]: Back clip, [17]: FOV
                let values: Vec<f64> = vec![
                    r.data[0] as f64, r.data[1] as f64, r.data[2] as f64,
                    r.data[4] as f64, r.data[5] as f64, r.data[6] as f64,
                    r.data[8] as f64, r.data[9] as f64, r.data[10] as f64,
                    view.position.x as f64, view.position.y as f64, view.position.z as f64,
                    view.origin.x as f64, view.origin.y as f64, view.origin.z as f64,
                    view.clip_front as f64, view.clip_back as f64, view.fov as f64,
                ];

                Some(IpcResponse::Value {
                    id: *id,
                    value: serde_json::json!(values),
                })
            }
        }
    }

    /// Evaluate all visible selections for all molecules.
    ///
    /// Delegates to `SelectionManager::evaluate_visible` which builds full
    /// evaluation contexts (implicit object-name selections, pre-evaluated
    /// named selections) for each molecule.
    fn evaluate_visible_selections(&self) -> Vec<(String, SelectionResult)> {
        self.state.selections.evaluate_visible(
            &self.state.registry,
            Default::default(),
        )
    }

    /// Sync the current "sele" selection to the sequence panel highlights
    fn sync_selection_to_sequence(&mut self) {
        use std::collections::HashSet;

        let mut highlighted = HashSet::new();

        // Look for the "sele" selection (the default selection name)
        if let Some(entry) = self.state.selections.get("sele") {
            let expression = entry.expression.clone();

            for name in self.state.registry.names() {
                if let Some(mol_obj) = self.state.registry.get_molecule(name) {
                    let mol = mol_obj.molecule();

                    let result = match pymol_select::select(mol, &expression) {
                        Ok(r) => r,
                        Err(_) => continue,
                    };

                    if !result.any() {
                        continue;
                    }

                    // Map selected atoms to residues
                    for residue in mol.residues() {
                        let any_selected = residue.atom_range.clone().any(|idx| {
                            result.contains_index(idx)
                        });
                        if any_selected {
                            highlighted.insert(ResidueRef {
                                object_name: name.to_string(),
                                chain_id: residue.chain().to_string(),
                                resv: residue.resv(),
                            });
                        }
                    }
                }
            }
        }

        self.sequence_panel.update_highlights(highlighted);
    }

    /// Handle object list actions
    fn handle_object_action(&mut self, action: ObjectAction) {
        match action {
            ObjectAction::None => {}
            ObjectAction::ToggleEnabled(name) => {
                if let Some(obj) = self.state.registry.get_mut(&name) {
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
                let names: Vec<String> = self.state.registry.names().map(|s| s.to_string()).collect();
                for name in names {
                    if let Some(obj) = self.state.registry.get_mut(&name) {
                        obj.enable();
                    }
                }
                self.needs_redraw = true;
            }
            ObjectAction::DisableAll => {
                let names: Vec<String> = self.state.registry.names().map(|s| s.to_string()).collect();
                for name in names {
                    if let Some(obj) = self.state.registry.get_mut(&name) {
                        obj.disable();
                    }
                }
                self.needs_redraw = true;
            }
            ObjectAction::ShowAll(name) => {
                if let Some(mol) = self.state.registry.get_molecule_mut(&name) {
                    mol.show(RepMask::ALL);
                    self.needs_redraw = true;
                }
            }
            ObjectAction::HideAll(name) => {
                if let Some(mol) = self.state.registry.get_molecule_mut(&name) {
                    mol.hide_all();
                    self.needs_redraw = true;
                }
            }
            ObjectAction::ShowRep(name, rep) => {
                if let Some(mol) = self.state.registry.get_molecule_mut(&name) {
                    mol.show(rep);
                    self.needs_redraw = true;
                }
            }
            ObjectAction::HideRep(name, rep) => {
                if let Some(mol) = self.state.registry.get_molecule_mut(&name) {
                    mol.hide(rep);
                    self.needs_redraw = true;
                }
            }
            ObjectAction::SetColor(name, color) => {
                if let Some(color_idx) = self.state.named_colors.get_by_name(&color).map(|(idx, _)| idx) {
                    if let Some(mol_obj) = self.state.registry.get_molecule_mut(&name) {
                        for atom in mol_obj.molecule_mut().atoms_mut() {
                            atom.repr.colors.base = color_idx as i32;
                        }
                        self.needs_redraw = true;
                    }
                }
            }
            ObjectAction::Delete(name) => {
                self.state.registry.remove(&name);
                self.needs_redraw = true;
            }
            ObjectAction::ZoomTo(name) => {
                let _ = self.execute_command(&format!("zoom {}", name));
            }
            ObjectAction::CenterOn(name) => {
                let _ = self.execute_command(&format!("center {}", name));
            }
            ObjectAction::DeleteSelection(name) => {
                self.state.selections.remove(&name);
                self.output.print_info(format!("Deleted selection \"{}\"", name));
                self.needs_redraw = true;
            }
            ObjectAction::ToggleSelectionEnabled(name) => {
                if let Some(entry) = self.state.selections.get_mut(&name) {
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
        let mouse_pos = self.input.mouse_position();
        let over_viewport = self.view.is_over_viewport(mouse_pos);

        // Also skip camera updates when egui is using the pointer (e.g. dragging a
        // floating window that overlaps the viewport).
        let egui_using_pointer = self.view.egui_ctx.is_using_pointer();

        if !over_viewport || egui_using_pointer {
            // Still need to consume the deltas to avoid accumulation
            self.input.take_camera_deltas();
            return;
        }

        let deltas = self.input.take_camera_deltas();
        for delta in deltas {
            // Clear raytraced overlay on any camera change
            if self.state.raytraced_image.is_some() {
                self.state.raytraced_image = None;
            }
            match delta {
                CameraDelta::Rotate { x, y } => {
                    self.state.camera.rotate_x(x);
                    self.state.camera.rotate_y(y);
                    self.needs_redraw = true;
                }
                CameraDelta::Translate(v) => {
                    let vh = self.view.viewport_rect.map(|r| r.height()).unwrap_or(768.0);
                    let scale = self.state.camera.screen_vertex_scale(vh);
                    self.state.camera.translate(v * scale);
                    self.needs_redraw = true;
                }
                CameraDelta::Zoom(factor) => {
                    self.state.camera.zoom(factor);
                    self.needs_redraw = true;
                }
                CameraDelta::Clip { front, back } => {
                    let view = self.state.camera.view_mut();
                    view.clip_front = (view.clip_front + front).max(0.01);
                    view.clip_back = (view.clip_back + back).max(view.clip_front + 0.01);
                    self.needs_redraw = true;
                }
                CameraDelta::SlabScale(raw_delta) => {
                    let mws = self
                        .state
                        .settings
                        .get_float(pymol_settings::id::mouse_wheel_scale);
                    let scale = 1.0 + 0.04 * mws * raw_delta;

                    let view = self.state.camera.view_mut();
                    let avg = (view.clip_front + view.clip_back) * 0.5;
                    let half_width = (view.clip_back - avg).max(0.1);
                    let new_half = (half_width * scale).max(0.1);

                    view.clip_front = (avg - new_half).max(0.01);
                    view.clip_back = (avg + new_half).max(view.clip_front + 0.1);
                    self.needs_redraw = true;
                }
            }
        }
    }

    /// Ray-cast at a screen position and return the closest hit.
    ///
    /// Converts screen coordinates to a picking ray via the camera matrices,
    /// then tests against all enabled objects in the registry.
    fn pick_at(&mut self, screen_pos: (f32, f32)) -> Option<PickHit> {
        let vp = *self.view.viewport_rect.as_ref()?;
        let vp_x = screen_pos.0 - vp.min.x;
        let vp_y = screen_pos.1 - vp.min.y;

        let view_mat = self.state.camera.view_matrix();
        let mut view_inv = view_mat
            .inverse()
            .unwrap_or(lin_alg::f32::Mat4::new_identity());
        normalize_matrix(&mut view_inv);

        let ray = Ray::from_screen(
            vp_x, vp_y,
            vp.width(), vp.height(),
            &self.state.camera.projection_matrix().data,
            &view_inv.data,
        );

        self.picker.pick_ray(&ray, &self.state.registry)
    }

    /// Process hover picking — detect atom under cursor and update hover indicator.
    ///
    /// Called once per frame when the cursor moves over the viewport.
    /// Skipped while dragging (any mouse button held) to avoid interference with camera control.
    fn process_hover(&mut self) {
        if self.input.any_button_pressed() {
            return;
        }

        let mouse_pos = self.input.mouse_position();
        if !self.view.is_over_viewport(mouse_pos) {
            if self.hover_hit.is_some() {
                self.hover_hit = None;
                self.needs_redraw = true;
            }
            return;
        }

        let new_hit = self.pick_at(mouse_pos);

        let changed = match (&self.hover_hit, &new_hit) {
            (None, None) => false,
            (Some(_), None) | (None, Some(_)) => true,
            (Some(a), Some(b)) => {
                a.object_name != b.object_name || a.atom_index != b.atom_index
            }
        };

        if changed {
            self.hover_hit = new_hit;
            self.needs_redraw = true;
        }
    }

    /// Process a mouse click for atom/object selection.
    ///
    /// Clicking on unselected atoms adds them to `sele`.
    /// Clicking on already-selected atoms removes them from `sele`.
    /// Clicking on empty space clears `sele`.
    fn process_click(&mut self) {
        let mouse_pos = self.input.mouse_position();
        let hit = self.pick_at(mouse_pos);
        let mode = self.state.settings.get_int(pymol_settings::id::mouse_selection_mode);

        if let Some(ref hit) = hit {
            if let Some(mol_obj) = self.state.registry.get_molecule(&hit.object_name) {
                let mol = mol_obj.molecule();
                if let Some(expr) = pick_expression_for_hit(hit, mode, mol) {
                    // Check if hovered atoms overlap with current sele
                    let overlaps_sele = self.state.selections.get("sele").is_some_and(|entry| {
                        let sel = expand_pick_to_selection(hit, mode, mol);
                        pymol_select::select(mol, &entry.expression)
                            .map(|sele| sel.intersection(&sele).any())
                            .unwrap_or(false)
                    });

                    let has_sele = self.state.selections.contains("sele");
                    if let Some(cmd) = build_sele_command(&expr, overlaps_sele, has_sele) {
                        let _ = self.execute_command(&cmd);
                    }
                }
            }
        } else if self.state.selections.contains("sele") {
            let _ = self.execute_command("deselect sele");
        }

        // Sync selection highlights to sequence viewer immediately
        self.sync_selection_to_sequence();
    }

    /// Handle keyboard input
    fn handle_key(&mut self, event: KeyEvent) {
        if event.state != ElementState::Pressed {
            return;
        }

        // Don't process app key bindings when command input has focus
        // Let egui handle the keyboard events for text input
        if self.command_line.has_focus {
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
        if self.view.window.is_some() {
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
        match pollster::block_on(self.init_gpu(window.clone())) {
            Ok((device_name, backend)) => {
                self.output.print_info(format!("DEVICE: {}", device_name));
                self.output.print_info(format!("ENGINE: {}", backend));
            }
            Err(e) => {
                log::error!("GPU initialization failed: {}", e);
                event_loop.exit();
                return;
            }
        }

        // Load pending file if any (use command executor for consistency)
        if let Some(path) = self.pending_load_file.take() {
            // Quote the path to handle spaces and special characters
            // Errors are displayed in the GUI output, no need to propagate
            let _ = self.execute_command(&format!("load \"{}\"", path));
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
            WindowEvent::PinchGesture { delta, .. } => {
                self.input.handle_pinch_zoom(*delta);
            }
            WindowEvent::ModifiersChanged(modifiers) => {
                self.input.handle_modifiers(modifiers.state());
            }
            _ => {}
        }

        // Pass events to egui
        let egui_wants_input = if let Some(egui_state) = &mut self.view.egui_state {
            if let Some(window) = &self.view.window {
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

                // Detect atom under cursor for hover highlight
                self.process_hover();

                // Update movie (before camera animation)
                let movie_frame_changed = self.state.movie.update();
                if movie_frame_changed {
                    // Apply movie frame (view interpolation)
                    if let Some(view) = self.state.movie.interpolated_view() {
                        self.state.camera.set_view(view);
                    }
                    self.needs_redraw = true;
                }

                // Update rock animation (Y-axis oscillation)
                if self.state.movie.is_rock_enabled() {
                    // Rock parameters: ~15 degrees amplitude, smooth oscillation
                    let amplitude = 45.0_f32.to_radians();
                    let speed = 5.0; // radians per second for the sine wave phase
                    let rock_delta = self.state.movie.update_rock(dt, amplitude, speed);
                    // Apply the delta rotation directly (update_rock returns the instantaneous angle)
                    self.state.camera.rotate_y(rock_delta * dt);
                    // Clear raytraced overlay on camera change
                    self.state.raytraced_image = None;
                    self.needs_redraw = true;
                }

                // Update camera animation
                if self.state.camera.update(dt) {
                    self.needs_redraw = true;
                }

                // Always render on RedrawRequested to ensure UI is up to date
                match self.render() {
                    Ok(()) => {}
                    Err(wgpu::SurfaceError::Lost | wgpu::SurfaceError::Outdated) => {
                        if let Some(config) = &self.view.surface_config {
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
                    || self.state.camera.is_animating()
                    || self.state.movie.is_playing()
                    || self.state.movie.is_rock_enabled()
                    || self.input.any_button_pressed()
                {
                    self.request_redraw();
                }
            }

            WindowEvent::MouseInput { state, button, .. } => {
                // Click detection for picking (use viewport check, not egui_wants_input,
                // because egui reports press events as consumed even over the 3D viewport)
                if button == winit::event::MouseButton::Left {
                    let mouse_pos = self.input.mouse_position();
                    if state == ElementState::Pressed {
                        if self.view.is_over_viewport(mouse_pos) {
                            self.click_start_pos = Some(mouse_pos);
                        }
                    } else if state == ElementState::Released {
                        if let Some(start) = self.click_start_pos.take() {
                            let dx = (mouse_pos.0 - start.0).abs();
                            let dy = (mouse_pos.1 - start.1).abs();
                            if dx + dy < 5.0 {
                                self.process_click();
                            }
                        }
                    }
                }
                self.request_redraw();
            }

            WindowEvent::CursorMoved { .. } => {
                // Input already handled above
                // Always request redraw on cursor move so egui hover states update
                self.request_redraw();
                // Mark 3D scene as needing redraw when dragging or hovering over viewport
                if self.input.any_button_pressed() || self.view.is_over_viewport(self.input.mouse_position()) {
                    self.needs_redraw = true;
                }
            }

            WindowEvent::MouseWheel { .. } => {
                // Input already handled above
                // Use viewport_rect.contains() directly instead of the removed viewport_hovered field
                let mouse_pos = self.input.mouse_position();
                if self.view.is_over_viewport(mouse_pos) && !egui_wants_input {
                    self.needs_redraw = true;
                }
                self.request_redraw();
            }

            WindowEvent::PinchGesture { .. } => {
                let mouse_pos = self.input.mouse_position();
                if self.view.is_over_viewport(mouse_pos) && !egui_wants_input {
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

            // File dragged over window — show drop hint
            WindowEvent::HoveredFile(path) => {
                self.drag_hover_path = Some(path);
                self.request_redraw();
            }

            // Drag cancelled without dropping
            WindowEvent::HoveredFileCancelled => {
                self.drag_hover_path = None;
                self.request_redraw();
            }

            // File dropped — load it
            WindowEvent::DroppedFile(path) => {
                self.drag_hover_path = None;
                let path_str = path.to_string_lossy();
                log::info!("Drag-and-drop: loading \"{}\"", path_str);
                match self.execute_command(&format!("load \"{}\"", path_str)) {
                    Ok(_) => {}
                    Err(e) => self.output.print_error(format!(
                        "Failed to load \"{}\": {}",
                        path_str, e
                    )),
                }
                self.request_redraw();
            }

            _ => {}
        }
    }

    fn about_to_wait(&mut self, event_loop: &ActiveEventLoop) {
        // Process IPC requests - this is essential for both headless mode and
        // GUI mode when controlled externally (e.g., from Python).
        // We must ALWAYS call process_ipc() to accept new connections,
        // not just when a client is already connected.
        if self.ipc_server.is_some() {
            self.process_ipc();

            // When IPC is enabled, use a timeout-based wake-up to periodically
            // check for connections and requests. This is needed for:
            // - Headless mode: no user interaction generates events
            // - GUI mode: window may be unfocused, no events coming in
            // This avoids 100% CPU usage while still being responsive.
            use std::time::{Duration, Instant};
            use winit::event_loop::ControlFlow;

            // Wake up every 50ms to check for IPC requests
            event_loop.set_control_flow(ControlFlow::WaitUntil(
                Instant::now() + Duration::from_millis(50)
            ));

            // Check for quit request
            if self.quit_requested {
                event_loop.exit();
            }

            // In GUI mode, request a redraw if IPC commands modified state
            if !self.headless && self.needs_redraw {
                self.request_redraw();
            }
        }

        // Process async tasks
        self.process_async_tasks();
    }
}
