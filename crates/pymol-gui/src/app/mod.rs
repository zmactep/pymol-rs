//! Main Application
//!
//! The App struct is the main entry point that combines the Viewer, CommandExecutor,
//! and egui UI into a single application. Components communicate through a unified
//! [`MessageBus`](pymol_framework::message::MessageBus).

mod commands;
mod event_loop;
mod input;
mod render;
mod ui;

use std::sync::Arc;
use std::time::Instant;

use pymol_cmd::{CommandExecutor, DynamicCommand};
use pymol_render::ShadingManager;
use pymol_scene::{KeyBinding, KeyBindings, MoleculeObject, Session, SessionAdapter};
// Re-export SelectionEntry for use in UI
pub use pymol_scene::SelectionEntry;

use crate::component_store::ComponentStore;
use pymol_framework::message::MessageBus;

use crate::async_tasks::{TaskContext, TaskRunner};
use crate::layout::{Layout, pymol_classic};
use pymol_framework::model::ViewportModel;
use crate::plugin_manager::PluginManager;
use crate::view::AppView;

/// Type alias for key action callbacks
pub type KeyAction = Arc<dyn Fn(&mut App) + Send + Sync>;

/// Frame-level bookkeeping: timing, frame counter, quit flag.
pub(crate) struct FrameState {
    pub last_frame: Instant,
    pub frame_count: u32,
    pub quit_requested: bool,
}

impl FrameState {
    fn new() -> Self {
        Self {
            last_frame: Instant::now(),
            frame_count: 0,
            quit_requested: false,
        }
    }
}

/// Main application state
pub struct App {
    // Core scene + rendering
    pub state: Session,
    pub view: AppView,
    executor: CommandExecutor,
    shading: ShadingManager,

    // Component system
    pub(crate) components: ComponentStore,
    pub(crate) layout: Layout,
    bus: MessageBus,

    // 3D viewport (not a Component — special wgpu handling)
    pub(crate) viewport: ViewportModel,

    // Frame bookkeeping
    pub(crate) frame: FrameState,
    /// Whether the scene changed and a redraw should be scheduled.
    pub(crate) scene_dirty: bool,
    /// Whether the viewport image data changed (even if size stayed the same).
    pub(crate) image_dirty: bool,

    // Input
    key_bindings: KeyBindings<KeyAction>,
    /// File currently being dragged over the window (for hover UI hint).
    pub(crate) drag_hover_path: Option<std::path::PathBuf>,

    // Background tasks
    task_runner: TaskRunner,

    // Plugin system
    pub(crate) plugin_manager: PluginManager,

    // Init / mode
    headless: bool,
    pending_load_file: Option<String>,
}

// ============================================================================
// TaskContext implementation for App
// ============================================================================

impl TaskContext for App {
    fn add_molecule(&mut self, name: &str, mut mol: pymol_mol::ObjectMolecule) {
        // Apply DSS (Define Secondary Structure) if auto_dss is enabled
        let auto_dss = self.state.settings.behavior.auto_dss;
        if auto_dss {
            use pymol_mol::dss::{assign_secondary_structure, DssSettings};
            let settings = DssSettings::default();
            assign_secondary_structure(&mut mol, 0, &settings);
        }

        self.state.registry.add(MoleculeObject::with_name(mol, name));

        // Update effective frame count for movie navigation
        let max_states = self.state.registry.iter()
            .map(|obj| obj.n_states())
            .max()
            .unwrap_or(1);
        self.state.movie.set_n_object_states(max_states);
    }

    fn execute_command(&mut self, cmd: &str) {
        // Ignore the result - TaskContext doesn't propagate errors
        let _ = self.execute_command(cmd, false);
    }

    fn print_info(&mut self, msg: String) {
        self.bus.print_info(msg);
    }

    fn print_warning(&mut self, msg: String) {
        self.bus.print_warning(msg);
    }

    fn print_error(&mut self, msg: String) {
        self.bus.print_error(msg);
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
            shading: ShadingManager::new(),
            components: {
                let mut store = ComponentStore::new();
                for component in crate::components::default_components() {
                    store.add_boxed(component);
                }
                store
            },
            layout: pymol_classic(),
            bus: MessageBus::new(),
            viewport: ViewportModel::new(),
            frame: FrameState::new(),
            scene_dirty: true,
            image_dirty: false,
            key_bindings: KeyBindings::new(),
            drag_hover_path: None,
            task_runner: TaskRunner::new(),
            plugin_manager: PluginManager::new(),
            headless,
            pending_load_file: None,
        };

        // Set up default key bindings
        app.setup_default_key_bindings();

        app
    }

    /// Check if the application is running in headless mode
    pub fn is_headless(&self) -> bool {
        self.headless
    }

    /// Show the window (makes it visible)
    pub fn show_window(&self) {
        self.view.show_window();
    }

    /// Hide the window (makes it invisible)
    pub fn hide_window(&self) {
        self.view.hide_window();
    }

    /// Returns whether the window is currently visible
    pub fn is_window_visible(&self) -> Option<bool> {
        self.view.is_window_visible()
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

    /// Mark the scene as dirty (needs re-rendering).
    ///
    /// This schedules a window redraw via the event loop. Call this whenever
    /// state that affects rendering has changed (camera, objects, selections, etc.).
    pub(crate) fn mark_dirty(&mut self) {
        self.scene_dirty = true;
        self.view.request_redraw();
    }

    /// Set or clear the viewport image overlay and mark scene dirty.
    pub(crate) fn set_viewport_image(&mut self, image: Option<pymol_scene::ViewportImage>) {
        self.state.viewport_image = image;
        self.scene_dirty = true;
        self.image_dirty = true;
    }

    /// Clear the viewport image overlay and mark scene dirty.
    pub(crate) fn clear_viewport_image(&mut self) {
        self.state.viewport_image = None;
        self.scene_dirty = true;
    }

    /// Load plugins from a directory.
    pub fn load_plugins(&mut self, dir: &std::path::Path) {
        self.plugin_manager.load_dir(
            dir,
            &mut self.executor,
            &mut self.components,
            &mut self.layout,
        );
    }

    /// Poll plugins and process their queued operations.
    ///
    /// Three sequential phases, each borrowing non-overlapping fields:
    ///
    /// 1. **Poll** — calls `poll()` on handlers that need it, providing
    ///    read-only state and the message bus.
    /// 2. **Execute** — runs queued commands via `execute_builtin_command_internal`,
    ///    stores results for delivery in the next cycle.
    /// 3. **Register/Unregister** — processes dynamic command registrations
    ///    and unregistrations in the `CommandRegistry`.
    pub(crate) fn poll_plugins(&mut self) {
        if !self.plugin_manager.any_needs_poll() && !self.plugin_manager.has_triggered_hotkeys() {
            return;
        }

        // Phase 1: Poll all handlers
        {
            let all_names: Vec<String> = self
                .executor
                .registry()
                .all_names()
                .map(|s| s.to_string())
                .collect();
            let setting_names = pymol_settings::setting_names();
            let mut setting_names_refs: Vec<&str> = setting_names.to_vec();
            setting_names_refs.extend(self.executor.dynamic_settings().names().iter().map(String::as_str));

            let (gpu_device, gpu_queue) = self.view.gpu.render_context.as_ref()
                .map(|c| (c.device(), c.queue()))
                .unzip();

            let shared = pymol_framework::component::SharedContext {
                registry: &self.state.registry,
                camera: &self.state.camera,
                selections: &self.state.selections,
                named_colors: &self.state.named_colors,
                movie: &self.state.movie,
                settings: &self.state.settings,
                clear_color: self.state.clear_color,
                gpu_device,
                gpu_queue,
                viewport_image: self.state.viewport_image.as_ref(),
                command_names: &all_names,
                command_registry: self.executor.registry(),
                setting_names: &setting_names_refs,
                dynamic_settings: Some(self.executor.dynamic_settings()),
            };

            self.plugin_manager.poll_all(&shared, &mut self.bus);
        }

        // Phase 2: Execute queued commands, collect results
        {
            let requests = self.plugin_manager.take_pending_executions();
            if !requests.is_empty() {
                let mut results = Vec::with_capacity(requests.len());
                for req in requests {
                    let result = self.execute_builtin_command_internal(&req.command, req.silent);
                    results.push(pymol_plugin::registrar::CommandResult {
                        id: req.id,
                        result: result.map_err(|e: String| e),
                    });
                }
                self.plugin_manager.store_command_results(results);
            }
        }

        // Phase 2b: Execute queued viewer mutations
        {
            let mutations = self.plugin_manager.take_pending_mutations();
            if !mutations.is_empty() {
                let default_size = self.view.viewport_rect
                    .map(|r| (r.width().max(1.0) as u32, r.height().max(1.0) as u32))
                    .unwrap_or((1024, 768));
                let mut dirty_proxy = self.scene_dirty;
                let mut adapter = SessionAdapter {
                    session: &mut self.state,
                    render_context: self.view.gpu.render_context.as_ref(),
                    default_size,
                    needs_redraw: &mut dirty_proxy,
                    async_fetch_fn: None,
                };
                for mutation in mutations {
                    mutation(&mut adapter);
                }
                drop(adapter);
                if dirty_proxy {
                    self.scene_dirty = true;
                }
            }
        }

        // Phase 3: Process dynamic command registrations/unregistrations
        {
            let registrations = self.plugin_manager.take_pending_registrations();
            let unregistrations = self.plugin_manager.take_pending_unregistrations();

            for reg in registrations {
                let invocations = self.plugin_manager.invocations_handle();
                let dyn_cmd = DynamicCommand::new(
                    reg.name.clone(),
                    reg.help,
                    invocations,
                );
                self.executor.registry_mut().register_boxed(Box::new(dyn_cmd));
                log::info!("Registered dynamic command: {}", reg.name);
            }

            for name in unregistrations {
                if self.executor.registry_mut().unregister(&name) {
                    log::info!("Unregistered dynamic command: {}", name);
                }
            }

            // Process hotkey registration/unregistration changes
            self.plugin_manager.apply_hotkey_changes();
        }
    }
}
