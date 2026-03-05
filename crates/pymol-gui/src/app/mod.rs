//! Main Application
//!
//! The App struct is the main entry point that combines the Viewer, CommandExecutor,
//! and egui UI into a single application. Components communicate through a unified
//! [`MessageBus`](pymol_framework::message::MessageBus).

mod commands;
mod event_loop;
mod input;
mod ipc;
mod render;
mod ui;

use std::sync::Arc;
use std::time::Instant;

use pymol_cmd::CommandExecutor;
use pymol_render::ShadingManager;
use pymol_scene::{KeyBinding, KeyBindings, MoleculeObject, Session};
// Re-export SelectionEntry for use in UI
pub use pymol_scene::SelectionEntry;

use pymol_framework::component_store::ComponentStore;
use pymol_framework::message::MessageBus;

use crate::async_tasks::{TaskContext, TaskRunner};
use crate::ipc::{ExternalCommandRegistry, IpcServer};
use crate::layout::{Layout, pymol_classic};
use crate::model::ViewportModel;
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

/// IPC server and externally-registered command names.
pub(crate) struct IpcContext {
    pub server: Option<IpcServer>,
    pub external_commands: ExternalCommandRegistry,
}

impl IpcContext {
    fn new() -> Self {
        Self {
            server: None,
            external_commands: ExternalCommandRegistry::new(),
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

    // Input
    key_bindings: KeyBindings<KeyAction>,
    /// File currently being dragged over the window (for hover UI hint).
    pub(crate) drag_hover_path: Option<std::path::PathBuf>,

    // Background tasks
    task_runner: TaskRunner,

    // IPC
    pub(crate) ipc: IpcContext,

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
        let auto_dss = self.state.settings.get_bool(pymol_settings::id::auto_dss);
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
            key_bindings: KeyBindings::new(),
            drag_hover_path: None,
            task_runner: TaskRunner::new(),
            ipc: IpcContext::new(),
            plugin_manager: PluginManager::new(),
            headless,
            pending_load_file: None,
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

        app.ipc.server = Some(server);
        app.bus.print_info("IPC server enabled");

        Ok(app)
    }

    /// Check if IPC is enabled
    pub fn ipc_enabled(&self) -> bool {
        self.ipc.server.is_some()
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

    /// Load plugins from a directory.
    pub fn load_plugins(&mut self, dir: &std::path::Path) {
        self.plugin_manager.load_dir(
            dir,
            self.executor.registry_mut(),
            &mut self.components,
            &mut self.layout,
        );
    }
}
