//! Main Application
//!
//! The App struct is the main entry point that combines the Viewer, CommandExecutor,
//! and egui UI into a single application.

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
use pymol_scene::{
    InputState, KeyBinding, KeyBindings, MoleculeObject, Picker, PickHit, Session,
};
// Re-export SelectionEntry for use in UI
pub use pymol_scene::SelectionEntry;

use crate::async_tasks::{TaskContext, TaskRunner};
use crate::ipc::{ExternalCommandRegistry, IpcServer};
use crate::state::{CommandLineState, OutputBufferState};
use crate::ui::{ObjectListPanel, SequencePanel};
use crate::ui::sequence::ResidueRef;
use crate::view::AppView;

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
        let _ = self.execute_command(cmd, false);
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

    /// Request a redraw of the window
    fn request_redraw(&mut self) {
        self.needs_redraw = true;
        self.view.request_redraw();
    }
}
