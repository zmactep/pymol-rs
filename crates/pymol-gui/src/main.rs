//! PyMOL-RS GUI Application Entry Point
//!
//! Run with:
//! ```bash
//! cargo run -p pymol-gui
//! cargo run -p pymol-gui -- protein.pdb
//! cargo run -p pymol-gui -- --headless
//! ```

use std::path::PathBuf;

use clap::Parser;
use pymol_gui::App;
use pymol_gui::menu;
use winit::event_loop::{ControlFlow, EventLoop};

/// PyMOL-RS GUI Application
#[derive(Parser, Debug)]
#[command(name = "pymol-rs")]
#[command(author, version, about = "PyMOL-RS molecular visualization", long_about = None)]
struct Args {
    /// Load initial state from JSON file
    #[arg(long, value_name = "STATE_FILE")]
    state: Option<PathBuf>,

    /// Run in headless mode (no window displayed)
    #[arg(long)]
    headless: bool,

    /// Directory to load plugins from
    #[arg(long, value_name = "DIR")]
    plugin_dir: Option<PathBuf>,

    /// Files to load at startup
    #[arg(value_name = "FILE")]
    files: Vec<PathBuf>,
}

fn main() {
    // Initialize logging
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    log::info!("Starting PyMOL-RS GUI");

    // Parse command-line arguments
    let args = Args::parse();

    // Log headless mode if enabled
    if args.headless {
        log::info!("Headless mode enabled");
    }

    // When running inside a macOS .app bundle, the working directory is
    // typically "/" or the bundle location — switch to the user's home.
    if pymol_gui::bundle::bundle_contents_dir().is_some() {
        if let Some(home) = dirs::home_dir() {
            let _ = std::env::set_current_dir(&home);
        }
    }

    // Create the application
    let mut app = App::new(args.headless);

    // Load initial state if provided
    if let Some(ref state_path) = args.state {
        log::info!("Loading initial state from: {:?}", state_path);
        // TODO: Implement state loading
    }

    // Load plugins: bundled first (inside .app), then user directory
    {
        if let Some(contents) = pymol_gui::bundle::bundle_contents_dir() {
            let bundle_plugins = contents.join("PlugIns");
            if bundle_plugins.is_dir() {
                log::info!("Loading bundled plugins from {:?}", bundle_plugins);
                app.load_plugins(&bundle_plugins);
            }
        }

        let user_plugin_dir = args.plugin_dir
            .unwrap_or_else(pymol_settings::paths::plugin_dir);
        if user_plugin_dir.is_dir() {
            app.load_plugins(&user_plugin_dir);
        }
    }

    // Queue files to load
    for file_path in &args.files {
        if file_path.exists() || file_path.extension().is_some() {
            app.queue_load_file(file_path.to_string_lossy().to_string());
        }
    }

    // Create event loop (initializes NSApplication on macOS)
    let event_loop = EventLoop::new().expect("Failed to create event loop");
    event_loop.set_control_flow(ControlFlow::Wait);

    // Build native menu bar (platform init deferred to App::resumed)
    let app_menu = menu::build_menu();
    app.set_native_menu(app_menu);

    // Run the application
    if let Err(e) = event_loop.run_app(&mut app) {
        log::error!("Application error: {}", e);
    }
}
