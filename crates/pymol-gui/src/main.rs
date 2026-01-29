//! PyMOL-RS GUI Application Entry Point
//!
//! Run with:
//! ```bash
//! cargo run -p pymol-gui
//! cargo run -p pymol-gui -- protein.pdb
//! cargo run -p pymol-gui -- --ipc /tmp/pymol-rs.sock
//! cargo run -p pymol-gui -- --headless --ipc /tmp/pymol-rs.sock
//! ```

use std::path::PathBuf;

use clap::Parser;
use pymol_gui::App;
use winit::event_loop::{ControlFlow, EventLoop};

/// PyMOL-RS GUI Application
#[derive(Parser, Debug)]
#[command(name = "pymol-rs")]
#[command(author, version, about = "PyMOL-RS molecular visualization", long_about = None)]
struct Args {
    /// Enable IPC server for external control (e.g., from Python)
    /// Provide the path to the Unix domain socket
    #[arg(long, value_name = "SOCKET_PATH")]
    ipc: Option<PathBuf>,

    /// Load initial state from JSON file
    #[arg(long, value_name = "STATE_FILE")]
    state: Option<PathBuf>,

    /// Run in headless mode (no window displayed)
    #[arg(long)]
    headless: bool,

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

    // Create the application
    let mut app = if let Some(ref socket_path) = args.ipc {
        log::info!("IPC mode enabled: {:?}", socket_path);
        match App::with_ipc(socket_path, args.headless) {
            Ok(app) => app,
            Err(e) => {
                log::error!("Failed to create app with IPC: {}", e);
                std::process::exit(1);
            }
        }
    } else {
        App::new(args.headless)
    };

    // Load initial state if provided
    if let Some(ref state_path) = args.state {
        log::info!("Loading initial state from: {:?}", state_path);
        // TODO: Implement state loading
    }

    // Queue files to load
    for file_path in &args.files {
        if file_path.exists() || file_path.extension().is_some() {
            app.queue_load_file(file_path.to_string_lossy().to_string());
        }
    }

    // Create event loop
    let event_loop = EventLoop::new().expect("Failed to create event loop");
    event_loop.set_control_flow(ControlFlow::Wait);

    // Run the application
    if let Err(e) = event_loop.run_app(&mut app) {
        log::error!("Application error: {}", e);
    }
}
