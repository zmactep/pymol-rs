//! PyMOL-RS GUI Application Entry Point
//!
//! Run with:
//! ```bash
//! cargo run -p pymol-gui
//! cargo run -p pymol-gui -- protein.pdb
//! ```

use std::path::Path;

use pymol_gui::App;
use winit::event_loop::{ControlFlow, EventLoop};

fn main() {
    // Initialize logging
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    log::info!("Starting PyMOL-RS GUI");

    // Create the application
    let mut app = App::new();

    // Check for file argument
    let args: Vec<String> = std::env::args().collect();
    let file_arg = args.iter().skip(1).find(|arg| {
        let path = Path::new(arg);
        path.extension().is_some() || path.exists()
    });

    if let Some(file_path) = file_arg {
        // Queue file to load after initialization
        app.queue_load_file(file_path.clone());
    }

    // Create event loop
    let event_loop = EventLoop::new().expect("Failed to create event loop");
    event_loop.set_control_flow(ControlFlow::Wait);

    // Run the application
    if let Err(e) = event_loop.run_app(&mut app) {
        log::error!("Application error: {}", e);
    }
}
