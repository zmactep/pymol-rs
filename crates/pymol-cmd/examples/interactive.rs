//! Interactive PyMOL command viewer
//!
//! This example provides both a 3D molecular visualization window and a
//! command-line interface for testing PyMOL commands interactively.
//!
//! ## Usage
//!
//! ```bash
//! # Run with no file (shows help)
//! cargo run --example interactive
//!
//! # Run with a PDB file
//! cargo run --example interactive -- protein.pdb
//! ```
//!
//! ## Commands
//!
//! Type PyMOL commands at the `PyMOL> ` prompt:
//!
//! - `load protein.pdb` - Load a molecule file
//! - `show cartoon` - Show cartoon representation
//! - `hide lines` - Hide lines representation
//! - `color green, chain A` - Color chain A green
//! - `zoom` - Zoom to fit all objects
//! - `bg_color white` - Set background color
//! - `orient` - Orient to principal axes
//! - `quit` or `exit` - Exit the viewer
//!
//! ## Mouse Controls
//!
//! - Left drag: Rotate
//! - Middle drag: Pan
//! - Right drag / Scroll: Zoom
//!
//! ## Keyboard Shortcuts
//!
//! - 1-8: Toggle representations (lines, sticks, spheres, cartoon, surface, mesh, dots, ribbon)
//! - R: Reset view
//! - O: Toggle orthographic/perspective
//! - H: Hide all representations
//! - A: Show default (lines + sticks)
//!
//! ## Command History
//!
//! - Up arrow: Previous command
//! - Down arrow: Next command

use std::path::Path;
use std::thread;

use pymol_cmd::CommandExecutor;
use rustyline::error::ReadlineError;
use rustyline::DefaultEditor;
use pymol_mol::RepMask;
use pymol_scene::{KeyBinding, KeyCode, Viewer};
use winit::application::ApplicationHandler;
use winit::event::WindowEvent;
use winit::event_loop::{ActiveEventLoop, ControlFlow, EventLoop, EventLoopProxy};
use winit::window::WindowId;

/// User events sent from the stdin reader thread
#[derive(Debug, Clone)]
enum UserEvent {
    /// A command string to execute
    Command(String),
    /// Quit the application
    Quit,
}

/// Interactive viewer that combines a Viewer with a CommandExecutor
struct InteractiveViewer {
    /// The molecular viewer
    viewer: Viewer,
    /// The command executor
    executor: CommandExecutor,
}

impl InteractiveViewer {
    /// Create a new interactive viewer
    fn new() -> Self {
        let mut viewer = Viewer::new();

        // Register default key bindings
        // Camera controls
        viewer.bind_key(KeyCode::KeyR, |v| v.reset_view());
        viewer.bind_key(KeyCode::KeyO, |v| {
            v.camera_mut().toggle_projection();
            v.request_redraw();
        });

        // Representation hotkeys (number keys)
        viewer.bind_key(KeyCode::Digit1, |v| v.toggle_representation(RepMask::LINES));
        viewer.bind_key(KeyCode::Digit2, |v| v.toggle_representation(RepMask::STICKS));
        viewer.bind_key(KeyCode::Digit3, |v| v.toggle_representation(RepMask::SPHERES));
        viewer.bind_key(KeyCode::Digit4, |v| v.toggle_representation(RepMask::CARTOON));
        viewer.bind_key(KeyCode::Digit5, |v| v.toggle_representation(RepMask::SURFACE));
        viewer.bind_key(KeyCode::Digit6, |v| v.toggle_representation(RepMask::MESH));
        viewer.bind_key(KeyCode::Digit7, |v| v.toggle_representation(RepMask::DOTS));
        viewer.bind_key(KeyCode::Digit8, |v| v.toggle_representation(RepMask::RIBBON));

        // Representation management
        viewer.bind_key(KeyCode::KeyH, |v| v.hide_all_representations());
        viewer.bind_key(KeyCode::KeyA, |v| v.show_default_representations());

        // Surface quality controls
        viewer.bind_key(KeyCode::Equal, |v| v.increase_surface_quality());
        viewer.bind_key(KeyCode::Minus, |v| v.decrease_surface_quality());

        // Ctrl+R: Reset view (alternative binding)
        viewer.bind_key(KeyBinding::new(KeyCode::KeyR).ctrl(), |v| {
            log::info!("Ctrl+R: Reset view");
            v.reset_view();
        });

        Self {
            viewer,
            executor: CommandExecutor::new(),
        }
    }

    /// Execute a command and handle errors
    fn execute_command(&mut self, cmd: &str) {
        let cmd = cmd.trim();
        if cmd.is_empty() {
            return;
        }

        // Execute the command
        match self.executor.do_(&mut self.viewer, cmd) {
            Ok(()) => {
                log::debug!("Command executed: {}", cmd);
            }
            Err(e) => {
                eprintln!("Error: {}", e);
            }
        }

        // Request redraw after command execution
        self.viewer.request_redraw();
    }
}

impl ApplicationHandler<UserEvent> for InteractiveViewer {
    fn resumed(&mut self, event_loop: &ActiveEventLoop) {
        // Forward to Viewer's ApplicationHandler implementation
        <Viewer as ApplicationHandler>::resumed(&mut self.viewer, event_loop);
    }

    fn window_event(
        &mut self,
        event_loop: &ActiveEventLoop,
        window_id: WindowId,
        event: WindowEvent,
    ) {
        // Forward to Viewer's ApplicationHandler implementation
        <Viewer as ApplicationHandler>::window_event(&mut self.viewer, event_loop, window_id, event);
    }

    fn user_event(&mut self, event_loop: &ActiveEventLoop, event: UserEvent) {
        match event {
            UserEvent::Command(cmd) => {
                self.execute_command(&cmd);
            }
            UserEvent::Quit => {
                log::info!("Quit requested");
                event_loop.exit();
            }
        }
    }
}

/// Spawn a thread that reads commands from stdin and sends them via the event loop proxy
fn spawn_stdin_reader(proxy: EventLoopProxy<UserEvent>) {
    thread::spawn(move || {
        let mut rl = DefaultEditor::new().expect("Failed to create line editor");

        loop {
            match rl.readline("PyMOL> ") {
                Ok(line) => {
                    let cmd = line.trim();

                    // Handle quit commands specially
                    if cmd == "quit" || cmd == "exit" {
                        let _ = proxy.send_event(UserEvent::Quit);
                        break;
                    }

                    // Add non-empty commands to history
                    if !cmd.is_empty() {
                        let _ = rl.add_history_entry(cmd);
                    }

                    // Send command to main thread
                    if proxy.send_event(UserEvent::Command(cmd.to_string())).is_err() {
                        // Event loop has exited
                        break;
                    }
                }
                Err(ReadlineError::Interrupted | ReadlineError::Eof) => {
                    // Ctrl+C or Ctrl+D
                    log::info!("EOF/interrupt received, exiting...");
                    let _ = proxy.send_event(UserEvent::Quit);
                    break;
                }
                Err(e) => {
                    eprintln!("Error reading input: {}", e);
                    break;
                }
            }
        }
    });
}

fn main() {
    // Initialize logging
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    // Create the interactive viewer
    let mut app = InteractiveViewer::new();

    // Check for file argument
    let args: Vec<String> = std::env::args().collect();
    let file_arg = args.iter().skip(1).find(|arg| {
        let path = Path::new(arg);
        path.extension().is_some() || path.exists()
    });

    if let Some(file_path) = file_arg {
        // Load molecule from file via command executor
        let load_cmd = format!("load {}", file_path);
        app.execute_command(&load_cmd);
    }

    // Set a nice dark blue background
    app.viewer.set_background_color(0.0, 0.0, 0.1);

    // Center the view on all objects
    app.viewer.center_all();

    // Create event loop with user events
    let event_loop = EventLoop::<UserEvent>::with_user_event()
        .build()
        .expect("Failed to create event loop");

    // Get proxy for sending events from stdin thread
    let proxy = event_loop.create_proxy();

    // Spawn stdin reader thread
    spawn_stdin_reader(proxy);

    // Print help
    println!();
    println!("=== PyMOL-RS Interactive Command Viewer ===");
    println!();
    println!("Type PyMOL commands at the prompt. Examples:");
    println!("  load protein.pdb    - Load a molecule file");
    println!("  show cartoon        - Show cartoon representation");
    println!("  hide lines          - Hide lines representation");
    println!("  color green, chain A - Color chain A green");
    println!("  zoom                - Zoom to fit all objects");
    println!("  bg_color white      - Set background color");
    println!("  quit / exit         - Exit the viewer");
    println!();
    println!("Mouse: Left=Rotate, Middle=Pan, Right/Scroll=Zoom");
    println!("Keys: 1-8=Reps, R=Reset, O=Ortho, H=Hide, A=Default");
    println!("History: Up/Down arrows to navigate previous commands");
    println!();

    // Run the event loop
    event_loop.set_control_flow(ControlFlow::Wait);
    if let Err(e) = event_loop.run_app(&mut app) {
        log::error!("Event loop error: {}", e);
    }
}
