//! Event Loop
//!
//! winit ApplicationHandler implementation — window creation, event dispatch,
//! and per-frame orchestration.

use std::sync::Arc;
use std::time::Instant;

use winit::application::ApplicationHandler;
use winit::dpi::PhysicalSize;
use winit::event::{ElementState, WindowEvent};
use winit::event_loop::ActiveEventLoop;
use winit::window::{Window, WindowId};

use super::App;

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
            let _ = self.execute_command(&format!("load \"{}\"", path), false);
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
                match self.execute_command(&format!("load \"{}\"", path_str), false) {
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
            use std::time::Duration;
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
