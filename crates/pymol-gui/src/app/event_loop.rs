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
use pymol_framework::message::AppMessage;

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
                self.bus.print_info(format!("DEVICE: {}", device_name));
                self.bus.print_info(format!("ENGINE: {}", backend));
            }
            Err(e) => {
                log::error!("GPU initialization failed: {}", e);
                event_loop.exit();
                return;
            }
        }

        // Load pending file if any (use command executor for consistency)
        if let Some(path) = self.pending_load_file.take() {
            let _ = self.execute_command(&format!("load \"{}\"", path), false);
        }

        window.request_redraw();
    }

    fn window_event(&mut self, event_loop: &ActiveEventLoop, _window_id: WindowId, event: WindowEvent) {
        // Track raw input state before egui gets it
        self.preprocess_input(&event);

        // Forward to egui — returns true if egui consumed the event
        let egui_consumed = self.forward_to_egui(&event);

        match event {
            WindowEvent::CloseRequested => {
                event_loop.exit();
                return; // skip mark_dirty
            }

            WindowEvent::Resized(size) => {
                self.resize(size);
            }

            WindowEvent::RedrawRequested => {
                self.on_redraw(event_loop);
                return; // RedrawRequested manages its own redraw scheduling
            }

            WindowEvent::MouseInput { state, button, .. } => {
                self.handle_click_detection(state, button);
            }

            WindowEvent::CursorMoved { .. } => {
                if self.viewport.input.any_button_pressed()
                    || self.view.is_over_viewport(self.viewport.input.mouse_position())
                {
                    self.scene_dirty = true;
                }
            }

            WindowEvent::MouseWheel { .. } | WindowEvent::PinchGesture { .. } => {
                let mouse_pos = self.viewport.input.mouse_position();
                if self.view.is_over_viewport(mouse_pos) && !egui_consumed {
                    self.scene_dirty = true;
                }
            }

            WindowEvent::KeyboardInput { event, .. } => {
                if !egui_consumed {
                    self.handle_key(event);
                }
            }

            WindowEvent::HoveredFile(path) => {
                self.drag_hover_path = Some(path);
            }

            WindowEvent::HoveredFileCancelled => {
                self.drag_hover_path = None;
            }

            WindowEvent::DroppedFile(path) => {
                self.handle_file_drop(path);
            }

            _ => return, // no redraw needed for unhandled events
        }

        // Single redraw request for all handled events
        self.mark_dirty();
    }

    fn about_to_wait(&mut self, event_loop: &ActiveEventLoop) {
        // Poll plugins that need periodic processing
        self.poll_plugins();

        // Timeout-based wake-up when plugins need polling
        if self.plugin_manager.any_needs_poll() {
            use std::time::Duration;
            use winit::event_loop::ControlFlow;
            event_loop.set_control_flow(ControlFlow::WaitUntil(
                Instant::now() + Duration::from_millis(50),
            ));
        }

        if self.frame.quit_requested {
            event_loop.exit();
        }

        if !self.headless && self.scene_dirty {
            self.mark_dirty();
        }

        // Process async tasks
        self.process_async_tasks();
    }
}

// =============================================================================
// Event handling helpers
// =============================================================================

impl App {
    /// Track raw mouse/modifier state before egui processing.
    fn preprocess_input(&mut self, event: &WindowEvent) {
        match event {
            WindowEvent::MouseInput { state, button, .. } => {
                self.viewport.input.handle_mouse_button(*state, *button);
            }
            WindowEvent::CursorMoved { position, .. } => {
                self.viewport.input.handle_mouse_motion((position.x, position.y));
            }
            WindowEvent::MouseWheel { delta, .. } => {
                self.viewport.input.handle_scroll(*delta);
            }
            WindowEvent::PinchGesture { delta, .. } => {
                self.viewport.input.handle_pinch_zoom(*delta);
            }
            WindowEvent::ModifiersChanged(modifiers) => {
                self.viewport.input.handle_modifiers(modifiers.state());
            }
            _ => {}
        }
    }

    /// Forward event to egui; returns true if egui consumed it.
    fn forward_to_egui(&mut self, event: &WindowEvent) -> bool {
        if let (Some(egui_state), Some(window)) = (&mut self.view.egui.state, &self.view.window) {
            egui_state.on_window_event(window, event).consumed
        } else {
            false
        }
    }

    /// Full frame update: animations, state sync, render, and repaint scheduling.
    fn on_redraw(&mut self, event_loop: &ActiveEventLoop) {
        let now = Instant::now();
        let dt = (now - self.frame.last_frame).as_secs_f32();
        self.frame.last_frame = now;
        self.frame.frame_count = self.frame.frame_count.saturating_add(1);

        // Process external events (async tasks, plugins)
        self.process_async_tasks();
        self.poll_plugins();

        // Process input and update camera
        self.process_input();
        self.process_hover();

        // Update movie and animations
        self.update_animations(dt);

        // Update selection/hover indicators before rendering
        let names: Vec<String> = self.state.registry.names().map(|s: &str| s.to_string()).collect();
        self.update_indicators(&names);

        // Render frame
        match self.render() {
            Ok(()) => {}
            Err(wgpu::SurfaceError::Lost | wgpu::SurfaceError::Outdated) => {
                if let Some(config) = &self.view.gpu.surface_config {
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
        self.scene_dirty = false;

        if self.frame.quit_requested {
            event_loop.exit();
        }

        // Schedule continuous redraws when animation is active
        if self.needs_continuous_redraw() {
            self.mark_dirty();
        }
    }

    /// Update movie playback, rock animation, and camera interpolation.
    fn update_animations(&mut self, dt: f32) {
        // Movie frame advance
        if self.state.movie.update() {
            // Sync display state of all objects to the new frame
            let current_frame = self.state.movie.current_frame();
            let state_index = self.state.movie.frame_to_state(current_frame);
            for name in self.state.registry.names().map(|s| s.to_string()).collect::<Vec<_>>() {
                if let Some(obj) = self.state.registry.get_molecule_mut(&name) {
                    obj.set_display_state(state_index);
                }
            }
            if let Some(view) = self.state.movie.interpolated_view() {
                self.state.camera.set_view(view);
            }
            self.scene_dirty = true;
        }

        // Rock animation (Y-axis oscillation)
        if self.state.movie.is_rock_enabled() {
            let amplitude = 45.0_f32.to_radians();
            let speed = 5.0;
            let rock_delta = self.state.movie.update_rock(dt, amplitude, speed);
            self.state.camera.rotate_y(rock_delta * dt);
            self.state.raytraced_image = None;
            self.scene_dirty = true;
        }

        // Camera animation (zoom/pan/rotate lerp)
        if self.state.camera.update(dt) {
            self.scene_dirty = true;
        }
    }

    /// Whether the event loop should keep requesting redraws.
    fn needs_continuous_redraw(&self) -> bool {
        self.frame.frame_count < 5
            || self.state.camera.is_animating()
            || self.state.movie.is_playing()
            || self.state.movie.is_rock_enabled()
            || self.viewport.input.any_button_pressed()
    }

    /// Detect click vs. drag for atom picking.
    fn handle_click_detection(&mut self, state: ElementState, button: winit::event::MouseButton) {
        if button != winit::event::MouseButton::Left {
            return;
        }

        let mouse_pos = self.viewport.input.mouse_position();
        if state == ElementState::Pressed {
            if self.view.is_over_viewport(mouse_pos) {
                self.viewport.click_start_pos = Some(mouse_pos);
            }
        } else if state == ElementState::Released {
            if let Some(start) = self.viewport.click_start_pos.take() {
                let dx = (mouse_pos.0 - start.0).abs();
                let dy = (mouse_pos.1 - start.1).abs();
                if dx + dy < 5.0 {
                    self.process_click();
                }
            }
        }
    }

    /// Handle a file dropped onto the window.
    fn handle_file_drop(&mut self, path: std::path::PathBuf) {
        self.drag_hover_path = None;
        let path_str = path.to_string_lossy();
        log::info!("Drag-and-drop: loading \"{}\"", path_str);
        match self.execute_command(&format!("load \"{}\"", path_str), false) {
            Ok(_) => {}
            Err(e) => self.bus.send(AppMessage::PrintError(format!(
                "Failed to load \"{}\": {}",
                path_str, e,
            ))),
        }
    }
}
