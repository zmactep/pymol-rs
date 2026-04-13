//! Event Loop
//!
//! winit ApplicationHandler implementation — window creation, event dispatch,
//! and per-frame orchestration.

use std::sync::Arc;
use std::time::Instant;

use winit::application::ApplicationHandler;
use winit::dpi::PhysicalSize;
use winit::event::WindowEvent;
use pymol_scene::ButtonState;
use winit::event_loop::ActiveEventLoop;
use winit::window::{Window, WindowId};

use super::App;
use pymol_framework::message::AppMessage;

/// Convert a path to a command-safe string.
///
/// On Windows the native separator is `\`, which the command parser treats as
/// an escape character inside quoted strings.  Replacing with `/` is harmless
/// on Windows (the OS accepts both separators) and avoids parse failures.
fn command_safe_path(path: &std::path::Path) -> String {
    path.to_string_lossy().replace('\\', "/")
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

        // macOS: disable automatic window tabbing (removes "Show Tab Bar" from View menu)
        #[cfg(target_os = "macos")]
        {
            use objc2::MainThreadMarker;
            use objc2_app_kit::NSWindow;
            // Safety: resumed() is always called on the main thread by winit
            let mtm = unsafe { MainThreadMarker::new_unchecked() };
            NSWindow::setAllowsAutomaticWindowTabbing(false, mtm);
        }

        // Activate native menu bar (must happen after window creation).
        // The AppMenu is kept alive in self.native_menu — dropping the
        // muda::Menu would deallocate the underlying NSMenu / Win32 menu.
        if let Some(ref app_menu) = self.native_menu {
            #[cfg(target_os = "macos")]
            app_menu.menu.init_for_nsapp();

            #[cfg(target_os = "windows")]
            {
                use winit::raw_window_handle::{HasWindowHandle, RawWindowHandle};
                if let Ok(handle) = window.window_handle() {
                    if let RawWindowHandle::Win32(win32) = handle.as_raw() {
                        // Safety: the HWND is valid — we just created the window above.
                        unsafe {
                            let _ = app_menu.menu.init_for_hwnd(win32.hwnd.get() as isize);
                        }
                    }
                }
            }
        }

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

        // Execute startup script if it exists
        {
            let rc_path = pymol_settings::paths::pymolrc_path();
            if rc_path.is_file() {
                log::info!("Loading startup script: {}", rc_path.display());
                if let Err(e) = self.execute_command(
                    &format!("run \"{}\"", command_safe_path(&rc_path)),
                    false,
                ) {
                    log::warn!("Error in pymolrc: {}", e);
                }
            }
        }

        // Load pending file if any (use command executor for consistency)
        if let Some(path) = self.pending_load_file.take() {
            let path = path.replace('\\', "/");
            let _ = self.execute_command(&format!("load \"{}\"", &path), false);
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
                let egui_using_pointer = self.view.egui.ctx.is_using_pointer();
                if !egui_using_pointer {
                    self.handle_click_detection(state, button);
                } else {
                    // Clear stale click state when egui handles the event
                    // (e.g. press was on viewport but release lands on a popup)
                    self.viewport.click_start_pos = None;
                }
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
        // Process native menu events
        self.process_menu_events();

        // Sync menu check states with actual settings
        self.sync_menu_check_states();

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

        // Request redraw when plugins have notifications to display.
        // Without this, notification overlays and component state changes
        // (like busy→idle transitions) won't render until a user event.
        if !self.plugin_manager.notification_messages().is_empty() {
            self.mark_dirty();
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
// Native menu event handling
// =============================================================================

impl App {
    /// Poll and dispatch native menu events from `muda`.
    pub(crate) fn process_menu_events(&mut self) {
        // Collect pending menu event IDs to avoid holding an immutable borrow
        // on self.native_menu while calling &mut self methods.
        let ids = match &self.native_menu {
            Some(m) => &m.ids,
            None => return,
        };

        let mut pending = Vec::new();
        while let Ok(event) = muda::MenuEvent::receiver().try_recv() {
            pending.push(event.id);
        }
        if pending.is_empty() {
            return;
        }

        // Snapshot the IDs we need to compare against (all are small, cloneable).
        let run = ids.run.clone();
        let open = ids.open.clone();
        let fetch = ids.fetch.clone();
        let save = ids.save.clone();
        let export_png = ids.export_png.clone();
        let export_movie = ids.export_movie.clone();
        let select_all = ids.select_all.clone();
        let deselect_all = ids.deselect_all.clone();
        let reset_view = ids.reset_view.clone();
        let zoom_all = ids.zoom_all.clone();
        let orient = ids.orient.clone();
        let center = ids.center.clone();
        let opaque_background = ids.opaque_background.clone();
        let bg_white = ids.bg_white.clone();
        let bg_black = ids.bg_black.clone();
        let transparent_panels = ids.transparent_panels.clone();
        let help_commands = ids.help_commands.clone();

        // Snapshot panel component IDs for matching
        let panel_ids: Vec<String> = self
            .layout
            .panels
            .iter()
            .map(|p| p.component_id.clone())
            .collect();

        for id in pending {
            if id == run {
                self.menu_run_script();
            } else if id == open {
                self.menu_open_file();
            } else if id == fetch {
                self.fetch_dialog.open();
            } else if id == save {
                self.menu_save_file();
            } else if id == export_png {
                self.menu_export_png();
            } else if id == export_movie {
                self.menu_export_movie();
            } else if id == select_all {
                let _ = self.execute_command("select sele, all", false);
            } else if id == deselect_all {
                let _ = self.execute_command("deselect", false);
            } else if id == reset_view {
                let _ = self.execute_command("reset", false);
            } else if id == zoom_all {
                let _ = self.execute_command("zoom", false);
            } else if id == orient {
                let _ = self.execute_command("orient", false);
            } else if id == center {
                let _ = self.execute_command("center", false);
            } else if id == opaque_background {
                let current = self.state.settings.ui.opaque_background;
                let cmd = if current { "set opaque_background, off" } else { "set opaque_background, on" };
                let _ = self.execute_command(cmd, false);
            } else if id == bg_white {
                let _ = self.execute_command("bg_color white", false);
            } else if id == bg_black {
                let _ = self.execute_command("bg_color black", false);
            } else if id == transparent_panels {
                let current = self.state.settings.ui.transparent_panels;
                let cmd = if current { "set transparent_panels, off" } else { "set transparent_panels, on" };
                let _ = self.execute_command(cmd, false);
            } else if id == help_commands {
                let _ = self.execute_command("help", false);
            } else if let Some(component_id) = panel_ids.iter().find(|p| id == muda::MenuId::new(p.as_str())) {
                let visible = self
                    .layout
                    .panels
                    .iter()
                    .find(|p| p.component_id == *component_id)
                    .is_some_and(|p| p.visible);
                if visible {
                    self.bus.send(AppMessage::HidePanel(component_id.clone()));
                } else {
                    self.bus.send(AppMessage::ShowPanel(component_id.clone()));
                }
            }
        }

        self.mark_dirty();
        self.view.egui.ctx.request_repaint();
    }

    /// Sync native menu check states with actual application settings.
    pub(crate) fn sync_menu_check_states(&mut self) {
        if let Some(ref mut app_menu) = self.native_menu {
            app_menu.state.transparent_panels = self.state.settings.ui.transparent_panels;
            app_menu.state.opaque_background = self.state.settings.ui.opaque_background;
            app_menu.state.bg_color = self.state.clear_color;
            app_menu.state.panels = self
                .layout
                .panels
                .iter()
                .map(|p| {
                    let title = self
                        .components
                        .get_title(&p.component_id)
                        .unwrap_or_else(|| p.component_id.clone());
                    (p.component_id.clone(), title, p.visible)
                })
                .collect();
            app_menu.sync();
        }
    }

    /// Collect extra extensions registered by plugins for file format handlers.
    fn plugin_format_extensions(&self) -> Vec<String> {
        self.executor
            .format_handlers()
            .keys()
            .cloned()
            .collect()
    }

    /// Collect extra extensions registered by plugins for script handlers.
    fn plugin_script_extensions(&self) -> Vec<String> {
        self.executor
            .script_handlers()
            .keys()
            .cloned()
            .collect()
    }

    // ── File dialogs ────────────────────────────────────────────────────

    /// Run Script… — open a script file.
    fn menu_run_script(&mut self) {
        let mut builtin: Vec<&str> = vec!["pml"];
        let plugin_exts = self.plugin_script_extensions();
        let plugin_refs: Vec<&str> = plugin_exts.iter().map(|s| s.as_str()).collect();
        builtin.extend(&plugin_refs);

        let dialog = rfd::FileDialog::new()
            .set_title("Run Script")
            .add_filter("Script files", &builtin)
            .add_filter("All files", &["*"]);

        if let Some(path) = dialog.pick_file() {
            let _ = self.execute_command(&format!("run \"{}\"", command_safe_path(&path)), false);
        }
    }

    /// Trajectory-only file extensions.
    const TRAJECTORY_EXTS: &'static [&'static str] = &["xtc", "trr"];

    /// Open… — load a molecular structure or trajectory file.
    fn menu_open_file(&mut self) {
        let mut structure_exts: Vec<&str> = vec![
            "pdb", "ent", "cif", "mmcif", "bcif",
            "sdf", "mol", "sd", "mol2", "ml2",
            "xyz", "gro", "ccp4", "map", "mrc",
        ];
        let plugin_exts = self.plugin_format_extensions();
        let plugin_refs: Vec<&str> = plugin_exts.iter().map(|s| s.as_str()).collect();
        structure_exts.extend(&plugin_refs);
        let session_exts: Vec<&str> = vec![
            "pse", "pze", "prs",
        ];

        // Combined filter: structures + trajectories
        let mut all_exts = structure_exts.clone();
        all_exts.extend(Self::TRAJECTORY_EXTS);
        all_exts.extend(&session_exts);

        let dialog = rfd::FileDialog::new()
            .set_title("Open")
            .add_filter("All supported files", &all_exts)
            .add_filter("Session files", &session_exts)
            .add_filter("Molecular files", &structure_exts)
            .add_filter("Trajectory files", Self::TRAJECTORY_EXTS)
            .add_filter("All files", &["*"]);

        if let Some(path) = dialog.pick_file() {
            let ext = path.extension()
                .and_then(|e| e.to_str())
                .unwrap_or("")
                .to_lowercase();
            let path_str = command_safe_path(&path);

            if Self::TRAJECTORY_EXTS.contains(&ext.as_str()) {
                let _ = self.execute_command(
                    &format!("load_traj \"{}\"", path_str),
                    false,
                );
            } else {
                let _ = self.execute_command(
                    &format!("load \"{}\"", path_str),
                    false,
                );
            }
        }
    }

    /// Save… — save structure or session.
    fn menu_save_file(&mut self) {
        let mut builtin: Vec<&str> = vec![
            "pdb", "cif", "mmcif",
            "sdf", "mol", "mol2", "ml2",
            "xyz", "gro"
        ];
        let plugin_exts = self.plugin_format_extensions();
        let plugin_refs: Vec<&str> = plugin_exts.iter().map(|s| s.as_str()).collect();
        builtin.extend(&plugin_refs);

        let dialog = rfd::FileDialog::new()
            .set_title("Save")
            .add_filter("Molecular files", &builtin)
            .add_filter("Session", &["prs"])
            .add_filter("Scripts", &["pml"])
            .add_filter("All files", &["*"]);

        if let Some(path) = dialog.save_file() {
            let _ = self.execute_command(&format!("save \"{}\"", command_safe_path(&path)), false);
        }
    }

    /// Export PNG… — save a PNG screenshot.
    fn menu_export_png(&mut self) {
        let dialog = rfd::FileDialog::new()
            .set_title("Export PNG")
            .add_filter("PNG image", &["png"])
            .set_file_name("screenshot.png");

        if let Some(path) = dialog.save_file() {
            let _ = self.execute_command(&format!("png \"{}\"", command_safe_path(&path)), false);
        }
    }

    /// Export Movie… — render movie to video file.
    fn menu_export_movie(&mut self) {
        let dialog = rfd::FileDialog::new()
            .set_title("Export Movie")
            .add_filter("Video files", &["mp4", "mov", "webm"])
            .set_file_name("movie.mp4");

        if let Some(path) = dialog.save_file() {
            let _ = self.execute_command(&format!("mproduce \"{}\"", command_safe_path(&path)), false);
        }
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
                self.viewport.input.handle_mouse_button((*state).into(), (*button).into());
            }
            WindowEvent::CursorMoved { position, .. } => {
                self.viewport.input.handle_mouse_motion((position.x, position.y));
            }
            WindowEvent::MouseWheel { delta, .. } => {
                self.viewport.input.handle_scroll((*delta).into());
            }
            WindowEvent::PinchGesture { delta, .. } => {
                self.viewport.input.handle_pinch_zoom(*delta);
            }
            WindowEvent::ModifiersChanged(modifiers) => {
                self.viewport.input.handle_modifiers(modifiers.state().into());
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
            // Sync display state, object transforms, and camera to the new frame
            let current_frame = self.state.movie.current_frame();
            let state_index = self.state.movie.frame_to_state(current_frame);
            for name in self.state.registry.names().map(|s| s.to_string()).collect::<Vec<_>>() {
                if let Some(obj) = self.state.registry.get_molecule_mut(&name) {
                    obj.set_display_state(state_index);
                }
            }
            self.state.apply_movie_object_transforms();
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
            self.clear_viewport_image();
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
    fn handle_click_detection(&mut self, state: winit::event::ElementState, button: winit::event::MouseButton) {
        if button != winit::event::MouseButton::Left {
            return;
        }

        let state: ButtonState = state.into();
        let mouse_pos = self.viewport.input.mouse_position();
        if state == ButtonState::Pressed {
            if self.view.is_over_viewport(mouse_pos) {
                self.viewport.click_start_pos = Some(mouse_pos);
            }
        } else if state == ButtonState::Released {
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
        let path_str = command_safe_path(&path);
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
