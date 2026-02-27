//! Input & Interaction
//!
//! Mouse/keyboard input processing, camera control, atom picking,
//! hover detection, click-based selection, and keyboard shortcuts.

use pymol_render::picking::Ray;
use pymol_scene::{
    expand_pick_to_selection, normalize_matrix, pick_expression_for_hit,
    CameraDelta, KeyBinding, PickHit,
};
use pymol_select::build_sele_command;
use winit::event::{ElementState, KeyEvent};
use winit::keyboard::PhysicalKey;

use super::App;

impl App {
    /// Setup default keyboard shortcuts
    pub(crate) fn setup_default_key_bindings(&mut self) {
        // TODO: Implement default key bindings
    }

    /// Process accumulated input and update camera
    ///
    /// This processes mouse deltas accumulated in InputState and applies them
    /// to the camera. Called once per frame in update().
    pub(crate) fn process_input(&mut self) {
        // Check if mouse is over viewport using stored rect and current mouse position
        let mouse_pos = self.input.mouse_position();
        let over_viewport = self.view.is_over_viewport(mouse_pos);

        // Also skip camera updates when egui is using the pointer (e.g. dragging a
        // floating window that overlaps the viewport).
        let egui_using_pointer = self.view.egui_ctx.is_using_pointer();

        if !over_viewport || egui_using_pointer {
            // Still need to consume the deltas to avoid accumulation
            self.input.take_camera_deltas();
            return;
        }

        let deltas = self.input.take_camera_deltas();
        for delta in deltas {
            // Clear raytraced overlay on any camera change
            if self.state.raytraced_image.is_some() {
                self.state.raytraced_image = None;
            }
            match delta {
                CameraDelta::Rotate { x, y } => {
                    self.state.camera.rotate_x(x);
                    self.state.camera.rotate_y(y);
                    self.needs_redraw = true;
                }
                CameraDelta::Translate(v) => {
                    let vh = self.view.viewport_rect.map(|r| r.height()).unwrap_or(768.0);
                    let scale = self.state.camera.screen_vertex_scale(vh);
                    self.state.camera.translate(v * scale);
                    self.needs_redraw = true;
                }
                CameraDelta::Zoom(factor) => {
                    self.state.camera.zoom(factor);
                    self.needs_redraw = true;
                }
                CameraDelta::Clip { front, back } => {
                    let view = self.state.camera.view_mut();
                    view.clip_front = (view.clip_front + front).max(0.01);
                    view.clip_back = (view.clip_back + back).max(view.clip_front + 0.01);
                    self.needs_redraw = true;
                }
                CameraDelta::SlabScale(raw_delta) => {
                    let mws = self
                        .state
                        .settings
                        .get_float(pymol_settings::id::mouse_wheel_scale);
                    let scale = 1.0 + 0.04 * mws * raw_delta;

                    let view = self.state.camera.view_mut();
                    let avg = (view.clip_front + view.clip_back) * 0.5;
                    let half_width = (view.clip_back - avg).max(0.1);
                    let new_half = (half_width * scale).max(0.1);

                    view.clip_front = (avg - new_half).max(0.01);
                    view.clip_back = (avg + new_half).max(view.clip_front + 0.1);
                    self.needs_redraw = true;
                }
            }
        }
    }

    /// Ray-cast at a screen position and return the closest hit.
    ///
    /// Converts screen coordinates to a picking ray via the camera matrices,
    /// then tests against all enabled objects in the registry.
    pub(crate) fn pick_at(&mut self, screen_pos: (f32, f32)) -> Option<PickHit> {
        let vp = *self.view.viewport_rect.as_ref()?;
        let vp_x = screen_pos.0 - vp.min.x;
        let vp_y = screen_pos.1 - vp.min.y;

        let view_mat = self.state.camera.view_matrix();
        let mut view_inv = view_mat
            .inverse()
            .unwrap_or(lin_alg::f32::Mat4::new_identity());
        normalize_matrix(&mut view_inv);

        let ray = Ray::from_screen(
            vp_x, vp_y,
            vp.width(), vp.height(),
            &self.state.camera.projection_matrix().data,
            &view_inv.data,
        );

        self.picker.pick_ray(&ray, &self.state.registry)
    }

    /// Process hover picking — detect atom under cursor and update hover indicator.
    ///
    /// Called once per frame when the cursor moves over the viewport.
    /// Skipped while dragging (any mouse button held) to avoid interference with camera control.
    pub(crate) fn process_hover(&mut self) {
        if self.input.any_button_pressed() {
            return;
        }

        let mouse_pos = self.input.mouse_position();
        if !self.view.is_over_viewport(mouse_pos) {
            if self.hover_hit.is_some() {
                self.hover_hit = None;
                self.needs_redraw = true;
            }
            return;
        }

        let new_hit = self.pick_at(mouse_pos);

        let changed = match (&self.hover_hit, &new_hit) {
            (None, None) => false,
            (Some(_), None) | (None, Some(_)) => true,
            (Some(a), Some(b)) => {
                a.object_name != b.object_name || a.atom_index != b.atom_index
            }
        };

        if changed {
            self.hover_hit = new_hit;
            self.needs_redraw = true;
        }
    }

    /// Process a mouse click for atom/object selection.
    ///
    /// Clicking on unselected atoms adds them to `sele`.
    /// Clicking on already-selected atoms removes them from `sele`.
    /// Clicking on empty space clears `sele`.
    pub(crate) fn process_click(&mut self) {
        let mouse_pos = self.input.mouse_position();
        let hit = self.pick_at(mouse_pos);
        let mode = self.state.settings.get_int(pymol_settings::id::mouse_selection_mode);

        if let Some(ref hit) = hit {
            if let Some(mol_obj) = self.state.registry.get_molecule(&hit.object_name) {
                let mol = mol_obj.molecule();
                if let Some(expr) = pick_expression_for_hit(hit, mode, mol) {
                    // Check if hovered atoms overlap with current sele
                    let overlaps_sele = self.state.selections.get("sele").is_some_and(|entry| {
                        let sel = expand_pick_to_selection(hit, mode, mol);
                        pymol_select::select(mol, &entry.expression)
                            .map(|sele| sel.intersection(&sele).any())
                            .unwrap_or(false)
                    });

                    let has_sele = self.state.selections.contains("sele");
                    if let Some(cmd) = build_sele_command(&expr, overlaps_sele, has_sele) {
                        let _ = self.execute_command(&cmd, false);
                    }
                }
            }
        } else if self.state.selections.contains("sele") {
            let _ = self.execute_command("deselect sele", false);
        }

        // Sync selection highlights to sequence viewer immediately
        self.sync_selection_to_sequence();
    }

    /// Handle keyboard input
    pub(crate) fn handle_key(&mut self, event: KeyEvent) {
        if event.state != ElementState::Pressed {
            return;
        }

        // Don't process app key bindings when command input has focus
        // Let egui handle the keyboard events for text input
        if self.command_line.has_focus {
            return;
        }

        if let PhysicalKey::Code(key) = event.physical_key {
            let binding = KeyBinding {
                key,
                ctrl: self.input.ctrl_held(),
                shift: self.input.shift_held(),
                alt: self.input.alt_held(),
            };

            if let Some(action) = self.key_bindings.get_cloned(&binding) {
                // Need to work around the borrow checker here
                // Clone the action Arc and execute
                let action = action.clone();
                action(self);
            }
        }
    }
}
