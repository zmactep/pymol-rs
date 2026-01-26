//! Input handling for camera control
//!
//! Handles mouse and keyboard input for interactive camera manipulation.

use lin_alg::f32::Vec3;
use winit::event::{ElementState, MouseButton, MouseScrollDelta};
use winit::keyboard::ModifiersState;

/// Camera movement delta
#[derive(Debug, Clone)]
pub enum CameraDelta {
    /// Rotation around axes (in radians)
    Rotate { x: f32, y: f32 },
    /// Translation in view space
    Translate(Vec3),
    /// Zoom factor (positive = zoom in)
    Zoom(f32),
    /// Clipping plane adjustment
    Clip { front: f32, back: f32 },
}

/// Mouse button indices
const LEFT: usize = 0;
const RIGHT: usize = 1;
const MIDDLE: usize = 2;

/// Input state for camera control
///
/// Tracks mouse position, button state, and modifier keys to generate
/// camera movement deltas.
#[derive(Debug, Clone)]
pub struct InputState {
    /// Mouse button states (left, right, middle)
    mouse_buttons: [bool; 3],
    /// Current mouse position (physical pixels)
    mouse_pos: (f32, f32),
    /// Previous mouse position (for delta calculation)
    prev_mouse_pos: (f32, f32),
    /// Accumulated mouse delta since last update
    mouse_delta: (f32, f32),
    /// Accumulated scroll delta since last update
    scroll_delta: f32,
    /// Current modifier keys
    modifiers: ModifiersState,
    /// Rotation sensitivity (radians per pixel)
    pub rotate_sensitivity: f32,
    /// Translation sensitivity (units per pixel)
    pub translate_sensitivity: f32,
    /// Zoom sensitivity (factor per scroll unit)
    pub zoom_sensitivity: f32,
    /// Clip sensitivity (units per pixel)
    pub clip_sensitivity: f32,
}

impl Default for InputState {
    fn default() -> Self {
        Self {
            mouse_buttons: [false; 3],
            mouse_pos: (0.0, 0.0),
            prev_mouse_pos: (0.0, 0.0),
            mouse_delta: (0.0, 0.0),
            scroll_delta: 0.0,
            modifiers: ModifiersState::empty(),
            rotate_sensitivity: 0.005,
            translate_sensitivity: 0.02,
            zoom_sensitivity: 0.1,
            clip_sensitivity: 0.1,
        }
    }
}

impl InputState {
    /// Create a new input state with default settings
    pub fn new() -> Self {
        Self::default()
    }

    /// Handle a mouse button event
    pub fn handle_mouse_button(&mut self, state: ElementState, button: MouseButton) {
        let pressed = state == ElementState::Pressed;

        // Reset delta when starting a new drag to prevent accumulated
        // passive mouse movement from being applied
        if pressed {
            self.mouse_delta = (0.0, 0.0);
        }

        match button {
            MouseButton::Left => self.mouse_buttons[LEFT] = pressed,
            MouseButton::Right => self.mouse_buttons[RIGHT] = pressed,
            MouseButton::Middle => self.mouse_buttons[MIDDLE] = pressed,
            _ => {}
        }
    }

    /// Handle mouse movement
    pub fn handle_mouse_motion(&mut self, position: (f64, f64)) {
        let new_pos = (position.0 as f32, position.1 as f32);

        // Accumulate delta (only if we have a previous valid position)
        if self.mouse_pos != (0.0, 0.0) || self.prev_mouse_pos != (0.0, 0.0) {
            self.mouse_delta.0 += new_pos.0 - self.mouse_pos.0;
            self.mouse_delta.1 += new_pos.1 - self.mouse_pos.1;
        }

        self.prev_mouse_pos = self.mouse_pos;
        self.mouse_pos = new_pos;
    }

    /// Handle mouse scroll
    pub fn handle_scroll(&mut self, delta: MouseScrollDelta) {
        let scroll = match delta {
            MouseScrollDelta::LineDelta(_, y) => y,
            MouseScrollDelta::PixelDelta(pos) => pos.y as f32 / 100.0,
        };
        self.scroll_delta += scroll;
    }

    /// Handle modifier key changes
    pub fn handle_modifiers(&mut self, modifiers: ModifiersState) {
        self.modifiers = modifiers;
    }

    /// Check if left mouse button is pressed
    pub fn left_pressed(&self) -> bool {
        self.mouse_buttons[LEFT]
    }

    /// Check if right mouse button is pressed
    pub fn right_pressed(&self) -> bool {
        self.mouse_buttons[RIGHT]
    }

    /// Check if middle mouse button is pressed
    pub fn middle_pressed(&self) -> bool {
        self.mouse_buttons[MIDDLE]
    }

    /// Check if shift is held
    pub fn shift_held(&self) -> bool {
        self.modifiers.shift_key()
    }

    /// Check if ctrl is held
    pub fn ctrl_held(&self) -> bool {
        self.modifiers.control_key()
    }

    /// Check if alt is held
    pub fn alt_held(&self) -> bool {
        self.modifiers.alt_key()
    }

    /// Get current mouse position
    pub fn mouse_position(&self) -> (f32, f32) {
        self.mouse_pos
    }

    /// Process accumulated input and return camera deltas
    ///
    /// This consumes the accumulated deltas and returns the corresponding
    /// camera movements. Call this once per frame.
    ///
    /// The default mapping follows PyMOL conventions:
    /// - Left drag: Rotate (X/Y rotation)
    /// - Middle drag: Translate (X/Y panning)
    /// - Right drag: Zoom (Y), Clip (Shift+Y), Slab (Ctrl+Y)
    /// - Scroll: Zoom
    /// - Shift+Left drag: Translate (X/Y panning)
    /// - Ctrl+Left drag: Zoom (Y), Rotate-Z (X)
    pub fn take_camera_deltas(&mut self) -> Vec<CameraDelta> {
        let mut deltas = Vec::new();
        let (dx, dy) = self.mouse_delta;

        // Handle scroll for zoom
        if self.scroll_delta.abs() > 0.001 {
            deltas.push(CameraDelta::Zoom(self.scroll_delta * self.zoom_sensitivity));
        }

        // Handle mouse drag
        if dx.abs() > 0.001 || dy.abs() > 0.001 {
            if self.mouse_buttons[LEFT] {
                if self.modifiers.shift_key() {
                    // Shift+Left: Pan
                    deltas.push(CameraDelta::Translate(Vec3::new(
                        dx * self.translate_sensitivity,
                        -dy * self.translate_sensitivity,
                        0.0,
                    )));
                } else if self.modifiers.control_key() {
                    // Ctrl+Left: Zoom (Y) + Rotate-Z (X)
                    deltas.push(CameraDelta::Zoom(dy * self.zoom_sensitivity * 0.1));
                    deltas.push(CameraDelta::Rotate {
                        x: 0.0,
                        y: dx * self.rotate_sensitivity,
                    });
                } else {
                    // Left: Rotate
                    deltas.push(CameraDelta::Rotate {
                        x: dy * self.rotate_sensitivity,
                        y: dx * self.rotate_sensitivity,
                    });
                }
            } else if self.mouse_buttons[MIDDLE] {
                // Middle: Pan
                deltas.push(CameraDelta::Translate(Vec3::new(
                    dx * self.translate_sensitivity,
                    -dy * self.translate_sensitivity,
                    0.0,
                )));
            } else if self.mouse_buttons[RIGHT] {
                if self.modifiers.shift_key() {
                    // Shift+Right: Clip planes
                    deltas.push(CameraDelta::Clip {
                        front: -dy * self.clip_sensitivity,
                        back: dy * self.clip_sensitivity,
                    });
                } else if self.modifiers.control_key() {
                    // Ctrl+Right: Slab (move both clip planes together)
                    deltas.push(CameraDelta::Clip {
                        front: -dy * self.clip_sensitivity,
                        back: -dy * self.clip_sensitivity,
                    });
                } else {
                    // Right: Zoom
                    deltas.push(CameraDelta::Zoom(dy * self.zoom_sensitivity * 0.1));
                }
            }
        }

        // Reset accumulated deltas
        self.mouse_delta = (0.0, 0.0);
        self.scroll_delta = 0.0;

        deltas
    }

    /// Reset all state
    pub fn reset(&mut self) {
        self.mouse_buttons = [false; 3];
        self.mouse_delta = (0.0, 0.0);
        self.scroll_delta = 0.0;
        self.modifiers = ModifiersState::empty();
    }

    /// Check if any mouse button is pressed
    pub fn any_button_pressed(&self) -> bool {
        self.mouse_buttons.iter().any(|&b| b)
    }
}

/// PyMOL mouse button modes
///
/// PyMOL supports different mouse button configurations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(dead_code)]
pub enum MouseMode {
    /// Default 3-button mouse mode
    ThreeButton,
    /// Two-button mouse mode (right button does zoom)
    TwoButton,
    /// Single-button mouse mode
    OneButton,
}

impl Default for MouseMode {
    fn default() -> Self {
        MouseMode::ThreeButton
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_input_state_default() {
        let state = InputState::new();
        assert!(!state.left_pressed());
        assert!(!state.right_pressed());
        assert!(!state.middle_pressed());
    }

    #[test]
    fn test_mouse_button_handling() {
        let mut state = InputState::new();

        state.handle_mouse_button(ElementState::Pressed, MouseButton::Left);
        assert!(state.left_pressed());

        state.handle_mouse_button(ElementState::Released, MouseButton::Left);
        assert!(!state.left_pressed());
    }

    #[test]
    fn test_mouse_motion() {
        let mut state = InputState::new();

        state.handle_mouse_motion((100.0, 200.0));
        assert_eq!(state.mouse_position(), (100.0, 200.0));

        state.handle_mouse_motion((110.0, 205.0));
        assert_eq!(state.mouse_position(), (110.0, 205.0));
        assert!((state.mouse_delta.0 - 10.0).abs() < 0.001);
        assert!((state.mouse_delta.1 - 5.0).abs() < 0.001);
    }

    #[test]
    fn test_camera_deltas_rotation() {
        let mut state = InputState::new();

        // Press left button and move
        state.handle_mouse_button(ElementState::Pressed, MouseButton::Left);
        state.handle_mouse_motion((100.0, 100.0));
        state.handle_mouse_motion((110.0, 105.0));

        let deltas = state.take_camera_deltas();
        assert_eq!(deltas.len(), 1);

        if let CameraDelta::Rotate { x, y } = &deltas[0] {
            assert!(x.abs() > 0.0);
            assert!(y.abs() > 0.0);
        } else {
            panic!("Expected Rotate delta");
        }
    }

    #[test]
    fn test_camera_deltas_scroll_zoom() {
        let mut state = InputState::new();

        state.handle_scroll(MouseScrollDelta::LineDelta(0.0, 1.0));

        let deltas = state.take_camera_deltas();
        assert_eq!(deltas.len(), 1);

        if let CameraDelta::Zoom(z) = &deltas[0] {
            assert!(*z > 0.0);
        } else {
            panic!("Expected Zoom delta");
        }
    }
}
