//! DOM event to InputState mapping.
//!
//! Converts JS mouse/keyboard event parameters into the framework-agnostic
//! input types used by `pymol_scene::InputState`.

use pymol_scene::{ButtonState, InputState, Modifiers, MouseButton, ScrollDelta};

/// Modifier flag bits matching the JS `MouseEvent.buttons` / modifier booleans.
pub const MOD_SHIFT: u32 = 1;
pub const MOD_CTRL: u32 = 2;
pub const MOD_ALT: u32 = 4;
pub const MOD_META: u32 = 8;

/// Convert a JS modifier bitmask to `Modifiers`.
pub fn modifiers_from_bits(bits: u32) -> Modifiers {
    Modifiers {
        shift: bits & MOD_SHIFT != 0,
        ctrl: bits & MOD_CTRL != 0,
        alt: bits & MOD_ALT != 0,
        logo: bits & MOD_META != 0,
    }
}

/// Convert a JS mouse button index (0=left, 1=middle, 2=right) to `MouseButton`.
pub fn mouse_button_from_js(button: u32) -> MouseButton {
    match button {
        0 => MouseButton::Left,
        1 => MouseButton::Middle,
        2 => MouseButton::Right,
        n => MouseButton::Other(n as u16),
    }
}

/// Apply a mouse-down event to the input state.
pub fn handle_mouse_down(input: &mut InputState, x: f32, y: f32, button: u32, mods: u32) {
    input.handle_modifiers(modifiers_from_bits(mods));
    input.handle_mouse_motion((x as f64, y as f64));
    input.handle_mouse_button(ButtonState::Pressed, mouse_button_from_js(button));
}

/// Apply a mouse-move event to the input state.
pub fn handle_mouse_move(input: &mut InputState, x: f32, y: f32, mods: u32) {
    input.handle_modifiers(modifiers_from_bits(mods));
    input.handle_mouse_motion((x as f64, y as f64));
}

/// Apply a mouse-up event to the input state.
pub fn handle_mouse_up(input: &mut InputState, x: f32, y: f32, button: u32) {
    input.handle_mouse_motion((x as f64, y as f64));
    input.handle_mouse_button(ButtonState::Released, mouse_button_from_js(button));
}

/// Apply a wheel event to the input state.
pub fn handle_wheel(input: &mut InputState, delta_y: f32, mods: u32) {
    input.handle_modifiers(modifiers_from_bits(mods));
    // Normalize: browsers typically send delta in pixels; treat as line-based for consistency.
    input.handle_scroll(ScrollDelta::LineDelta(0.0, -delta_y / 100.0));
}
