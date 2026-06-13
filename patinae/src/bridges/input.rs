//! Input bridge: Slint ViewportState globals → patinae-scene InputState.
//!
//! Slint exposes boolean "is pressed" state for mouse buttons. InputState
//! needs press/release events. This bridge tracks previous frame state
//! and emits transitions.

use patinae_scene::{ButtonState, InputState, Modifiers, MouseButton, ScrollDelta};

use crate::ViewportState;

/// Click detection threshold in logical pixels (scale-independent).
const CLICK_THRESHOLD_LP: f32 = 5.0;

#[derive(Debug, Clone, Copy)]
pub struct PointerSnapshot {
    pub mouse_logical: (f32, f32),
    pub press_logical: (f32, f32),
    pub left_pressed: bool,
    pub middle_pressed: bool,
    pub right_pressed: bool,
    pub suppress_click: bool,
}

/// Translates Slint viewport input globals into patinae-scene [`InputState`] events.
pub struct InputBridge {
    prev_left: bool,
    prev_middle: bool,
    prev_right: bool,
    /// Mouse position at left button press (logical pixels), for click vs drag detection.
    click_start_pos: Option<(f32, f32)>,
    /// Last Slint click event serial consumed. Slint keeps this as a durable
    /// event counter so quick press+release pairs are not lost between frames.
    last_click_serial: Option<i32>,
    /// A click that passed the threshold test, ready to be consumed by the picker.
    pending_click: Option<(f32, f32)>,
}

impl InputBridge {
    pub fn new() -> Self {
        Self {
            prev_left: false,
            prev_middle: false,
            prev_right: false,
            click_start_pos: None,
            last_click_serial: None,
            pending_click: None,
        }
    }

    /// Consume and return a pending click position (physical pixels), if any.
    pub fn take_pending_click(&mut self) -> Option<(f32, f32)> {
        self.pending_click.take()
    }

    /// Read current state from Slint ViewportState and push events into InputState.
    ///
    /// `winit_modifiers` is `(shift, ctrl, alt, super)` from winit's
    /// `ModifiersChanged` event — always up-to-date, unlike Slint's
    /// `pointer-event` which only fires on mouse button changes.
    pub fn sync(
        &mut self,
        input: &mut InputState,
        vp: &ViewportState,
        scale_factor: f32,
        winit_modifiers: (bool, bool, bool, bool),
        pointer: Option<PointerSnapshot>,
    ) {
        let (mouse_logical, press_logical, left_now, middle_now, right_now, suppress_click) =
            if let Some(pointer) = pointer {
                (
                    pointer.mouse_logical,
                    pointer.press_logical,
                    pointer.left_pressed,
                    pointer.middle_pressed,
                    pointer.right_pressed,
                    pointer.suppress_click,
                )
            } else {
                (
                    (vp.get_mouse_x(), vp.get_mouse_y()),
                    (vp.get_press_x(), vp.get_press_y()),
                    vp.get_left_pressed(),
                    vp.get_middle_pressed(),
                    vp.get_right_pressed(),
                    vp.get_suppress_click(),
                )
            };
        let mouse_phys = (
            mouse_logical.0 * scale_factor,
            mouse_logical.1 * scale_factor,
        );
        let press_phys = (
            press_logical.0 * scale_factor,
            press_logical.1 * scale_factor,
        );
        let click_serial = vp.get_click_serial();
        let slint_click = if self
            .last_click_serial
            .replace(click_serial)
            .is_some_and(|prev| prev != click_serial)
        {
            Some((vp.get_click_x(), vp.get_click_y()))
        } else {
            None
        };

        // Button transitions
        let left_started = !self.prev_left && left_now;
        let left_released = self.prev_left && !left_now;
        let any_button_started =
            left_started || (!self.prev_middle && middle_now) || (!self.prev_right && right_now);

        if any_button_started {
            input.handle_mouse_motion((press_phys.0 as f64, press_phys.1 as f64));
        }

        Self::detect_transition(&mut self.prev_left, left_now, input, MouseButton::Left);
        Self::detect_transition(
            &mut self.prev_middle,
            middle_now,
            input,
            MouseButton::Middle,
        );
        Self::detect_transition(&mut self.prev_right, right_now, input, MouseButton::Right);

        input.handle_mouse_motion((mouse_phys.0 as f64, mouse_phys.1 as f64));

        // Click detection: record press position, detect short-distance release.
        // Threshold is in logical pixels (scale-independent).
        if let Some(click_logical) = slint_click {
            self.click_start_pos = None;
            let dx = click_logical.0 - press_logical.0;
            let dy = click_logical.1 - press_logical.1;
            if (dx * dx + dy * dy).sqrt() < CLICK_THRESHOLD_LP && !suppress_click {
                self.pending_click = Some((
                    click_logical.0 * scale_factor,
                    click_logical.1 * scale_factor,
                ));
            }
            if suppress_click {
                vp.set_suppress_click(false);
            }
        } else if left_started {
            self.click_start_pos = Some(press_logical);
        } else if left_released {
            if let Some(start) = self.click_start_pos.take() {
                let dx = mouse_logical.0 - start.0;
                let dy = mouse_logical.1 - start.1;
                if (dx * dx + dy * dy).sqrt() < CLICK_THRESHOLD_LP && !suppress_click {
                    self.pending_click = Some(mouse_phys);
                }
            }
            if suppress_click {
                vp.set_suppress_click(false);
            }
        } else if !left_now {
            self.click_start_pos = None;
            if suppress_click {
                vp.set_suppress_click(false);
            }
        }

        // Modifiers: use winit state (always current, even mid-drag)
        let (shift, ctrl, alt, logo) = winit_modifiers;
        input.handle_modifiers(Modifiers {
            shift,
            ctrl,
            alt,
            logo,
        });

        // Scroll (accumulated by Slint, consumed here)
        let scroll = vp.get_scroll_delta();
        if scroll != 0.0 {
            input.handle_scroll(ScrollDelta::PixelDelta(0.0, scroll as f64));
            vp.set_scroll_delta(0.0);
        }
    }

    fn detect_transition(
        prev: &mut bool,
        current: bool,
        input: &mut InputState,
        button: MouseButton,
    ) {
        if current != *prev {
            let state = if current {
                ButtonState::Pressed
            } else {
                ButtonState::Released
            };
            input.handle_mouse_button(state, button);
            *prev = current;
        }
    }
}
