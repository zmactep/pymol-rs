//! Viewport Model
//!
//! 3D viewport interaction state: input, hover detection, click tracking.
//! No egui dependency.

use patinae_scene::{InputState, PickHit};

use super::sequence::ResidueRef;

/// Viewport model: 3D molecular viewer interaction state.
pub struct ViewportModel {
    /// Mouse/keyboard input handler (camera control with sensitivity)
    pub input: InputState,
    /// Current hover hit (atom under cursor)
    pub hover_hit: Option<PickHit>,
    /// Current sequence viewer hover (for 3D highlight)
    pub sequence_hover: Option<ResidueRef>,
    /// Mouse position at left button press (for click vs drag detection)
    pub click_start_pos: Option<(f32, f32)>,
}

impl Default for ViewportModel {
    fn default() -> Self {
        Self::new()
    }
}

impl ViewportModel {
    pub fn new() -> Self {
        Self {
            input: InputState::new(),
            hover_hit: None,
            sequence_hover: None,
            click_start_pos: None,
        }
    }
}
