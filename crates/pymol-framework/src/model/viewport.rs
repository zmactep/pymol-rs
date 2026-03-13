//! Viewport Model
//!
//! 3D viewport interaction state: input, picker, hover detection, click tracking.
//! No egui dependency.

use pymol_scene::{InputState, Picker, PickHit};

use super::sequence::ResidueRef;

/// Viewport model: 3D molecular viewer interaction state.
pub struct ViewportModel {
    /// Mouse/keyboard input handler (camera control with sensitivity)
    pub input: InputState,
    /// CPU ray-casting picker for hover detection
    pub picker: Picker,
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
            picker: Picker::new(),
            hover_hit: None,
            sequence_hover: None,
            click_start_pos: None,
        }
    }
}
