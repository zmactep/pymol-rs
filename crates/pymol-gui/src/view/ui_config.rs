//! UI Configuration
//!
//! Contains UI layout configuration such as panel visibility and dimensions.

/// UI layout configuration
#[derive(Debug, Clone)]
pub struct UiConfig {
    /// Show the output panel at top
    pub show_output_panel: bool,
    /// Show the right-side control panel
    pub show_control_panel: bool,
    /// Show the sequence viewer panel at bottom
    pub show_sequence_panel: bool,
    /// Whether the sequence panel is floating (true) or docked at bottom (false)
    pub sequence_panel_floating: bool,
    /// Default size for the floating sequence window
    pub sequence_window_size: [f32; 2],
    /// Output panel height in pixels
    pub output_panel_height: f32,
    /// Right panel width in pixels
    pub right_panel_width: f32,
}

impl Default for UiConfig {
    fn default() -> Self {
        Self::new()
    }
}

impl UiConfig {
    /// Create a new UI config with default values
    pub fn new() -> Self {
        Self {
            show_output_panel: true,
            show_control_panel: true,
            show_sequence_panel: false,
            sequence_panel_floating: false,
            sequence_window_size: [600.0, 150.0],
            output_panel_height: 150.0,
            right_panel_width: 200.0,
        }
    }
}
