//! Toolbar Data Model
//!
//! Types describing toolbar structure: groups, buttons, and actions.

// ---------------------------------------------------------------------------
// Data model
// ---------------------------------------------------------------------------

/// Action triggered by clicking a toolbar button.
pub enum ToolbarAction {
    /// Execute one or more CLI commands in sequence.
    Commands(Vec<String>),
}

/// A single toolbar button with icon, label, and action.
pub struct ToolbarButton {
    /// Unique identifier used as texture cache key.
    pub icon_id: &'static str,
    /// Embedded PNG image bytes (`include_bytes!`).
    pub icon_bytes: &'static [u8],
    /// Short label displayed below the icon.
    pub label: &'static str,
    /// Tooltip shown on hover.
    pub tooltip: &'static str,
    /// Action to perform when clicked.
    pub action: ToolbarAction,
}

/// A named group of toolbar buttons, rendered with a group label.
pub struct ToolbarGroup {
    /// Group label shown below the row of buttons.
    pub name: &'static str,
    /// Buttons in this group.
    pub buttons: Vec<ToolbarButton>,
}

// ---------------------------------------------------------------------------
// Default toolbar configuration
// ---------------------------------------------------------------------------

/// Build the default toolbar groups.
///
/// Each button uses an empty `icon_bytes` slice, which triggers placeholder
/// icon generation. Replace with `include_bytes!("../assets/icons/foo.png")`
/// when real icons are available.
pub fn default_toolbar_groups() -> Vec<ToolbarGroup> {
    vec![
        ToolbarGroup {
            name: "File",
            buttons: vec![
                ToolbarButton {
                    icon_id: "tb_open",
                    icon_bytes: &[],
                    label: "Open",
                    tooltip: "Open file (load)",
                    action: ToolbarAction::Commands(vec![]),
                },
                ToolbarButton {
                    icon_id: "tb_save",
                    icon_bytes: &[],
                    label: "Save",
                    tooltip: "Save session",
                    action: ToolbarAction::Commands(vec!["save".into()]),
                },
            ],
        },
        ToolbarGroup {
            name: "Images",
            buttons: vec![
                ToolbarButton {
                    icon_id: "tb_snapshot",
                    icon_bytes: &[],
                    label: "Snapshot",
                    tooltip: "Save PNG screenshot",
                    action: ToolbarAction::Commands(vec!["png".into()]),
                },
                ToolbarButton {
                    icon_id: "tb_ray",
                    icon_bytes: &[],
                    label: "Ray",
                    tooltip: "Ray trace current view",
                    action: ToolbarAction::Commands(vec!["ray".into()]),
                },
            ],
        },
        ToolbarGroup {
            name: "Atoms",
            buttons: vec![
                ToolbarButton {
                    icon_id: "tb_show_lines",
                    icon_bytes: &[],
                    label: "Show",
                    tooltip: "Show lines representation",
                    action: ToolbarAction::Commands(vec!["show lines".into()]),
                },
                ToolbarButton {
                    icon_id: "tb_hide_lines",
                    icon_bytes: &[],
                    label: "Hide",
                    tooltip: "Hide lines representation",
                    action: ToolbarAction::Commands(vec!["hide lines".into()]),
                },
            ],
        },
        ToolbarGroup {
            name: "Cartoons",
            buttons: vec![
                ToolbarButton {
                    icon_id: "tb_show_cartoon",
                    icon_bytes: &[],
                    label: "Show",
                    tooltip: "Show cartoon representation",
                    action: ToolbarAction::Commands(vec!["show cartoon".into()]),
                },
                ToolbarButton {
                    icon_id: "tb_hide_cartoon",
                    icon_bytes: &[],
                    label: "Hide",
                    tooltip: "Hide cartoon representation",
                    action: ToolbarAction::Commands(vec!["hide cartoon".into()]),
                },
            ],
        },
        ToolbarGroup {
            name: "Styles",
            buttons: vec![
                ToolbarButton {
                    icon_id: "tb_sticks",
                    icon_bytes: &[],
                    label: "Stick",
                    tooltip: "Show sticks representation",
                    action: ToolbarAction::Commands(vec!["show sticks".into()]),
                },
                ToolbarButton {
                    icon_id: "tb_spheres",
                    icon_bytes: &[],
                    label: "Sphere",
                    tooltip: "Show spheres representation",
                    action: ToolbarAction::Commands(vec!["show spheres".into()]),
                },
                ToolbarButton {
                    icon_id: "tb_ball_stick",
                    icon_bytes: &[],
                    label: "Ball\nstick",
                    tooltip: "Show ball-and-stick representation",
                    action: ToolbarAction::Commands(vec!["show nb_spheres".into()]),
                },
            ],
        },
        ToolbarGroup {
            name: "Background",
            buttons: vec![
                ToolbarButton {
                    icon_id: "tb_bg_white",
                    icon_bytes: &[],
                    label: "White",
                    tooltip: "Set white background",
                    action: ToolbarAction::Commands(vec!["bg_color white".into()]),
                },
                ToolbarButton {
                    icon_id: "tb_bg_black",
                    icon_bytes: &[],
                    label: "Black",
                    tooltip: "Set black background",
                    action: ToolbarAction::Commands(vec!["bg_color black".into()]),
                },
            ],
        },
        ToolbarGroup {
            name: "Lighting",
            buttons: vec![
                ToolbarButton {
                    icon_id: "tb_light_simple",
                    icon_bytes: &[],
                    label: "Simple",
                    tooltip: "Simple flat lighting",
                    action: ToolbarAction::Commands(vec![
                        "set ambient, 0.4".into(),
                        "set spec_reflect, 0.0".into(),
                    ]),
                },
                ToolbarButton {
                    icon_id: "tb_light_soft",
                    icon_bytes: &[],
                    label: "Soft",
                    tooltip: "Soft diffuse lighting",
                    action: ToolbarAction::Commands(vec![
                        "set ambient, 0.6".into(),
                        "set spec_reflect, 0.3".into(),
                    ]),
                },
                ToolbarButton {
                    icon_id: "tb_light_full",
                    icon_bytes: &[],
                    label: "Full",
                    tooltip: "Full specular lighting",
                    action: ToolbarAction::Commands(vec![
                        "set ambient, 0.2".into(),
                        "set spec_reflect, 1.0".into(),
                    ]),
                },
            ],
        },
        ToolbarGroup {
            name: "Selection",
            buttons: vec![
                ToolbarButton {
                    icon_id: "tb_select_all",
                    icon_bytes: &[],
                    label: "Inspect",
                    tooltip: "Select all atoms",
                    action: ToolbarAction::Commands(vec!["select sele, all".into()]),
                },
            ],
        },
    ]
}
