//! PyMOL-RS Toolbar Plugin
//!
//! Quick-action toolbar with configurable groups of icon buttons.
//! Each button triggers CLI commands or custom actions.

mod component;
mod model;
mod panel;

use pymol_plugin::prelude::PanelConfig;
use pymol_plugin::pymol_plugin;

use component::ToolbarComponent;

pymol_plugin! {
    name: "toolbar",
    description: "Quick-action toolbar with configurable icon button groups",
    commands: [],
    components: [
        (ToolbarComponent::new(), PanelConfig::top(72.0))
    ],
}
