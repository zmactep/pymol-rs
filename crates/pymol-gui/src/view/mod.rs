//! View Management
//!
//! This module contains all rendering and display resources:
//! - `AppView` - GPU resources, window, surface, egui integration
//! - `UiConfig` - UI layout configuration (panel visibility, dimensions)

mod app_view;
mod ui_config;

pub use app_view::AppView;
pub use ui_config::UiConfig;
