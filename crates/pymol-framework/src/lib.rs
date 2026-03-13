//! PyMOL-RS Framework
//!
//! Shared application infrastructure extracted from `pymol-gui` so that both
//! the host application and dynamically loaded plugins share the same types.
//!
//! ## Modules
//!
//! - [`message`] — `AppMessage` enum and `MessageBus` (per-frame two-phase queue)
//! - [`component`] — `Component` trait and `SharedContext`
//! - [`layout`] — `Slot`, `PanelSlot`, `PanelConfig`, `Layout` (data + state, NO egui rendering)
//! - [`topics`] — typed pub/sub helpers over `AppMessage::Custom`
//! - [`model`] — pure domain models (no UI dependency) for panels and viewport
//! - [`completion`] — command-line tab completion engine

pub mod component;
pub mod completion;
pub mod layout;
pub mod message;
pub mod model;
pub mod topics;
