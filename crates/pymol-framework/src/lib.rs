//! PyMOL-RS Framework
//!
//! Shared application infrastructure extracted from `pymol-gui` so that both
//! the host application and dynamically loaded plugins share the same types.
//!
//! ## Modules
//!
//! - [`message`] — `AppMessage` enum and `MessageBus` (per-frame two-phase queue)
//! - [`component`] — `Component` trait and `SharedContext`
//! - [`component_store`] — `ComponentStore` (typed + dynamic access)
//! - [`layout`] — `Slot`, `PanelSlot`, `PanelConfig`, `Layout` (data + state, NO egui rendering)
//! - [`topics`] — typed pub/sub helpers over `AppMessage::Custom`

pub mod component;
pub mod component_store;
pub mod layout;
pub mod message;
pub mod topics;
