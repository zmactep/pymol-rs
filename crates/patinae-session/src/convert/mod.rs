//! PSE → Session conversion
//!
//! Converts the intermediate [`crate::pse::PseSession`] representation into the
//! live [`patinae_scene::Session`] type used by patinae for rendering and
//! interaction.

mod pymol_colors;
mod to_session;

pub use to_session::pse_to_session;
