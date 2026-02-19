//! PSE → Session conversion
//!
//! Converts the intermediate [`PseSession`] representation into the
//! live [`Session`] type used by pymol-rs for rendering and interaction.

mod pymol_colors;
mod to_session;

pub use to_session::pse_to_session;
