//! PyMOL-RS Color System
//!
//! This crate provides color management for PyMOL-RS, including:
//! - Named colors (PyMOL's built-in color palette)
//! - Color ramps for continuous coloring
//! - Per-atom coloring schemes (by element, chain, etc.)

mod color;
mod named;
mod ramp;
mod scheme;
mod error;

pub use color::{Color, ColorIndex, SCHEME_NAMES};
pub use named::NamedColors;
pub use ramp::{ColorRamp, RampType};
pub use scheme::{ChainColors, ColorScheme, ElementColors, ss_color};
pub use error::ColorError;

/// Re-export commonly used types
pub mod prelude {
    pub use crate::{ChainColors, Color, ColorIndex, ColorRamp, ColorScheme, ElementColors, NamedColors};
}
