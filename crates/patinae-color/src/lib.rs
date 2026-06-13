//! Color system.
//!
//! Color management for molecular visualization:
//! - Named colors (runtime registry)
//! - Gradients for continuous coloring (spectrum, b-factor)
//! - Per-scheme palettes (element, chain, SS, residue type)
//! - Themed palettes composing all schemes into one swappable unit

mod color;
mod constants;
mod error;
mod gradient;
mod named;
mod resolver;
mod scheme;
mod theme;

pub use color::{Color, ColorIndex, SCHEME_NAMES};
pub use constants::*;
pub use error::ColorError;
pub use gradient::Gradient;
pub use named::NamedPalette;
pub use resolver::ColorResolver;
pub use scheme::{ChainPalette, ColorScheme, ElementPalette, ResiduePalette, SsPalette};
pub use theme::ThemedPalette;
