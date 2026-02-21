//! Core color types

use serde::{Deserialize, Serialize};

/// An RGB color with values in the range [0.0, 1.0]
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Color {
    pub r: f32,
    pub g: f32,
    pub b: f32,
}

impl Color {
    /// Create a new color from RGB values (0.0-1.0)
    pub const fn new(r: f32, g: f32, b: f32) -> Self {
        Color { r, g, b }
    }

    /// Create a color from RGB bytes (0-255)
    pub fn from_rgb8(r: u8, g: u8, b: u8) -> Self {
        Color {
            r: r as f32 / 255.0,
            g: g as f32 / 255.0,
            b: b as f32 / 255.0,
        }
    }

    /// Create a color from a packed 0x00RRGGBB integer
    pub fn from_packed_rgb(packed: i32) -> Self {
        Self::from_rgb8(
            ((packed >> 16) & 0xFF) as u8,
            ((packed >> 8) & 0xFF) as u8,
            (packed & 0xFF) as u8,
        )
    }

    /// Convert to a packed 0x00RRGGBB integer
    pub fn to_packed_rgb(&self) -> i32 {
        let [r, g, b] = self.to_rgb8();
        ((r as i32) << 16) | ((g as i32) << 8) | (b as i32)
    }

    /// Create a color from a hex string (e.g., "#FF0000", "0xFF0000", or "FF0000")
    pub fn from_hex(hex: &str) -> Option<Self> {
        let hex = hex.trim_start_matches('#').trim_start_matches("0x");
        if hex.len() != 6 {
            return None;
        }

        let r = u8::from_str_radix(&hex[0..2], 16).ok()?;
        let g = u8::from_str_radix(&hex[2..4], 16).ok()?;
        let b = u8::from_str_radix(&hex[4..6], 16).ok()?;

        Some(Self::from_rgb8(r, g, b))
    }

    /// Convert to RGB bytes (0-255)
    pub fn to_rgb8(&self) -> [u8; 3] {
        [
            (self.r * 255.0).clamp(0.0, 255.0) as u8,
            (self.g * 255.0).clamp(0.0, 255.0) as u8,
            (self.b * 255.0).clamp(0.0, 255.0) as u8,
        ]
    }

    /// Convert to array
    pub fn to_array(&self) -> [f32; 3] {
        [self.r, self.g, self.b]
    }

    /// Convert to array with alpha
    pub fn to_rgba(&self, alpha: f32) -> [f32; 4] {
        [self.r, self.g, self.b, alpha]
    }

    /// Linear interpolation between two colors
    pub fn lerp(&self, other: &Color, t: f32) -> Color {
        let t = t.clamp(0.0, 1.0);
        Color {
            r: self.r + (other.r - self.r) * t,
            g: self.g + (other.g - self.g) * t,
            b: self.b + (other.b - self.b) * t,
        }
    }

    /// Common colors
    pub const WHITE: Color = Color::new(1.0, 1.0, 1.0);
    pub const BLACK: Color = Color::new(0.0, 0.0, 0.0);
    pub const RED: Color = Color::new(1.0, 0.0, 0.0);
    pub const GREEN: Color = Color::new(0.0, 1.0, 0.0);
    pub const BLUE: Color = Color::new(0.0, 0.0, 1.0);
    pub const YELLOW: Color = Color::new(1.0, 1.0, 0.0);
    pub const CYAN: Color = Color::new(0.0, 1.0, 1.0);
    pub const MAGENTA: Color = Color::new(1.0, 0.0, 1.0);
    pub const GRAY: Color = Color::new(0.5, 0.5, 0.5);
}

impl Default for Color {
    fn default() -> Self {
        Color::WHITE
    }
}

impl From<[f32; 3]> for Color {
    fn from(arr: [f32; 3]) -> Self {
        Color::new(arr[0], arr[1], arr[2])
    }
}

impl From<Color> for [f32; 3] {
    fn from(c: Color) -> Self {
        c.to_array()
    }
}

/// A color index referencing either a named color or a custom color
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum ColorIndex {
    /// Index into the named color table
    Named(u32),
    /// Special color: by element
    ByElement,
    /// Special color: by chain
    ByChain,
    /// Special color: by secondary structure
    BySS,
    /// Special color: by B-factor (temperature factor)
    ByBFactor,
    /// Special color: by residue type
    ByResidueType,
    /// Atomic colors (per-atom custom colors)
    Atomic,
}

/// Primary color scheme names for autocompletion (one per variant).
/// Keep in sync with [`ColorIndex::from_scheme_name`].
pub const SCHEME_NAMES: &[&str] = &["atomic", "chain", "ss", "b"];

impl ColorIndex {
    /// Parse a color scheme name alias into a `ColorIndex`.
    ///
    /// Recognizes common PyMOL scheme names like "atomic", "chain", "ss", "b_factor", etc.
    pub fn from_scheme_name(name: &str) -> Option<Self> {
        match name.to_lowercase().as_str() {
            "atomic" | "cpk" | "element" | "by_element" => Some(ColorIndex::ByElement),
            "chain" | "by_chain" | "chainbow" => Some(ColorIndex::ByChain),
            "ss" | "secondary_structure" | "by_ss" | "dssp" => Some(ColorIndex::BySS),
            "b" | "b_factor" | "bfactor" | "by_b" => Some(ColorIndex::ByBFactor),
            _ => None,
        }
    }
}

impl Default for ColorIndex {
    fn default() -> Self {
        ColorIndex::Named(0) // White
    }
}

impl From<ColorIndex> for i32 {
    fn from(index: ColorIndex) -> i32 {
        match index {
            ColorIndex::Named(idx) => idx as i32,
            ColorIndex::ByElement | ColorIndex::Atomic => -1,
            ColorIndex::ByChain => -2,
            ColorIndex::BySS => -3,
            ColorIndex::ByBFactor => -4,
            ColorIndex::ByResidueType => -5,
        }
    }
}

impl From<i32> for ColorIndex {
    fn from(value: i32) -> Self {
        match value {
            -1 => ColorIndex::ByElement,
            -2 => ColorIndex::ByChain,
            -3 => ColorIndex::BySS,
            -4 => ColorIndex::ByBFactor,
            -5 => ColorIndex::ByResidueType,
            c if c >= 0 => ColorIndex::Named(c as u32),
            _ => ColorIndex::default(),
        }
    }
}
