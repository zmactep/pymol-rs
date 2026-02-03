//! Color resolution for molecular visualization
//!
//! Resolves `ColorIndex` values from atoms to final RGBA colors.

use pymol_color::{ChainColors, Color, ColorRamp, ElementColors, NamedColors};
use pymol_mol::{Atom, ObjectMolecule};
use std::collections::HashMap;

/// Resolves color indices to final RGBA values
///
/// The `ColorResolver` takes references to the various color tables and
/// resolves `ColorIndex` values to concrete RGBA colors based on atom
/// properties like element, chain, B-factor, etc.
pub struct ColorResolver<'a> {
    /// Named color table
    named_colors: &'a NamedColors,
    /// Element-based colors (CPK)
    element_colors: &'a ElementColors,
    /// Chain-based colors
    #[allow(dead_code)]
    chain_colors: &'a ChainColors,
    /// Color ramps for continuous coloring
    color_ramps: HashMap<String, &'a ColorRamp>,
    /// Default transparency (alpha)
    default_alpha: f32,
    /// B-factor range for coloring (min, max)
    b_factor_range: (f32, f32),
}

impl<'a> ColorResolver<'a> {
    /// Create a new color resolver with the given color tables
    pub fn new(
        named_colors: &'a NamedColors,
        element_colors: &'a ElementColors,
        chain_colors: &'a ChainColors,
    ) -> Self {
        Self {
            named_colors,
            element_colors,
            chain_colors,
            color_ramps: HashMap::new(),
            default_alpha: 1.0,
            b_factor_range: (0.0, 100.0),
        }
    }

    /// Set the default alpha (transparency) value
    pub fn with_alpha(mut self, alpha: f32) -> Self {
        self.default_alpha = alpha;
        self
    }

    /// Add a color ramp for continuous coloring
    pub fn with_ramp(mut self, name: impl Into<String>, ramp: &'a ColorRamp) -> Self {
        self.color_ramps.insert(name.into(), ramp);
        self
    }

    /// Set the B-factor range for coloring
    pub fn with_b_factor_range(mut self, min: f32, max: f32) -> Self {
        self.b_factor_range = (min, max);
        self
    }

    /// Resolve an atom's color to RGBA
    ///
    /// Uses the atom's `color` field to determine the coloring scheme:
    /// - Negative values indicate special schemes (by element, chain, etc.)
    /// - Non-negative values are indices into the named color table
    pub fn resolve_atom(&self, atom: &Atom, _molecule: &ObjectMolecule) -> [f32; 4] {
        let color = self.resolve_color_index(atom);
        color.to_rgba(self.default_alpha)
    }

    /// Resolve an atom's color with custom transparency
    pub fn resolve_atom_with_alpha(&self, atom: &Atom, alpha: f32) -> [f32; 4] {
        let color = self.resolve_color_index(atom);
        color.to_rgba(alpha)
    }

    /// Resolve a color index to a Color
    fn resolve_color_index(&self, atom: &Atom) -> Color {
        self.resolve_color_index_value(atom.colors.base, atom)
    }

    /// Resolve a color index value to a Color
    ///
    /// Takes an explicit color index rather than reading from atom.colors.base,
    /// allowing per-representation colors to override the atom color.
    fn resolve_color_index_value(&self, color_idx: i32, atom: &Atom) -> Color {
        // Interpret the color index:
        // - Positive: index into named color table
        // - -1: color by element
        // - Other negatives: various special schemes

        match color_idx {
            -1 => {
                // Color by element (CPK)
                self.element_colors.get(atom.element as u8)
            }
            -2 => {
                // Color by chain
                self.resolve_by_chain(atom)
            }
            -3 => {
                // Color by secondary structure
                self.resolve_by_secondary_structure(atom)
            }
            -4 => {
                // Color by B-factor
                self.resolve_by_b_factor(atom)
            }
            idx if idx >= 0 => {
                // Named color index
                self.named_colors
                    .get_by_index(idx as u32)
                    .unwrap_or(Color::WHITE)
            }
            _ => {
                // Unknown scheme, default to element color
                self.element_colors.get(atom.element as u8)
            }
        }
    }

    /// Resolve color by chain
    fn resolve_by_chain(&self, atom: &Atom) -> Color {
        // Use chain name (ChainColors::get is a static method)
        ChainColors::get(&atom.chain)
    }

    /// Resolve color by secondary structure
    fn resolve_by_secondary_structure(&self, atom: &Atom) -> Color {
        Self::ss_color(atom.ss_type)
    }

    /// Get color for a secondary structure type (static helper)
    ///
    /// Uses PyMOL's default colors from util.cbss():
    /// - Helix: red
    /// - Sheet: yellow
    /// - Loop: green
    fn ss_color(ss_type: pymol_mol::SecondaryStructure) -> Color {
        use pymol_mol::SecondaryStructure;

        match ss_type {
            SecondaryStructure::Helix | SecondaryStructure::Helix310 | SecondaryStructure::HelixPi => {
                // Helix: red (PyMOL default)
                Color::from_rgb8(255, 0, 0)
            }
            SecondaryStructure::Sheet => {
                // Sheet: yellow (PyMOL default)
                Color::from_rgb8(255, 255, 0)
            }
            SecondaryStructure::Loop | SecondaryStructure::Turn | SecondaryStructure::Bend => {
                // Loop/coil: green (PyMOL default)
                Color::from_rgb8(0, 255, 0)
            }
        }
    }

    /// Resolve color for cartoon representation
    ///
    /// Uses colors.cartoon if set, otherwise falls back to colors.base.
    /// When colors.cartoon is None and colors.base is -1 (by element),
    /// uses green as the default cartoon color.
    pub fn resolve_cartoon(&self, atom: &Atom, _molecule: &ObjectMolecule) -> [f32; 4] {
        // Use cartoon color if set, otherwise fall back to base color
        let color_idx = atom.colors.cartoon_or_base();
        let color = match color_idx {
            -1 => {
                // Default: green for cartoon (when explicitly "by element" or default)
                Color::new(0.0, 1.0, 0.0)
            }
            _ => self.resolve_color_index_value(color_idx, atom),
        };
        color.to_rgba(self.default_alpha)
    }

    /// Resolve color for ribbon representation
    ///
    /// Uses colors.ribbon if set, otherwise falls back to colors.base.
    pub fn resolve_ribbon(&self, atom: &Atom, _molecule: &ObjectMolecule) -> [f32; 4] {
        let color_idx = atom.colors.ribbon_or_base();
        let color = self.resolve_color_index_value(color_idx, atom);
        color.to_rgba(self.default_alpha)
    }

    /// Resolve color for stick representation
    ///
    /// Uses colors.stick if set, otherwise falls back to colors.base.
    pub fn resolve_stick(&self, atom: &Atom, _molecule: &ObjectMolecule) -> [f32; 4] {
        let color_idx = atom.colors.stick_or_base();
        let color = self.resolve_color_index_value(color_idx, atom);
        color.to_rgba(self.default_alpha)
    }

    /// Resolve color for line representation
    ///
    /// Uses colors.line if set, otherwise falls back to colors.base.
    pub fn resolve_line(&self, atom: &Atom, _molecule: &ObjectMolecule) -> [f32; 4] {
        let color_idx = atom.colors.line_or_base();
        let color = self.resolve_color_index_value(color_idx, atom);
        color.to_rgba(self.default_alpha)
    }

    /// Resolve color for sphere representation
    ///
    /// Uses colors.sphere if set, otherwise falls back to colors.base.
    pub fn resolve_sphere(&self, atom: &Atom, _molecule: &ObjectMolecule) -> [f32; 4] {
        let color_idx = atom.colors.sphere_or_base();
        let color = self.resolve_color_index_value(color_idx, atom);
        color.to_rgba(self.default_alpha)
    }

    /// Resolve color for surface representation
    ///
    /// Uses colors.surface if set, otherwise falls back to colors.base.
    pub fn resolve_surface(&self, atom: &Atom, _molecule: &ObjectMolecule) -> [f32; 4] {
        let color_idx = atom.colors.surface_or_base();
        let color = self.resolve_color_index_value(color_idx, atom);
        color.to_rgba(self.default_alpha)
    }

    /// Resolve color for mesh representation
    ///
    /// Uses colors.mesh if set, otherwise falls back to colors.base.
    pub fn resolve_mesh(&self, atom: &Atom, _molecule: &ObjectMolecule) -> [f32; 4] {
        let color_idx = atom.colors.mesh_or_base();
        let color = self.resolve_color_index_value(color_idx, atom);
        color.to_rgba(self.default_alpha)
    }

    /// Resolve color by B-factor using the blue-white-red ramp
    fn resolve_by_b_factor(&self, atom: &Atom) -> Color {
        let (min, max) = self.b_factor_range;
        let range = max - min;
        if range <= 0.0 {
            return Color::WHITE;
        }

        // Normalize B-factor to [0, 1]
        let t = ((atom.b_factor - min) / range).clamp(0.0, 1.0);

        // Use blue-white-red ramp if available, otherwise interpolate manually
        if let Some(ramp) = self.color_ramps.get("b_factor") {
            ramp.get_color(t)
        } else {
            // Default blue-white-red interpolation
            if t < 0.5 {
                // Blue to white
                let s = t * 2.0;
                Color::new(s, s, 1.0)
            } else {
                // White to red
                let s = (t - 0.5) * 2.0;
                Color::new(1.0, 1.0 - s, 1.0 - s)
            }
        }
    }

    /// Get a color by name
    pub fn get_named_color(&self, name: &str) -> Option<Color> {
        self.named_colors.get_by_name(name).map(|(_, c)| c)
    }

    /// Get a color by index
    pub fn get_color_by_index(&self, index: u32) -> Option<Color> {
        self.named_colors.get_by_index(index)
    }

    /// Get the element color for an atomic number
    pub fn get_element_color(&self, atomic_number: u8) -> Color {
        self.element_colors.get(atomic_number)
    }

    /// Get the chain color for a chain ID
    pub fn get_chain_color(&self, chain: &str) -> Color {
        ChainColors::get(chain)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use pymol_mol::Element;

    // ChainColors is a unit struct with static methods, no need to instantiate
    static CHAIN_COLORS: ChainColors = ChainColors;

    #[test]
    fn test_element_coloring() {
        let named = NamedColors::default();
        let elements = ElementColors::default();
        let resolver = ColorResolver::new(&named, &elements, &CHAIN_COLORS);

        // Create a carbon atom
        let mut atom = Atom::default();
        atom.element = Element::Carbon;
        atom.colors.base = -1; // By element

        let molecule = ObjectMolecule::new("test");
        let color = resolver.resolve_atom(&atom, &molecule);

        // Carbon should be dark gray
        assert!(color[0] > 0.4 && color[0] < 0.7);
        assert!(color[1] > 0.4 && color[1] < 0.7);
        assert!(color[2] > 0.4 && color[2] < 0.7);
        assert_eq!(color[3], 1.0); // Full opacity
    }

    #[test]
    fn test_chain_coloring() {
        let named = NamedColors::default();
        let elements = ElementColors::default();
        let resolver = ColorResolver::new(&named, &elements, &CHAIN_COLORS);

        let mut atom = Atom::default();
        atom.chain = "A".to_string();
        atom.colors.base = -2; // By chain

        let molecule = ObjectMolecule::new("test");
        let color_a = resolver.resolve_atom(&atom, &molecule);

        atom.chain = "B".to_string();
        let color_b = resolver.resolve_atom(&atom, &molecule);

        // Different chains should have different colors
        assert!(color_a[0] != color_b[0] || color_a[1] != color_b[1] || color_a[2] != color_b[2]);
    }
}
