//! Color resolution for molecular visualization
//!
//! Resolves `ColorIndex` values from atoms to final RGBA colors.

use pymol_color::{ChainColors, Color, ColorIndex, ColorRamp, ElementColors, NamedColors};
use pymol_mol::Atom;
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
    ) -> Self {
        Self {
            named_colors,
            element_colors,
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
    pub fn resolve_atom(&self, atom: &Atom) -> [f32; 4] {
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
        self.resolve_color_index_value(atom.repr.colors.base, atom)
    }

    /// Resolve a color index value to a Color
    ///
    /// Takes an explicit color index rather than reading from atom.colors.base,
    /// allowing per-representation colors to override the atom color.
    fn resolve_color_index_value(&self, color_idx: i32, atom: &Atom) -> Color {
        match ColorIndex::from(color_idx) {
            ColorIndex::ByElement | ColorIndex::Atomic => {
                self.element_colors.get(atom.element as u8)
            }
            ColorIndex::ByChain => self.resolve_by_chain(atom),
            ColorIndex::BySS => self.resolve_by_secondary_structure(atom),
            ColorIndex::ByBFactor => self.resolve_by_b_factor(atom),
            ColorIndex::Named(idx) => {
                self.named_colors
                    .get_by_index(idx)
                    .unwrap_or(Color::WHITE)
            }
            _ => self.element_colors.get(atom.element as u8),
        }
    }

    /// Resolve color by chain
    fn resolve_by_chain(&self, atom: &Atom) -> Color {
        // Use chain name (ChainColors::get is a static method)
        ChainColors::get(&atom.residue.chain)
    }

    /// Resolve color by secondary structure
    fn resolve_by_secondary_structure(&self, atom: &Atom) -> Color {
        pymol_color::ss_color(atom.ss_type as u8)
    }

    /// Resolve color for cartoon representation
    ///
    /// Uses colors.cartoon if set, otherwise falls back to colors.base.
    pub fn resolve_cartoon(&self, atom: &Atom) -> [f32; 4] {
        let color_idx = atom.repr.colors.cartoon_or_base();
        let color = self.resolve_color_index_value(color_idx, atom);
        color.to_rgba(self.default_alpha)
    }

    /// Generic rep-color resolution with 3-level fallback:
    /// 1. Per-atom rep color (if != COLOR_UNSET)
    /// 2. Settings default from SettingResolver (if >= 0)
    /// 3. Atom's base color
    ///
    /// A `default_color` of -1 means "no settings-level override, use atom color".
    pub fn resolve_rep_color(&self, atom: &Atom, per_atom_color: i32, default_color: i32) -> [f32; 4] {
        let color_idx = if per_atom_color != pymol_mol::COLOR_UNSET {
            per_atom_color
        } else if default_color >= 0 {
            default_color
        } else {
            atom.repr.colors.base
        };
        let color = self.resolve_color_index_value(color_idx, atom);
        color.to_rgba(self.default_alpha)
    }

    /// Resolve color for cartoon with an object-level default override
    pub fn resolve_cartoon_with_default(&self, atom: &Atom, default_color: i32) -> [f32; 4] {
        self.resolve_rep_color(atom, atom.repr.colors.cartoon, default_color)
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

}

#[cfg(test)]
mod tests {
    use super::*;
    use pymol_mol::Element;

    #[test]
    fn test_element_coloring() {
        let named = NamedColors::default();
        let elements = ElementColors::default();
        let resolver = ColorResolver::new(&named, &elements);

        // Create a carbon atom
        let mut atom = Atom::default();
        atom.element = Element::Carbon;
        atom.repr.colors.base = -1; // By element

        let color = resolver.resolve_atom(&atom);

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
        let resolver = ColorResolver::new(&named, &elements);

        let mut atom = Atom::default();
        atom.set_residue("ALA", 1, "A");
        atom.repr.colors.base = -2; // By chain

        let color_a = resolver.resolve_atom(&atom);

        atom.set_residue("ALA", 1, "B");
        let color_b = resolver.resolve_atom(&atom);

        // Different chains should have different colors
        assert!(color_a[0] != color_b[0] || color_a[1] != color_b[1] || color_a[2] != color_b[2]);
    }

    #[test]
    fn test_cartoon_with_default_priority() {
        use pymol_mol::COLOR_UNSET;

        let named = NamedColors::default();
        let elements = ElementColors::default();
        let resolver = ColorResolver::new(&named, &elements);

        let mut atom = Atom::default();
        atom.element = Element::Carbon;
        atom.repr.colors.base = -1; // By element

        // Priority 3: no override, no per-atom cartoon → uses base color (by element)
        let color = resolver.resolve_cartoon_with_default(&atom, -1);
        // Carbon element color (gray)
        assert!(color[0] > 0.4 && color[0] < 0.7, "Should be element color");

        // Priority 2: object-level default overrides base
        // Named color index 2 = red in default NamedColors
        let color = resolver.resolve_cartoon_with_default(&atom, 2);
        assert!(color[0] > 0.9, "Object default should override to red");
        assert!(color[1] < 0.1);

        // Priority 1: per-atom cartoon color overrides object default
        atom.repr.colors.cartoon = 3; // green
        let color = resolver.resolve_cartoon_with_default(&atom, 2);
        assert!(color[1] > 0.9, "Per-atom cartoon should override to green");
        assert!(color[0] < 0.1);

        // Reset: COLOR_UNSET means "not set", falls through to object default
        atom.repr.colors.cartoon = COLOR_UNSET;
        let color = resolver.resolve_cartoon_with_default(&atom, 2);
        assert!(color[0] > 0.9, "COLOR_UNSET should fall through to object default");
    }
}
