//! Coloring schemes for molecular visualization

use serde::{Deserialize, Serialize};

use crate::Color;

/// Color scheme for coloring molecules
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ColorScheme {
    /// Color by element type (CPK coloring)
    Element,
    /// Color by chain ID
    Chain,
    /// Color by secondary structure
    SecondaryStructure,
    /// Color by B-factor (temperature factor)
    BFactor,
    /// Color by residue type (hydrophobic, polar, charged)
    ResidueType,
    /// Color by residue index (rainbow along chain)
    ResidueIndex,
    /// Color by atom index
    AtomIndex,
    /// Uniform color (all atoms same color)
    Uniform,
}

/// Element-based coloring (CPK colors)
#[derive(Debug, Serialize, Deserialize)]
pub struct ElementColors {
    #[serde(
        serialize_with = "serialize_color_array",
        deserialize_with = "deserialize_color_array"
    )]
    colors: [Color; 119], // Elements 0-118
}

fn serialize_color_array<S: serde::Serializer>(
    colors: &[Color; 119],
    serializer: S,
) -> Result<S::Ok, S::Error> {
    use serde::ser::SerializeSeq;
    let mut seq = serializer.serialize_seq(Some(119))?;
    for c in colors.iter() {
        seq.serialize_element(c)?;
    }
    seq.end()
}

fn deserialize_color_array<'de, D: serde::Deserializer<'de>>(
    deserializer: D,
) -> Result<[Color; 119], D::Error> {
    let v: Vec<Color> = Vec::deserialize(deserializer)?;
    v.try_into()
        .map_err(|v: Vec<Color>| serde::de::Error::custom(format!("expected 119 colors, got {}", v.len())))
}

impl ElementColors {
    /// Create a new element color table with standard CPK colors
    pub fn new() -> Self {
        let mut colors = [Color::MAGENTA; 119]; // Default for unknown

        // Standard CPK/Jmol colors
        colors[0] = Color::MAGENTA; // Unknown/placeholder
        colors[1] = Color::WHITE; // H - Hydrogen
        colors[2] = Color::new(0.85, 1.0, 1.0); // He - Helium
        colors[3] = Color::new(0.8, 0.5, 1.0); // Li - Lithium
        colors[4] = Color::new(0.76, 1.0, 0.0); // Be - Beryllium
        colors[5] = Color::new(1.0, 0.71, 0.71); // B - Boron
        colors[6] = Color::new(0.56, 0.56, 0.56); // C - Carbon (gray)
        colors[7] = Color::new(0.19, 0.31, 0.97); // N - Nitrogen (blue)
        colors[8] = Color::new(1.0, 0.05, 0.05); // O - Oxygen (red)
        colors[9] = Color::new(0.56, 0.88, 0.31); // F - Fluorine
        colors[10] = Color::new(0.7, 0.89, 0.96); // Ne - Neon
        colors[11] = Color::new(0.67, 0.36, 0.95); // Na - Sodium
        colors[12] = Color::new(0.54, 1.0, 0.0); // Mg - Magnesium
        colors[13] = Color::new(0.75, 0.65, 0.65); // Al - Aluminum
        colors[14] = Color::new(0.94, 0.78, 0.63); // Si - Silicon
        colors[15] = Color::new(1.0, 0.5, 0.0); // P - Phosphorus (orange)
        colors[16] = Color::new(1.0, 1.0, 0.19); // S - Sulfur (yellow)
        colors[17] = Color::new(0.12, 0.94, 0.12); // Cl - Chlorine
        colors[18] = Color::new(0.5, 0.82, 0.89); // Ar - Argon
        colors[19] = Color::new(0.56, 0.25, 0.83); // K - Potassium
        colors[20] = Color::new(0.24, 1.0, 0.0); // Ca - Calcium
        colors[26] = Color::new(0.88, 0.4, 0.2); // Fe - Iron
        colors[29] = Color::new(0.78, 0.5, 0.2); // Cu - Copper
        colors[30] = Color::new(0.49, 0.5, 0.69); // Zn - Zinc
        colors[35] = Color::new(0.65, 0.16, 0.16); // Br - Bromine
        colors[53] = Color::new(0.58, 0.0, 0.58); // I - Iodine

        ElementColors { colors }
    }

    /// Get the color for an element by atomic number
    pub fn get(&self, atomic_number: u8) -> Color {
        self.colors
            .get(atomic_number as usize)
            .copied()
            .unwrap_or(Color::MAGENTA)
    }

    /// Get the color for an element by symbol
    pub fn get_by_symbol(&self, symbol: &str) -> Color {
        let atomic_number = match symbol.to_uppercase().as_str() {
            "H" => 1,
            "HE" => 2,
            "LI" => 3,
            "BE" => 4,
            "B" => 5,
            "C" => 6,
            "N" => 7,
            "O" => 8,
            "F" => 9,
            "NE" => 10,
            "NA" => 11,
            "MG" => 12,
            "AL" => 13,
            "SI" => 14,
            "P" => 15,
            "S" => 16,
            "CL" => 17,
            "AR" => 18,
            "K" => 19,
            "CA" => 20,
            "FE" => 26,
            "CU" => 29,
            "ZN" => 30,
            "BR" => 35,
            "I" => 53,
            _ => 0,
        };
        self.get(atomic_number)
    }

    /// Set a custom color for an element
    pub fn set(&mut self, atomic_number: u8, color: Color) {
        if (atomic_number as usize) < self.colors.len() {
            self.colors[atomic_number as usize] = color;
        }
    }
}

impl Default for ElementColors {
    fn default() -> Self {
        Self::new()
    }
}

/// Chain-based coloring
#[derive(Debug, Serialize, Deserialize)]
pub struct ChainColors;

impl ChainColors {
    /// Standard chain colors (cyclic)
    const COLORS: [Color; 10] = [
        Color::new(0.0, 1.0, 0.0),   // A - Green
        Color::new(0.0, 1.0, 1.0),   // B - Cyan
        Color::new(1.0, 0.5, 0.0),   // C - Orange
        Color::new(1.0, 1.0, 0.0),   // D - Yellow
        Color::new(1.0, 0.0, 1.0),   // E - Magenta
        Color::new(0.5, 0.5, 1.0),   // F - Light blue
        Color::new(1.0, 0.5, 0.5),   // G - Salmon
        Color::new(0.5, 1.0, 0.5),   // H - Light green
        Color::new(1.0, 0.75, 0.8),  // I - Pink
        Color::new(0.6, 0.6, 0.6),   // J - Gray
    ];

    /// Get color for a chain ID
    pub fn get(chain_id: &str) -> Color {
        if chain_id.is_empty() {
            return Color::WHITE;
        }

        // Use first character of chain ID
        let c = chain_id.chars().next().unwrap();
        let idx = if c.is_ascii_alphabetic() {
            (c.to_ascii_uppercase() as usize - 'A' as usize) % Self::COLORS.len()
        } else if c.is_ascii_digit() {
            (c as usize - '0' as usize) % Self::COLORS.len()
        } else {
            0
        };

        Self::COLORS[idx]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_element_colors() {
        let colors = ElementColors::new();

        // Carbon should be gray
        let carbon = colors.get_by_symbol("C");
        assert!(carbon.r < 0.6 && carbon.g < 0.6 && carbon.b < 0.6);

        // Nitrogen should be blue
        let nitrogen = colors.get_by_symbol("N");
        assert!(nitrogen.b > nitrogen.r && nitrogen.b > nitrogen.g);

        // Oxygen should be red
        let oxygen = colors.get_by_symbol("O");
        assert!(oxygen.r > oxygen.g && oxygen.r > oxygen.b);
    }

    #[test]
    fn test_chain_colors() {
        let color_a = ChainColors::get("A");
        let color_b = ChainColors::get("B");

        // Different chains should have different colors
        assert_ne!(color_a, color_b);
    }
}
