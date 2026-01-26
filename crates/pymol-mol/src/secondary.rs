//! Secondary structure types
//!
//! Provides secondary structure classification for protein residues.

use std::fmt;

/// Secondary structure type for protein residues
///
/// Based on DSSP classification used in PDB files and PyMOL.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[repr(u8)]
pub enum SecondaryStructure {
    /// Loop, turn, or coil (no regular secondary structure)
    #[default]
    Loop = 0,
    /// Alpha helix
    Helix = 1,
    /// Beta sheet/strand
    Sheet = 2,
    /// 3-10 helix
    Helix310 = 3,
    /// Pi helix
    HelixPi = 4,
    /// Turn
    Turn = 5,
    /// Bend
    Bend = 6,
}

impl SecondaryStructure {
    /// Create from single character code (as used in PDB/PyMOL)
    ///
    /// - 'H' = Helix (alpha helix)
    /// - 'S' or 'E' = Sheet (extended beta strand)
    /// - 'L', ' ', or other = Loop
    /// - 'G' = 3-10 helix
    /// - 'I' = Pi helix
    /// - 'T' = Turn
    /// - 'B' = Bend
    pub fn from_char(c: char) -> Self {
        match c {
            'H' | 'h' => SecondaryStructure::Helix,
            'S' | 's' | 'E' | 'e' => SecondaryStructure::Sheet,
            'G' | 'g' => SecondaryStructure::Helix310,
            'I' | 'i' => SecondaryStructure::HelixPi,
            'T' | 't' => SecondaryStructure::Turn,
            'B' | 'b' => SecondaryStructure::Bend,
            _ => SecondaryStructure::Loop,
        }
    }

    /// Convert to single character code
    pub fn to_char(&self) -> char {
        match self {
            SecondaryStructure::Loop => 'L',
            SecondaryStructure::Helix => 'H',
            SecondaryStructure::Sheet => 'S',
            SecondaryStructure::Helix310 => 'G',
            SecondaryStructure::HelixPi => 'I',
            SecondaryStructure::Turn => 'T',
            SecondaryStructure::Bend => 'B',
        }
    }

    /// Convert to PyMOL-style two-character code
    pub fn to_pymol_code(&self) -> &'static str {
        match self {
            SecondaryStructure::Loop => "L ",
            SecondaryStructure::Helix => "H ",
            SecondaryStructure::Sheet => "S ",
            SecondaryStructure::Helix310 => "G ",
            SecondaryStructure::HelixPi => "I ",
            SecondaryStructure::Turn => "T ",
            SecondaryStructure::Bend => "B ",
        }
    }

    /// Check if this is any type of helix
    #[inline]
    pub fn is_helix(&self) -> bool {
        matches!(
            self,
            SecondaryStructure::Helix
                | SecondaryStructure::Helix310
                | SecondaryStructure::HelixPi
        )
    }

    /// Check if this is a beta sheet/strand
    #[inline]
    pub fn is_sheet(&self) -> bool {
        *self == SecondaryStructure::Sheet
    }

    /// Check if this is a loop or coil
    #[inline]
    pub fn is_loop(&self) -> bool {
        matches!(
            self,
            SecondaryStructure::Loop | SecondaryStructure::Turn | SecondaryStructure::Bend
        )
    }

    /// Check if this is regular secondary structure (helix or sheet)
    #[inline]
    pub fn is_regular(&self) -> bool {
        self.is_helix() || self.is_sheet()
    }

    /// Get a human-readable name
    pub fn name(&self) -> &'static str {
        match self {
            SecondaryStructure::Loop => "Loop",
            SecondaryStructure::Helix => "Alpha Helix",
            SecondaryStructure::Sheet => "Beta Sheet",
            SecondaryStructure::Helix310 => "3-10 Helix",
            SecondaryStructure::HelixPi => "Pi Helix",
            SecondaryStructure::Turn => "Turn",
            SecondaryStructure::Bend => "Bend",
        }
    }
}

impl fmt::Display for SecondaryStructure {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_char())
    }
}

impl From<char> for SecondaryStructure {
    fn from(c: char) -> Self {
        SecondaryStructure::from_char(c)
    }
}

impl From<SecondaryStructure> for char {
    fn from(ss: SecondaryStructure) -> Self {
        ss.to_char()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_secondary_structure_from_char() {
        assert_eq!(SecondaryStructure::from_char('H'), SecondaryStructure::Helix);
        assert_eq!(SecondaryStructure::from_char('S'), SecondaryStructure::Sheet);
        assert_eq!(SecondaryStructure::from_char('E'), SecondaryStructure::Sheet);
        assert_eq!(SecondaryStructure::from_char('L'), SecondaryStructure::Loop);
        assert_eq!(SecondaryStructure::from_char(' '), SecondaryStructure::Loop);
        assert_eq!(SecondaryStructure::from_char('G'), SecondaryStructure::Helix310);
    }

    #[test]
    fn test_secondary_structure_to_char() {
        assert_eq!(SecondaryStructure::Helix.to_char(), 'H');
        assert_eq!(SecondaryStructure::Sheet.to_char(), 'S');
        assert_eq!(SecondaryStructure::Loop.to_char(), 'L');
    }

    #[test]
    fn test_secondary_structure_classification() {
        assert!(SecondaryStructure::Helix.is_helix());
        assert!(SecondaryStructure::Helix310.is_helix());
        assert!(SecondaryStructure::Sheet.is_sheet());
        assert!(SecondaryStructure::Loop.is_loop());
        assert!(SecondaryStructure::Turn.is_loop());

        assert!(SecondaryStructure::Helix.is_regular());
        assert!(SecondaryStructure::Sheet.is_regular());
        assert!(!SecondaryStructure::Loop.is_regular());
    }

    #[test]
    fn test_secondary_structure_default() {
        assert_eq!(SecondaryStructure::default(), SecondaryStructure::Loop);
    }

    #[test]
    fn test_secondary_structure_display() {
        assert_eq!(format!("{}", SecondaryStructure::Helix), "H");
    }
}
