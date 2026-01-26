//! Atom flags and bitfields
//!
//! Provides bitflags for atom properties and state.

use bitflags::bitflags;

bitflags! {
    /// Atom flags indicating various atom properties and states
    ///
    /// These flags are used for selection, display, and manipulation of atoms.
    #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
    pub struct AtomFlags: u32 {
        /// Atom is in focus (for sculpting)
        const FOCUS = 0x00000001;
        /// Atom is free to move (for sculpting)
        const FREE = 0x00000002;
        /// Atom is restrained (for sculpting)
        const RESTRAIN = 0x00000004;
        /// Atom is fixed (for sculpting)
        const FIX = 0x00000008;
        /// Atom is excluded from calculations
        const EXCLUDE = 0x00000010;
        /// Atom is part of a study
        const STUDY = 0x00000020;
        /// Atom is part of a protein
        const PROTEIN = 0x00000040;
        /// Atom is part of a nucleic acid
        const NUCLEIC = 0x00000080;
        /// Atom should be exfoliated (removed from surface)
        const EXFOLIATE = 0x01000000;
        /// Atom should be ignored
        const IGNORE = 0x02000000;
        /// Atom should not be smoothed
        const NO_SMOOTH = 0x04000000;
        /// Atom is part of a polymer
        const POLYMER = 0x08000000;
        /// Atom is a solvent molecule
        const SOLVENT = 0x10000000;
        /// Atom is an organic molecule
        const ORGANIC = 0x20000000;
        /// Atom is an inorganic molecule
        const INORGANIC = 0x40000000;
        /// Atom is a guide atom (for cartoon representation)
        const GUIDE = 0x80000000;

        /// Classification mask (polymer, solvent, organic, inorganic, guide)
        const CLASS = 0xF8000000;
        /// Mask for non-class flags
        const CLASS_MASK = 0x07FFFFFF;
    }
}

impl AtomFlags {
    /// Check if this atom is part of a biomolecule (protein or nucleic acid)
    #[inline]
    pub fn is_biomolecule(&self) -> bool {
        self.intersects(AtomFlags::PROTEIN | AtomFlags::NUCLEIC)
    }

    /// Check if this atom is a solvent
    #[inline]
    pub fn is_solvent(&self) -> bool {
        self.contains(AtomFlags::SOLVENT)
    }

    /// Check if this atom is organic
    #[inline]
    pub fn is_organic(&self) -> bool {
        self.contains(AtomFlags::ORGANIC)
    }

    /// Check if this atom is inorganic
    #[inline]
    pub fn is_inorganic(&self) -> bool {
        self.contains(AtomFlags::INORGANIC)
    }

    /// Get the classification bits only
    #[inline]
    pub fn classification(&self) -> AtomFlags {
        *self & AtomFlags::CLASS
    }
}

/// Atom geometry describing the bonding arrangement
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[repr(i8)]
pub enum AtomGeometry {
    /// No specific geometry
    #[default]
    None = 5,
    /// Single bond / terminal
    Single = 1,
    /// Linear geometry (sp hybridization)
    Linear = 2,
    /// Planar/trigonal geometry (sp2 hybridization)
    Planar = 3,
    /// Tetrahedral geometry (sp3 hybridization)
    Tetrahedral = 4,
}

impl AtomGeometry {
    /// Create geometry from raw value (as stored in files)
    pub fn from_raw(value: i8) -> Self {
        match value {
            1 => AtomGeometry::Single,
            2 => AtomGeometry::Linear,
            3 => AtomGeometry::Planar,
            4 => AtomGeometry::Tetrahedral,
            _ => AtomGeometry::None,
        }
    }

    /// Get the raw value
    pub fn as_raw(&self) -> i8 {
        *self as i8
    }

    /// Get expected bond count for this geometry
    pub fn expected_bonds(&self) -> Option<u8> {
        match self {
            AtomGeometry::None => None,
            AtomGeometry::Single => Some(1),
            AtomGeometry::Linear => Some(2),
            AtomGeometry::Planar => Some(3),
            AtomGeometry::Tetrahedral => Some(4),
        }
    }
}

/// Atom stereochemistry
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[repr(u8)]
pub enum Stereo {
    /// No stereochemistry specified
    #[default]
    None = 0,
    /// Odd parity (used in SDF format)
    Odd = 1,
    /// Even parity (used in SDF format)
    Even = 2,
    /// Either/unmarked
    Either = 3,
}

impl Stereo {
    /// Create from raw SDF stereo value
    pub fn from_sdf(value: u8) -> Self {
        match value {
            1 => Stereo::Odd,
            2 => Stereo::Even,
            3 => Stereo::Either,
            _ => Stereo::None,
        }
    }
}

/// Atom chirality (R/S nomenclature)
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[repr(u8)]
pub enum Chirality {
    /// No chirality specified
    #[default]
    None = 0,
    /// R configuration
    R = 1,
    /// S configuration
    S = 2,
    /// Unknown/unspecified chirality
    Unknown = 3,
}

impl Chirality {
    /// Create from character ('R', 'S', '?', or other)
    pub fn from_char(c: char) -> Self {
        match c {
            'R' | 'r' => Chirality::R,
            'S' | 's' => Chirality::S,
            '?' => Chirality::Unknown,
            _ => Chirality::None,
        }
    }

    /// Convert to character
    pub fn to_char(&self) -> char {
        match self {
            Chirality::None => ' ',
            Chirality::R => 'R',
            Chirality::S => 'S',
            Chirality::Unknown => '?',
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_atom_flags() {
        let mut flags = AtomFlags::empty();
        assert!(!flags.is_solvent());

        flags |= AtomFlags::SOLVENT;
        assert!(flags.is_solvent());
        assert!(flags.contains(AtomFlags::SOLVENT));

        flags |= AtomFlags::PROTEIN;
        assert!(flags.is_biomolecule());
    }

    #[test]
    fn test_atom_geometry() {
        assert_eq!(AtomGeometry::from_raw(4), AtomGeometry::Tetrahedral);
        assert_eq!(AtomGeometry::Tetrahedral.expected_bonds(), Some(4));
        assert_eq!(AtomGeometry::Planar.as_raw(), 3);
    }

    #[test]
    fn test_chirality() {
        assert_eq!(Chirality::from_char('R'), Chirality::R);
        assert_eq!(Chirality::R.to_char(), 'R');
    }
}
