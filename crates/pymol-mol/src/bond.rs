//! Bond data structure
//!
//! Provides the `Bond` struct and `BondOrder` enum for representing molecular bonds.

use crate::index::AtomIndex;

/// Bond order enumeration
///
/// Represents the type/order of a chemical bond.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[repr(u8)]
pub enum BondOrder {
    /// Unknown or unspecified bond order
    #[default]
    Unknown = 0,
    /// Single bond
    Single = 1,
    /// Double bond
    Double = 2,
    /// Triple bond
    Triple = 3,
    /// Aromatic/delocalized bond (1.5 order)
    Aromatic = 4,
}

impl BondOrder {
    /// Create from raw numeric value
    pub fn from_raw(value: i8) -> Self {
        match value {
            1 => BondOrder::Single,
            2 => BondOrder::Double,
            3 => BondOrder::Triple,
            4 => BondOrder::Aromatic,
            _ => BondOrder::Unknown,
        }
    }

    /// Get the raw numeric value
    #[inline]
    pub fn as_raw(&self) -> i8 {
        *self as i8
    }

    /// Get the effective bond order as a float (aromatic = 1.5)
    #[inline]
    pub fn as_float(&self) -> f32 {
        match self {
            BondOrder::Unknown => 0.0,
            BondOrder::Single => 1.0,
            BondOrder::Double => 2.0,
            BondOrder::Triple => 3.0,
            BondOrder::Aromatic => 1.5,
        }
    }

    /// Check if this is a multiple bond (double, triple, or aromatic)
    #[inline]
    pub fn is_multiple(&self) -> bool {
        matches!(
            self,
            BondOrder::Double | BondOrder::Triple | BondOrder::Aromatic
        )
    }

    /// Check if this is an aromatic bond
    #[inline]
    pub fn is_aromatic(&self) -> bool {
        *self == BondOrder::Aromatic
    }
}

impl From<i8> for BondOrder {
    fn from(value: i8) -> Self {
        BondOrder::from_raw(value)
    }
}

impl From<BondOrder> for i8 {
    fn from(order: BondOrder) -> Self {
        order.as_raw()
    }
}

impl std::fmt::Display for BondOrder {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            BondOrder::Unknown => write!(f, "?"),
            BondOrder::Single => write!(f, "-"),
            BondOrder::Double => write!(f, "="),
            BondOrder::Triple => write!(f, "#"),
            BondOrder::Aromatic => write!(f, ":"),
        }
    }
}

/// Bond stereochemistry
///
/// Used for representing E/Z isomerism and other stereochemical properties.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[repr(u8)]
pub enum BondStereo {
    /// No stereo specified
    #[default]
    None = 0,
    /// Up (wedge) bond - for 3D stereo
    Up = 1,
    /// Down (dash) bond - for 3D stereo
    Down = 2,
    /// Either up or down (wavy)
    Either = 3,
    /// E (trans) isomer
    E = 4,
    /// Z (cis) isomer
    Z = 5,
}

impl BondStereo {
    /// Create from SDF stereo value
    pub fn from_sdf(value: u8) -> Self {
        match value {
            1 => BondStereo::Up,
            4 => BondStereo::Either,
            6 => BondStereo::Down,
            _ => BondStereo::None,
        }
    }

    /// Convert to SDF stereo value
    pub fn to_sdf(&self) -> u8 {
        match self {
            BondStereo::None => 0,
            BondStereo::Up => 1,
            BondStereo::Down => 6,
            BondStereo::Either => 4,
            BondStereo::E => 0,
            BondStereo::Z => 0,
        }
    }

    /// Check if this is a wedge bond (up or down)
    #[inline]
    pub fn is_wedge(&self) -> bool {
        matches!(self, BondStereo::Up | BondStereo::Down | BondStereo::Either)
    }

    /// Check if this is a double bond stereo (E or Z)
    #[inline]
    pub fn is_double_bond_stereo(&self) -> bool {
        matches!(self, BondStereo::E | BondStereo::Z)
    }
}

/// Symmetry operation for crystallographic bonds
///
/// Represents the symmetry operation applied to the second atom of a bond
/// (first atom implicitly uses identity operation 1_555).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub struct SymOp {
    /// Symmetry operation index (1-based, 0 = no symmetry)
    pub op: u8,
    /// Translation in x
    pub tx: i8,
    /// Translation in y
    pub ty: i8,
    /// Translation in z
    pub tz: i8,
}

impl SymOp {
    /// Identity operation (no symmetry transformation)
    pub const IDENTITY: SymOp = SymOp {
        op: 0,
        tx: 0,
        ty: 0,
        tz: 0,
    };

    /// Check if this is the identity operation
    #[inline]
    pub fn is_identity(&self) -> bool {
        self.op == 0 && self.tx == 0 && self.ty == 0 && self.tz == 0
    }

    /// Create from PDB-style symmetry string (e.g., "1_555")
    pub fn from_pdb_string(s: &str) -> Option<Self> {
        let parts: Vec<&str> = s.split('_').collect();
        if parts.len() != 2 {
            return None;
        }

        let op = parts[0].parse::<u8>().ok()?;
        let trans = parts[1];
        if trans.len() != 3 {
            return None;
        }

        let chars: Vec<char> = trans.chars().collect();
        let tx = (chars[0].to_digit(10)? as i8) - 5;
        let ty = (chars[1].to_digit(10)? as i8) - 5;
        let tz = (chars[2].to_digit(10)? as i8) - 5;

        Some(SymOp { op, tx, ty, tz })
    }

    /// Convert to PDB-style symmetry string
    pub fn to_pdb_string(&self) -> String {
        format!(
            "{}_{}{}{}",
            self.op,
            (self.tx + 5) as u8,
            (self.ty + 5) as u8,
            (self.tz + 5) as u8
        )
    }
}

/// Bond data structure
///
/// Represents a chemical bond between two atoms.
/// By convention, `atom1 < atom2` (indices are ordered).
#[derive(Debug, Clone)]
pub struct Bond {
    /// Index of the first atom (always atom1 < atom2)
    pub atom1: AtomIndex,
    /// Index of the second atom
    pub atom2: AtomIndex,
    /// Bond order (single, double, triple, aromatic)
    pub order: BondOrder,
    /// Stereochemistry
    pub stereo: BondStereo,
    /// Symmetry operation for second atom (for crystallographic bonds)
    pub symop: SymOp,
    /// Unique ID for per-bond settings
    pub unique_id: Option<i32>,
    /// Whether this bond has per-bond settings
    pub has_setting: bool,
}

impl Default for Bond {
    fn default() -> Self {
        Bond {
            atom1: AtomIndex::invalid(),
            atom2: AtomIndex::invalid(),
            order: BondOrder::Single,
            stereo: BondStereo::None,
            symop: SymOp::IDENTITY,
            unique_id: None,
            has_setting: false,
        }
    }
}

impl Bond {
    /// Create a new bond between two atoms
    ///
    /// The atom indices are automatically ordered so that atom1 < atom2.
    pub fn new(a1: AtomIndex, a2: AtomIndex, order: BondOrder) -> Self {
        let (atom1, atom2) = if a1.0 <= a2.0 { (a1, a2) } else { (a2, a1) };
        Bond {
            atom1,
            atom2,
            order,
            ..Default::default()
        }
    }

    /// Create a single bond
    pub fn single(a1: AtomIndex, a2: AtomIndex) -> Self {
        Self::new(a1, a2, BondOrder::Single)
    }

    /// Create a double bond
    pub fn double(a1: AtomIndex, a2: AtomIndex) -> Self {
        Self::new(a1, a2, BondOrder::Double)
    }

    /// Create a triple bond
    pub fn triple(a1: AtomIndex, a2: AtomIndex) -> Self {
        Self::new(a1, a2, BondOrder::Triple)
    }

    /// Create an aromatic bond
    pub fn aromatic(a1: AtomIndex, a2: AtomIndex) -> Self {
        Self::new(a1, a2, BondOrder::Aromatic)
    }

    /// Check if this bond involves the given atom
    #[inline]
    pub fn involves(&self, atom: AtomIndex) -> bool {
        self.atom1 == atom || self.atom2 == atom
    }

    /// Get the other atom in the bond
    ///
    /// Returns `None` if the given atom is not part of this bond.
    #[inline]
    pub fn other(&self, atom: AtomIndex) -> Option<AtomIndex> {
        if self.atom1 == atom {
            Some(self.atom2)
        } else if self.atom2 == atom {
            Some(self.atom1)
        } else {
            None
        }
    }

    /// Check if this bond connects to a symmetry-related atom
    #[inline]
    pub fn has_symop(&self) -> bool {
        !self.symop.is_identity()
    }

    /// Check if this bond connects the two given atoms (in any order)
    #[inline]
    pub fn connects(&self, a1: AtomIndex, a2: AtomIndex) -> bool {
        (self.atom1 == a1 && self.atom2 == a2) || (self.atom1 == a2 && self.atom2 == a1)
    }
}

impl std::fmt::Display for Bond {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}{}{}", self.atom1, self.order, self.atom2)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bond_order() {
        assert_eq!(BondOrder::from_raw(1), BondOrder::Single);
        assert_eq!(BondOrder::Double.as_float(), 2.0);
        assert_eq!(BondOrder::Aromatic.as_float(), 1.5);
        assert!(BondOrder::Double.is_multiple());
        assert!(!BondOrder::Single.is_multiple());
    }

    #[test]
    fn test_bond_creation() {
        let bond = Bond::new(AtomIndex(5), AtomIndex(3), BondOrder::Double);
        // Should be ordered
        assert_eq!(bond.atom1, AtomIndex(3));
        assert_eq!(bond.atom2, AtomIndex(5));
        assert_eq!(bond.order, BondOrder::Double);
    }

    #[test]
    fn test_bond_involves() {
        let bond = Bond::single(AtomIndex(1), AtomIndex(2));
        assert!(bond.involves(AtomIndex(1)));
        assert!(bond.involves(AtomIndex(2)));
        assert!(!bond.involves(AtomIndex(3)));
    }

    #[test]
    fn test_bond_other() {
        let bond = Bond::single(AtomIndex(1), AtomIndex(2));
        assert_eq!(bond.other(AtomIndex(1)), Some(AtomIndex(2)));
        assert_eq!(bond.other(AtomIndex(2)), Some(AtomIndex(1)));
        assert_eq!(bond.other(AtomIndex(3)), None);
    }

    #[test]
    fn test_symop() {
        let symop = SymOp::from_pdb_string("1_555").unwrap();
        assert_eq!(symop.op, 1);
        assert_eq!(symop.tx, 0);
        assert_eq!(symop.ty, 0);
        assert_eq!(symop.tz, 0);

        let symop2 = SymOp::from_pdb_string("2_645").unwrap();
        assert_eq!(symop2.op, 2);
        assert_eq!(symop2.tx, 1);
        assert_eq!(symop2.ty, -1);
        assert_eq!(symop2.tz, 0);

        assert_eq!(symop2.to_pdb_string(), "2_645");
    }

    #[test]
    fn test_bond_display() {
        let bond = Bond::double(AtomIndex(1), AtomIndex(2));
        assert_eq!(format!("{}", bond), "1=2");
    }
}
