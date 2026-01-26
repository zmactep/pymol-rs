//! Type-safe index wrappers
//!
//! Provides newtype wrappers around indices to prevent accidentally mixing
//! atom indices with bond indices or coordinate set indices.

use std::fmt;

/// Invalid index marker value
pub const INVALID_INDEX: u32 = u32::MAX;

/// Type-safe index into an atom array
///
/// Prevents accidentally using a bond index where an atom index is expected.
#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Default)]
#[repr(transparent)]
pub struct AtomIndex(pub u32);

impl AtomIndex {
    /// Create a new atom index
    #[inline]
    pub const fn new(index: u32) -> Self {
        AtomIndex(index)
    }

    /// Get the raw index value
    #[inline]
    pub const fn as_usize(&self) -> usize {
        self.0 as usize
    }

    /// Get the raw u32 value
    #[inline]
    pub const fn as_u32(&self) -> u32 {
        self.0
    }

    /// Check if this is a valid index
    #[inline]
    pub const fn is_valid(&self) -> bool {
        self.0 != INVALID_INDEX
    }

    /// Create an invalid index
    #[inline]
    pub const fn invalid() -> Self {
        AtomIndex(INVALID_INDEX)
    }
}

impl fmt::Debug for AtomIndex {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_valid() {
            write!(f, "AtomIndex({})", self.0)
        } else {
            write!(f, "AtomIndex(INVALID)")
        }
    }
}

impl fmt::Display for AtomIndex {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_valid() {
            write!(f, "{}", self.0)
        } else {
            write!(f, "INVALID")
        }
    }
}

impl From<u32> for AtomIndex {
    #[inline]
    fn from(index: u32) -> Self {
        AtomIndex(index)
    }
}

impl From<usize> for AtomIndex {
    #[inline]
    fn from(index: usize) -> Self {
        AtomIndex(index as u32)
    }
}

impl From<AtomIndex> for u32 {
    #[inline]
    fn from(index: AtomIndex) -> Self {
        index.0
    }
}

impl From<AtomIndex> for usize {
    #[inline]
    fn from(index: AtomIndex) -> Self {
        index.0 as usize
    }
}

/// Type-safe index into a bond array
///
/// Prevents accidentally using an atom index where a bond index is expected.
#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Default)]
#[repr(transparent)]
pub struct BondIndex(pub u32);

impl BondIndex {
    /// Create a new bond index
    #[inline]
    pub const fn new(index: u32) -> Self {
        BondIndex(index)
    }

    /// Get the raw index value
    #[inline]
    pub const fn as_usize(&self) -> usize {
        self.0 as usize
    }

    /// Get the raw u32 value
    #[inline]
    pub const fn as_u32(&self) -> u32 {
        self.0
    }

    /// Check if this is a valid index
    #[inline]
    pub const fn is_valid(&self) -> bool {
        self.0 != INVALID_INDEX
    }

    /// Create an invalid index
    #[inline]
    pub const fn invalid() -> Self {
        BondIndex(INVALID_INDEX)
    }
}

impl fmt::Debug for BondIndex {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_valid() {
            write!(f, "BondIndex({})", self.0)
        } else {
            write!(f, "BondIndex(INVALID)")
        }
    }
}

impl fmt::Display for BondIndex {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_valid() {
            write!(f, "{}", self.0)
        } else {
            write!(f, "INVALID")
        }
    }
}

impl From<u32> for BondIndex {
    #[inline]
    fn from(index: u32) -> Self {
        BondIndex(index)
    }
}

impl From<usize> for BondIndex {
    #[inline]
    fn from(index: usize) -> Self {
        BondIndex(index as u32)
    }
}

impl From<BondIndex> for u32 {
    #[inline]
    fn from(index: BondIndex) -> Self {
        index.0
    }
}

impl From<BondIndex> for usize {
    #[inline]
    fn from(index: BondIndex) -> Self {
        index.0 as usize
    }
}

/// Type-safe index into a coordinate set array
///
/// Used to identify a specific state/frame in a multi-state molecule.
#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Default)]
#[repr(transparent)]
pub struct StateIndex(pub u32);

impl StateIndex {
    /// Create a new state index
    #[inline]
    pub const fn new(index: u32) -> Self {
        StateIndex(index)
    }

    /// Get the raw index value
    #[inline]
    pub const fn as_usize(&self) -> usize {
        self.0 as usize
    }

    /// Get the raw u32 value
    #[inline]
    pub const fn as_u32(&self) -> u32 {
        self.0
    }

    /// Check if this is a valid index
    #[inline]
    pub const fn is_valid(&self) -> bool {
        self.0 != INVALID_INDEX
    }

    /// Create an invalid index
    #[inline]
    pub const fn invalid() -> Self {
        StateIndex(INVALID_INDEX)
    }
}

impl fmt::Debug for StateIndex {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_valid() {
            write!(f, "StateIndex({})", self.0)
        } else {
            write!(f, "StateIndex(INVALID)")
        }
    }
}

impl fmt::Display for StateIndex {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_valid() {
            write!(f, "{}", self.0)
        } else {
            write!(f, "INVALID")
        }
    }
}

impl From<u32> for StateIndex {
    #[inline]
    fn from(index: u32) -> Self {
        StateIndex(index)
    }
}

impl From<usize> for StateIndex {
    #[inline]
    fn from(index: usize) -> Self {
        StateIndex(index as u32)
    }
}

impl From<StateIndex> for u32 {
    #[inline]
    fn from(index: StateIndex) -> Self {
        index.0
    }
}

impl From<StateIndex> for usize {
    #[inline]
    fn from(index: StateIndex) -> Self {
        index.0 as usize
    }
}

/// Index within a coordinate set (maps to atom index via lookup table)
#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Default)]
#[repr(transparent)]
pub struct CoordIndex(pub u32);

impl CoordIndex {
    /// Create a new coord index
    #[inline]
    pub const fn new(index: u32) -> Self {
        CoordIndex(index)
    }

    /// Get the raw index value
    #[inline]
    pub const fn as_usize(&self) -> usize {
        self.0 as usize
    }

    /// Get the raw u32 value
    #[inline]
    pub const fn as_u32(&self) -> u32 {
        self.0
    }

    /// Check if this is a valid index
    #[inline]
    pub const fn is_valid(&self) -> bool {
        self.0 != INVALID_INDEX
    }

    /// Create an invalid index
    #[inline]
    pub const fn invalid() -> Self {
        CoordIndex(INVALID_INDEX)
    }
}

impl fmt::Debug for CoordIndex {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_valid() {
            write!(f, "CoordIndex({})", self.0)
        } else {
            write!(f, "CoordIndex(INVALID)")
        }
    }
}

impl fmt::Display for CoordIndex {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_valid() {
            write!(f, "{}", self.0)
        } else {
            write!(f, "INVALID")
        }
    }
}

impl From<u32> for CoordIndex {
    #[inline]
    fn from(index: u32) -> Self {
        CoordIndex(index)
    }
}

impl From<usize> for CoordIndex {
    #[inline]
    fn from(index: usize) -> Self {
        CoordIndex(index as u32)
    }
}

impl From<CoordIndex> for u32 {
    #[inline]
    fn from(index: CoordIndex) -> Self {
        index.0
    }
}

impl From<CoordIndex> for usize {
    #[inline]
    fn from(index: CoordIndex) -> Self {
        index.0 as usize
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_atom_index() {
        let idx = AtomIndex::new(42);
        assert_eq!(idx.as_usize(), 42);
        assert_eq!(idx.as_u32(), 42);
        assert!(idx.is_valid());

        let invalid = AtomIndex::invalid();
        assert!(!invalid.is_valid());
    }

    #[test]
    fn test_bond_index() {
        let idx = BondIndex::new(10);
        assert_eq!(idx.as_usize(), 10);
        assert!(idx.is_valid());
    }

    #[test]
    fn test_state_index() {
        let idx = StateIndex::new(5);
        assert_eq!(idx.as_usize(), 5);
        assert!(idx.is_valid());
    }

    #[test]
    fn test_index_conversions() {
        let atom_idx: AtomIndex = 100u32.into();
        assert_eq!(u32::from(atom_idx), 100);

        let atom_idx: AtomIndex = 50usize.into();
        assert_eq!(usize::from(atom_idx), 50);
    }

    #[test]
    fn test_index_display() {
        assert_eq!(format!("{}", AtomIndex::new(5)), "5");
        assert_eq!(format!("{}", AtomIndex::invalid()), "INVALID");
        assert_eq!(format!("{:?}", AtomIndex::new(5)), "AtomIndex(5)");
        assert_eq!(format!("{:?}", AtomIndex::invalid()), "AtomIndex(INVALID)");
    }

    #[test]
    fn test_index_ordering() {
        let a = AtomIndex::new(1);
        let b = AtomIndex::new(2);
        assert!(a < b);
        assert!(b > a);
    }
}
