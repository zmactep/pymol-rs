//! Error types for molecular operations
//!
//! Provides error types for operations on molecular data structures.

use thiserror::Error;

/// Errors that can occur when working with molecular data
#[derive(Error, Debug, Clone)]
pub enum MolError {
    /// Atom index is out of bounds
    #[error("Atom index {0} is out of bounds (max: {1})")]
    AtomIndexOutOfBounds(u32, usize),

    /// Bond index is out of bounds
    #[error("Bond index {0} is out of bounds (max: {1})")]
    BondIndexOutOfBounds(u32, usize),

    /// State/coordinate set index is out of bounds
    #[error("State index {0} is out of bounds (max: {1})")]
    StateIndexOutOfBounds(usize, usize),

    /// Invalid element symbol
    #[error("Invalid element symbol: {0}")]
    InvalidElement(String),

    /// Coordinate count doesn't match atom count
    #[error("Coordinate count mismatch: expected {expected}, got {actual}")]
    CoordinateMismatch { expected: usize, actual: usize },

    /// Attempting to add a duplicate bond
    #[error("Duplicate bond between atoms {0} and {1}")]
    DuplicateBond(u32, u32),

    /// Invalid bond (self-loop or invalid atoms)
    #[error("Invalid bond: atom1={0}, atom2={1}")]
    InvalidBond(u32, u32),

    /// Atom has no coordinates in the specified state
    #[error("Atom {atom} has no coordinates in state {state}")]
    NoCoordinates { atom: u32, state: usize },

    /// Operation requires a non-discrete molecule
    #[error("Operation not supported for discrete molecules")]
    DiscreteNotSupported,

    /// Operation requires a discrete molecule
    #[error("Operation requires a discrete molecule")]
    RequiresDiscrete,

    /// Invalid unique ID
    #[error("Invalid unique ID: {0}")]
    InvalidUniqueId(i32),

    /// General I/O error
    #[error("I/O error: {0}")]
    Io(String),

    /// Parse error
    #[error("Parse error: {0}")]
    Parse(String),
}

impl MolError {
    /// Create an atom out of bounds error
    pub fn atom_out_of_bounds(index: u32, max: usize) -> Self {
        MolError::AtomIndexOutOfBounds(index, max)
    }

    /// Create a bond out of bounds error
    pub fn bond_out_of_bounds(index: u32, max: usize) -> Self {
        MolError::BondIndexOutOfBounds(index, max)
    }

    /// Create a state out of bounds error
    pub fn state_out_of_bounds(index: usize, max: usize) -> Self {
        MolError::StateIndexOutOfBounds(index, max)
    }
}

/// Result type for molecular operations
pub type MolResult<T> = Result<T, MolError>;
