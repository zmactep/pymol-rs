//! Computational algorithms for PyMOL-RS
//!
//! This crate provides general computational algorithms used across pymol-rs,
//! organized into three submodules:
//!
//! - [`linalg`] — Linear algebra: SVD, matrix operations
//! - [`align`] — Structural alignment: Kabsch, sequence alignment, CE
//! - [`symmetry`] — Crystallographic symmetry: unit cell math, space group operations

pub mod linalg;
pub mod align;
pub mod symmetry;

// Re-export linalg submodules at crate root for backward compatibility
pub use linalg::svd3;
pub use linalg::{
    Svd3,
    mat3x3_to_mat4, left_multiply_mat4, transform_mat4, is_identity_mat4,
    transform_3x3, invert_3x3,
};

// Re-export alignment submodules at crate root for backward compatibility
pub use align::substitution_matrix;

// Re-export alignment types/functions at crate root
pub use align::{
    kabsch, rmsd, apply_transform, KabschResult,
    global_align, AlignedPair, AlignmentResult, AlignmentScoring,
    superpose, SuperposeParams, SuperposeResult,
    ce_align, CeParams, CeResult,
};

// Re-export symmetry submodules at crate root for convenience
pub use symmetry::crystal;
pub use symmetry::space_groups;

/// Errors from alignment algorithms
#[derive(Debug, thiserror::Error)]
pub enum AlignError {
    #[error("Coordinate arrays have different lengths: {0} vs {1}")]
    LengthMismatch(usize, usize),

    #[error("Not enough atoms for alignment (need at least 3, got {0})")]
    TooFewAtoms(usize),

    #[error("SVD failed to converge")]
    SvdFailed,

    #[error("No matching residues found between structures")]
    NoMatches,

    #[error("All pairs rejected as outliers")]
    AllRejected,
}
