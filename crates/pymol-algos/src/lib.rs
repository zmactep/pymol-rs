//! Computational algorithms for PyMOL-RS
//!
//! This crate provides general computational algorithms used across pymol-rs:
//! - Analytical 3×3 SVD decomposition
//! - Kabsch algorithm for optimal rigid-body superposition
//! - Needleman-Wunsch global sequence alignment
//! - Iterative superposition with outlier rejection

mod svd3;
mod kabsch;
mod sequence_align;
mod superpose;

pub use svd3::Svd3;
pub use kabsch::{kabsch, rmsd, apply_transform, KabschResult};
pub use sequence_align::{global_align, AlignedPair, AlignmentResult, AlignmentScoring};
pub use superpose::{superpose, SuperposeParams, SuperposeResult};

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
