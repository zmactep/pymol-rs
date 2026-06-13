//! Computational molecular algorithms.
//!
//! This crate provides general computational algorithms used across the workspace,
//! organized into three submodules:
//!
//! - [`linalg`] — Linear algebra: SVD, matrix operations
//! - [`align`] — Structural alignment: Kabsch, sequence alignment, CE
//! - [`symmetry`] — Crystallographic symmetry: unit cell math, space group operations

pub mod align;
pub mod dss;
pub mod linalg;
pub mod surface;
pub mod symmetry;

// Re-export linalg submodules at crate root for backward compatibility
pub use linalg::svd3;
pub use linalg::Svd3;

// Re-export alignment submodules at crate root for backward compatibility
pub use align::substitution_matrix;

// Re-export alignment types/functions at crate root
pub use align::{
    apply_transform, ce_align, global_align, kabsch, rmsd, superpose, AlignedPair, AlignmentResult,
    AlignmentScoring, CeParams, CeResult, KabschResult, SuperposeParams, SuperposeResult,
};

// Re-export symmetry submodules at crate root for convenience
pub use symmetry::crystal;
pub use symmetry::space_groups;

// Re-export DSS types at crate root
pub use dss::{
    dss as assign_dss, AngleWindow, BackboneResidue, DssParams, Dssp, DsspParams, PyMolDss,
    SecondaryStructureAssigner, SsType,
};

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
