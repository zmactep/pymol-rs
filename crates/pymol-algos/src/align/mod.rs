//! Structural alignment algorithms
//!
//! - Kabsch algorithm for optimal rigid-body superposition
//! - Needleman-Wunsch global sequence alignment
//! - Iterative superposition with outlier rejection
//! - Combinatorial Extension (CE) structural alignment

pub mod kabsch;
mod ce;
mod sequence_align;
mod superpose;
pub mod substitution_matrix;

pub use kabsch::{kabsch, rmsd, apply_transform, KabschResult};
pub use sequence_align::{global_align, AlignedPair, AlignmentResult, AlignmentScoring};
pub use superpose::{superpose, SuperposeParams, SuperposeResult};
pub use ce::{ce_align, CeParams, CeResult};
