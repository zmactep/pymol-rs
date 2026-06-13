//! Structural alignment algorithms
//!
//! - Kabsch algorithm for optimal rigid-body superposition
//! - Needleman-Wunsch global sequence alignment
//! - Iterative superposition with outlier rejection
//! - Combinatorial Extension (CE) structural alignment

mod ce;
pub mod kabsch;
mod sequence_align;
pub mod substitution_matrix;
mod superpose;

pub use ce::{ce_align, CeParams, CeResult};
pub use kabsch::{apply_transform, kabsch, rmsd, KabschResult};
pub use sequence_align::{global_align, AlignedPair, AlignmentResult, AlignmentScoring};
pub use superpose::{superpose, SuperposeParams, SuperposeResult};
