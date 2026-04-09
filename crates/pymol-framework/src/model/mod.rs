//! Domain Models
//!
//! Pure data types without egui dependency. Testable, serializable, headless-compatible.

pub mod command_line;
pub mod output;
pub mod sequence;
pub mod viewport;

pub use command_line::CommandLineModel;
pub use output::{OutputModel, OutputKind, OutputMessage};
pub use sequence::{SequenceModel, SequenceUiState, SequenceColorContext, SeqObject, SeqChain, SeqResidue, ResidueRef, ResidueKind};
pub use viewport::ViewportModel;
