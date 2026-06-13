//! Domain Models
//!
//! Pure data types without egui dependency. Testable, serializable, headless-compatible.

pub mod command_line;
pub mod movie;
pub mod output;
pub mod scene;
pub mod sequence;
pub mod viewport;

pub use command_line::CommandLineModel;
pub use movie::{MovieBlock, MovieLane, MovieMarker, MovieTick, MovieTimeline, MovieTimelineModel};
pub use output::{OutputKind, OutputMessage, OutputModel};
pub use scene::{
    SceneColorContext, SceneEntry, SceneGroup, SceneMapVisualKind, SceneModel, SceneObject,
    SceneObjectKind, SceneSubchain, SidebarColor,
};
pub use sequence::{
    ResidueKind, ResidueRef, SeqChain, SeqObject, SeqResidue, SequenceColorContext, SequenceModel,
    SequenceUiState,
};
pub use viewport::ViewportModel;
