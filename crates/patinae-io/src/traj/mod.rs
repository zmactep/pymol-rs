//! Trajectory file format readers (coordinate-only formats)
//!
//! Trajectory formats like XTC and TRR contain only per-frame coordinates
//! without topology (atoms, bonds). They are loaded onto existing molecular
//! objects via the `load_traj` command.

pub mod trr;
pub mod xtc;

pub use trr::TrrReader;
pub use xtc::XtcReader;
