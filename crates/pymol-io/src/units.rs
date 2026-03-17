//! Unit conversion constants for GROMACS formats (GRO, XTC, TRR)
//!
//! GROMACS stores coordinates in nanometers; PyMOL uses Angstroms.

/// Nanometers to Angstroms conversion factor
pub(crate) const NM_TO_ANGSTROM: f32 = 10.0;

/// Angstroms to nanometers conversion factor
pub(crate) const ANGSTROM_TO_NM: f32 = 0.1;
