//! Secondary structure assignment algorithms
//!
//! Provides molecule-independent algorithms for assigning secondary structure
//! (helix, sheet, loop) to protein backbone residues.
//!
//! # Input
//!
//! All algorithms operate on [`BackboneResidue`] slices — pre-extracted backbone
//! atom positions with chain/residue metadata. The caller (typically `pymol-mol`)
//! is responsible for extracting this data from a molecular structure.
//!
//! # Algorithms
//!
//! - [`dss`] — PyMOL's DSS algorithm (H-bonds + phi/psi dihedral angles)

mod pymol;

pub use pymol::{dss, AngleWindow, DssParams};

use lin_alg::f32::Vec3;

/// Backbone residue data for secondary structure algorithms
///
/// Contains the minimum information needed by SS assignment algorithms:
/// backbone atom positions, chain/residue identifiers, N-H bond direction,
/// and peptide bond connectivity.
#[derive(Debug, Clone)]
pub struct BackboneResidue {
    /// Alpha carbon (CA) position
    pub ca: Vec3,
    /// Amide nitrogen (N) position
    pub n: Vec3,
    /// Carbonyl carbon (C) position
    pub c: Vec3,
    /// Carbonyl oxygen (O) position
    pub o: Vec3,
    /// Chain identifier
    pub chain: String,
    /// Residue sequence number
    pub resv: i32,
    /// Direction of N-H bond (unit vector), if available.
    /// Pre-computed from molecular topology by the caller.
    pub nh_direction: Option<Vec3>,
    /// Whether this residue is peptide-bonded to the previous entry in the slice.
    /// Used for chain break detection within algorithms.
    pub bonded_to_prev: bool,
}

/// Secondary structure type assigned by an algorithm
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SsType {
    /// Loop, turn, or coil
    #[default]
    Loop,
    /// Alpha helix
    Helix,
    /// Beta sheet/strand
    Sheet,
}
