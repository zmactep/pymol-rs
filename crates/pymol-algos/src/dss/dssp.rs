//! DSSP secondary structure assignment algorithm
//!
//! Implements the Kabsch & Sander DSSP algorithm for secondary structure
//! assignment based on hydrogen bond energy patterns.
//!
//! Currently a stub — returns Loop for all residues.
//!
//! # References
//!
//! - Kabsch W, Sander C (1983). "Dictionary of protein secondary structure:
//!   pattern recognition of hydrogen-bonded and geometrical features."
//!   Biopolymers 22(12):2577-637.

use super::{BackboneResidue, SecondaryStructureAssigner, SsType};

/// Parameters for the DSSP algorithm
#[derive(Debug, Clone)]
pub struct DsspParams {
    /// Hydrogen bond energy cutoff in kcal/mol (Kabsch & Sander default: -0.5)
    pub hbond_energy_cutoff: f32,
}

impl Default for DsspParams {
    fn default() -> Self {
        Self {
            hbond_energy_cutoff: -0.5,
        }
    }
}

/// DSSP secondary structure assigner.
///
/// Implements the Kabsch & Sander algorithm which classifies residues based on
/// hydrogen bond energy patterns into helices (alpha, 3-10, pi), sheets, turns,
/// and bends — then collapses to the three canonical types (Helix, Sheet, Loop).
///
/// Currently a stub — returns Loop for all residues.
#[derive(Default)]
pub struct Dssp {
    pub params: DsspParams,
}

impl Dssp {
    pub fn new(params: DsspParams) -> Self {
        Self { params }
    }
}

impl SecondaryStructureAssigner for Dssp {
    fn assign(&self, residues: &[BackboneResidue]) -> Vec<SsType> {
        // TODO: Implement full DSSP algorithm
        vec![SsType::Loop; residues.len()]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dssp_stub_returns_all_loop() {
        let dssp = Dssp::default();
        let residues = vec![];
        assert!(dssp.assign(&residues).is_empty());
    }
}
