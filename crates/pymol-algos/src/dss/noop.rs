//! No-op secondary structure assignment — returns Loop for every residue.

use super::{BackboneResidue, SecondaryStructureAssigner, SsType};

/// No-op assigner that assigns [`SsType::Loop`] to all residues.
#[derive(Default)]
pub struct NoOp;

impl SecondaryStructureAssigner for NoOp {
    fn assign(&self, residues: &[BackboneResidue]) -> Vec<SsType> {
        vec![SsType::Loop; residues.len()]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn noop_returns_all_loop() {
        let noop = NoOp;
        assert!(noop.assign(&[]).is_empty());
    }
}
