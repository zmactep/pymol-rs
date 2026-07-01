//! Free-energy scoring for design: binding (affinity) and folding (stability).
//!
//! Both metrics are built from single-point MM + implicit-solvent (GB/SA)
//! energies of atom subsets ([`crate::energy::energy_for_mask`]):
//!
//! - **Affinity:** `ΔG_bind = E(complex) − E(receptor) − E(ligand)`. The caller
//!   computes `ΔΔG_aff = ΔG_bind(mutant) − ΔG_bind(WT)`.
//! - **Stability:** the folded monomer energy referenced against an unfolded
//!   proxy (the mutated residue's energy in isolation). The caller computes
//!   `ΔΔG_stab = [E_fold(mut) − E_fold(WT)] − [E_iso(mut) − E_iso(WT)]`.
//!
//! This module only evaluates one state; the mutate/scan tools combine WT and
//! mutant states into the ΔΔG values.

use std::collections::HashSet;

use lin_alg::f32::Vec3;

use crate::energy::energy_for_mask;
use crate::topology::{EnergyBreakdown, EnergySettings, ParameterizedSystem};

/// Single-point energy of an atom subset.
pub fn total_energy(
    system: &ParameterizedSystem,
    coords: &[Vec3],
    mask: &HashSet<usize>,
    settings: &EnergySettings,
) -> EnergyBreakdown {
    energy_for_mask(system, coords, mask, settings)
}

/// MM/GBSA binding score for one state.
#[derive(Debug, Clone, Copy)]
pub struct BindingScore {
    pub complex: EnergyBreakdown,
    pub receptor: EnergyBreakdown,
    pub ligand: EnergyBreakdown,
}

impl BindingScore {
    /// `ΔG_bind = E(complex) − E(receptor) − E(ligand)`, in kJ/mol.
    pub fn delta_g(&self) -> f64 {
        self.complex.total() - self.receptor.total() - self.ligand.total()
    }
}

/// Computes the binding score for one conformation. `receptor` and `ligand` are
/// topology-index masks; the complex is their union.
pub fn binding_score(
    system: &ParameterizedSystem,
    coords: &[Vec3],
    receptor: &HashSet<usize>,
    ligand: &HashSet<usize>,
    settings: &EnergySettings,
) -> BindingScore {
    let complex: HashSet<usize> = receptor.union(ligand).copied().collect();
    BindingScore {
        complex: energy_for_mask(system, coords, &complex, settings),
        receptor: energy_for_mask(system, coords, receptor, settings),
        ligand: energy_for_mask(system, coords, ligand, settings),
    }
}

/// Folding-stability score for one state: the folded monomer energy and the
/// isolated mutated-residue energy (the unfolded reference proxy).
#[derive(Debug, Clone, Copy)]
pub struct StabilityScore {
    pub folded: EnergyBreakdown,
    pub isolated: EnergyBreakdown,
}

impl StabilityScore {
    /// Folded energy referenced against the unfolded (isolated-residue) proxy.
    pub fn reduced(&self) -> f64 {
        self.folded.total() - self.isolated.total()
    }
}

/// Computes the stability score for one conformation. `monomer` is the whole
/// single-chain mask; `residue` is the mutated residue's atoms (its isolated
/// energy stands in for the unfolded reference).
pub fn stability_score(
    system: &ParameterizedSystem,
    coords: &[Vec3],
    monomer: &HashSet<usize>,
    residue: &HashSet<usize>,
    settings: &EnergySettings,
) -> StabilityScore {
    StabilityScore {
        folded: energy_for_mask(system, coords, monomer, settings),
        isolated: energy_for_mask(system, coords, residue, settings),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn delta_g_is_complex_minus_parts() {
        let score = BindingScore {
            complex: EnergyBreakdown {
                coulomb: -100.0,
                ..Default::default()
            },
            receptor: EnergyBreakdown {
                coulomb: -40.0,
                ..Default::default()
            },
            ligand: EnergyBreakdown {
                coulomb: -10.0,
                ..Default::default()
            },
        };
        assert!((score.delta_g() - (-50.0)).abs() < 1e-9);
    }
}
