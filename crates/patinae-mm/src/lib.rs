//! Molecular-mechanics engine for Patinae.
//!
//! This crate is the compute core shared by molecular-design tooling: it parses
//! GROMACS force-field directories, parametrizes loaded structures with
//! RTP-strict residue templates, and evaluates MM + implicit-solvent (GB/SA)
//! energies. Higher milestones add analytic gradients, energy minimization, and
//! structure-building (hydrogens, side chains, caps, mutations).
//!
//! The crate is deliberately free of any plugin/FFI or command-layer
//! dependencies so the engine can be unit-tested in isolation.

pub mod build;
pub mod cap;
pub mod energy;
pub mod hbond;
pub mod forcefields;
pub mod gradient;
pub mod gromacs;
pub mod minimize;
pub mod mutate;
pub mod parametrize;
pub mod protonation;
pub mod residue_geometry;
pub mod rotamers;
pub mod score;
pub mod sidechain;
pub mod topology;

pub use build::{hydrogens_to_add, HydrogenAddition};
pub use cap::{build_ace, build_nme, CapAtom, CapBuild};
pub use energy::{run as run_energy, RunInput};
pub use forcefields::{resolve_force_field_path, DEFAULT_FORCE_FIELD, FORCE_FIELD_HINT};
pub use hbond::{optimize_polar_hydrogens, PolarH};
pub use minimize::{minimize, Algorithm, MinimizeOptions, MinimizeResult};
pub use mutate::build_mutant;
pub use parametrize::{parameterize, parameterize_with, ParameterizeInput};
pub use protonation::{protonation_variants, variant_for};
pub use residue_geometry::{ideal_side_chain, ideal_side_chain_on, is_supported, IdealSideChain};
pub use rotamers::{optimize_rotamer, SideAtom};
pub use score::{binding_score, stability_score, total_energy, BindingScore, StabilityScore};
pub use sidechain::{plan_side_chains, ResidueBuild, SideChainPlan};
pub use topology::{
    AtomKey, EnergySettings, EnergySummary, MoleculeSnapshot, ParameterizedSystem, RunSelection,
    SelectedMolecule,
};
