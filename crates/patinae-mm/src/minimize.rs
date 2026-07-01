//! Cartesian energy minimization in implicit solvent.
//!
//! Drives [`crate::gradient`] with steepest-descent or nonlinear conjugate
//! gradient and a backtracking line search. Backbone atoms can be position
//! restrained. Convergence is on the maximum per-atom force.

use std::collections::HashSet;

use lin_alg::f32::Vec3;

use crate::gradient::{energy_and_forces, max_force, PositionRestraint, Vec3d};
use crate::topology::{EnergySettings, ParameterizedSystem};

/// Minimization algorithm.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Algorithm {
    SteepestDescent,
    ConjugateGradient,
}

impl Algorithm {
    /// Parses a panel/command token into an algorithm (default: CG).
    pub fn from_token(token: &str) -> Self {
        match token.trim().to_ascii_lowercase().as_str() {
            "sd" | "steepest" | "steepest_descent" | "steepest descent" => Self::SteepestDescent,
            _ => Self::ConjugateGradient,
        }
    }
}

/// Minimization controls.
#[derive(Debug, Clone)]
pub struct MinimizeOptions {
    pub algorithm: Algorithm,
    pub max_steps: usize,
    /// Convergence threshold on the maximum force, in kJ/mol/nm.
    pub force_tol_kj_mol_nm: f64,
    /// When set, backbone heavy atoms (N, CA, C, O) are restrained to their
    /// starting positions with this force constant (kJ/mol/Å²).
    pub restrain_backbone_k: Option<f64>,
    pub settings: EnergySettings,
}

impl Default for MinimizeOptions {
    fn default() -> Self {
        Self {
            algorithm: Algorithm::ConjugateGradient,
            max_steps: 5000,
            force_tol_kj_mol_nm: 100.0,
            restrain_backbone_k: None,
            settings: EnergySettings::default(),
        }
    }
}

/// Outcome of a minimization run.
#[derive(Debug, Clone)]
pub struct MinimizeResult {
    pub steps: usize,
    pub energy_start: f64,
    pub energy_final: f64,
    /// Maximum per-atom force at the final step, in kJ/mol/nm.
    pub max_force_kj_mol_nm: f64,
    pub rmsd_a: f64,
    pub converged: bool,
}

const ANGSTROM_PER_NM: f64 = 10.0;

/// Minimizes the masked atoms of `system` starting from `coords`.
///
/// `mask` is both the set of atoms allowed to move *and* the interaction set
/// over which energy and forces are summed. Use this for whole-selection
/// minimization where the selection is the entire system of interest. For a
/// shell relaxation that must feel a fixed environment, use
/// [`minimize_with_environment`].
///
/// Returns the minimized coordinates (length == `system.atoms.len()`; atoms
/// outside `mask` are returned unchanged) together with a summary.
pub fn minimize(
    system: &ParameterizedSystem,
    coords: &[Vec3],
    mask: &HashSet<usize>,
    options: &MinimizeOptions,
) -> (Vec<Vec3>, MinimizeResult) {
    minimize_with_environment(system, coords, mask, mask, options)
}

/// Minimizes only the `movable` atoms while summing energy and forces over the
/// (larger) `energy_mask`. Atoms in `energy_mask \ movable` are held fixed but
/// still exert forces on the movable atoms — this is what lets a relaxing
/// side-chain shell feel clashes against the rigid protein around it.
///
/// `movable` should be a subset of `energy_mask`; any movable atom not in
/// `energy_mask` would move without contributing to (or feeling) the energy.
/// Convergence and RMSD are measured over `movable` only.
pub fn minimize_with_environment(
    system: &ParameterizedSystem,
    coords: &[Vec3],
    energy_mask: &HashSet<usize>,
    movable: &HashSet<usize>,
    options: &MinimizeOptions,
) -> (Vec<Vec3>, MinimizeResult) {
    let mut pos: Vec<Vec3d> = coords
        .iter()
        .map(|c| [f64::from(c.x), f64::from(c.y), f64::from(c.z)])
        .collect();
    let start = pos.clone();

    let restraints = build_restraints(system, &pos, movable, options.restrain_backbone_k);
    let settings = options.settings;
    let force_tol_a = options.force_tol_kj_mol_nm / ANGSTROM_PER_NM;

    let mask = energy_mask;
    let movable: Vec<usize> = movable.iter().copied().collect();
    let movable_set: HashSet<usize> = movable.iter().copied().collect();

    let (mut energy, mut forces) = energy_and_forces(system, &pos, mask, &settings, &restraints);
    let energy_start = energy;
    let mut gmax = max_force(&forces, &movable_set);

    // Search direction (CG); initialized to steepest descent.
    let mut dir = forces.clone();
    let mut prev_forces = forces.clone();
    let mut step_len = 0.01; // Å scaling for the trial step

    let mut steps = 0;
    let mut converged = gmax <= force_tol_a;
    while steps < options.max_steps && !converged {
        // Direction for this iteration.
        if options.algorithm == Algorithm::SteepestDescent || steps == 0 {
            dir = forces.clone();
        } else {
            let beta = polak_ribiere(&forces, &prev_forces, &movable).max(0.0);
            for &i in &movable {
                for d in 0..3 {
                    dir[i][d] = forces[i][d] + beta * dir[i][d];
                }
            }
        }
        // Ensure the direction is a descent direction; otherwise reset.
        if dot_masked(&dir, &forces, &movable) <= 0.0 {
            dir = forces.clone();
        }

        let dir_max = movable
            .iter()
            .map(|&i| norm(dir[i]))
            .fold(0.0_f64, f64::max)
            .max(1e-12);

        // Backtracking line search along `dir`.
        let mut alpha = step_len / dir_max;
        let mut accepted = false;
        for _ in 0..20 {
            let trial = step(&pos, &dir, &movable, alpha);
            let (e_trial, f_trial) =
                energy_and_forces(system, &trial, mask, &settings, &restraints);
            if e_trial < energy {
                pos = trial;
                energy = e_trial;
                prev_forces = std::mem::replace(&mut forces, f_trial);
                accepted = true;
                step_len = (alpha * dir_max * 1.2).min(0.5);
                break;
            }
            alpha *= 0.5;
        }
        if !accepted {
            // Stuck: shrink the trial step and stop if it underflows.
            step_len *= 0.5;
            if step_len < 1e-6 {
                break;
            }
            steps += 1;
            continue;
        }

        gmax = max_force(&forces, &movable_set);
        converged = gmax <= force_tol_a;
        steps += 1;
    }

    // Write minimized positions back into a full-length coordinate set.
    let mut out_coords = coords.to_vec();
    for &i in &movable {
        out_coords[i] = Vec3::new(pos[i][0] as f32, pos[i][1] as f32, pos[i][2] as f32);
    }

    let result = MinimizeResult {
        steps,
        energy_start,
        energy_final: energy,
        max_force_kj_mol_nm: gmax * ANGSTROM_PER_NM,
        rmsd_a: rmsd(&start, &pos, &movable_set),
        converged,
    };
    (out_coords, result)
}

fn build_restraints(
    system: &ParameterizedSystem,
    pos: &[Vec3d],
    mask: &HashSet<usize>,
    k: Option<f64>,
) -> Vec<PositionRestraint> {
    let Some(k) = k else {
        return Vec::new();
    };
    let mut restraints = Vec::new();
    for &i in mask {
        let name = system.atoms[i].atom_name.as_str();
        if matches!(name, "N" | "CA" | "C" | "O") {
            restraints.push(PositionRestraint {
                atom: i,
                reference: pos[i],
                k,
            });
        }
    }
    restraints
}

fn step(pos: &[Vec3d], dir: &[Vec3d], movable: &[usize], alpha: f64) -> Vec<Vec3d> {
    let mut out = pos.to_vec();
    for &i in movable {
        for d in 0..3 {
            out[i][d] += alpha * dir[i][d];
        }
    }
    out
}

fn polak_ribiere(forces: &[Vec3d], prev: &[Vec3d], movable: &[usize]) -> f64 {
    let mut num = 0.0;
    let mut den = 0.0;
    for &i in movable {
        for d in 0..3 {
            num += forces[i][d] * (forces[i][d] - prev[i][d]);
            den += prev[i][d] * prev[i][d];
        }
    }
    if den < 1e-12 {
        0.0
    } else {
        num / den
    }
}

fn dot_masked(a: &[Vec3d], b: &[Vec3d], movable: &[usize]) -> f64 {
    let mut s = 0.0;
    for &i in movable {
        for d in 0..3 {
            s += a[i][d] * b[i][d];
        }
    }
    s
}

fn norm(a: Vec3d) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

fn rmsd(a: &[Vec3d], b: &[Vec3d], mask: &HashSet<usize>) -> f64 {
    if mask.is_empty() {
        return 0.0;
    }
    let mut sum = 0.0;
    for &i in mask {
        let dx = a[i][0] - b[i][0];
        let dy = a[i][1] - b[i][1];
        let dz = a[i][2] - b[i][2];
        sum += dx * dx + dy * dy + dz * dz;
    }
    (sum / mask.len() as f64).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::topology::{
        AtomKey, CoordinateSource, ForceFieldDefaults, ParameterizedAtom, ParameterizedBond,
        ParameterizedSystem,
    };
    use std::collections::HashMap;

    fn two_atom_bonded_system() -> ParameterizedSystem {
        let mk = |name: &str| ParameterizedAtom {
            key: AtomKey {
                object: "m".into(),
                atom_index: 0,
            },
            source: CoordinateSource::Loaded(AtomKey {
                object: "m".into(),
                atom_index: 0,
            }),
            synthetic_parent: None,
            residue_key: "r".into(),
            residue_name: "ALA".into(),
            atom_name: name.into(),
            atom_type: "t".into(),
            charge: 0.0,
            mass: 12.0,
            sigma_a: 0.0,
            epsilon_kj: 0.0,
            gb_radius_a: 1.5,
        };
        ParameterizedSystem {
            ff_name: "test".into(),
            source: "test".into(),
            selection: "all".into(),
            defaults: ForceFieldDefaults::default(),
            atoms: vec![mk("CA"), mk("CB")],
            bonds: vec![ParameterizedBond {
                a: 0,
                b: 1,
                length_nm: Some(0.15),
                force_kj_mol_nm2: Some(250_000.0),
            }],
            angles: Vec::new(),
            exclusions: HashSet::new(),
            one_four_pairs: HashSet::new(),
            key_to_index: HashMap::new(),
            report: String::new(),
        }
    }

    #[test]
    fn environment_atom_stays_fixed_but_exerts_force() {
        // Two bonded atoms stretched to 2.0 Å (equilibrium 1.5 Å). Only atom 1
        // is movable; atom 0 is in the energy set but frozen. The bond should
        // still relax — atom 1 moves toward atom 0, which itself does not move.
        let sys = two_atom_bonded_system();
        let energy_mask: HashSet<usize> = [0, 1].into_iter().collect();
        let movable: HashSet<usize> = [1].into_iter().collect();
        let coords = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(2.0, 0.0, 0.0)];
        let opts = MinimizeOptions {
            algorithm: Algorithm::ConjugateGradient,
            max_steps: 2000,
            force_tol_kj_mol_nm: 10.0,
            restrain_backbone_k: None,
            settings: EnergySettings::default(),
        };
        let (out, result) = minimize_with_environment(&sys, &coords, &energy_mask, &movable, &opts);
        assert_eq!(out[0], coords[0], "frozen environment atom must not move");
        let r = (f64::from(out[1].x) - f64::from(out[0].x)).abs();
        assert!(result.energy_final < result.energy_start);
        assert!(
            (r - 1.5).abs() < 0.05,
            "bond length {r} should relax to ~1.5 Å"
        );
    }

    #[test]
    fn relaxes_a_stretched_bond_to_equilibrium() {
        let sys = two_atom_bonded_system();
        let mask: HashSet<usize> = [0, 1].into_iter().collect();
        // Start stretched to 2.0 Å; equilibrium is 1.5 Å (0.15 nm).
        let coords = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(2.0, 0.0, 0.0)];
        let opts = MinimizeOptions {
            algorithm: Algorithm::ConjugateGradient,
            max_steps: 2000,
            force_tol_kj_mol_nm: 10.0,
            restrain_backbone_k: None,
            settings: EnergySettings::default(),
        };
        let (out, result) = minimize(&sys, &coords, &mask, &opts);
        let r = (f64::from(out[1].x) - f64::from(out[0].x)).abs();
        assert!(result.energy_final < result.energy_start);
        assert!(
            (r - 1.5).abs() < 0.05,
            "bond length {r} should relax to ~1.5 Å"
        );
    }
}
