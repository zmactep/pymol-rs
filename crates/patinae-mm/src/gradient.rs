//! Analytic energy and forces for the molecular-mechanics potential.
//!
//! This mirrors the energy terms in [`crate::energy`] but additionally returns
//! the negative gradient (force) on every atom, so the potential can drive an
//! energy minimizer ([`crate::minimize`]).
//!
//! Terms covered: harmonic bonds, harmonic angles, Coulomb, Lennard-Jones, and
//! the Generalized-Born polar solvation energy. The implicit-solvent
//! surface-area (SA) term is a soft heuristic and is intentionally excluded
//! from the force so that the analytic gradient is exactly the derivative of
//! [`potential_energy`] — this is what the finite-difference tests check.
//!
//! Positions are in ångström (matching loaded coordinates); energies are in
//! kJ/mol; forces are in kJ/mol/Å.

use std::collections::HashSet;

use crate::topology::{ordered_pair, CombinationRule, EnergySettings, ParameterizedSystem};

const COULOMB_KJ_NM_PER_MOL: f64 = 138.935_458;
const ANGSTROM_PER_NM: f64 = 10.0;

/// A 3-component force/position accumulator in f64.
pub type Vec3d = [f64; 3];

fn sub(a: Vec3d, b: Vec3d) -> Vec3d {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn norm(a: Vec3d) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}
fn dot(a: Vec3d, b: Vec3d) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn scale(a: Vec3d, s: f64) -> Vec3d {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// A harmonic position restraint pinning one atom toward a reference point.
#[derive(Debug, Clone, Copy)]
pub struct PositionRestraint {
    pub atom: usize,
    pub reference: Vec3d,
    /// Force constant in kJ/mol/Å².
    pub k: f64,
}

/// Evaluates the potential energy of the masked subsystem (no forces).
///
/// Sums exactly the same terms as [`energy_and_forces`], which makes it the
/// reference for finite-difference gradient checks.
pub fn potential_energy(
    system: &ParameterizedSystem,
    coords: &[Vec3d],
    mask: &HashSet<usize>,
    settings: &EnergySettings,
    restraints: &[PositionRestraint],
) -> f64 {
    let (energy, _) = energy_and_forces(system, coords, mask, settings, restraints);
    energy
}

/// Evaluates the potential energy and the force (negative gradient) on every
/// atom of the masked subsystem.
///
/// Forces on atoms outside `mask` are left at zero. Interactions are summed
/// only between atoms that are both in `mask`, matching
/// [`crate::energy::run`]'s `energy_for_mask` semantics.
pub fn energy_and_forces(
    system: &ParameterizedSystem,
    coords: &[Vec3d],
    mask: &HashSet<usize>,
    settings: &EnergySettings,
    restraints: &[PositionRestraint],
) -> (f64, Vec<Vec3d>) {
    let mut energy = 0.0;
    let mut forces = vec![[0.0; 3]; system.atoms.len()];

    energy += bonds(system, coords, mask, &mut forces);
    energy += angles(system, coords, mask, &mut forces);
    energy += nonbonded(system, coords, mask, &mut forces);
    energy += generalized_born(system, coords, mask, settings, &mut forces);
    energy += restraint_energy(coords, restraints, &mut forces);

    (energy, forces)
}

fn bonds(
    system: &ParameterizedSystem,
    coords: &[Vec3d],
    mask: &HashSet<usize>,
    forces: &mut [Vec3d],
) -> f64 {
    let mut energy = 0.0;
    for bond in &system.bonds {
        if !mask.contains(&bond.a) || !mask.contains(&bond.b) {
            continue;
        }
        let (Some(length), Some(force_k)) = (bond.length_nm, bond.force_kj_mol_nm2) else {
            continue;
        };
        let u = sub(coords[bond.a], coords[bond.b]);
        let r = norm(u).max(1e-6);
        let r_nm = r / ANGSTROM_PER_NM;
        let diff = r_nm - length;
        energy += 0.5 * force_k * diff * diff;
        // dE/dr (per Å) = force_k * diff * (1/10)
        let de_dr = force_k * diff / ANGSTROM_PER_NM;
        let dir = scale(u, 1.0 / r);
        accumulate_pair(forces, bond.a, bond.b, dir, de_dr);
    }
    energy
}

fn angles(
    system: &ParameterizedSystem,
    coords: &[Vec3d],
    mask: &HashSet<usize>,
    forces: &mut [Vec3d],
) -> f64 {
    let mut energy = 0.0;
    for angle in &system.angles {
        if !mask.contains(&angle.a) || !mask.contains(&angle.b) || !mask.contains(&angle.c) {
            continue;
        }
        let (Some(theta0_deg), Some(force_k)) = (angle.angle_deg, angle.force_kj_mol_rad2) else {
            continue;
        };
        let p = sub(coords[angle.a], coords[angle.b]);
        let q = sub(coords[angle.c], coords[angle.b]);
        let lp = norm(p).max(1e-6);
        let lq = norm(q).max(1e-6);
        let cos_t = (dot(p, q) / (lp * lq)).clamp(-1.0, 1.0);
        let theta = cos_t.acos();
        let diff = theta - theta0_deg.to_radians();
        energy += 0.5 * force_k * diff * diff;

        let sin_t = (1.0 - cos_t * cos_t).sqrt();
        if sin_t < 1e-6 {
            continue;
        }
        let de_dtheta = force_k * diff;
        let p_hat = scale(p, 1.0 / lp);
        let q_hat = scale(q, 1.0 / lq);
        // dtheta/dp = -(q_hat - cos_t p_hat)/(lp sin_t); force_a = -dE/dtheta * dtheta/dp
        let dtheta_dp = scale(sub(q_hat, scale(p_hat, cos_t)), -1.0 / (lp * sin_t));
        let dtheta_dq = scale(sub(p_hat, scale(q_hat, cos_t)), -1.0 / (lq * sin_t));
        let fa = scale(dtheta_dp, -de_dtheta);
        let fc = scale(dtheta_dq, -de_dtheta);
        add_force(forces, angle.a, fa);
        add_force(forces, angle.c, fc);
        add_force(forces, angle.b, scale(add(fa, fc), -1.0));
    }
    energy
}

fn nonbonded(
    system: &ParameterizedSystem,
    coords: &[Vec3d],
    mask: &HashSet<usize>,
    forces: &mut [Vec3d],
) -> f64 {
    let atoms: Vec<usize> = mask.iter().copied().collect();
    let mut energy = 0.0;
    for i in 0..atoms.len() {
        for j in (i + 1)..atoms.len() {
            let a = atoms[i];
            let b = atoms[j];
            let pair = ordered_pair(a, b);
            if system.exclusions.contains(&pair) {
                continue;
            }
            let atom_a = &system.atoms[a];
            let atom_b = &system.atoms[b];
            let u = sub(coords[a], coords[b]);
            let r = norm(u).max(0.001);
            let dir = scale(u, 1.0 / r);
            let is_14 = system.one_four_pairs.contains(&pair);

            // Coulomb: E = qq_scale * C * qa qb / (r/10)
            let qq_scale = if is_14 { system.defaults.fudge_qq } else { 1.0 };
            let q_prod = atom_a.charge * atom_b.charge;
            if q_prod != 0.0 {
                let r_nm = r / ANGSTROM_PER_NM;
                energy += qq_scale * COULOMB_KJ_NM_PER_MOL * q_prod / r_nm;
                // dE/dr (per Å) = -qq_scale C q q / (r_nm^2) * (1/10)
                let de_dr =
                    -qq_scale * COULOMB_KJ_NM_PER_MOL * q_prod / (r_nm * r_nm) / ANGSTROM_PER_NM;
                accumulate_pair(forces, a, b, dir, de_dr);
            }

            // Lennard-Jones
            let (sigma_a, epsilon) = combine_lj(
                system.defaults.combination_rule,
                atom_a.sigma_a,
                atom_a.epsilon_kj,
                atom_b.sigma_a,
                atom_b.epsilon_kj,
            );
            if sigma_a > 0.0 && epsilon > 0.0 {
                let lj_scale = if is_14 { system.defaults.fudge_lj } else { 1.0 };
                let sr6 = (sigma_a / r).powi(6);
                let sr12 = sr6 * sr6;
                energy += lj_scale * 4.0 * epsilon * (sr12 - sr6);
                // dE/dr = lj_scale * 4 eps * (-12 sr12 + 6 sr6) / r
                let de_dr = lj_scale * 4.0 * epsilon * (-12.0 * sr12 + 6.0 * sr6) / r;
                accumulate_pair(forces, a, b, dir, de_dr);
            }
        }
    }
    energy
}

fn generalized_born(
    system: &ParameterizedSystem,
    coords: &[Vec3d],
    mask: &HashSet<usize>,
    settings: &EnergySettings,
    forces: &mut [Vec3d],
) -> f64 {
    let atoms: Vec<usize> = mask.iter().copied().collect();
    let dielectric = 1.0 / settings.solute_dielectric - 1.0 / settings.solvent_dielectric;
    let salt_factor = (1.0 + settings.salt_m.max(0.0)).sqrt();
    let prefactor = -0.5 * COULOMB_KJ_NM_PER_MOL * dielectric / salt_factor;
    let mut energy = 0.0;

    // Self terms (i == i): r = 0, fgb = born = radius_i. Force-free.
    for &i in &atoms {
        let atom_i = &system.atoms[i];
        let ri_nm = atom_i.gb_radius_a.max(0.5) / ANGSTROM_PER_NM;
        energy += prefactor * atom_i.charge * atom_i.charge / ri_nm;
    }

    // Cross terms: each unordered pair is counted twice in the energy sum.
    for i in 0..atoms.len() {
        for j in (i + 1)..atoms.len() {
            let a = atoms[i];
            let b = atoms[j];
            let atom_a = &system.atoms[a];
            let atom_b = &system.atoms[b];
            let q_prod = atom_a.charge * atom_b.charge;
            let ra_nm = atom_a.gb_radius_a.max(0.5) / ANGSTROM_PER_NM;
            let rb_nm = atom_b.gb_radius_a.max(0.5) / ANGSTROM_PER_NM;
            let born = (ra_nm * rb_nm).sqrt();

            let u = sub(coords[a], coords[b]);
            let r_a = norm(u).max(1e-6);
            let r_nm = r_a / ANGSTROM_PER_NM;
            let r2 = r_nm * r_nm;
            let b2 = born * born;
            let exp_term = (-r2 / (4.0 * b2)).exp();
            let d = r2 + b2 * exp_term;
            let fgb = d.sqrt().max(0.001);

            // Factor 2: both ordered (a,b) and (b,a) appear in the energy sum.
            energy += 2.0 * prefactor * q_prod / fgb;

            // dD/dr_nm = 2r - (r/2) exp_term ; dfgb/dr_nm = dD/dr_nm / (2 fgb)
            let dd_dr = 2.0 * r_nm - 0.5 * r_nm * exp_term;
            let dfgb_dr = dd_dr / (2.0 * fgb);
            // dE/dr_nm = 2 prefactor q_prod * (-1/fgb^2) * dfgb_dr ; convert to per Å
            let de_dr_nm = 2.0 * prefactor * q_prod * (-1.0 / (fgb * fgb)) * dfgb_dr;
            let de_dr = de_dr_nm / ANGSTROM_PER_NM;
            let dir = scale(u, 1.0 / r_a);
            accumulate_pair(forces, a, b, dir, de_dr);
        }
    }
    energy
}

fn restraint_energy(
    coords: &[Vec3d],
    restraints: &[PositionRestraint],
    forces: &mut [Vec3d],
) -> f64 {
    let mut energy = 0.0;
    for r in restraints {
        let u = sub(coords[r.atom], r.reference);
        energy += 0.5 * r.k * dot(u, u);
        // F = -k u
        add_force(forces, r.atom, scale(u, -r.k));
    }
    energy
}

fn combine_lj(
    rule: CombinationRule,
    sigma_a: f64,
    eps_a: f64,
    sigma_b: f64,
    eps_b: f64,
) -> (f64, f64) {
    let epsilon = (eps_a * eps_b).sqrt();
    let sigma = match rule {
        CombinationRule::Geometric => (sigma_a * sigma_b).sqrt(),
        CombinationRule::C6C12 | CombinationRule::LorentzBerthelot => (sigma_a + sigma_b) * 0.5,
    };
    (sigma, epsilon)
}

fn add(a: Vec3d, b: Vec3d) -> Vec3d {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn add_force(forces: &mut [Vec3d], idx: usize, f: Vec3d) {
    forces[idx][0] += f[0];
    forces[idx][1] += f[1];
    forces[idx][2] += f[2];
}

/// Applies `F = -dE/dr * dir` to atom `a` and the opposite to atom `b`, where
/// `dir` is the unit vector from `b` to `a`.
fn accumulate_pair(forces: &mut [Vec3d], a: usize, b: usize, dir: Vec3d, de_dr: f64) {
    let fa = scale(dir, -de_dr);
    add_force(forces, a, fa);
    add_force(forces, b, scale(fa, -1.0));
}

/// Largest per-atom force magnitude over the masked atoms, in kJ/mol/Å.
pub fn max_force(forces: &[Vec3d], mask: &HashSet<usize>) -> f64 {
    mask.iter()
        .map(|&i| norm(forces[i]))
        .fold(0.0_f64, f64::max)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::topology::{
        AtomKey, CoordinateSource, ForceFieldDefaults, ParameterizedAngle, ParameterizedAtom,
        ParameterizedBond,
    };
    use std::collections::HashMap;

    fn atom(charge: f64, sigma_a: f64, epsilon_kj: f64, gb_radius_a: f64) -> ParameterizedAtom {
        ParameterizedAtom {
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
            atom_name: "X".into(),
            atom_type: "t".into(),
            charge,
            mass: 12.0,
            sigma_a,
            epsilon_kj,
            gb_radius_a,
        }
    }

    fn system(atoms: Vec<ParameterizedAtom>) -> ParameterizedSystem {
        ParameterizedSystem {
            ff_name: "test".into(),
            source: "test".into(),
            selection: "all".into(),
            defaults: ForceFieldDefaults::default(),
            atoms,
            bonds: Vec::new(),
            angles: Vec::new(),
            exclusions: HashSet::new(),
            one_four_pairs: HashSet::new(),
            key_to_index: HashMap::new(),
            report: String::new(),
        }
    }

    fn full_mask(n: usize) -> HashSet<usize> {
        (0..n).collect()
    }

    /// Central finite-difference check: analytic force == -dE/dx numerically.
    fn check_gradient(sys: &ParameterizedSystem, coords: &[Vec3d], settings: &EnergySettings) {
        let mask = full_mask(sys.atoms.len());
        let (_, forces) = energy_and_forces(sys, coords, &mask, settings, &[]);
        let h = 1e-5;
        for i in 0..sys.atoms.len() {
            for d in 0..3 {
                let mut plus = coords.to_vec();
                let mut minus = coords.to_vec();
                plus[i][d] += h;
                minus[i][d] -= h;
                let e_plus = potential_energy(sys, &plus, &mask, settings, &[]);
                let e_minus = potential_energy(sys, &minus, &mask, settings, &[]);
                let numeric = -(e_plus - e_minus) / (2.0 * h);
                let analytic = forces[i][d];
                let tol = 1e-3 * (1.0 + analytic.abs());
                assert!(
                    (numeric - analytic).abs() < tol,
                    "atom {i} dim {d}: analytic {analytic} vs numeric {numeric}"
                );
            }
        }
    }

    #[test]
    fn bond_force_matches_finite_difference() {
        let mut sys = system(vec![atom(0.0, 0.0, 0.0, 1.5), atom(0.0, 0.0, 0.0, 1.5)]);
        sys.bonds.push(ParameterizedBond {
            a: 0,
            b: 1,
            length_nm: Some(0.15),
            force_kj_mol_nm2: Some(250_000.0),
        });
        let coords = vec![[0.0, 0.0, 0.0], [1.7, 0.0, 0.0]];
        check_gradient(&sys, &coords, &EnergySettings::default());
    }

    #[test]
    fn angle_force_matches_finite_difference() {
        let mut sys = system(vec![
            atom(0.0, 0.0, 0.0, 1.5),
            atom(0.0, 0.0, 0.0, 1.5),
            atom(0.0, 0.0, 0.0, 1.5),
        ]);
        sys.angles.push(ParameterizedAngle {
            a: 0,
            b: 1,
            c: 2,
            angle_deg: Some(109.5),
            force_kj_mol_rad2: Some(400.0),
        });
        let coords = vec![[1.0, 0.2, 0.0], [0.0, 0.0, 0.0], [-0.3, 1.0, 0.1]];
        check_gradient(&sys, &coords, &EnergySettings::default());
    }

    #[test]
    fn nonbonded_force_matches_finite_difference() {
        let sys = system(vec![atom(0.6, 3.2, 0.6, 1.7), atom(-0.5, 3.0, 0.5, 1.6)]);
        let coords = vec![[0.0, 0.0, 0.0], [3.4, 0.5, 0.2]];
        check_gradient(&sys, &coords, &EnergySettings::default());
    }

    #[test]
    fn generalized_born_force_matches_finite_difference() {
        let sys = system(vec![
            atom(0.8, 0.0, 0.0, 1.7),
            atom(-0.7, 0.0, 0.0, 1.5),
            atom(0.3, 0.0, 0.0, 2.0),
        ]);
        let coords = vec![[0.0, 0.0, 0.0], [2.5, 0.0, 0.0], [1.0, 2.2, 0.3]];
        check_gradient(&sys, &coords, &EnergySettings::default());
    }
}
