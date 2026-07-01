use std::collections::{HashMap, HashSet};

use lin_alg::f32::Vec3;

use crate::topology::{
    ordered_pair, CombinationRule, CoordinateSource, EnergyBreakdown, EnergySettings,
    EnergySummary, FrameEnergy, MoleculeSnapshot, ParameterizedSystem, RebuiltHydrogenSource,
    RunSelection,
};

const COULOMB_KJ_NM_PER_MOL: f64 = 138.935_458;
const ANGSTROM_PER_NM: f64 = 10.0;
const SASA_PROBE_A: f64 = 1.4;

#[derive(Debug, Clone)]
pub struct RunInput {
    pub topology: ParameterizedSystem,
    pub molecules: Vec<MoleculeSnapshot>,
    pub selection: RunSelection,
    pub start_frame: usize,
    pub end_frame: usize,
    pub interval: usize,
    pub settings: EnergySettings,
}

pub fn run(input: &RunInput) -> Result<EnergySummary, String> {
    validate_run_input(input)?;
    let molecule_map: HashMap<&str, &MoleculeSnapshot> = input
        .molecules
        .iter()
        .map(|mol| (mol.object.as_str(), mol))
        .collect();

    let complex_mask = input.selection.complex_mask();
    let mut frames = Vec::new();
    for frame in (input.start_frame..input.end_frame).step_by(input.interval.max(1)) {
        let coords = coordinates_for_frame(&input.topology, &molecule_map, frame)?;
        let complex = energy_for_mask(&input.topology, &coords, &complex_mask, &input.settings);
        let receptor = energy_for_mask(
            &input.topology,
            &coords,
            &input.selection.receptor,
            &input.settings,
        );
        let ligand = energy_for_mask(
            &input.topology,
            &coords,
            &input.selection.ligand,
            &input.settings,
        );
        frames.push(FrameEnergy {
            frame: frame + 1,
            complex,
            receptor,
            ligand,
            delta: complex - receptor - ligand,
        });
    }

    let mean_delta = mean_breakdown(frames.iter().map(|frame| frame.delta));
    let std_delta = std_breakdown(frames.iter().map(|frame| frame.delta), mean_delta);
    Ok(EnergySummary {
        frames,
        mean_delta,
        std_delta,
    })
}

fn validate_run_input(input: &RunInput) -> Result<(), String> {
    if input.selection.receptor.is_empty() {
        return Err("receptor selection has no parameterized atoms".into());
    }
    if input.selection.ligand.is_empty() {
        return Err("ligand selection has no parameterized atoms".into());
    }
    if !input
        .selection
        .receptor
        .is_disjoint(&input.selection.ligand)
    {
        return Err("receptor and ligand selections overlap".into());
    }
    if input.start_frame >= input.end_frame {
        return Err("frame range is empty".into());
    }
    if input.settings.solute_dielectric <= 0.0 || input.settings.solvent_dielectric <= 0.0 {
        return Err("dielectric constants must be positive".into());
    }
    Ok(())
}

/// Resolves coordinates for every topology atom at `frame`, rebuilding any
/// synthetic hydrogens from their control atoms. Shared with the minimizer.
pub fn coordinates_for_frame(
    topology: &ParameterizedSystem,
    molecule_map: &HashMap<&str, &MoleculeSnapshot>,
    frame: usize,
) -> Result<Vec<Vec3>, String> {
    let mut coords = vec![None; topology.atoms.len()];
    for (idx, atom) in topology.atoms.iter().enumerate() {
        let CoordinateSource::Loaded(key) = &atom.source else {
            continue;
        };
        let Some(molecule) = molecule_map.get(key.object.as_str()) else {
            return Err(format!("missing molecule object {}", key.object));
        };
        let Some(coord) = molecule.selected_coord(key.atom_index, frame) else {
            return Err(format!(
                "missing coordinates for {} atom {} at frame {}",
                key.object,
                key.atom_index + 1,
                frame + 1
            ));
        };
        coords[idx] = Some(coord);
    }

    let mut progress = true;
    while progress && coords.iter().any(Option::is_none) {
        progress = false;
        for (idx, atom) in topology.atoms.iter().enumerate() {
            if coords[idx].is_some() {
                continue;
            }
            let CoordinateSource::RebuiltHydrogen(source) = &atom.source else {
                continue;
            };
            let control_coords = source
                .controls
                .iter()
                .map(|&control| coords.get(control).copied().flatten())
                .collect::<Option<Vec<_>>>();
            let Some(control_coords) = control_coords else {
                continue;
            };
            coords[idx] = Some(rebuilt_hydrogen_coord(source, &control_coords));
            progress = true;
        }
    }

    coords
        .into_iter()
        .enumerate()
        .map(|(idx, coord)| {
            coord.ok_or_else(|| {
                let atom = &topology.atoms[idx];
                format!(
                    "could not rebuild coordinates for {}:{} at frame {}",
                    atom.residue_key,
                    atom.atom_name,
                    frame + 1
                )
            })
        })
        .collect()
}

/// Single-point MM + GB/SA energy over the atoms in `mask`. Shared with scoring.
pub fn energy_for_mask(
    topology: &ParameterizedSystem,
    coords: &[Vec3],
    mask: &HashSet<usize>,
    settings: &EnergySettings,
) -> EnergyBreakdown {
    let bonded = bonded_energy(topology, coords, mask);
    let (coulomb, lj) = nonbonded_energy(topology, coords, mask);
    let gb = gb_energy(topology, coords, mask, settings);
    let sa = sa_energy(topology, coords, mask, settings);
    EnergyBreakdown {
        bonded,
        coulomb,
        lj,
        gb,
        sa,
    }
}

fn bonded_energy(topology: &ParameterizedSystem, coords: &[Vec3], mask: &HashSet<usize>) -> f64 {
    let bond_energy = topology
        .bonds
        .iter()
        .filter(|bond| mask.contains(&bond.a) && mask.contains(&bond.b))
        .filter_map(|bond| {
            let (Some(length), Some(force)) = (bond.length_nm, bond.force_kj_mol_nm2) else {
                return None;
            };
            let r_nm = distance_a(coords[bond.a], coords[bond.b]) / ANGSTROM_PER_NM;
            Some(0.5 * force * (r_nm - length).powi(2))
        })
        .sum::<f64>();

    let angle_energy = topology
        .angles
        .iter()
        .filter(|angle| {
            mask.contains(&angle.a) && mask.contains(&angle.b) && mask.contains(&angle.c)
        })
        .filter_map(|angle| {
            let (Some(theta0), Some(force)) = (angle.angle_deg, angle.force_kj_mol_rad2) else {
                return None;
            };
            let theta = angle_rad(coords[angle.a], coords[angle.b], coords[angle.c]);
            Some(0.5 * force * (theta - theta0.to_radians()).powi(2))
        })
        .sum::<f64>();

    bond_energy + angle_energy
}

fn nonbonded_energy(
    topology: &ParameterizedSystem,
    coords: &[Vec3],
    mask: &HashSet<usize>,
) -> (f64, f64) {
    let atoms: Vec<usize> = mask.iter().copied().collect();
    let mut coulomb = 0.0;
    let mut lj = 0.0;

    for i in 0..atoms.len() {
        for j in (i + 1)..atoms.len() {
            let a = atoms[i];
            let b = atoms[j];
            let pair = ordered_pair(a, b);
            if topology.exclusions.contains(&pair) {
                continue;
            }
            let atom_a = &topology.atoms[a];
            let atom_b = &topology.atoms[b];
            let r_a = distance_a(coords[a], coords[b]).max(0.001);
            let r_nm = r_a / ANGSTROM_PER_NM;
            let scale = if topology.one_four_pairs.contains(&pair) {
                topology.defaults.fudge_qq
            } else {
                1.0
            };
            coulomb += scale * COULOMB_KJ_NM_PER_MOL * atom_a.charge * atom_b.charge / r_nm;

            let (sigma_a, epsilon) = combine_lj(
                topology.defaults.combination_rule,
                atom_a.sigma_a,
                atom_a.epsilon_kj,
                atom_b.sigma_a,
                atom_b.epsilon_kj,
            );
            let lj_scale = if topology.one_four_pairs.contains(&pair) {
                topology.defaults.fudge_lj
            } else {
                1.0
            };
            if sigma_a > 0.0 && epsilon > 0.0 {
                let sr6 = (sigma_a / r_a).powi(6);
                lj += lj_scale * 4.0 * epsilon * (sr6 * sr6 - sr6);
            }
        }
    }

    (coulomb, lj)
}

fn gb_energy(
    topology: &ParameterizedSystem,
    coords: &[Vec3],
    mask: &HashSet<usize>,
    settings: &EnergySettings,
) -> f64 {
    let atoms: Vec<usize> = mask.iter().copied().collect();
    let dielectric = 1.0 / settings.solute_dielectric - 1.0 / settings.solvent_dielectric;
    let salt_factor = (1.0 + settings.salt_m.max(0.0)).sqrt();
    let mut sum = 0.0;

    for &i in &atoms {
        for &j in &atoms {
            let atom_i = &topology.atoms[i];
            let atom_j = &topology.atoms[j];
            let ri_nm = atom_i.gb_radius_a.max(0.5) / ANGSTROM_PER_NM;
            let rj_nm = atom_j.gb_radius_a.max(0.5) / ANGSTROM_PER_NM;
            let r_nm = if i == j {
                0.0
            } else {
                distance_a(coords[i], coords[j]) / ANGSTROM_PER_NM
            };
            let born = (ri_nm * rj_nm).sqrt();
            let fgb = (r_nm * r_nm + born * born * (-r_nm * r_nm / (4.0 * born * born)).exp())
                .sqrt()
                .max(0.001);
            sum += atom_i.charge * atom_j.charge / fgb;
        }
    }

    -0.5 * COULOMB_KJ_NM_PER_MOL * dielectric * sum / salt_factor
}

fn sa_energy(
    topology: &ParameterizedSystem,
    coords: &[Vec3],
    mask: &HashSet<usize>,
    settings: &EnergySettings,
) -> f64 {
    let mut total_area = 0.0;
    for &i in mask {
        let radius = topology.atoms[i].gb_radius_a + SASA_PROBE_A;
        let mut exposure = 1.0;
        for &j in mask {
            if i == j {
                continue;
            }
            let neighbor_radius = topology.atoms[j].gb_radius_a + SASA_PROBE_A;
            let dist = distance_a(coords[i], coords[j]);
            let overlap =
                ((radius + neighbor_radius - dist) / (2.0 * radius)).clamp(0.0, 1.0) * 0.35;
            exposure *= 1.0 - overlap;
        }
        total_area += 4.0 * std::f64::consts::PI * radius * radius * exposure.max(0.05);
    }
    settings.sa_gamma * total_area + settings.sa_offset
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

fn distance_a(a: Vec3, b: Vec3) -> f64 {
    let dx = f64::from(a.x - b.x);
    let dy = f64::from(a.y - b.y);
    let dz = f64::from(a.z - b.z);
    (dx * dx + dy * dy + dz * dz).sqrt()
}

fn angle_rad(a: Vec3, b: Vec3, c: Vec3) -> f64 {
    let ab = [
        f64::from(a.x - b.x),
        f64::from(a.y - b.y),
        f64::from(a.z - b.z),
    ];
    let cb = [
        f64::from(c.x - b.x),
        f64::from(c.y - b.y),
        f64::from(c.z - b.z),
    ];
    let dot = ab[0] * cb[0] + ab[1] * cb[1] + ab[2] * cb[2];
    let len_ab = (ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2]).sqrt();
    let len_cb = (cb[0] * cb[0] + cb[1] * cb[1] + cb[2] * cb[2]).sqrt();
    (dot / (len_ab * len_cb).max(0.000_001))
        .clamp(-1.0, 1.0)
        .acos()
}

fn rebuilt_hydrogen_coord(source: &RebuiltHydrogenSource, controls: &[Vec3]) -> Vec3 {
    let parent = controls[0];
    let base = rebuilt_base_direction(parent, &controls[1..]);
    let (u, v) = orthonormal_basis(base);
    let dir = rebuilt_direction(base, u, v, source);
    vec3_add(parent, vec3_scale(dir, source.bond_length_a))
}

fn rebuilt_base_direction(parent: Vec3, controls: &[Vec3]) -> Vec3 {
    if controls.is_empty() {
        return Vec3::new(1.0, 0.0, 0.0);
    }
    let mut centroid = Vec3::new(0.0, 0.0, 0.0);
    for &coord in controls {
        centroid = vec3_add(centroid, coord);
    }
    centroid = vec3_scale(centroid, 1.0 / controls.len() as f32);
    normalize_or(vec3_sub(parent, centroid), Vec3::new(1.0, 0.0, 0.0))
}

fn rebuilt_direction(base: Vec3, u: Vec3, v: Vec3, source: &RebuiltHydrogenSource) -> Vec3 {
    if source.count <= 1 {
        return base;
    }
    let slot = source.slot.min(source.count - 1);
    let phi = 2.0 * std::f32::consts::PI * slot as f32 / source.count as f32;
    let cone_cos = match (source.method, source.count) {
        (_, 2) => 0.577_350_26,
        (4 | 10 | 11, _) => 0.333_333_34,
        (_, 3..) => 0.333_333_34,
        _ => 1.0,
    };
    let cone_sin = (1.0_f32 - cone_cos * cone_cos).max(0.0).sqrt();
    let around = vec3_add(vec3_scale(u, phi.cos()), vec3_scale(v, phi.sin()));
    normalize_or(
        vec3_add(vec3_scale(base, cone_cos), vec3_scale(around, cone_sin)),
        base,
    )
}

fn orthonormal_basis(axis: Vec3) -> (Vec3, Vec3) {
    let helper = if axis.x.abs() < 0.8 {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let u = normalize_or(cross(axis, helper), Vec3::new(0.0, 0.0, 1.0));
    let v = normalize_or(cross(axis, u), Vec3::new(0.0, 1.0, 0.0));
    (u, v)
}

fn vec3_add(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(a.x + b.x, a.y + b.y, a.z + b.z)
}

fn vec3_sub(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(a.x - b.x, a.y - b.y, a.z - b.z)
}

fn vec3_scale(v: Vec3, scale: f32) -> Vec3 {
    Vec3::new(v.x * scale, v.y * scale, v.z * scale)
}

fn cross(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
    )
}

fn normalize_or(v: Vec3, fallback: Vec3) -> Vec3 {
    let len = (v.x * v.x + v.y * v.y + v.z * v.z).sqrt();
    if len > 0.000_001 {
        return vec3_scale(v, 1.0 / len);
    }
    fallback
}

fn mean_breakdown(values: impl Iterator<Item = EnergyBreakdown>) -> EnergyBreakdown {
    let mut count = 0.0;
    let mut sum = EnergyBreakdown::default();
    for value in values {
        count += 1.0;
        sum.bonded += value.bonded;
        sum.coulomb += value.coulomb;
        sum.lj += value.lj;
        sum.gb += value.gb;
        sum.sa += value.sa;
    }
    if count == 0.0 {
        return sum;
    }
    EnergyBreakdown {
        bonded: sum.bonded / count,
        coulomb: sum.coulomb / count,
        lj: sum.lj / count,
        gb: sum.gb / count,
        sa: sum.sa / count,
    }
}

fn std_breakdown(
    values: impl Iterator<Item = EnergyBreakdown>,
    mean: EnergyBreakdown,
) -> EnergyBreakdown {
    let values: Vec<_> = values.collect();
    let count = values.len() as f64;
    if count <= 1.0 {
        return EnergyBreakdown::default();
    }
    let mut sum = EnergyBreakdown::default();
    for value in values {
        sum.bonded += (value.bonded - mean.bonded).powi(2);
        sum.coulomb += (value.coulomb - mean.coulomb).powi(2);
        sum.lj += (value.lj - mean.lj).powi(2);
        sum.gb += (value.gb - mean.gb).powi(2);
        sum.sa += (value.sa - mean.sa).powi(2);
    }
    let denom = count - 1.0;
    EnergyBreakdown {
        bonded: (sum.bonded / denom).sqrt(),
        coulomb: (sum.coulomb / denom).sqrt(),
        lj: (sum.lj / denom).sqrt(),
        gb: (sum.gb / denom).sqrt(),
        sa: (sum.sa / denom).sqrt(),
    }
}

#[cfg(test)]
mod tests {
    use std::collections::{HashMap, HashSet};

    use lin_alg::f32::Vec3;

    use super::*;
    use crate::topology::{
        AtomKey, CoordinateSource, EnergySettings, ForceFieldDefaults, ParameterizedAtom,
        ParameterizedSystem, RebuiltHydrogenSource,
    };

    fn two_atom_topology() -> ParameterizedSystem {
        let atoms = vec![
            ParameterizedAtom {
                key: AtomKey {
                    object: "mol".into(),
                    atom_index: 0,
                },
                source: CoordinateSource::Loaded(AtomKey {
                    object: "mol".into(),
                    atom_index: 0,
                }),
                synthetic_parent: None,
                residue_key: "A".into(),
                residue_name: "REC".into(),
                atom_name: "A".into(),
                atom_type: "A".into(),
                charge: 1.0,
                mass: 1.0,
                sigma_a: 1.0,
                epsilon_kj: 0.0,
                gb_radius_a: 1.5,
            },
            ParameterizedAtom {
                key: AtomKey {
                    object: "mol".into(),
                    atom_index: 1,
                },
                source: CoordinateSource::Loaded(AtomKey {
                    object: "mol".into(),
                    atom_index: 1,
                }),
                synthetic_parent: None,
                residue_key: "B".into(),
                residue_name: "LIG".into(),
                atom_name: "B".into(),
                atom_type: "B".into(),
                charge: -1.0,
                mass: 1.0,
                sigma_a: 1.0,
                epsilon_kj: 0.0,
                gb_radius_a: 1.5,
            },
        ];
        ParameterizedSystem {
            ff_name: "tiny.ff".into(),
            source: "tiny.ff".into(),
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

    #[test]
    fn coulomb_delta_matches_two_atom_cross_term() {
        let topology = two_atom_topology();
        let coords = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(10.0, 0.0, 0.0)];
        let receptor = HashSet::from([0]);
        let ligand = HashSet::from([1]);
        let complex = HashSet::from([0, 1]);
        let settings = EnergySettings {
            sa_gamma: 0.0,
            ..EnergySettings::default()
        };

        let complex_e = energy_for_mask(&topology, &coords, &complex, &settings);
        let receptor_e = energy_for_mask(&topology, &coords, &receptor, &settings);
        let ligand_e = energy_for_mask(&topology, &coords, &ligand, &settings);
        let delta = complex_e - receptor_e - ligand_e;

        assert!((delta.coulomb + COULOMB_KJ_NM_PER_MOL).abs() < 1e-6);
    }

    #[test]
    fn rebuilt_hydrogen_coordinate_uses_requested_bond_length() {
        let source = RebuiltHydrogenSource {
            controls: vec![0],
            bond_length_a: 1.09,
            slot: 0,
            count: 1,
            method: 1,
        };
        let coord = rebuilt_hydrogen_coord(&source, &[Vec3::new(2.0, 0.0, 0.0)]);

        assert!((distance_a(coord, Vec3::new(2.0, 0.0, 0.0)) - 1.09).abs() < 1e-6);
    }
}
