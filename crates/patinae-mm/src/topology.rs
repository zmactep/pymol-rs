use std::collections::{HashMap, HashSet};
use std::fmt;

use lin_alg::f32::Vec3;
use patinae_mol::{AtomIndex, ObjectMolecule};
use patinae_select::SelectionResult;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CombinationRule {
    C6C12,
    LorentzBerthelot,
    Geometric,
}

impl CombinationRule {
    pub fn from_i32(value: i32) -> Self {
        match value {
            1 => Self::C6C12,
            3 => Self::Geometric,
            _ => Self::LorentzBerthelot,
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct ForceFieldDefaults {
    pub combination_rule: CombinationRule,
    pub fudge_lj: f64,
    pub fudge_qq: f64,
}

impl Default for ForceFieldDefaults {
    fn default() -> Self {
        Self {
            combination_rule: CombinationRule::LorentzBerthelot,
            fudge_lj: 1.0,
            fudge_qq: 1.0,
        }
    }
}

#[derive(Debug, Clone)]
pub struct AtomType {
    pub mass: f64,
    pub sigma_a: f64,
    pub epsilon_kj: f64,
}

#[derive(Debug, Clone)]
pub struct BondType {
    pub atoms: [String; 2],
    pub length_nm: f64,
    pub force_kj_mol_nm2: f64,
}

#[derive(Debug, Clone)]
pub struct AngleType {
    pub atoms: [String; 3],
    pub angle_deg: f64,
    pub force_kj_mol_rad2: f64,
}

#[derive(Debug, Clone)]
pub struct ResidueAtom {
    pub name: String,
    pub atom_type: String,
    pub charge: f64,
}

#[derive(Debug, Clone)]
pub struct ResidueTemplate {
    pub name: String,
    pub atoms: Vec<ResidueAtom>,
    pub bonds: Vec<[String; 2]>,
    pub exclusions: Vec<[String; 2]>,
    pub impropers: Vec<[String; 4]>,
}

#[derive(Debug, Clone)]
pub struct HydrogenRule {
    pub count: usize,
    pub method: i32,
    pub name: String,
    pub controls: Vec<String>,
}

impl HydrogenRule {
    pub fn expanded_names(&self) -> Vec<String> {
        if self.count <= 1 {
            return vec![self.name.clone()];
        }
        (1..=self.count)
            .map(|idx| format!("{}{idx}", self.name))
            .collect()
    }
}

impl ResidueTemplate {
    pub fn atom_names(&self) -> HashSet<&str> {
        self.atoms.iter().map(|atom| atom.name.as_str()).collect()
    }
}

#[derive(Debug, Clone, Default)]
pub struct TerminalPatchSummary {
    pub n_terminal_files: usize,
    pub c_terminal_files: usize,
    pub entries: Vec<String>,
}

#[derive(Debug, Clone)]
pub struct ForceField {
    pub name: String,
    pub defaults: ForceFieldDefaults,
    pub atom_types: HashMap<String, AtomType>,
    pub bond_types: Vec<BondType>,
    pub angle_types: Vec<AngleType>,
    pub residues: HashMap<String, ResidueTemplate>,
    pub hydrogen_rules: HashMap<String, Vec<HydrogenRule>>,
    pub terminal_patches: TerminalPatchSummary,
}

impl ForceField {
    pub fn atom_type(&self, name: &str) -> Option<&AtomType> {
        self.atom_types.get(name)
    }

    pub fn residue(&self, name: &str) -> Option<&ResidueTemplate> {
        self.residues.get(name)
    }

    pub fn bond_type(&self, a: &str, b: &str) -> Option<&BondType> {
        self.bond_types.iter().find(|ty| {
            (ty.atoms[0] == a && ty.atoms[1] == b) || (ty.atoms[0] == b && ty.atoms[1] == a)
        })
    }

    pub fn angle_type(&self, a: &str, b: &str, c: &str) -> Option<&AngleType> {
        self.angle_types.iter().find(|ty| {
            (ty.atoms[0] == a && ty.atoms[1] == b && ty.atoms[2] == c)
                || (ty.atoms[0] == c && ty.atoms[1] == b && ty.atoms[2] == a)
        })
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct AtomKey {
    pub object: String,
    pub atom_index: usize,
}

#[derive(Debug, Clone)]
pub enum CoordinateSource {
    Loaded(AtomKey),
    RebuiltHydrogen(RebuiltHydrogenSource),
}

#[derive(Debug, Clone)]
pub struct RebuiltHydrogenSource {
    pub controls: Vec<usize>,
    pub bond_length_a: f32,
    pub slot: usize,
    pub count: usize,
    pub method: i32,
}

#[derive(Debug, Clone)]
pub struct ParameterizedAtom {
    pub key: AtomKey,
    pub source: CoordinateSource,
    pub synthetic_parent: Option<usize>,
    pub residue_key: String,
    pub residue_name: String,
    pub atom_name: String,
    pub atom_type: String,
    pub charge: f64,
    pub mass: f64,
    pub sigma_a: f64,
    pub epsilon_kj: f64,
    pub gb_radius_a: f64,
}

#[derive(Debug, Clone)]
pub struct ParameterizedBond {
    pub a: usize,
    pub b: usize,
    pub length_nm: Option<f64>,
    pub force_kj_mol_nm2: Option<f64>,
}

#[derive(Debug, Clone)]
pub struct ParameterizedAngle {
    pub a: usize,
    pub b: usize,
    pub c: usize,
    pub angle_deg: Option<f64>,
    pub force_kj_mol_rad2: Option<f64>,
}

#[derive(Debug, Clone)]
pub struct ParameterizedSystem {
    pub ff_name: String,
    pub source: String,
    pub selection: String,
    pub defaults: ForceFieldDefaults,
    pub atoms: Vec<ParameterizedAtom>,
    pub bonds: Vec<ParameterizedBond>,
    pub angles: Vec<ParameterizedAngle>,
    pub exclusions: HashSet<(usize, usize)>,
    pub one_four_pairs: HashSet<(usize, usize)>,
    pub key_to_index: HashMap<AtomKey, usize>,
    pub report: String,
}

impl ParameterizedSystem {
    pub fn index_for(&self, key: &AtomKey) -> Option<usize> {
        self.key_to_index.get(key).copied()
    }
}

#[derive(Debug, Clone)]
pub struct SelectedMolecule {
    pub object: String,
    pub molecule: ObjectMolecule,
    pub selection: SelectionResult,
}

#[derive(Debug, Clone)]
pub struct MoleculeSnapshot {
    pub object: String,
    pub molecule: ObjectMolecule,
}

impl MoleculeSnapshot {
    pub fn selected_coord(&self, atom_index: usize, state: usize) -> Option<Vec3> {
        self.molecule.get_coord(AtomIndex(atom_index as u32), state)
    }
}

#[derive(Debug, Clone)]
pub struct RunSelection {
    pub receptor: HashSet<usize>,
    pub ligand: HashSet<usize>,
}

impl RunSelection {
    pub fn complex_mask(&self) -> HashSet<usize> {
        self.receptor.union(&self.ligand).copied().collect()
    }
}

#[derive(Debug, Clone, Copy)]
pub struct EnergySettings {
    pub solute_dielectric: f64,
    pub solvent_dielectric: f64,
    pub salt_m: f64,
    pub sa_gamma: f64,
    pub sa_offset: f64,
}

impl Default for EnergySettings {
    fn default() -> Self {
        Self {
            solute_dielectric: 1.0,
            solvent_dielectric: 80.0,
            salt_m: 0.150,
            sa_gamma: 0.022_677_8,
            sa_offset: 0.0,
        }
    }
}

#[derive(Debug, Clone)]
pub struct FrameEnergy {
    pub frame: usize,
    pub complex: EnergyBreakdown,
    pub receptor: EnergyBreakdown,
    pub ligand: EnergyBreakdown,
    pub delta: EnergyBreakdown,
}

#[derive(Debug, Clone, Copy, Default)]
pub struct EnergyBreakdown {
    pub bonded: f64,
    pub coulomb: f64,
    pub lj: f64,
    pub gb: f64,
    pub sa: f64,
}

impl EnergyBreakdown {
    pub fn total(self) -> f64 {
        self.bonded + self.coulomb + self.lj + self.gb + self.sa
    }
}

impl std::ops::Sub for EnergyBreakdown {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self {
            bonded: self.bonded - rhs.bonded,
            coulomb: self.coulomb - rhs.coulomb,
            lj: self.lj - rhs.lj,
            gb: self.gb - rhs.gb,
            sa: self.sa - rhs.sa,
        }
    }
}

#[derive(Debug, Clone)]
pub struct EnergySummary {
    pub frames: Vec<FrameEnergy>,
    pub mean_delta: EnergyBreakdown,
    pub std_delta: EnergyBreakdown,
}

impl fmt::Display for EnergySummary {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "MM/GBSA results over {} frame(s)", self.frames.len())?;
        writeln!(
            f,
            "Delta mean kJ/mol: total={:.3}, bonded={:.3}, coulomb={:.3}, lj={:.3}, gb={:.3}, sa={:.3}",
            self.mean_delta.total(),
            self.mean_delta.bonded,
            self.mean_delta.coulomb,
            self.mean_delta.lj,
            self.mean_delta.gb,
            self.mean_delta.sa
        )?;
        writeln!(
            f,
            "Delta std  kJ/mol: total={:.3}, bonded={:.3}, coulomb={:.3}, lj={:.3}, gb={:.3}, sa={:.3}",
            self.std_delta.total(),
            self.std_delta.bonded,
            self.std_delta.coulomb,
            self.std_delta.lj,
            self.std_delta.gb,
            self.std_delta.sa
        )
    }
}

pub fn ordered_pair(a: usize, b: usize) -> (usize, usize) {
    if a <= b {
        (a, b)
    } else {
        (b, a)
    }
}
