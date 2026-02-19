pub mod reader;

use std::collections::HashMap;

/// A fully deserialized PSE session (intermediate representation).
#[derive(Debug, Clone)]
pub struct PseSession {
    pub version: i64,
    pub settings: Vec<PseSetting>,
    pub view: [f64; 18],
    pub names: Vec<Option<PseNameEntry>>,
    pub colors: Vec<PseColor>,
    pub view_dict: HashMap<String, [f64; 18]>,
    pub scene_order: Vec<String>,
    pub unique_settings: Vec<PseUniqueSetting>,
}

#[derive(Debug, Clone)]
pub struct PseSetting {
    pub id: u16,
    pub type_code: u8,
    pub value: PseSettingValue,
}

#[derive(Debug, Clone, PartialEq)]
pub enum PseSettingValue {
    Bool(bool),
    Int(i64),
    Float(f64),
    Float3([f64; 3]),
    Color(i64),
    String(String),
}

#[derive(Debug, Clone)]
pub struct PseColor {
    pub name: String,
    pub rgb: [f32; 3],
}

#[derive(Debug, Clone)]
pub struct PseUniqueSetting {
    pub unique_id: u64,
    pub setting_id: u16,
    pub type_code: u8,
    pub value: PseSettingValue,
}

#[derive(Debug, Clone)]
pub enum PseNameEntry {
    Object(PseObject),
    Selection(PseSelection),
}

#[derive(Debug, Clone)]
pub struct PseObject {
    pub name: String,
    /// 1=molecule, 2=map, 3=mesh, 4=surface, etc.
    pub type_code: i64,
    pub visible: bool,
    pub rep_mask: u32,
    pub color: i64,
    pub data: PseObjectData,
}

#[derive(Debug, Clone)]
pub enum PseObjectData {
    Molecule(PseMolecule),
    /// Maps, CGOs, surfaces, etc. — not yet parsed.
    Unsupported,
}

#[derive(Debug, Clone)]
pub struct PseMolecule {
    pub atoms: Vec<PseAtom>,
    pub bonds: Vec<PseBond>,
    pub coord_sets: Vec<Option<PseCoordSet>>,
    pub discrete: bool,
    pub n_discrete: i64,
    pub symmetry: Option<PseSymmetry>,
}

#[derive(Debug, Clone)]
pub struct PseAtom {
    pub resv: i64,
    pub chain: String,
    pub alt: String,
    pub resi: String,
    pub segi: String,
    pub resn: String,
    pub name: String,
    pub elem: String,
    pub text_type: String,
    pub label: String,
    pub ss: String,
    pub b: f64,
    pub q: f64,
    pub vdw: f64,
    pub partial_charge: f64,
    pub formal_charge: i64,
    pub color: i64,
    pub id: i64,
    pub unique_id: u64,
    /// Per-atom visible representations bitmask (from PSE field [20]).
    pub visible_reps: u32,
    /// Whether this is a HETATM record (from PSE field [22]).
    pub hetatm: bool,
}

#[derive(Debug, Clone)]
pub struct PseBond {
    pub index0: usize,
    pub index1: usize,
    pub order: i64,
    pub id: i64,
    pub stereo: i64,
    pub unique_id: u64,
}

#[derive(Debug, Clone)]
pub struct PseCoordSet {
    pub n_atom: usize,
    /// Flattened `[x, y, z, x, y, z, ...]` in index order.
    pub coords: Vec<f64>,
    /// Mapping from coordinate index to atom index.
    /// If empty, assumes identity mapping (coord[i] → atom[i]).
    pub idx_to_atm: Vec<usize>,
}

#[derive(Debug, Clone)]
pub struct PseSymmetry {
    /// `[a, b, c, alpha, beta, gamma]`
    pub cell: [f64; 6],
    pub space_group: String,
}

#[derive(Debug, Clone)]
pub struct PseSelection {
    pub name: String,
    pub visible: bool,
    /// `(object_index_in_names, atom_indices)`
    pub members: Vec<(usize, Vec<usize>)>,
}
