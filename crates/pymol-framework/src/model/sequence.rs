//! Sequence Model
//!
//! Pure domain types and logic for the sequence viewer.
//! No egui dependency — testable, serializable, headless-compatible.
//!
//! Contains: sequence cache, domain types (ResidueRef, SeqObject, etc.),
//! and pure transformation functions (build_seq_object, compress_resi_list, etc.).

use std::collections::HashSet;
use std::ops::Range;

use pymol_mol::{
    is_capping_group, is_ion, is_standard_amino_acid, is_standard_nucleotide, is_water,
    nucleotide_to_char, residue_to_char, three_to_one, ObjectMolecule,
};
use pymol_scene::{Object, ObjectRegistry};

/// Reference to a specific residue in a loaded object
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct ResidueRef {
    pub object_name: String,
    pub chain_id: String,
    pub resv: i32,
}

/// Semantic classification of a residue for sequence viewer coloring.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ResidueKind {
    /// One of the 20 standard amino acids — PyMOL atom color
    AminoAcidCanonical,
    /// Known AA variant (HIP, CYX, MSE, SEP, ACE, NME, …)
    AminoAcidNonCanonical,
    /// Standard unmodified nucleotide
    NucleotideCanonical,
    /// Modified or unknown nucleotide
    NucleotideNonCanonical,
    /// Non-polymer HETATM: ligand, lipid, cofactor
    Ligand,
}

impl ResidueKind {
    pub fn is_polymer(self) -> bool {
        !matches!(self, ResidueKind::Ligand)
    }
}

/// Single residue entry in the sequence viewer
#[derive(Debug, Clone)]
pub struct SeqResidue {
    pub chain: String,
    pub resn: String,
    pub resv: i32,
    pub display_char: char,
    /// Full display label: single char for polymer residues, "[ATP]" for ligands
    pub display_label: String,
    /// Number of character cells this residue occupies in the sequence row
    pub char_width: usize,
    pub atom_range: Range<usize>,
    /// Base color index of the CA atom (or first atom), for display coloring
    pub color_index: i32,
    /// Drives color and font size in the sequence viewer
    pub kind: ResidueKind,
}

/// Sequence data for one chain
#[derive(Debug, Clone)]
pub struct SeqChain {
    pub chain_id: String,
    pub residues: Vec<SeqResidue>,
}

/// Sequence data for one object
#[derive(Debug, Clone)]
pub struct SeqObject {
    pub object_name: String,
    pub chains: Vec<SeqChain>,
}

/// Sequence panel interaction state (no egui dependency).
///
/// Tracks drag selection and hover state for the sequence viewer.
pub struct SequenceUiState {
    /// Active drag start: (object_index, chain_index, start_residue_index)
    pub drag_start: Option<(usize, usize, usize)>,
    /// Current drag end residue index (updated each frame while dragging)
    pub drag_end: Option<usize>,
    /// Current sequence hover state for change detection
    pub current_hover: Option<ResidueRef>,
}

impl Default for SequenceUiState {
    fn default() -> Self {
        Self::new()
    }
}

impl SequenceUiState {
    pub fn new() -> Self {
        Self {
            drag_start: None,
            drag_end: None,
            current_hover: None,
        }
    }

    /// Clear any active drag state (called on dock/float transition).
    pub fn clear_drag(&mut self) {
        self.drag_start = None;
        self.drag_end = None;
    }
}

/// Sequence model: cached sequence data and highlighted residues.
pub struct SequenceModel {
    /// Cached sequence data per object
    pub sequences: Vec<SeqObject>,
    /// Number of objects when cache was last built
    cached_object_count: usize,
    /// Number of enabled objects when cache was last built
    cached_enabled_count: usize,
    /// Registry generation when cache was last built
    cached_generation: u64,
    /// Currently highlighted residues from 3D selection
    pub highlighted: HashSet<ResidueRef>,
}

impl Default for SequenceModel {
    fn default() -> Self {
        Self::new()
    }
}

impl SequenceModel {
    pub fn new() -> Self {
        Self {
            sequences: Vec::new(),
            cached_object_count: 0,
            cached_enabled_count: 0,
            cached_generation: 0,
            highlighted: HashSet::new(),
        }
    }

    /// Update the highlighted residues from the current 3D selection.
    pub fn update_highlights(&mut self, highlighted: HashSet<ResidueRef>) {
        self.highlighted = highlighted;
    }

    /// Rebuild the sequence cache from the registry
    pub fn rebuild_cache(&mut self, registry: &ObjectRegistry) {
        self.sequences.clear();

        for name in registry.names() {
            if let Some(mol_obj) = registry.get_molecule(name) {
                if !mol_obj.is_enabled() {
                    continue;
                }
                let mol = mol_obj.molecule();
                let seq_obj = build_seq_object(name, mol);
                if !seq_obj.chains.is_empty() {
                    self.sequences.push(seq_obj);
                }
            }
        }

        self.cached_object_count = registry.names().count();
        self.cached_enabled_count = registry.enabled_objects().count();
        self.cached_generation = registry.generation();
    }

    /// Check if cache needs rebuilding
    pub fn needs_rebuild(&self, registry: &ObjectRegistry) -> bool {
        registry.names().count() != self.cached_object_count
            || registry.enabled_objects().count() != self.cached_enabled_count
            || registry.generation() != self.cached_generation
    }

    /// Compute the desired panel height based on the number of objects and chains.
    pub fn desired_height(&mut self, registry: &ObjectRegistry) -> f32 {
        if self.needs_rebuild(registry) {
            self.rebuild_cache(registry);
        }
        if self.sequences.is_empty() {
            return 30.0;
        }
        let row_height = 18.0;
        let ruler_height = 10.0;
        let spacing = 2.0;
        let chain_height = ruler_height + row_height;
        let total_chains: usize = self.sequences.iter().map(|s| s.chains.len()).sum();
        let visible_chains = total_chains.min(4) as f32;
        let obj_spacing = self.sequences.len() as f32 * spacing;
        visible_chains * chain_height + obj_spacing + 16.0
    }
}

// =============================================================================
// Pure functions (no egui dependency)
// =============================================================================

/// Compress a sorted list of residue numbers into range notation.
/// E.g., `[74, 75, 76, 77, 80, 85, 86, 87]` → `"74-77+80+85-87"`
pub fn compress_resi_list(values: &[i32]) -> String {
    if values.is_empty() {
        return String::new();
    }
    let mut parts = Vec::new();
    let mut start = values[0];
    let mut end = values[0];
    for &v in &values[1..] {
        if v == end + 1 {
            end = v;
        } else {
            if start == end {
                parts.push(start.to_string());
            } else {
                parts.push(format!("{}-{}", start, end));
            }
            start = v;
            end = v;
        }
    }
    if start == end {
        parts.push(start.to_string());
    } else {
        parts.push(format!("{}-{}", start, end));
    }
    parts.join("+")
}

/// Build sequence data for one molecule object
///
/// ChainIterator groups *consecutive* atoms by chain ID. In mmCIF/PDB files,
/// HETATM records often follow all ATOM records, so the same chain ID (e.g. "A")
/// can appear multiple times. We merge these into a single SeqChain.
pub fn build_seq_object(object_name: &str, mol: &ObjectMolecule) -> SeqObject {
    let mut chain_index: std::collections::HashMap<String, usize> = Default::default();
    let mut chains: Vec<SeqChain> = Vec::new();

    for chain_view in mol.chains() {
        for residue_view in chain_view.residues() {
            let resn = residue_view.resn();

            let (display_label, kind) = if is_water(resn) || is_ion(resn) {
                continue;
            } else if is_capping_group(resn) {
                (format!("[{}]", resn), ResidueKind::AminoAcidNonCanonical)
            } else if residue_view.is_protein() {
                if is_standard_amino_acid(resn) {
                    let ch = three_to_one(resn).expect("standard AA must have 1-letter code");
                    (ch.to_string(), ResidueKind::AminoAcidCanonical)
                } else if let Some(ch) = three_to_one(resn) {
                    (ch.to_string(), ResidueKind::AminoAcidNonCanonical)
                } else {
                    (format!("[{}]", resn), ResidueKind::AminoAcidNonCanonical)
                }
            } else if residue_view.is_nucleic() {
                if is_standard_nucleotide(resn) {
                    let ch = nucleotide_to_char(resn).expect("standard nucleotide must have 1-letter code");
                    (ch.to_string(), ResidueKind::NucleotideCanonical)
                } else if let Some(ch) = nucleotide_to_char(resn) {
                    (ch.to_string(), ResidueKind::NucleotideNonCanonical)
                } else {
                    (format!("[{}]", resn), ResidueKind::NucleotideNonCanonical)
                }
            } else {
                (format!("[{}]", resn), ResidueKind::Ligand)
            };

            let display_char = if kind.is_polymer() {
                residue_to_char(resn)
            } else {
                '?'
            };
            let char_width = display_label.len();

            let color_index = residue_view
                .ca()
                .map(|(_, atom)| atom.repr.colors.base)
                .unwrap_or_else(|| {
                    residue_view
                        .iter()
                        .next()
                        .map(|a| a.repr.colors.base)
                        .unwrap_or(-1)
                });

            let seq_residue = SeqResidue {
                chain: residue_view.chain().to_string(),
                resn: resn.to_string(),
                resv: residue_view.resv(),
                display_char,
                display_label,
                char_width,
                atom_range: residue_view.atom_range.clone(),
                color_index,
                kind,
            };

            let chain_id = chain_view.id().to_string();
            if let Some(&idx) = chain_index.get(&chain_id) {
                chains[idx].residues.push(seq_residue);
            } else {
                let idx = chains.len();
                chain_index.insert(chain_id.clone(), idx);
                chains.push(SeqChain {
                    chain_id,
                    residues: vec![seq_residue],
                });
            }
        }
    }

    SeqObject {
        object_name: object_name.to_string(),
        chains,
    }
}

/// Extract a one-letter-per-residue sequence string from a chain, skipping ligands.
pub fn chain_to_sequence(chain: &SeqChain) -> String {
    chain
        .residues
        .iter()
        .filter(|r| r.kind.is_polymer())
        .map(|r| if r.display_char == '?' { 'X' } else { r.display_char })
        .collect()
}

/// Extract sequence of only the highlighted residues from a chain.
pub fn chain_highlighted_sequence(
    chain: &SeqChain,
    highlighted_resvs: &HashSet<i32>,
) -> String {
    chain
        .residues
        .iter()
        .filter(|r| r.kind.is_polymer())
        .filter(|r| highlighted_resvs.contains(&r.resv))
        .map(|r| if r.display_char == '?' { 'X' } else { r.display_char })
        .collect()
}

/// Format `(header, sequence)` pairs as FASTA with 60-character line wrapping.
pub fn format_fasta(entries: &[(String, String)]) -> String {
    let mut result = String::new();
    for (header, seq) in entries {
        result.push('>');
        result.push_str(header);
        result.push('\n');
        for chunk in seq.as_bytes().chunks(60) {
            result.push_str(std::str::from_utf8(chunk).unwrap());
            result.push('\n');
        }
    }
    if result.ends_with('\n') {
        result.pop();
    }
    result
}

/// Collect all polymer sequences from all objects/chains.
/// Returns `(text, residue_count)` — single chain = plain, multiple = FASTA.
pub fn collect_all_sequences(sequences: &[SeqObject]) -> (String, usize) {
    let mut entries: Vec<(String, String)> = Vec::new();
    let mut total = 0;
    for seq_obj in sequences {
        for chain in &seq_obj.chains {
            let seq = chain_to_sequence(chain);
            if !seq.is_empty() {
                total += seq.len();
                entries.push((
                    format!("{}:{}", seq_obj.object_name, chain.chain_id),
                    seq,
                ));
            }
        }
    }
    if entries.is_empty() {
        return (String::new(), 0);
    }
    let text = if entries.len() == 1 {
        entries.into_iter().next().unwrap().1
    } else {
        format_fasta(&entries)
    };
    (text, total)
}
