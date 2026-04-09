//! Sequence Model
//!
//! Pure domain types and logic for the sequence viewer.
//! No egui dependency — testable, serializable, headless-compatible.
//!
//! Contains: sequence cache, domain types (ResidueRef, SeqObject, etc.),
//! and pure transformation functions (build_seq_object, compress_resi_list, etc.).

use std::collections::HashSet;
use std::ops::Range;

use pymol_color::{
    ChainColors, Color, ColorIndex, ElementColors, NamedColors, residue_type_color, spectrum_color,
    ss_color,
};

/// Everything the model needs from the outside world to rebuild its cache.
///
/// Bundles color-lookup tables so that `rebuild_cache` / `build_seq_object`
/// take a single context reference instead of a growing parameter list.
pub struct SequenceColorContext<'a> {
    pub named_colors: &'a NamedColors,
    pub element_colors: &'a ElementColors,
}
use pymol_mol::{
    is_capping_group, is_ion, is_standard_amino_acid, is_standard_nucleotide, is_water,
    nucleotide_to_char, residue_to_char, three_to_one, Atom, ObjectMolecule,
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
    /// Pre-resolved color for the default residue coloring mode (RGB).
    /// For canonical amino acids this is the named-color of the representative atom;
    /// for other kinds the render path uses fixed constants instead.
    pub named_color: [u8; 3],
    /// Pre-resolved color matching the 3D viewer (RGB)
    pub viewer_color: [u8; 3],
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
    /// When true, residue symbols are colored to match the 3D viewer
    pub viewer_colors: bool,
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
            viewer_colors: false,
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
    /// Checksum of all atoms' base color indices (for fast change detection)
    cached_color_checksum: u64,
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
            cached_color_checksum: 0,
            highlighted: HashSet::new(),
        }
    }

    /// Update the highlighted residues from the current 3D selection.
    pub fn update_highlights(&mut self, highlighted: HashSet<ResidueRef>) {
        self.highlighted = highlighted;
    }

    /// Rebuild the sequence cache from the registry.
    pub fn rebuild_cache(&mut self, registry: &ObjectRegistry, color_ctx: &SequenceColorContext) {
        self.sequences.clear();

        for name in registry.names() {
            if let Some(mol_obj) = registry.get_molecule(name) {
                if !mol_obj.is_enabled() {
                    continue;
                }
                let mol = mol_obj.molecule();
                let resv_range = mol_resv_range(mol);
                let b_range = mol_b_factor_range(mol);
                let seq_obj =
                    build_seq_object(name, mol, color_ctx, resv_range, b_range);
                if !seq_obj.chains.is_empty() {
                    self.sequences.push(seq_obj);
                }
            }
        }

        self.cached_object_count = registry.names().count();
        self.cached_enabled_count = registry.enabled_objects().count();
        self.cached_generation = registry.generation();
        self.cached_color_checksum = color_checksum(registry);
    }

    /// Check if cache needs rebuilding.
    ///
    /// When `viewer_colors` is active, also detects color-only changes
    /// (e.g. `color rainbow`) which don't bump the registry generation counter.
    pub fn needs_rebuild(&self, registry: &ObjectRegistry, viewer_colors: bool) -> bool {
        registry.names().count() != self.cached_object_count
            || registry.enabled_objects().count() != self.cached_enabled_count
            || registry.generation() != self.cached_generation
            || (viewer_colors && color_checksum(registry) != self.cached_color_checksum)
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

/// Compute a fast checksum of all atoms' base color indices across enabled molecules.
///
/// Uses FNV-style hashing — O(atoms) but no allocations. Detects any color change
/// including partial selections like `color red, chain A`.
fn color_checksum(registry: &ObjectRegistry) -> u64 {
    let mut hash: u64 = 0xcbf29ce484222325; // FNV offset basis
    for name in registry.names() {
        if let Some(mol_obj) = registry.get_molecule(name) {
            if !mol_obj.is_enabled() {
                continue;
            }
            for atom in mol_obj.molecule().atoms() {
                hash ^= atom.repr.colors.base as u64;
                hash = hash.wrapping_mul(0x100000001b3); // FNV prime
            }
        }
    }
    hash
}

/// Compute (min, max) residue sequence number range for a molecule.
fn mol_resv_range(mol: &ObjectMolecule) -> (i32, i32) {
    let mut min = i32::MAX;
    let mut max = i32::MIN;
    for atom in mol.atoms() {
        let v = atom.residue.resv;
        if v < min { min = v; }
        if v > max { max = v; }
    }
    if min > max { (1, 100) } else { (min, max) }
}

/// Compute (min, max) B-factor range for a molecule.
fn mol_b_factor_range(mol: &ObjectMolecule) -> (f32, f32) {
    let mut min = f32::MAX;
    let mut max = f32::MIN;
    for atom in mol.atoms() {
        if atom.b_factor < min { min = atom.b_factor; }
        if atom.b_factor > max { max = atom.b_factor; }
    }
    if min > max { (0.0, 100.0) } else { (min, max) }
}

/// Resolve an atom's named color index to RGB, if it has a named color.
fn resolve_named_color(atom: &Atom, named_colors: &NamedColors) -> Option<[u8; 3]> {
    if atom.repr.colors.base >= 0 {
        if let Some(color) = named_colors.get_by_index(atom.repr.colors.base as u32) {
            return Some([
                (color.r * 255.0) as u8,
                (color.g * 255.0) as u8,
                (color.b * 255.0) as u8,
            ]);
        }
    }
    None
}

/// Resolve an atom's base color index to an RGB triple, matching 3D viewer logic.
fn resolve_viewer_color(
    atom: &Atom,
    named_colors: &NamedColors,
    element_colors: &ElementColors,
    resv_range: (i32, i32),
    b_factor_range: (f32, f32),
) -> [u8; 3] {
    let color = match ColorIndex::from(atom.repr.colors.base) {
        ColorIndex::ByElement | ColorIndex::Atomic => {
            element_colors.get(atom.element as u8)
        }
        ColorIndex::ByChain => ChainColors::get(&atom.residue.chain),
        ColorIndex::BySS => ss_color(atom.ss_type as u8),
        ColorIndex::ByResidueType => {
            residue_type_color(&atom.residue.resn)
                .unwrap_or_else(|| element_colors.get(atom.element as u8))
        }
        ColorIndex::ByResidueIndex => {
            let (min, max) = resv_range;
            let range = max - min;
            if range <= 0 {
                spectrum_color(0.5)
            } else {
                let t = ((atom.residue.resv - min) as f32) / (range as f32);
                spectrum_color(t)
            }
        }
        ColorIndex::ByBFactor => {
            let (min, max) = b_factor_range;
            let range = max - min;
            if range <= 0.0 {
                Color::WHITE
            } else {
                let t = ((atom.b_factor - min) / range).clamp(0.0, 1.0);
                if t < 0.5 {
                    let s = t * 2.0;
                    Color::new(s, s, 1.0)
                } else {
                    let s = (t - 0.5) * 2.0;
                    Color::new(1.0, 1.0 - s, 1.0 - s)
                }
            }
        }
        ColorIndex::Named(idx) => {
            named_colors.get_by_index(idx).unwrap_or(Color::WHITE)
        }
    };
    [
        (color.r * 255.0) as u8,
        (color.g * 255.0) as u8,
        (color.b * 255.0) as u8,
    ]
}

/// Build sequence data for one molecule object
///
/// ChainIterator groups *consecutive* atoms by chain ID. In mmCIF/PDB files,
/// HETATM records often follow all ATOM records, so the same chain ID (e.g. "A")
/// can appear multiple times. We merge these into a single SeqChain.
pub fn build_seq_object(
    object_name: &str,
    mol: &ObjectMolecule,
    color_ctx: &SequenceColorContext,
    resv_range: (i32, i32),
    b_factor_range: (f32, f32),
) -> SeqObject {
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

            // Pick representative atom: CA for proteins, first atom otherwise
            let repr_atom = residue_view
                .ca()
                .map(|(_, atom)| atom)
                .or_else(|| residue_view.iter().next());

            let named_color = repr_atom
                .and_then(|a| resolve_named_color(a, color_ctx.named_colors))
                .unwrap_or([200, 200, 200]);

            let viewer_color = repr_atom
                .map(|a| {
                    resolve_viewer_color(
                        a,
                        color_ctx.named_colors,
                        color_ctx.element_colors,
                        resv_range,
                        b_factor_range,
                    )
                })
                .unwrap_or([200, 200, 200]);

            let seq_residue = SeqResidue {
                chain: residue_view.chain().to_string(),
                resn: resn.to_string(),
                resv: residue_view.resv(),
                display_char,
                display_label,
                char_width,
                atom_range: residue_view.atom_range.clone(),
                named_color,
                viewer_color,
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
