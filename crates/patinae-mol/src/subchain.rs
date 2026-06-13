//! Chain subchains
//!
//! A *subchain* is a coherent unit of atoms within a single chain — either a
//! contiguous biopolymer stretch, a single ligand residue (e.g. HEM), a run of
//! solvent waters, an ion, or a connected component of HET atoms in the
//! covalent-bond graph (e.g. a multi-residue N-glycan such as NAG–NAG–BMA–MAN).
//! Subchains are the natural unit for the Objects panel: each row in the
//! sidebar maps to one subchain.
//!
//! Two paths produce subchains:
//!
//! - [`ObjectMolecule::subchains`](crate::ObjectMolecule::subchains) — bond-aware.
//!   Walks a cached [`SubchainPartition`] that merges connected HET runs into
//!   a single subchain. This is what the Objects panel and any UI consumer
//!   should use.
//! - [`ChainView::subchains`](crate::ChainView::subchains) — predicate-based.
//!   Walks consecutive runs grouped by `(chain, hetatm, resn)`. Used for
//!   chain-scoped iteration where bond merging is irrelevant (cartoon /
//!   ribbon polymer subchains, tests).

use serde::{Deserialize, Serialize};
use std::borrow::Cow;
use std::ops::Range;

use crate::atom::Atom;
use crate::flags::AtomFlags;
use crate::index::AtomIndex;
use crate::residue::{atoms_same_residue, ResidueIterator};

// ============================================================================
// SubchainKind
// ============================================================================

/// Semantic category for a chain subchain.
///
/// Computed from the first atom of the subchain. The build pass guarantees
/// every atom in the subchain shares the same classification, so peeking at
/// the first atom is sufficient for downstream consumers.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum SubchainKind {
    /// Protein or nucleic acid (PROTEIN or NUCLEIC flags).
    Biopolymer,
    /// Small-molecule organic ligand (e.g. HEM, ATP) or connected glycan.
    Organic,
    /// Solvent (water).
    Solvent,
    /// Ions and other inorganic groups.
    Inorganic,
    /// Unclassified.
    Other,
}

impl SubchainKind {
    /// Classify a subchain from its first atom's flags.
    pub fn from_atom(atom: &Atom) -> Self {
        let flags = atom.state.flags;
        if flags.is_biomolecule() {
            SubchainKind::Biopolymer
        } else if flags.is_solvent() {
            SubchainKind::Solvent
        } else if flags.contains(AtomFlags::ORGANIC) {
            SubchainKind::Organic
        } else if flags.contains(AtomFlags::INORGANIC) {
            SubchainKind::Inorganic
        } else {
            SubchainKind::Other
        }
    }

    /// Stable string tag for serialization or UI binding.
    pub fn as_str(self) -> &'static str {
        match self {
            SubchainKind::Biopolymer => "biopolymer",
            SubchainKind::Organic => "organic",
            SubchainKind::Solvent => "solvent",
            SubchainKind::Inorganic => "inorganic",
            SubchainKind::Other => "other",
        }
    }

    /// Returns true for protein / nucleic-acid subchains.
    #[inline]
    pub fn is_biopolymer(self) -> bool {
        matches!(self, SubchainKind::Biopolymer)
    }
}

// ============================================================================
// Atom layout & label
// ============================================================================

/// How a subchain's atoms are laid out in the parent molecule's atom array.
///
/// `All` is the most popular case: a single subchain that covers the entire
/// owning molecule (e.g. a single-chain biopolymer with no HETs). `Range` is
/// the next-most-common (biopolymer runs, solvent runs, single-residue
/// ligands, contiguous glycans). `Indices` is only emitted by the partition
/// builder when a connected HET component spans non-contiguous atom indices.
#[derive(Debug, Clone)]
pub enum SubchainAtoms {
    /// Covers the entire owning molecule's atom slice. Equivalent to
    /// `Range(0..atom_count)` but saves bytes per entry and makes the
    /// "single subchain spans the whole molecule" case explicit. Emitted by
    /// the partition builder when there is exactly one entry covering all
    /// atoms.
    All { atom_count: u32 },
    /// Atoms occupy a contiguous half-open range of global atom indices.
    Range(Range<u32>),
    /// Atoms occupy a sorted list of global atom indices.
    Indices(Vec<u32>),
}

impl SubchainAtoms {
    /// Number of atoms in the subchain.
    #[inline]
    pub fn len(&self) -> usize {
        match self {
            SubchainAtoms::All { atom_count } => *atom_count as usize,
            SubchainAtoms::Range(r) => (r.end - r.start) as usize,
            SubchainAtoms::Indices(v) => v.len(),
        }
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Smallest atom index in the subchain (or `None` if empty).
    #[inline]
    pub fn min_atom_index(&self) -> Option<u32> {
        match self {
            SubchainAtoms::All { atom_count } => (*atom_count > 0).then_some(0),
            SubchainAtoms::Range(r) => (r.start < r.end).then_some(r.start),
            SubchainAtoms::Indices(v) => v.first().copied(),
        }
    }

    /// Returns `true` when atoms are contiguous (`All` or `Range`).
    #[inline]
    pub fn is_contiguous(&self) -> bool {
        matches!(self, SubchainAtoms::All { .. } | SubchainAtoms::Range(_))
    }
}

/// Display label for a subchain. Determined during partition build (or by
/// the predicate iterator for chain-scoped views).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SubchainLabel {
    /// Biopolymer subchain — display the chain identifier.
    Polymer,
    /// All atoms agree on residue name (HEM, ATP, single NAG, water …).
    Single(String),
    /// Multi-residue connected HET component — show
    /// `"{primary_resn}+{n-1}"` (e.g. `"NAG+8"` for the 1IGT N-glycan).
    Composite {
        residue_count: u32,
        primary_resn: String,
    },
}

// ============================================================================
// SubchainEntry / SubchainPartition
// ============================================================================

/// One concrete subchain inside a [`SubchainPartition`].
#[derive(Debug, Clone)]
pub struct SubchainEntry {
    pub chain_id: String,
    pub kind: SubchainKind,
    pub label: SubchainLabel,
    pub atoms: SubchainAtoms,
}

/// `atom_index → entry_index` lookup, packed for the common case.
///
/// `AllEntries` is the zero-allocation form: a single `All` partition entry
/// covers every atom in `0..atom_count`, so `entry_for_atom` can answer
/// directly without a `Vec` lookup. `Map` is the per-atom table emitted when
/// the partition has multiple entries (or non-contiguous coverage).
#[derive(Debug, Clone)]
pub(crate) enum AtomToEntryMap {
    /// Trivial: a single entry covers every atom.
    AllEntries { atom_count: u32 },
    /// Per-atom lookup table. `u32::MAX` for atoms not assigned to any entry.
    Map(Vec<u32>),
}

impl Default for AtomToEntryMap {
    fn default() -> Self {
        AtomToEntryMap::AllEntries { atom_count: 0 }
    }
}

/// Partition of a molecule's atoms into subchains.
///
/// Built lazily via
/// [`ObjectMolecule::subchain_partition`](crate::ObjectMolecule::subchain_partition)
/// and cached behind a `OnceLock`. Invalidated whenever bonds or atom
/// classification change.
#[derive(Debug, Clone, Default)]
pub struct SubchainPartition {
    pub(crate) entries: Vec<SubchainEntry>,
    pub(crate) atom_to_entry: AtomToEntryMap,
}

impl SubchainPartition {
    #[inline]
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    pub fn entries(&self) -> &[SubchainEntry] {
        &self.entries
    }

    /// Look up the partition entry index for a given atom.
    pub fn entry_for_atom(&self, atom: AtomIndex) -> Option<u32> {
        match &self.atom_to_entry {
            AtomToEntryMap::AllEntries { atom_count } => {
                (atom.as_usize() < *atom_count as usize).then_some(0)
            }
            AtomToEntryMap::Map(map) => {
                map.get(atom.as_usize()).copied().filter(|&e| e != u32::MAX)
            }
        }
    }

    /// Construct a [`SubchainView`] for entry `idx`. `atoms` must be the
    /// owning molecule's full atom slice.
    pub fn view_for<'a>(&'a self, idx: u32, atoms: &'a [Atom]) -> Option<SubchainView<'a>> {
        self.entries
            .get(idx as usize)
            .map(|e| SubchainView::from_entry(e, atoms))
    }

    /// Build a partition from atoms only — no bond-graph merging. Each raw
    /// run becomes one entry. Useful for tests on synthetic atoms.
    pub fn from_atoms(atoms: &[Atom]) -> Self {
        let raw = build_raw_runs(atoms);
        if raw.is_empty() {
            return SubchainPartition::default();
        }

        // Fast path: one raw run covers every atom — emit `All` directly and
        // skip the per-atom map allocation.
        if raw.len() == 1 {
            let r = &raw[0];
            let entry = SubchainEntry {
                chain_id: r.chain_id.clone(),
                kind: r.kind,
                label: match &r.het_resn {
                    Some(resn) => SubchainLabel::Single(resn.clone()),
                    None => SubchainLabel::Polymer,
                },
                atoms: SubchainAtoms::All {
                    atom_count: atoms.len() as u32,
                },
            };
            return SubchainPartition {
                entries: vec![entry],
                atom_to_entry: AtomToEntryMap::AllEntries {
                    atom_count: atoms.len() as u32,
                },
            };
        }

        let entries: Vec<SubchainEntry> = raw
            .iter()
            .map(|r| SubchainEntry {
                chain_id: r.chain_id.clone(),
                kind: r.kind,
                label: match &r.het_resn {
                    Some(resn) => SubchainLabel::Single(resn.clone()),
                    None => SubchainLabel::Polymer,
                },
                atoms: SubchainAtoms::Range(r.start..r.end),
            })
            .collect();
        let atom_to_entry = AtomToEntryMap::Map(build_atom_to_entry_vec(atoms.len(), &entries));
        SubchainPartition {
            entries,
            atom_to_entry,
        }
    }

    /// Build a partition with bond-graph merging. Used by
    /// [`ObjectMolecule::subchain_partition`](crate::ObjectMolecule::subchain_partition).
    pub fn from_molecule(atoms: &[Atom], bonds: &[crate::bond::Bond]) -> Self {
        build_partition_with_bonds(atoms, bonds)
    }
}

// ============================================================================
// Subchain view
// ============================================================================

/// A view into one subchain of a molecule.
///
/// Owns its `chain_id`, `kind`, and `label`; borrows the underlying atom
/// slice plus an atom-set descriptor. Iterators in this crate construct the
/// appropriate variant — predicate-based for chain-scoped iteration,
/// partition-based for the bond-aware `mol.subchains()` API.
#[derive(Debug, Clone)]
pub struct SubchainView<'a> {
    chain_id: String,
    pub kind: SubchainKind,
    pub label: SubchainLabel,
    /// Slice that `atom_set` indices (or range bounds) refer into.
    /// For the predicate path this is the chain slice; for the partition
    /// path this is the full molecule slice.
    pub(crate) atoms_all: &'a [Atom],
    /// Global atom-index of `atoms_all[0]`. Added to local Range bounds in
    /// `iter_indexed()` to recover global indices. For the partition path
    /// this is 0 (since `atoms_all` is the full molecule slice). For the
    /// predicate path this is the chain's atom_range.start.
    pub(crate) base: u32,
    pub(crate) atom_set: AtomSetRef<'a>,
}

/// Borrowed reference into a subchain's atom set.
#[derive(Debug, Clone)]
pub(crate) enum AtomSetRef<'a> {
    /// Local indices `[start, end)` into `atoms_all`. Global atom indices
    /// are `base + start .. base + end`.
    Range { start: u32, end: u32 },
    /// Global atom indices listed in this slice. The view's `base` is 0
    /// when this variant is used (i.e. `atoms_all` is the full molecule
    /// slice).
    Indices(&'a [u32]),
}

impl<'a> SubchainView<'a> {
    /// Construct a Range-layout view from a slice + base index.
    pub(crate) fn from_range(
        chain_id: String,
        kind: SubchainKind,
        label: SubchainLabel,
        atoms_all: &'a [Atom],
        base: u32,
        local_start: u32,
        local_end: u32,
    ) -> Self {
        SubchainView {
            chain_id,
            kind,
            label,
            atoms_all,
            base,
            atom_set: AtomSetRef::Range {
                start: local_start,
                end: local_end,
            },
        }
    }

    /// Construct a view from a partition entry. `atoms_all` is the full
    /// molecule slice; `base` is implicitly 0.
    pub(crate) fn from_entry(entry: &'a SubchainEntry, atoms_all: &'a [Atom]) -> Self {
        let atom_set = match &entry.atoms {
            SubchainAtoms::All { atom_count } => AtomSetRef::Range {
                start: 0,
                end: *atom_count,
            },
            SubchainAtoms::Range(r) => AtomSetRef::Range {
                start: r.start,
                end: r.end,
            },
            SubchainAtoms::Indices(v) => AtomSetRef::Indices(v.as_slice()),
        };
        SubchainView {
            chain_id: entry.chain_id.clone(),
            kind: entry.kind,
            label: entry.label.clone(),
            atoms_all,
            base: 0,
            atom_set,
        }
    }

    /// Chain identifier.
    #[inline]
    pub fn chain_id(&self) -> &str {
        &self.chain_id
    }

    /// Number of atoms in the subchain.
    #[inline]
    pub fn len(&self) -> usize {
        match &self.atom_set {
            AtomSetRef::Range { start, end } => (*end - *start) as usize,
            AtomSetRef::Indices(v) => v.len(),
        }
    }

    /// Returns `true` if the subchain has no atoms.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns `true` when atoms are stored as a contiguous range.
    #[inline]
    pub fn is_contiguous(&self) -> bool {
        matches!(self.atom_set, AtomSetRef::Range { .. })
    }

    /// Iterate over atoms in the subchain.
    pub fn iter(&self) -> SubchainAtomIter<'a> {
        SubchainAtomIter::new(self.atoms_all, &self.atom_set)
    }

    /// Iterate over `(global AtomIndex, &Atom)` pairs in the subchain.
    pub fn iter_indexed(&self) -> SubchainAtomIndexedIter<'a> {
        SubchainAtomIndexedIter::new(self.atoms_all, self.base, &self.atom_set)
    }

    /// Iterate over residues in the subchain.
    ///
    /// Only supported for contiguous (`Range`) subchains, which is always
    /// the case for biopolymer / solvent / inorganic / single-residue HET.
    /// Calling on a non-contiguous (`Indices`) subchain panics — there is no
    /// consumer that needs it today, and a residue iterator over scrambled
    /// indices would require allocating intermediate storage. Use
    /// [`iter`](Self::iter) plus manual residue grouping if needed.
    pub fn residues(&self) -> ResidueIterator<'a> {
        match &self.atom_set {
            AtomSetRef::Range { start, end } => {
                let s = *start as usize;
                let e = *end as usize;
                let global_base = self.base as usize + s;
                ResidueIterator::new(&self.atoms_all[s..e], global_base)
            }
            AtomSetRef::Indices(_) => {
                panic!("SubchainView::residues() called on non-contiguous subchain")
            }
        }
    }

    /// Display label — `Polymer` → chain id, `Single(s)` → `s`,
    /// `Composite { residue_count, primary_resn }` →
    /// `"{primary_resn}+{residue_count - 1}"`.
    pub fn label(&self) -> Cow<'_, str> {
        match &self.label {
            SubchainLabel::Polymer => Cow::Borrowed(self.chain_id.as_str()),
            SubchainLabel::Single(s) => Cow::Borrowed(s.as_str()),
            SubchainLabel::Composite {
                residue_count,
                primary_resn,
            } => Cow::Owned(format!(
                "{}+{}",
                primary_resn,
                residue_count.saturating_sub(1)
            )),
        }
    }

    /// First atom of the subchain (lowest atom index).
    pub fn first_atom(&self) -> Option<&'a Atom> {
        let mut it = self.iter();
        it.next()
    }
}

// ============================================================================
// Atom iterators
// ============================================================================

#[derive(Clone)]
pub struct SubchainAtomIter<'a> {
    atoms_all: &'a [Atom],
    layout: AtomSetRef<'a>,
    pos: u32,
}

impl<'a> SubchainAtomIter<'a> {
    fn new(atoms_all: &'a [Atom], atom_set: &AtomSetRef<'a>) -> Self {
        SubchainAtomIter {
            atoms_all,
            layout: atom_set.clone(),
            pos: 0,
        }
    }
}

impl<'a> Iterator for SubchainAtomIter<'a> {
    type Item = &'a Atom;

    fn next(&mut self) -> Option<Self::Item> {
        match &self.layout {
            AtomSetRef::Range { start, end } => {
                let i = *start + self.pos;
                if i >= *end {
                    return None;
                }
                self.pos += 1;
                self.atoms_all.get(i as usize)
            }
            AtomSetRef::Indices(idxs) => {
                let p = self.pos as usize;
                if p >= idxs.len() {
                    return None;
                }
                self.pos += 1;
                self.atoms_all.get(idxs[p] as usize)
            }
        }
    }
}

#[derive(Clone)]
pub struct SubchainAtomIndexedIter<'a> {
    atoms_all: &'a [Atom],
    base: u32,
    layout: AtomSetRef<'a>,
    pos: u32,
}

impl<'a> SubchainAtomIndexedIter<'a> {
    fn new(atoms_all: &'a [Atom], base: u32, atom_set: &AtomSetRef<'a>) -> Self {
        SubchainAtomIndexedIter {
            atoms_all,
            base,
            layout: atom_set.clone(),
            pos: 0,
        }
    }
}

impl<'a> Iterator for SubchainAtomIndexedIter<'a> {
    type Item = (AtomIndex, &'a Atom);

    fn next(&mut self) -> Option<Self::Item> {
        match &self.layout {
            AtomSetRef::Range { start, end } => {
                let local = *start + self.pos;
                if local >= *end {
                    return None;
                }
                self.pos += 1;
                let atom = self.atoms_all.get(local as usize)?;
                let global = self.base + local;
                Some((AtomIndex(global), atom))
            }
            AtomSetRef::Indices(idxs) => {
                let p = self.pos as usize;
                if p >= idxs.len() {
                    return None;
                }
                self.pos += 1;
                let global = idxs[p];
                let atom = self.atoms_all.get(global as usize)?;
                Some((AtomIndex(global), atom))
            }
        }
    }
}

// ============================================================================
// Predicate (used by the chain-scoped iterator and Pass 1 of partition build)
// ============================================================================

/// Predicate: do two consecutive atoms belong to the same *raw* subchain run?
///
/// Atoms in the same raw run share their chain identifier, their hetatm flag,
/// and — for non-polymer atoms — their residue name. Bond-graph merging on
/// top of this predicate is handled by [`SubchainPartition::from_molecule`].
#[inline]
pub fn atoms_same_subchain(a: &Atom, b: &Atom) -> bool {
    a.residue.chain == b.residue.chain
        && a.state.hetatm == b.state.hetatm
        && (!a.state.hetatm || a.residue.resn == b.residue.resn)
}

// ============================================================================
// Partition iterator
// ============================================================================

/// Iterator over a [`SubchainPartition`]'s entries, yielding `SubchainView`s
/// that borrow the underlying molecule atoms. This is what
/// [`ObjectMolecule::subchains`](crate::ObjectMolecule::subchains) returns.
pub struct PartitionSubchainIter<'a> {
    atoms_all: &'a [Atom],
    entries: std::slice::Iter<'a, SubchainEntry>,
}

impl<'a> PartitionSubchainIter<'a> {
    pub(crate) fn new(atoms_all: &'a [Atom], partition: &'a SubchainPartition) -> Self {
        PartitionSubchainIter {
            atoms_all,
            entries: partition.entries.iter(),
        }
    }
}

impl<'a> Iterator for PartitionSubchainIter<'a> {
    type Item = SubchainView<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        let entry = self.entries.next()?;
        Some(SubchainView::from_entry(entry, self.atoms_all))
    }
}

// ============================================================================
// Polymer subchain iterator
// ============================================================================

/// Iterate over only the biopolymer subchains of an underlying subchain
/// iterator. Used by representations like cartoon and ribbon that only
/// operate on protein / nucleic-acid backbones.
pub struct PolymerSubchainIterator<'a> {
    inner: crate::iterator::SubchainIterator<'a>,
}

impl<'a> PolymerSubchainIterator<'a> {
    pub(crate) fn new(inner: crate::iterator::SubchainIterator<'a>) -> Self {
        PolymerSubchainIterator { inner }
    }
}

impl<'a> Iterator for PolymerSubchainIterator<'a> {
    type Item = SubchainView<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        self.inner.by_ref().find(|sub| sub.kind.is_biopolymer())
    }
}

// ============================================================================
// Partition build
// ============================================================================

#[derive(Debug, Clone)]
pub(crate) struct RawRun {
    pub start: u32,
    pub end: u32,
    pub chain_id: String,
    pub kind: SubchainKind,
    /// `Some(resn)` for HET runs, `None` for biopolymer runs.
    pub het_resn: Option<String>,
}

pub(crate) fn build_raw_runs(atoms: &[Atom]) -> Vec<RawRun> {
    let mut runs = Vec::new();
    if atoms.is_empty() {
        return runs;
    }
    let mut start = 0usize;
    while start < atoms.len() {
        let first = &atoms[start];
        let mut end = start + 1;
        while end < atoms.len() && atoms_same_subchain(&atoms[end], first) {
            end += 1;
        }
        let kind = SubchainKind::from_atom(first);
        let het_resn = if first.state.hetatm {
            Some(first.residue.resn.to_string())
        } else {
            None
        };
        runs.push(RawRun {
            start: start as u32,
            end: end as u32,
            chain_id: first.residue.chain.clone(),
            kind,
            het_resn,
        });
        start = end;
    }
    runs
}

fn build_partition_with_bonds(atoms: &[Atom], bonds: &[crate::bond::Bond]) -> SubchainPartition {
    let raw = build_raw_runs(atoms);
    if raw.is_empty() {
        return SubchainPartition::default();
    }

    // Fast path: a single raw run already covers every atom — there is
    // nothing to union, no HashMap to bucket, and no per-atom map to fill.
    // Skip straight to the single-entry partition with the trivial map.
    if raw.len() == 1 {
        let entry = build_entry(atoms, &raw, &[0]);
        let atom_to_entry = match &entry.atoms {
            SubchainAtoms::All { atom_count } => AtomToEntryMap::AllEntries {
                atom_count: *atom_count,
            },
            _ => AtomToEntryMap::Map(build_atom_to_entry_vec(
                atoms.len(),
                std::slice::from_ref(&entry),
            )),
        };
        return SubchainPartition {
            entries: vec![entry],
            atom_to_entry,
        };
    }

    let mut atom_to_run = vec![u32::MAX; atoms.len()];
    for (ri, r) in raw.iter().enumerate() {
        for ai in r.start..r.end {
            atom_to_run[ai as usize] = ri as u32;
        }
    }

    // Pass 2 — union-find over HET Organic runs that are covalently
    // connected within the same chain.
    let mut uf = UnionFind::new(raw.len());
    for bond in bonds {
        let a1 = bond.atom1.as_usize();
        let a2 = bond.atom2.as_usize();
        if a1 >= atoms.len() || a2 >= atoms.len() {
            continue;
        }
        let r1 = atom_to_run[a1];
        let r2 = atom_to_run[a2];
        if r1 == u32::MAX || r2 == u32::MAX || r1 == r2 {
            continue;
        }
        let run1 = &raw[r1 as usize];
        let run2 = &raw[r2 as usize];
        if !is_mergeable_kind(run1.kind) || !is_mergeable_kind(run2.kind) {
            continue;
        }
        if run1.chain_id != run2.chain_id {
            continue;
        }
        // Defensive: a HET run shouldn't contain biopolymer atoms, but the
        // per-atom flags are the source of truth — never let a polymer atom
        // pull its run into a HET component.
        let f1 = atoms[a1].state.flags;
        let f2 = atoms[a2].state.flags;
        if f1.is_biomolecule() || f2.is_biomolecule() {
            continue;
        }
        uf.union(r1 as usize, r2 as usize);
    }

    // Pass 3 — gather members per component, build entries.
    let mut by_root: std::collections::HashMap<u32, Vec<u32>> = std::collections::HashMap::new();
    for ri in 0..raw.len() {
        let root = uf.find(ri) as u32;
        by_root.entry(root).or_default().push(ri as u32);
    }

    let mut entries: Vec<SubchainEntry> = by_root
        .into_values()
        .map(|members| build_entry(atoms, &raw, &members))
        .collect();

    entries.sort_by(|a, b| {
        compare_natural(a.chain_id.as_str(), b.chain_id.as_str()).then_with(|| {
            a.atoms
                .min_atom_index()
                .unwrap_or(u32::MAX)
                .cmp(&b.atoms.min_atom_index().unwrap_or(u32::MAX))
        })
    });

    let atom_to_entry = AtomToEntryMap::Map(build_atom_to_entry_vec(atoms.len(), &entries));
    SubchainPartition {
        entries,
        atom_to_entry,
    }
}

/// Compare two strings with natural (alphanumeric) ordering.
///
/// Splits each string into alternating non-digit and digit segments.
/// Non-digit segments are compared lexicographically; digit segments are
/// compared numerically. When one string is a prefix of the other, the
/// shorter one sorts first (e.g. `"A" < "A1"`).
fn compare_natural(a: &str, b: &str) -> std::cmp::Ordering {
    let mut achars = a.chars().peekable();
    let mut bchars = b.chars().peekable();

    loop {
        let a_next = achars.peek().map(|c| c.is_ascii_digit());
        let b_next = bchars.peek().map(|c| c.is_ascii_digit());

        match (a_next, b_next) {
            (None, None) => return std::cmp::Ordering::Equal,
            (None, Some(_)) => return std::cmp::Ordering::Less,
            (Some(_), None) => return std::cmp::Ordering::Greater,
            (Some(true), Some(true)) => {
                let a_num: String = (&mut achars).take_while(|c| c.is_ascii_digit()).collect();
                let b_num: String = (&mut bchars).take_while(|c| c.is_ascii_digit()).collect();
                match a_num
                    .parse::<u64>()
                    .unwrap_or(0)
                    .cmp(&b_num.parse::<u64>().unwrap_or(0))
                {
                    std::cmp::Ordering::Equal => {
                        let len_cmp = a_num.len().cmp(&b_num.len());
                        if len_cmp != std::cmp::Ordering::Equal {
                            return len_cmp;
                        }
                    }
                    non_eq => return non_eq,
                }
            }
            (Some(false), Some(false)) => {
                let a_str: String = (&mut achars).take_while(|c| !c.is_ascii_digit()).collect();
                let b_str: String = (&mut bchars).take_while(|c| !c.is_ascii_digit()).collect();
                match a_str.cmp(&b_str) {
                    std::cmp::Ordering::Equal => {}
                    non_eq => return non_eq,
                }
            }
            (Some(true), Some(false)) => return std::cmp::Ordering::Less,
            (Some(false), Some(true)) => return std::cmp::Ordering::Greater,
        }
    }
}

fn is_mergeable_kind(k: SubchainKind) -> bool {
    matches!(
        k,
        SubchainKind::Organic | SubchainKind::Inorganic | SubchainKind::Other
    )
}

/// Choose a representative kind for a connected HET component.
///
/// `members` is the union-find component, sorted by ascending atom-index
/// start (per `build_entry`). When members span multiple kinds we promote
/// the component to `Organic` if any member is `Organic` — otherwise the
/// first member's kind wins. Single-kind components keep their original
/// kind unchanged.
fn pick_component_kind(members: &[u32], raw: &[RawRun]) -> SubchainKind {
    let mut has_organic = false;
    for &ri in members {
        if raw[ri as usize].kind == SubchainKind::Organic {
            has_organic = true;
            break;
        }
    }
    if has_organic {
        SubchainKind::Organic
    } else {
        raw[members[0] as usize].kind
    }
}

fn build_entry(atoms: &[Atom], raw: &[RawRun], members: &[u32]) -> SubchainEntry {
    debug_assert!(!members.is_empty());
    let mut members = members.to_vec();
    members.sort_by_key(|&ri| raw[ri as usize].start);

    let chain_id = raw[members[0] as usize].chain_id.clone();
    let kind = pick_component_kind(&members, raw);

    let atoms_layout = if is_contiguous_runs(&members, raw) {
        let start = raw[members[0] as usize].start;
        let end = raw[*members.last().unwrap() as usize].end;
        if start == 0 && end as usize == atoms.len() {
            SubchainAtoms::All {
                atom_count: atoms.len() as u32,
            }
        } else {
            SubchainAtoms::Range(start..end)
        }
    } else {
        let mut idxs: Vec<u32> = Vec::new();
        for &ri in &members {
            let r = &raw[ri as usize];
            for ai in r.start..r.end {
                idxs.push(ai);
            }
        }
        idxs.sort_unstable();
        SubchainAtoms::Indices(idxs)
    };

    let label = build_label(atoms, &atoms_layout, kind);

    SubchainEntry {
        chain_id,
        kind,
        label,
        atoms: atoms_layout,
    }
}

fn is_contiguous_runs(members: &[u32], raw: &[RawRun]) -> bool {
    for w in members.windows(2) {
        let a = &raw[w[0] as usize];
        let b = &raw[w[1] as usize];
        if a.end != b.start {
            return false;
        }
    }
    true
}

fn build_label(atoms: &[Atom], layout: &SubchainAtoms, kind: SubchainKind) -> SubchainLabel {
    match kind {
        SubchainKind::Biopolymer => SubchainLabel::Polymer,
        SubchainKind::Solvent => {
            let resn = first_atom_resn(atoms, layout).unwrap_or_default();
            SubchainLabel::Single(resn)
        }
        SubchainKind::Organic | SubchainKind::Inorganic | SubchainKind::Other => {
            let (residue_count, primary_resn, all_same_resn) = scan_residues(atoms, layout);
            if residue_count <= 1 || all_same_resn {
                SubchainLabel::Single(primary_resn)
            } else {
                SubchainLabel::Composite {
                    residue_count,
                    primary_resn,
                }
            }
        }
    }
}

fn first_atom_resn(atoms: &[Atom], layout: &SubchainAtoms) -> Option<String> {
    let idx = layout.min_atom_index()?;
    atoms.get(idx as usize).map(|a| a.residue.resn.to_string())
}

/// Walk the layout, returning `(residue_count, primary_resn, all_same_resn)`.
/// `primary_resn` is the residue name of the atom with the smallest atom
/// index.
fn scan_residues(atoms: &[Atom], layout: &SubchainAtoms) -> (u32, String, bool) {
    let mut prev: Option<&Atom> = None;
    let mut count: u32 = 0;
    let mut primary_resn: String = String::new();
    let mut all_same = true;

    let walk: Box<dyn Iterator<Item = u32>> = match layout {
        SubchainAtoms::All { atom_count } => Box::new(0..*atom_count),
        SubchainAtoms::Range(r) => Box::new(r.clone()),
        SubchainAtoms::Indices(v) => Box::new(v.iter().copied()),
    };

    for ai in walk {
        let Some(atom) = atoms.get(ai as usize) else {
            continue;
        };
        let new_residue = match prev {
            None => true,
            Some(p) => !atoms_same_residue(p, atom),
        };
        if new_residue {
            count += 1;
            if count == 1 {
                primary_resn = atom.residue.resn.to_string();
            } else if *atom.residue.resn != primary_resn {
                all_same = false;
            }
        }
        prev = Some(atom);
    }
    (count, primary_resn, all_same)
}

fn build_atom_to_entry_vec(n_atoms: usize, entries: &[SubchainEntry]) -> Vec<u32> {
    let mut map = vec![u32::MAX; n_atoms];
    for (ei, entry) in entries.iter().enumerate() {
        match &entry.atoms {
            SubchainAtoms::All { atom_count } => {
                let upper = (*atom_count as usize).min(n_atoms);
                for slot in map.iter_mut().take(upper) {
                    *slot = ei as u32;
                }
            }
            SubchainAtoms::Range(r) => {
                for ai in r.start..r.end {
                    if let Some(slot) = map.get_mut(ai as usize) {
                        *slot = ei as u32;
                    }
                }
            }
            SubchainAtoms::Indices(v) => {
                for &ai in v {
                    if let Some(slot) = map.get_mut(ai as usize) {
                        *slot = ei as u32;
                    }
                }
            }
        }
    }
    map
}

// ============================================================================
// Union-find
// ============================================================================

struct UnionFind {
    parent: Vec<usize>,
    rank: Vec<u8>,
}

impl UnionFind {
    fn new(n: usize) -> Self {
        UnionFind {
            parent: (0..n).collect(),
            rank: vec![0; n],
        }
    }

    fn find(&mut self, x: usize) -> usize {
        let mut root = x;
        while self.parent[root] != root {
            root = self.parent[root];
        }
        let mut cur = x;
        while self.parent[cur] != root {
            let next = self.parent[cur];
            self.parent[cur] = root;
            cur = next;
        }
        root
    }

    fn union(&mut self, a: usize, b: usize) {
        let ra = self.find(a);
        let rb = self.find(b);
        if ra == rb {
            return;
        }
        let (smaller, larger) = if self.rank[ra] < self.rank[rb] {
            (ra, rb)
        } else {
            (rb, ra)
        };
        self.parent[smaller] = larger;
        if self.rank[ra] == self.rank[rb] {
            self.rank[larger] += 1;
        }
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::AtomResidue;
    use crate::bond::BondOrder;
    use crate::element::Element;
    use crate::flags::AtomFlags;
    use crate::index::AtomIndex;
    use crate::molecule::ObjectMolecule;
    use std::sync::Arc;

    fn add_residue_atoms(
        mol: &mut ObjectMolecule,
        chain: &str,
        resn: &str,
        resv: i32,
        names: &[&str],
        flags: AtomFlags,
        hetatm: bool,
    ) -> Vec<AtomIndex> {
        let res = Arc::new(AtomResidue::from_parts(chain, resn, resv, ' ', ""));
        let mut idxs = Vec::new();
        for name in names {
            let mut atom = Atom::new(*name, Element::Carbon);
            atom.residue = res.clone();
            atom.state.hetatm = hetatm;
            atom.state.flags |= flags;
            idxs.push(mol.add_atom(atom));
        }
        idxs
    }

    fn link(mol: &mut ObjectMolecule, a: AtomIndex, b: AtomIndex) {
        let _ = mol.add_bond(a, b, BondOrder::Single);
    }

    #[test]
    fn partition_groups_glycan_by_bonds() {
        // Three HET residues NAG–NAG–BMA on chain E, glycosidically bonded.
        let mut mol = ObjectMolecule::new("glycan");
        let nag1 = add_residue_atoms(
            &mut mol,
            "E",
            "NAG",
            1,
            &["C1", "C2", "O4"],
            AtomFlags::ORGANIC,
            true,
        );
        let nag2 = add_residue_atoms(
            &mut mol,
            "E",
            "NAG",
            2,
            &["C1", "C2", "O4"],
            AtomFlags::ORGANIC,
            true,
        );
        let bma = add_residue_atoms(
            &mut mol,
            "E",
            "BMA",
            3,
            &["C1", "C2", "O4"],
            AtomFlags::ORGANIC,
            true,
        );
        link(&mut mol, nag1[2], nag2[0]); // O4 -- C1
        link(&mut mol, nag2[2], bma[0]);

        let entries = mol.subchain_partition().entries();
        assert_eq!(entries.len(), 1, "glycan should merge into one entry");
        let e = &entries[0];
        assert_eq!(e.chain_id, "E");
        assert_eq!(e.kind, SubchainKind::Organic);
        match &e.label {
            SubchainLabel::Composite {
                residue_count,
                primary_resn,
            } => {
                assert_eq!(*residue_count, 3);
                assert_eq!(primary_resn, "NAG");
            }
            other => panic!("expected Composite, got {:?}", other),
        }
        assert_eq!(e.atoms.len(), 9);
        let view = mol.subchains().next().unwrap();
        assert_eq!(view.label(), "NAG+2");
    }

    #[test]
    fn partition_splits_unbonded_het() {
        let mut mol = ObjectMolecule::new("ligs");
        add_residue_atoms(
            &mut mol,
            "A",
            "HEM",
            1,
            &["FE", "NA"],
            AtomFlags::ORGANIC,
            true,
        );
        add_residue_atoms(
            &mut mol,
            "A",
            "PO4",
            2,
            &["P", "O1"],
            AtomFlags::ORGANIC,
            true,
        );
        let entries = mol.subchain_partition().entries();
        assert_eq!(entries.len(), 2);
        assert!(matches!(entries[0].label, SubchainLabel::Single(ref s) if s == "HEM"));
        assert!(matches!(entries[1].label, SubchainLabel::Single(ref s) if s == "PO4"));
    }

    #[test]
    fn partition_does_not_cross_chain_boundary() {
        // Two NAG residues, one on chain E one on chain F, bonded across chains.
        let mut mol = ObjectMolecule::new("xchain");
        let e1 = add_residue_atoms(
            &mut mol,
            "E",
            "NAG",
            1,
            &["C1", "O4"],
            AtomFlags::ORGANIC,
            true,
        );
        let f1 = add_residue_atoms(
            &mut mol,
            "F",
            "NAG",
            1,
            &["C1", "O4"],
            AtomFlags::ORGANIC,
            true,
        );
        link(&mut mol, e1[1], f1[0]);
        let entries = mol.subchain_partition().entries();
        assert_eq!(entries.len(), 2);
        assert_eq!(entries[0].chain_id, "E");
        assert_eq!(entries[1].chain_id, "F");
    }

    #[test]
    fn partition_does_not_absorb_polymer() {
        // Asn–NAG glycosidic bond: polymer ASN must not pull NAG into the
        // biopolymer subchain, and vice versa.
        let mut mol = ObjectMolecule::new("asn-nag");
        let asn = add_residue_atoms(
            &mut mol,
            "A",
            "ASN",
            1,
            &["N", "CA", "C", "O", "ND2"],
            AtomFlags::PROTEIN | AtomFlags::POLYMER,
            false,
        );
        let nag = add_residue_atoms(
            &mut mol,
            "A",
            "NAG",
            2,
            &["C1", "C2"],
            AtomFlags::ORGANIC,
            true,
        );
        link(&mut mol, asn[4], nag[0]); // ND2 -- C1
        let entries = mol.subchain_partition().entries();
        assert_eq!(entries.len(), 2);
        let kinds: Vec<_> = entries.iter().map(|e| e.kind).collect();
        assert!(kinds.contains(&SubchainKind::Biopolymer));
        assert!(kinds.contains(&SubchainKind::Organic));
    }

    #[test]
    fn partition_merges_bonded_inorganic() {
        // Two distinct Inorganic residues on the same chain joined by a
        // covalent bond — should fuse into one component.
        let mut mol = ObjectMolecule::new("inorg");
        let fe = add_residue_atoms(&mut mol, "A", "FE", 1, &["FE"], AtomFlags::INORGANIC, true);
        let cl = add_residue_atoms(&mut mol, "A", "CL", 2, &["CL"], AtomFlags::INORGANIC, true);
        link(&mut mol, fe[0], cl[0]);

        let entries = mol.subchain_partition().entries();
        assert_eq!(entries.len(), 1, "bonded inorganic ions should merge");
        let e = &entries[0];
        assert_eq!(e.chain_id, "A");
        assert_eq!(e.kind, SubchainKind::Inorganic);
        assert_eq!(e.atoms.len(), 2);
        match &e.label {
            SubchainLabel::Composite {
                residue_count,
                primary_resn,
            } => {
                assert_eq!(*residue_count, 2);
                assert_eq!(primary_resn, "FE");
            }
            other => panic!("expected Composite, got {:?}", other),
        }
    }

    #[test]
    fn partition_merges_organic_inorganic() {
        // Organic ligand (HEM) covalently bonded to a separate Inorganic
        // ion residue (FE) on the same chain. The component must merge and
        // be promoted to `Organic` kind so the label uses Composite logic.
        let mut mol = ObjectMolecule::new("metalorg");
        let hem = add_residue_atoms(
            &mut mol,
            "A",
            "HEM",
            1,
            &["NA", "NB", "NC"],
            AtomFlags::ORGANIC,
            true,
        );
        let fe = add_residue_atoms(&mut mol, "A", "FE2", 2, &["FE"], AtomFlags::INORGANIC, true);
        link(&mut mol, hem[0], fe[0]);

        let entries = mol.subchain_partition().entries();
        assert_eq!(entries.len(), 1, "organic+inorganic should merge");
        let e = &entries[0];
        assert_eq!(e.kind, SubchainKind::Organic);
        assert_eq!(e.atoms.len(), 4);
        match &e.label {
            SubchainLabel::Composite {
                residue_count,
                primary_resn,
            } => {
                assert_eq!(*residue_count, 2);
                assert_eq!(primary_resn, "HEM");
            }
            other => panic!("expected Composite, got {:?}", other),
        }
    }

    #[test]
    fn partition_does_not_merge_solvent() {
        // A bogus bond between a water residue and an adjacent ligand must
        // NOT pull the water into the ligand component. Solvent stays its
        // own entry.
        let mut mol = ObjectMolecule::new("solv-lig");
        let hoh = add_residue_atoms(&mut mol, "A", "HOH", 1, &["O"], AtomFlags::SOLVENT, true);
        let lig = add_residue_atoms(
            &mut mol,
            "A",
            "ATP",
            2,
            &["PA", "PB"],
            AtomFlags::ORGANIC,
            true,
        );
        link(&mut mol, hoh[0], lig[0]);

        let entries = mol.subchain_partition().entries();
        assert_eq!(entries.len(), 2, "solvent must not merge with ligand");
        let kinds: Vec<_> = entries.iter().map(|e| e.kind).collect();
        assert!(kinds.contains(&SubchainKind::Solvent));
        assert!(kinds.contains(&SubchainKind::Organic));
    }

    #[test]
    fn partition_groups_solvent_by_resn_not_bonds() {
        // 100 disconnected HOH residues should yield ONE solvent entry,
        // not 100.
        let mut mol = ObjectMolecule::new("waters");
        for i in 0..100 {
            add_residue_atoms(&mut mol, "S", "HOH", i, &["O"], AtomFlags::SOLVENT, true);
        }
        let entries = mol.subchain_partition().entries();
        assert_eq!(entries.len(), 1);
        assert_eq!(entries[0].kind, SubchainKind::Solvent);
        assert_eq!(entries[0].atoms.len(), 100);
    }

    #[test]
    fn partition_collapses_single_entry_to_all() {
        // Single-chain biopolymer with no HETs: the partition should
        // collapse to one entry whose atoms layout is `All`.
        let mut mol = ObjectMolecule::new("solo");
        add_residue_atoms(
            &mut mol,
            "A",
            "ALA",
            1,
            &["N", "CA", "C", "O"],
            AtomFlags::PROTEIN | AtomFlags::POLYMER,
            false,
        );
        add_residue_atoms(
            &mut mol,
            "A",
            "GLY",
            2,
            &["N", "CA", "C", "O"],
            AtomFlags::PROTEIN | AtomFlags::POLYMER,
            false,
        );
        let entries = mol.subchain_partition().entries();
        assert_eq!(entries.len(), 1);
        match &entries[0].atoms {
            SubchainAtoms::All { atom_count } => assert_eq!(*atom_count, 8),
            other => panic!("expected SubchainAtoms::All, got {:?}", other),
        }
    }

    #[test]
    fn multi_entry_partition_does_not_collapse() {
        // Polymer + ion: two entries, neither uses `All`.
        let mut mol = ObjectMolecule::new("multi");
        add_residue_atoms(
            &mut mol,
            "A",
            "ALA",
            1,
            &["N", "CA", "C", "O"],
            AtomFlags::PROTEIN | AtomFlags::POLYMER,
            false,
        );
        add_residue_atoms(
            &mut mol,
            "A",
            "MG",
            100,
            &["MG"],
            AtomFlags::INORGANIC,
            true,
        );
        let entries = mol.subchain_partition().entries();
        assert_eq!(entries.len(), 2);
        for e in entries {
            assert!(
                !matches!(e.atoms, SubchainAtoms::All { .. }),
                "multi-entry partition must not use All; got {:?}",
                e.atoms
            );
        }
    }

    #[test]
    fn from_atoms_collapses_single_entry_to_all() {
        // Same single-chain biopolymer shape, but built via
        // `SubchainPartition::from_atoms`.
        let mut mol = ObjectMolecule::new("from-atoms");
        add_residue_atoms(
            &mut mol,
            "A",
            "ALA",
            1,
            &["N", "CA", "C", "O"],
            AtomFlags::PROTEIN | AtomFlags::POLYMER,
            false,
        );
        let part = SubchainPartition::from_atoms(mol.atoms_slice());
        assert_eq!(part.entries().len(), 1);
        match &part.entries()[0].atoms {
            SubchainAtoms::All { atom_count } => assert_eq!(*atom_count, 4),
            other => panic!("expected All, got {:?}", other),
        }
    }

    #[test]
    fn view_for_all_iterates_all_atoms() {
        let mut mol = ObjectMolecule::new("walk");
        add_residue_atoms(
            &mut mol,
            "A",
            "ALA",
            1,
            &["N", "CA", "C", "O"],
            AtomFlags::PROTEIN | AtomFlags::POLYMER,
            false,
        );
        let n = mol.atoms_slice().len();
        let view = mol.subchains().next().unwrap();
        assert!(view.is_contiguous());
        assert_eq!(view.iter().count(), n);
        let globals: Vec<u32> = view.iter_indexed().map(|(i, _)| i.0).collect();
        assert_eq!(globals, (0..n as u32).collect::<Vec<_>>());
    }

    #[test]
    fn atom_to_entry_all_maps_every_atom() {
        let mut mol = ObjectMolecule::new("map-all");
        add_residue_atoms(
            &mut mol,
            "A",
            "ALA",
            1,
            &["N", "CA", "C", "O"],
            AtomFlags::PROTEIN | AtomFlags::POLYMER,
            false,
        );
        let part = mol.subchain_partition();
        let n = mol.atoms_slice().len();
        for i in 0..n {
            assert_eq!(
                part.entry_for_atom(AtomIndex(i as u32)),
                Some(0),
                "atom {} should map to entry 0",
                i
            );
        }
    }

    #[test]
    fn single_subchain_partition_uses_all_entries_map() {
        let mut mol_single = ObjectMolecule::new("solo-map");
        add_residue_atoms(
            &mut mol_single,
            "A",
            "ALA",
            1,
            &["N", "CA"],
            AtomFlags::PROTEIN | AtomFlags::POLYMER,
            false,
        );
        match &mol_single.subchain_partition().atom_to_entry {
            AtomToEntryMap::AllEntries { atom_count } => assert_eq!(*atom_count, 2),
            other => panic!("expected AllEntries, got {:?}", other),
        }

        let mut mol_multi = ObjectMolecule::new("multi-map");
        add_residue_atoms(
            &mut mol_multi,
            "A",
            "ALA",
            1,
            &["N", "CA"],
            AtomFlags::PROTEIN | AtomFlags::POLYMER,
            false,
        );
        add_residue_atoms(
            &mut mol_multi,
            "A",
            "MG",
            100,
            &["MG"],
            AtomFlags::INORGANIC,
            true,
        );
        match &mol_multi.subchain_partition().atom_to_entry {
            AtomToEntryMap::Map(_) => {}
            other => panic!("expected Map, got {:?}", other),
        }
    }

    #[test]
    fn entry_for_atom_out_of_range() {
        // AllEntries variant
        let mut mol_single = ObjectMolecule::new("solo-oor");
        add_residue_atoms(
            &mut mol_single,
            "A",
            "ALA",
            1,
            &["N", "CA"],
            AtomFlags::PROTEIN | AtomFlags::POLYMER,
            false,
        );
        let p1 = mol_single.subchain_partition();
        assert_eq!(p1.entry_for_atom(AtomIndex(99)), None);

        // Map variant
        let mut mol_multi = ObjectMolecule::new("multi-oor");
        add_residue_atoms(
            &mut mol_multi,
            "A",
            "ALA",
            1,
            &["N", "CA"],
            AtomFlags::PROTEIN | AtomFlags::POLYMER,
            false,
        );
        add_residue_atoms(
            &mut mol_multi,
            "A",
            "MG",
            100,
            &["MG"],
            AtomFlags::INORGANIC,
            true,
        );
        let p2 = mol_multi.subchain_partition();
        assert_eq!(p2.entry_for_atom(AtomIndex(99)), None);
    }

    #[test]
    fn partition_invalidates_on_rebond() {
        // Use two distinct resns so that the predicate splits them into two
        // raw runs; only bond-graph merging can union them.
        let mut mol = ObjectMolecule::new("invalidate");
        let a = add_residue_atoms(
            &mut mol,
            "A",
            "NAG",
            1,
            &["C1", "O4"],
            AtomFlags::ORGANIC,
            true,
        );
        let b = add_residue_atoms(
            &mut mol,
            "A",
            "BMA",
            2,
            &["C1", "O4"],
            AtomFlags::ORGANIC,
            true,
        );
        // Without an inter-residue bond -> two entries.
        assert_eq!(mol.subchain_partition().entries().len(), 2);

        // Add a bond — invalidation should kick in via add_bond.
        link(&mut mol, a[1], b[0]);
        assert_eq!(mol.subchain_partition().entries().len(), 1);
    }
}
