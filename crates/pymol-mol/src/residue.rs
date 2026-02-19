//! Residue and chain utilities
//!
//! Provides types and functions for working with residues and chains.
//! Since PyMOL stores residue/chain information inline in atoms rather than
//! as separate structures, this module provides view types for convenient access.

use serde::{Deserialize, Serialize};

use crate::atom::Atom;
use crate::index::AtomIndex;
use std::ops::Range;

// Re-export iterators from the dedicated iterator module
pub use crate::iterator::{ChainIterator, ResidueIterator};

/// Key for uniquely identifying a residue within a molecule
///
/// A residue is uniquely identified by its chain, name, number, and insertion code.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct ResidueKey {
    /// Chain identifier
    pub chain: String,
    /// Residue name (e.g., "ALA", "GLY")
    pub resn: String,
    /// Residue sequence number
    pub resv: i32,
    /// Insertion code
    pub inscode: char,
}

impl ResidueKey {
    /// Create a new residue key
    pub fn new(
        chain: impl Into<String>,
        resn: impl Into<String>,
        resv: i32,
        inscode: char,
    ) -> Self {
        ResidueKey {
            chain: chain.into(),
            resn: resn.into(),
            resv,
            inscode,
        }
    }

    /// Create a residue key from an atom
    pub fn from_atom(atom: &Atom) -> Self {
        atom.residue.key.clone()
    }

    /// Get a display string for the residue (e.g., "A/ALA`1")
    pub fn display(&self) -> String {
        if self.inscode != ' ' && self.inscode != '\0' {
            format!("{}/{}'{}{}", self.chain, self.resn, self.resv, self.inscode)
        } else {
            format!("{}/{}'{}",self.chain, self.resn, self.resv)
        }
    }
}

impl std::fmt::Display for ResidueKey {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.display())
    }
}

/// A view into the atoms of a single residue
///
/// This is a borrowed view that does not own the atom data.
/// It provides convenient access to residue-level operations.
#[derive(Debug)]
pub struct ResidueView<'a> {
    /// The residue key
    pub key: ResidueKey,
    /// Slice of atoms in this residue
    pub atoms: &'a [Atom],
    /// Range of atom indices in the parent molecule
    pub atom_range: Range<usize>,
}

impl<'a> ResidueView<'a> {
    /// Create a new residue view
    pub fn new(key: ResidueKey, atoms: &'a [Atom], atom_range: Range<usize>) -> Self {
        ResidueView {
            key,
            atoms,
            atom_range,
        }
    }

    /// Get the chain ID
    #[inline]
    pub fn chain(&self) -> &str {
        &self.key.chain
    }

    /// Get the residue name
    #[inline]
    pub fn resn(&self) -> &str {
        &self.key.resn
    }

    /// Get the residue number
    #[inline]
    pub fn resv(&self) -> i32 {
        self.key.resv
    }

    /// Get the number of atoms in this residue
    #[inline]
    pub fn len(&self) -> usize {
        self.atoms.len()
    }

    /// Check if the residue is empty
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.atoms.is_empty()
    }

    /// Get an atom by index within the residue
    pub fn get(&self, index: usize) -> Option<&'a Atom> {
        self.atoms.get(index)
    }

    /// Iterate over atoms in the residue
    pub fn iter(&self) -> impl Iterator<Item = &'a Atom> {
        self.atoms.iter()
    }

    /// Iterate over atom indices and atoms
    pub fn iter_indexed(&self) -> impl Iterator<Item = (AtomIndex, &'a Atom)> {
        self.atom_range
            .clone()
            .zip(self.atoms.iter())
            .map(|(idx, atom)| (AtomIndex(idx as u32), atom))
    }

    /// Find an atom by name
    pub fn find_by_name(&self, name: &str) -> Option<(AtomIndex, &'a Atom)> {
        for (idx, atom) in self.iter_indexed() {
            if &*atom.name == name {
                return Some((idx, atom));
            }
        }
        None
    }

    /// Get the C-alpha atom (for proteins)
    pub fn ca(&self) -> Option<(AtomIndex, &'a Atom)> {
        self.find_by_name("CA")
    }

    /// Check if this is a protein residue
    pub fn is_protein(&self) -> bool {
        self.atoms
            .first()
            .map(|a| a.state.flags.contains(crate::flags::AtomFlags::PROTEIN))
            .unwrap_or(false)
    }

    /// Check if this is a nucleic acid residue
    pub fn is_nucleic(&self) -> bool {
        self.atoms
            .first()
            .map(|a| a.state.flags.contains(crate::flags::AtomFlags::NUCLEIC))
            .unwrap_or(false)
    }
}

/// A view into the atoms of a single chain
#[derive(Debug)]
pub struct ChainView<'a> {
    /// Chain identifier
    pub chain_id: String,
    /// Slice of atoms in this chain
    pub atoms: &'a [Atom],
    /// Range of atom indices in the parent molecule
    pub atom_range: Range<usize>,
}

impl<'a> ChainView<'a> {
    /// Create a new chain view
    pub fn new(chain_id: impl Into<String>, atoms: &'a [Atom], atom_range: Range<usize>) -> Self {
        ChainView {
            chain_id: chain_id.into(),
            atoms,
            atom_range,
        }
    }

    /// Get the chain ID
    #[inline]
    pub fn id(&self) -> &str {
        &self.chain_id
    }

    /// Get the number of atoms in this chain
    #[inline]
    pub fn len(&self) -> usize {
        self.atoms.len()
    }

    /// Check if the chain is empty
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.atoms.is_empty()
    }

    /// Iterate over atoms in the chain
    pub fn iter(&self) -> impl Iterator<Item = &'a Atom> {
        self.atoms.iter()
    }

    /// Iterate over atom indices and atoms
    pub fn iter_indexed(&self) -> impl Iterator<Item = (AtomIndex, &'a Atom)> {
        self.atom_range
            .clone()
            .zip(self.atoms.iter())
            .map(|(idx, atom)| (AtomIndex(idx as u32), atom))
    }

    /// Iterate over residues in this chain
    pub fn residues(&self) -> ResidueIterator<'a> {
        ResidueIterator::new(self.atoms, self.atom_range.start)
    }
}

/// Check if two atoms are in the same residue
#[inline]
pub fn atoms_same_residue(a: &Atom, b: &Atom) -> bool {
    // Can use Arc pointer equality for efficiency if atoms share the same residue
    std::sync::Arc::ptr_eq(&a.residue, &b.residue)
        || (a.residue.chain == b.residue.chain
            && a.residue.resv == b.residue.resv
            && a.residue.inscode == b.residue.inscode
            && a.residue.resn == b.residue.resn)
}

/// Check if two atoms are in the same chain
#[inline]
pub fn atoms_same_chain(a: &Atom, b: &Atom) -> bool {
    a.residue.chain == b.residue.chain
}

/// Check if two atoms are in the same segment
#[inline]
pub fn atoms_same_segment(a: &Atom, b: &Atom) -> bool {
    a.residue.segi == b.residue.segi
}

/// Convert 3-letter amino acid code to 1-letter code
pub fn three_to_one(resn: &str) -> Option<char> {
    match resn {
        "ALA" => Some('A'), "ARG" => Some('R'), "ASN" => Some('N'),
        "ASP" => Some('D'), "CYS" => Some('C'), "GLN" => Some('Q'),
        "GLU" => Some('E'), "GLY" => Some('G'), "HIS" => Some('H'),
        "ILE" => Some('I'), "LEU" => Some('L'), "LYS" => Some('K'),
        "MET" => Some('M'), "PHE" => Some('F'), "PRO" => Some('P'),
        "SER" => Some('S'), "THR" => Some('T'), "TRP" => Some('W'),
        "TYR" => Some('Y'), "VAL" => Some('V'),
        // Histidine protonation variants
        "HID" | "HIE" | "HIP" => Some('H'),
        // Cysteine variants
        "CYX" => Some('C'),
        // Non-standard but common
        "MSE" => Some('M'), "SEC" => Some('U'), "PYL" => Some('O'),
        // Phosphorylated amino acids
        "SEP" => Some('S'), "TPO" => Some('T'), "PTR" => Some('Y'),
        // Methylated / acetylated / other common modifications
        "MLY" | "M3L" | "MLZ" | "KCX" => Some('K'),
        "OCS" | "CSO" | "CSS" | "CSD" | "CME" | "SCH" | "SMC" | "CSX" => Some('C'),
        "NEP" => Some('H'),
        "HSD" | "HSE" | "HSP" => Some('H'),
        "ASH" => Some('D'),
        "GLH" => Some('E'),
        "LYN" => Some('K'),
        "ARN" => Some('R'),
        "TYS" => Some('Y'),
        "FME" => Some('M'),
        "SAM" => Some('M'),
        "HYP" => Some('P'),
        "UNK" => Some('X'),
        _ => None,
    }
}

/// Convert nucleotide residue name to single character
pub fn nucleotide_to_char(resn: &str) -> Option<char> {
    match resn {
        // Standard bases
        "DA" | "A" | "ADE" => Some('A'),
        "DT" | "T" | "THY" => Some('T'),
        "DG" | "G" | "GUA" => Some('G'),
        "DC" | "C" | "CYT" => Some('C'),
        "U" | "URA" => Some('U'),
        "DI" | "I" => Some('I'),
        // Generic
        "N" => Some('n'),
        "DN" => Some('n'),
        // Modified RNA — map to lowercase parent base or 'n'
        "PSU" | "H2U" | "5MU" | "4SU" | "5BU" | "BRU" | "CBR" => Some('u'),
        "5MC" | "OMC" => Some('c'),
        "OMG" | "M2G" | "7MG" | "2MG" | "YYG" => Some('g'),
        "1MA" => Some('a'),
        // Modified DNA
        "5CM" => Some('c'),
        "8OG" => Some('g'),
        // More modified RNA/DNA
        "2MA" | "6MA" | "A23" | "DZ" => Some('a'),
        "IU" | "5IU" | "UFT" => Some('u'),
        "5IC" | "3MC" | "CFZ" => Some('c'),
        "G7M" | "O2G" | "RSQ" => Some('g'),
        "LCA" => Some('a'),
        "LCG" => Some('g'),
        "LCC" => Some('c'),
        "LCT" => Some('t'),
        _ => None,
    }
}

/// Get display character for any residue (amino acid, nucleotide, or '?' for unknown)
pub fn residue_to_char(resn: &str) -> char {
    three_to_one(resn)
        .or_else(|| nucleotide_to_char(resn))
        .unwrap_or('?')
}

/// Standard amino acid residue names (3-letter codes)
pub const AMINO_ACIDS: &[&str] = &[
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    // Histidine protonation variants
    "HID", "HIE", "HIP", "HSP", "HSD", "HSE",
    // Cysteine variants
    "CYX",
    // Non-standard but common
    "MSE", "SEC", "PYL",
    // Charged variants
    "ARGP", "ASPM", "GLUM", "LYSP",
];

/// Standard nucleotide residue names
pub const NUCLEOTIDES: &[&str] = &[
    // Standard DNA
    "DA", "DC", "DG", "DT", "DI",
    // Standard RNA
    "A", "C", "G", "U", "I",
    // Generic
    "N", "DN",
    // Alternative 3-letter names
    "ADE", "CYT", "GUA", "THY", "URA",
    // Common modified RNA (PDB/RCSB)
    "PSU", "5MC", "OMC", "OMG", "M2G", "5MU", "7MG", "2MG", "H2U",
    "YYG", "1MA", "4SU", "5BU", "BRU", "CBR",
    // Modified DNA
    "5CM", "8OG",
];

/// Check if a residue name is a standard amino acid
pub fn is_amino_acid(resn: &str) -> bool {
    AMINO_ACIDS.contains(&resn)
}

/// Check if a residue name is a nucleotide
pub fn is_nucleotide(resn: &str) -> bool {
    NUCLEOTIDES.contains(&resn)
}

/// Common water residue names
pub const WATER_NAMES: &[&str] = &["HOH", "WAT", "H2O", "DOD", "TIP", "TIP3", "SPC", "SOL"];

/// Check if a residue name is water
pub fn is_water(resn: &str) -> bool {
    WATER_NAMES.contains(&resn)
}

/// Common ion residue names (GROMACS / CHARMM / AMBER / OPLS conventions)
pub const ION_NAMES: &[&str] = &[
    // Monatomic cations
    "NA", "NA+", "K", "K+", "CA2", "MG2",
    "ZN", "ZN2", "FE", "FE2", "CU", "CU2",
    "MN", "MN2", "NI", "CO", "CD",
    "LI", "LI+", "RB", "CS", "BA", "SR",
    // Monatomic anions
    "CL", "CL-", "BR", "I-",
    // CHARMM / OPLS force field variants
    "SOD", "POT", "CLA", "CAL",
];

/// Check if a residue name is a known ion
pub fn is_ion(resn: &str) -> bool {
    ION_NAMES.contains(&resn)
}

/// Common lipid residue names (CHARMM36 / AMBER / Slipids / Martini conventions)
pub const LIPID_NAMES: &[&str] = &[
    // Phosphatidylcholines
    "POPC", "DPPC", "DMPC", "DSPC", "DOPC", "DLPC", "DAPC", "DEPC",
    // Phosphatidylethanolamines
    "POPE", "DPPE", "DMPE", "DSPE", "DOPE", "DLPE",
    // Phosphatidylglycerols / serines / inositols
    "POPG", "DPPG", "DMPG", "DOPG", "POPS", "DPPS", "DOPS", "POPI",
    // Sphingomyelins
    "PSM", "SSM", "NSM",
    // Cholesterol
    "CHOL", "CHL1", "CLR",
    // Ceramides / Cardiolipin
    "CER", "CERA", "CDL1", "CDL2",
];

/// Check if a residue name is a known lipid
pub fn is_lipid(resn: &str) -> bool {
    LIPID_NAMES.contains(&resn)
}

/// N/C-terminal capping groups that are chemically adjacent to the polymer
/// but have no standard 1-letter code. Always display in bracket notation.
pub const CAPPING_GROUPS: &[&str] = &[
    "ACE", // acetyl cap (N-terminus)
    "NME", // N-methylamide cap (C-terminus)
    "NHH", // amino cap (rare)
    "FOR", // formyl cap
    "NH2", // amide cap (C-terminus)
];

/// Check if a residue name is a capping group
pub fn is_capping_group(resn: &str) -> bool {
    CAPPING_GROUPS.contains(&resn)
}

/// Residue category for chain assignment
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ResidueCategory {
    Protein,
    Nucleic,
    Solvent,
    Ion,
    Lipid,
    Other,
}

/// Classify a residue by its name into a category
pub(crate) fn classify_residue(resn: &str) -> ResidueCategory {
    if is_amino_acid(resn) {
        ResidueCategory::Protein
    } else if is_nucleotide(resn) {
        ResidueCategory::Nucleic
    } else if is_water(resn) {
        ResidueCategory::Solvent
    } else if is_ion(resn) {
        ResidueCategory::Ion
    } else if is_lipid(resn) {
        ResidueCategory::Lipid
    } else {
        ResidueCategory::Other
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::AtomResidue;
    use crate::element::Element;
    use std::sync::Arc;

    fn make_test_atoms() -> Vec<Atom> {
        let mut atoms = Vec::new();

        // Chain A, residue 1 (ALA)
        let res_ala = Arc::new(AtomResidue::from_parts("A", "ALA", 1, ' ', ""));
        for name in &["N", "CA", "C", "O", "CB"] {
            let mut atom = Atom::new(*name, Element::Carbon);
            atom.residue = res_ala.clone();
            atoms.push(atom);
        }

        // Chain A, residue 2 (GLY)
        let res_gly = Arc::new(AtomResidue::from_parts("A", "GLY", 2, ' ', ""));
        for name in &["N", "CA", "C", "O"] {
            let mut atom = Atom::new(*name, Element::Carbon);
            atom.residue = res_gly.clone();
            atoms.push(atom);
        }

        // Chain B, residue 1 (SER)
        let res_ser = Arc::new(AtomResidue::from_parts("B", "SER", 1, ' ', ""));
        for name in &["N", "CA", "C", "O", "CB", "OG"] {
            let mut atom = Atom::new(*name, Element::Carbon);
            atom.residue = res_ser.clone();
            atoms.push(atom);
        }

        atoms
    }

    #[test]
    fn test_residue_key() {
        let key = ResidueKey::new("A", "ALA", 1, ' ');
        assert_eq!(key.chain, "A");
        assert_eq!(key.resn, "ALA");
        assert_eq!(key.resv, 1);
    }

    #[test]
    fn test_residue_iterator() {
        let atoms = make_test_atoms();
        let residues: Vec<_> = ResidueIterator::new(&atoms, 0).collect();

        assert_eq!(residues.len(), 3);
        assert_eq!(residues[0].resn(), "ALA");
        assert_eq!(residues[0].len(), 5);
        assert_eq!(residues[1].resn(), "GLY");
        assert_eq!(residues[1].len(), 4);
        assert_eq!(residues[2].resn(), "SER");
        assert_eq!(residues[2].len(), 6);
    }

    #[test]
    fn test_chain_iterator() {
        let atoms = make_test_atoms();
        let chains: Vec<_> = ChainIterator::new(&atoms, 0).collect();

        assert_eq!(chains.len(), 2);
        assert_eq!(chains[0].id(), "A");
        assert_eq!(chains[0].len(), 9);
        assert_eq!(chains[1].id(), "B");
        assert_eq!(chains[1].len(), 6);
    }

    #[test]
    fn test_residue_view_find() {
        let atoms = make_test_atoms();
        let residue = ResidueIterator::new(&atoms, 0).next().unwrap();

        let (idx, ca) = residue.ca().unwrap();
        assert_eq!(&*ca.name, "CA");
        assert_eq!(idx.as_usize(), 1);
    }

    #[test]
    fn test_atoms_same_residue() {
        let atoms = make_test_atoms();
        assert!(atoms_same_residue(&atoms[0], &atoms[1]));
        assert!(!atoms_same_residue(&atoms[0], &atoms[5]));
    }

    #[test]
    fn test_atoms_same_residue_shared() {
        let atoms = make_test_atoms();
        // Atoms 0 and 1 share the same Rc<AtomResidue>
        assert!(Arc::ptr_eq(&atoms[0].residue, &atoms[1].residue));
        assert!(atoms_same_residue(&atoms[0], &atoms[1]));
    }

    #[test]
    fn test_is_amino_acid() {
        assert!(is_amino_acid("ALA"));
        assert!(is_amino_acid("GLY"));
        assert!(!is_amino_acid("HOH"));
    }

    #[test]
    fn test_is_water() {
        assert!(is_water("HOH"));
        assert!(is_water("WAT"));
        assert!(!is_water("ALA"));
    }

    #[test]
    fn test_three_to_one_standard() {
        assert_eq!(three_to_one("ALA"), Some('A'));
        assert_eq!(three_to_one("ARG"), Some('R'));
        assert_eq!(three_to_one("ASN"), Some('N'));
        assert_eq!(three_to_one("ASP"), Some('D'));
        assert_eq!(three_to_one("CYS"), Some('C'));
        assert_eq!(three_to_one("GLN"), Some('Q'));
        assert_eq!(three_to_one("GLU"), Some('E'));
        assert_eq!(three_to_one("GLY"), Some('G'));
        assert_eq!(three_to_one("HIS"), Some('H'));
        assert_eq!(three_to_one("ILE"), Some('I'));
        assert_eq!(three_to_one("LEU"), Some('L'));
        assert_eq!(three_to_one("LYS"), Some('K'));
        assert_eq!(three_to_one("MET"), Some('M'));
        assert_eq!(three_to_one("PHE"), Some('F'));
        assert_eq!(three_to_one("PRO"), Some('P'));
        assert_eq!(three_to_one("SER"), Some('S'));
        assert_eq!(three_to_one("THR"), Some('T'));
        assert_eq!(three_to_one("TRP"), Some('W'));
        assert_eq!(three_to_one("TYR"), Some('Y'));
        assert_eq!(three_to_one("VAL"), Some('V'));
    }

    #[test]
    fn test_three_to_one_variants() {
        // Histidine protonation variants
        assert_eq!(three_to_one("HID"), Some('H'));
        assert_eq!(three_to_one("HIE"), Some('H'));
        assert_eq!(three_to_one("HIP"), Some('H'));
        // Cysteine variant
        assert_eq!(three_to_one("CYX"), Some('C'));
        // Non-standard
        assert_eq!(three_to_one("MSE"), Some('M'));
        assert_eq!(three_to_one("SEC"), Some('U'));
        assert_eq!(three_to_one("PYL"), Some('O'));
    }

    #[test]
    fn test_three_to_one_unknown() {
        assert_eq!(three_to_one("HOH"), None);
        assert_eq!(three_to_one("ZZZ"), None);
        assert_eq!(three_to_one(""), None);
    }

    #[test]
    fn test_nucleotide_to_char() {
        // DNA
        assert_eq!(nucleotide_to_char("DA"), Some('A'));
        assert_eq!(nucleotide_to_char("DT"), Some('T'));
        assert_eq!(nucleotide_to_char("DG"), Some('G'));
        assert_eq!(nucleotide_to_char("DC"), Some('C'));
        assert_eq!(nucleotide_to_char("DI"), Some('I'));
        // RNA
        assert_eq!(nucleotide_to_char("A"), Some('A'));
        assert_eq!(nucleotide_to_char("U"), Some('U'));
        assert_eq!(nucleotide_to_char("G"), Some('G'));
        assert_eq!(nucleotide_to_char("C"), Some('C'));
        assert_eq!(nucleotide_to_char("I"), Some('I'));
        // Alternative names
        assert_eq!(nucleotide_to_char("ADE"), Some('A'));
        assert_eq!(nucleotide_to_char("THY"), Some('T'));
        assert_eq!(nucleotide_to_char("GUA"), Some('G'));
        assert_eq!(nucleotide_to_char("CYT"), Some('C'));
        assert_eq!(nucleotide_to_char("URA"), Some('U'));
        // Generic
        assert_eq!(nucleotide_to_char("N"), Some('n'));
        assert_eq!(nucleotide_to_char("DN"), Some('n'));
        // Modified bases
        assert_eq!(nucleotide_to_char("PSU"), Some('u'));
        assert_eq!(nucleotide_to_char("5MC"), Some('c'));
        assert_eq!(nucleotide_to_char("OMG"), Some('g'));
        assert_eq!(nucleotide_to_char("1MA"), Some('a'));
        assert_eq!(nucleotide_to_char("8OG"), Some('g'));
        // Unknown
        assert_eq!(nucleotide_to_char("HOH"), None);
    }

    #[test]
    fn test_residue_to_char() {
        assert_eq!(residue_to_char("ALA"), 'A');
        assert_eq!(residue_to_char("DA"), 'A');
        assert_eq!(residue_to_char("HOH"), '?');
        assert_eq!(residue_to_char("UNK"), 'X');
    }

    #[test]
    fn test_is_ion() {
        assert!(is_ion("NA"));
        assert!(is_ion("CL"));
        assert!(is_ion("MG2"));
        assert!(is_ion("SOD"));
        assert!(is_ion("CLA"));
        assert!(!is_ion("ALA"));
        assert!(!is_ion("SOL"));
    }

    #[test]
    fn test_is_lipid() {
        assert!(is_lipid("POPC"));
        assert!(is_lipid("CHOL"));
        assert!(is_lipid("DPPC"));
        assert!(is_lipid("CHL1"));
        assert!(!is_lipid("ALA"));
        assert!(!is_lipid("SOL"));
        assert!(!is_lipid("NA"));
    }

    #[test]
    fn test_three_to_one_noncanonical() {
        assert_eq!(three_to_one("SEP"), Some('S'));
        assert_eq!(three_to_one("TPO"), Some('T'));
        assert_eq!(three_to_one("PTR"), Some('Y'));
        assert_eq!(three_to_one("MLY"), Some('K'));
        assert_eq!(three_to_one("HYP"), Some('P'));
        assert_eq!(three_to_one("UNK"), Some('X'));
    }

    #[test]
    fn test_is_capping_group() {
        assert!(is_capping_group("ACE"));
        assert!(is_capping_group("NME"));
        assert!(!is_capping_group("ALA"));
        assert!(!is_capping_group("HOH"));
    }

    #[test]
    fn test_classify_residue() {
        assert_eq!(classify_residue("ALA"), ResidueCategory::Protein);
        assert_eq!(classify_residue("DA"), ResidueCategory::Nucleic);
        assert_eq!(classify_residue("SOL"), ResidueCategory::Solvent);
        assert_eq!(classify_residue("NA"), ResidueCategory::Ion);
        assert_eq!(classify_residue("POPC"), ResidueCategory::Lipid);
        assert_eq!(classify_residue("UNK"), ResidueCategory::Other);
    }
}
