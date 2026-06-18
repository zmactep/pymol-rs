//! Bond generation and bond order assignment
//!
//! Extracted from molecule.rs — provides distance-based bond generation
//! using a spatial hash grid and protein-specific bond order assignment.

use std::collections::HashMap;
use std::ops::Range;

use lin_alg::f32::Vec3;

use crate::bond::{Bond, BondOrder};
use crate::ccd;
use crate::index::{AtomIndex, BondIndex};
use crate::molecule::ObjectMolecule;
use crate::residue::{atoms_same_residue, is_ion, is_water};
use crate::spatial::SpatialGrid;

/// Default additive bond tolerance in Angstroms.
///
/// A bond is generated when:  distance < cov_radius(a) + cov_radius(b) + tolerance
///
/// 0.45 Angstroms works well with Cordero 2008 covalent radii.
pub const DEFAULT_BOND_TOLERANCE: f32 = 0.45;

#[derive(Clone, Copy)]
struct ProteinTemplateBond {
    atom1: &'static str,
    atom2: &'static str,
    order: BondOrder,
}

macro_rules! aa_bonds {
    ($(($atom1:literal, $atom2:literal, $order:ident)),+ $(,)?) => {
        &[
            $(ProteinTemplateBond {
                atom1: $atom1,
                atom2: $atom2,
                order: BondOrder::$order,
            }),+
        ]
    };
}

static ALA_BONDS: &[ProteinTemplateBond] = aa_bonds![
    ("N", "CA", Single),
    ("CA", "C", Single),
    ("C", "O", Double),
    ("CA", "CB", Single),
];

static ARG_BONDS: &[ProteinTemplateBond] = aa_bonds![
    ("N", "CA", Single),
    ("CA", "C", Single),
    ("C", "O", Double),
    ("CA", "CB", Single),
    ("CB", "CG", Single),
    ("CG", "CD", Single),
    ("CD", "NE", Single),
    ("NE", "CZ", Single),
    ("CZ", "NH1", Double),
    ("CZ", "NH2", Single),
];

static ASN_BONDS: &[ProteinTemplateBond] = aa_bonds![
    ("N", "CA", Single),
    ("CA", "C", Single),
    ("C", "O", Double),
    ("CA", "CB", Single),
    ("CB", "CG", Single),
    ("CG", "OD1", Double),
    ("CG", "ND2", Single),
];

static ASP_BONDS: &[ProteinTemplateBond] = aa_bonds![
    ("N", "CA", Single),
    ("CA", "C", Single),
    ("C", "O", Double),
    ("CA", "CB", Single),
    ("CB", "CG", Single),
    ("CG", "OD1", Double),
    ("CG", "OD2", Single),
];

static CYS_BONDS: &[ProteinTemplateBond] = aa_bonds![
    ("N", "CA", Single),
    ("CA", "C", Single),
    ("C", "O", Double),
    ("CA", "CB", Single),
    ("CB", "SG", Single),
];

static GLN_BONDS: &[ProteinTemplateBond] = aa_bonds![
    ("N", "CA", Single),
    ("CA", "C", Single),
    ("C", "O", Double),
    ("CA", "CB", Single),
    ("CB", "CG", Single),
    ("CG", "CD", Single),
    ("CD", "OE1", Double),
    ("CD", "NE2", Single),
];

static GLU_BONDS: &[ProteinTemplateBond] = aa_bonds![
    ("N", "CA", Single),
    ("CA", "C", Single),
    ("C", "O", Double),
    ("CA", "CB", Single),
    ("CB", "CG", Single),
    ("CG", "CD", Single),
    ("CD", "OE1", Double),
    ("CD", "OE2", Single),
];

static GLY_BONDS: &[ProteinTemplateBond] =
    aa_bonds![("N", "CA", Single), ("CA", "C", Single), ("C", "O", Double),];

static HIS_BONDS: &[ProteinTemplateBond] = aa_bonds![
    ("N", "CA", Single),
    ("CA", "C", Single),
    ("C", "O", Double),
    ("CA", "CB", Single),
    ("CB", "CG", Single),
    ("CG", "ND1", Single),
    ("ND1", "CE1", Double),
    ("CE1", "NE2", Single),
    ("NE2", "CD2", Double),
    ("CD2", "CG", Single),
];

static ILE_BONDS: &[ProteinTemplateBond] = aa_bonds![
    ("N", "CA", Single),
    ("CA", "C", Single),
    ("C", "O", Double),
    ("CA", "CB", Single),
    ("CB", "CG1", Single),
    ("CG1", "CD1", Single),
    ("CB", "CG2", Single),
];

static LEU_BONDS: &[ProteinTemplateBond] = aa_bonds![
    ("N", "CA", Single),
    ("CA", "C", Single),
    ("C", "O", Double),
    ("CA", "CB", Single),
    ("CB", "CG", Single),
    ("CG", "CD1", Single),
    ("CG", "CD2", Single),
];

static LYS_BONDS: &[ProteinTemplateBond] = aa_bonds![
    ("N", "CA", Single),
    ("CA", "C", Single),
    ("C", "O", Double),
    ("CA", "CB", Single),
    ("CB", "CG", Single),
    ("CG", "CD", Single),
    ("CD", "CE", Single),
    ("CE", "NZ", Single),
];

static MET_BONDS: &[ProteinTemplateBond] = aa_bonds![
    ("N", "CA", Single),
    ("CA", "C", Single),
    ("C", "O", Double),
    ("CA", "CB", Single),
    ("CB", "CG", Single),
    ("CG", "SD", Single),
    ("SD", "CE", Single),
];

static PHE_BONDS: &[ProteinTemplateBond] = aa_bonds![
    ("N", "CA", Single),
    ("CA", "C", Single),
    ("C", "O", Double),
    ("CA", "CB", Single),
    ("CB", "CG", Single),
    ("CG", "CD1", Double),
    ("CD1", "CE1", Single),
    ("CE1", "CZ", Double),
    ("CZ", "CE2", Single),
    ("CE2", "CD2", Double),
    ("CD2", "CG", Single),
];

static PRO_BONDS: &[ProteinTemplateBond] = aa_bonds![
    ("N", "CA", Single),
    ("CA", "C", Single),
    ("C", "O", Double),
    ("CA", "CB", Single),
    ("CB", "CG", Single),
    ("CG", "CD", Single),
    ("CD", "N", Single),
];

static SER_BONDS: &[ProteinTemplateBond] = aa_bonds![
    ("N", "CA", Single),
    ("CA", "C", Single),
    ("C", "O", Double),
    ("CA", "CB", Single),
    ("CB", "OG", Single),
];

static THR_BONDS: &[ProteinTemplateBond] = aa_bonds![
    ("N", "CA", Single),
    ("CA", "C", Single),
    ("C", "O", Double),
    ("CA", "CB", Single),
    ("CB", "OG1", Single),
    ("CB", "CG2", Single),
];

static TRP_BONDS: &[ProteinTemplateBond] = aa_bonds![
    ("N", "CA", Single),
    ("CA", "C", Single),
    ("C", "O", Double),
    ("CA", "CB", Single),
    ("CB", "CG", Single),
    ("CG", "CD1", Double),
    ("CD1", "NE1", Single),
    ("NE1", "CE2", Single),
    ("CE2", "CD2", Single),
    ("CD2", "CG", Single),
    ("CD2", "CE3", Double),
    ("CE3", "CZ3", Single),
    ("CZ3", "CH2", Double),
    ("CH2", "CZ2", Single),
    ("CZ2", "CE2", Double),
];

static TYR_BONDS: &[ProteinTemplateBond] = aa_bonds![
    ("N", "CA", Single),
    ("CA", "C", Single),
    ("C", "O", Double),
    ("CA", "CB", Single),
    ("CB", "CG", Single),
    ("CG", "CD1", Double),
    ("CD1", "CE1", Single),
    ("CE1", "CZ", Double),
    ("CZ", "CE2", Single),
    ("CE2", "CD2", Double),
    ("CD2", "CG", Single),
    ("CZ", "OH", Single),
];

static VAL_BONDS: &[ProteinTemplateBond] = aa_bonds![
    ("N", "CA", Single),
    ("CA", "C", Single),
    ("C", "O", Double),
    ("CA", "CB", Single),
    ("CB", "CG1", Single),
    ("CB", "CG2", Single),
];

fn amino_acid_template(resn: &str) -> Option<&'static [ProteinTemplateBond]> {
    match resn {
        "ALA" => Some(ALA_BONDS),
        "ARG" => Some(ARG_BONDS),
        "ASN" => Some(ASN_BONDS),
        "ASP" => Some(ASP_BONDS),
        "CYS" => Some(CYS_BONDS),
        "GLN" => Some(GLN_BONDS),
        "GLU" => Some(GLU_BONDS),
        "GLY" => Some(GLY_BONDS),
        "HIS" => Some(HIS_BONDS),
        "ILE" => Some(ILE_BONDS),
        "LEU" => Some(LEU_BONDS),
        "LYS" => Some(LYS_BONDS),
        "MET" => Some(MET_BONDS),
        "PHE" => Some(PHE_BONDS),
        "PRO" => Some(PRO_BONDS),
        "SER" => Some(SER_BONDS),
        "THR" => Some(THR_BONDS),
        "TRP" => Some(TRP_BONDS),
        "TYR" => Some(TYR_BONDS),
        "VAL" => Some(VAL_BONDS),
        _ => None,
    }
}

fn push_bond(mol: &mut ObjectMolecule, a1: AtomIndex, a2: AtomIndex, order: BondOrder) {
    let bond = Bond::new(a1, a2, order);
    let bond_index = BondIndex(mol.bonds.len() as u32);
    mol.bonds.push(bond);
    mol.atom_bonds[a1.as_usize()].push(bond_index);
    mol.atom_bonds[a2.as_usize()].push(bond_index);
    mol.atoms[a1.as_usize()].state.bonded = true;
    mol.atoms[a2.as_usize()].state.bonded = true;
}

fn push_bond_if_absent(mol: &mut ObjectMolecule, a1: AtomIndex, a2: AtomIndex, order: BondOrder) {
    if mol.find_bond(a1, a2).is_none() {
        push_bond(mol, a1, a2, order);
    }
}

// ============================================================================
// Protein Double Bond Lookup Table
// ============================================================================

/// Double bonds in standard amino acid residues.
/// Format: (residue_pattern, atom1, atom2)
/// - "*" matches all residues (for backbone carbonyl)
/// - Multiple residue names for variants
static PROTEIN_DOUBLE_BONDS: &[(&[&str], &str, &str)] = &[
    // Backbone carbonyl C=O (all amino acids)
    (&["*"], "C", "O"),
    // ARG: guanidinium group CZ=NH1
    (&["ARG", "ARGP"], "CZ", "NH1"),
    // ASP: carboxylate CG=OD1
    (&["ASP", "ASPM"], "CG", "OD1"),
    // ASN: amide CG=OD1
    (&["ASN"], "CG", "OD1"),
    // GLU: carboxylate CD=OE1
    (&["GLU", "GLUM"], "CD", "OE1"),
    // GLN: amide CD=OE1
    (&["GLN"], "CD", "OE1"),
    // HIS and variants: imidazole ring (Kekulé alternating pattern)
    (&["HIS", "HID", "HIE", "HIP"], "ND1", "CE1"),
    (&["HIS", "HID", "HIE", "HIP"], "NE2", "CD2"),
    // PHE: benzene ring (Kekulé alternating pattern)
    (&["PHE"], "CG", "CD1"),
    (&["PHE"], "CE1", "CZ"),
    (&["PHE"], "CE2", "CD2"),
    // TYR: benzene ring (same Kekulé pattern as PHE)
    (&["TYR"], "CG", "CD1"),
    (&["TYR"], "CE1", "CZ"),
    (&["TYR"], "CE2", "CD2"),
    // TRP: indole ring system (Kekulé alternating pattern)
    (&["TRP"], "CG", "CD1"),
    (&["TRP"], "CE2", "CZ2"),
    (&["TRP"], "CH2", "CZ3"),
    (&["TRP"], "CE3", "CD2"),
];

/// Get bond order for a bond between two atoms in a known protein residue
fn get_protein_bond_order(resn: &str, name1: &str, name2: &str) -> BondOrder {
    let is_bond =
        |a: &str, b: &str| -> bool { (name1 == a && name2 == b) || (name1 == b && name2 == a) };

    for &(residues, atom1, atom2) in PROTEIN_DOUBLE_BONDS {
        if is_bond(atom1, atom2) {
            let matches = residues.iter().any(|&r| r == "*" || r == resn);
            if matches {
                return BondOrder::Double;
            }
        }
    }

    BondOrder::Single
}

/// Apply canonical amino-acid intra-residue templates.
fn apply_amino_acid_template_bonds(mol: &mut ObjectMolecule, template_ids: &mut [u32]) {
    let residue_info: Vec<(String, Range<usize>)> = mol
        .residues()
        .filter(|r| {
            let is_hetatm = r.atoms.first().is_some_and(|a| a.state.hetatm);
            !is_hetatm && amino_acid_template(&r.key.resn).is_some()
        })
        .map(|r| (r.key.resn.clone(), r.atom_range.clone()))
        .collect();

    let mut next_template_id = 1u32;
    for (resn, atom_range) in residue_info {
        let Some(template) = amino_acid_template(&resn) else {
            continue;
        };

        let mut alt_labels: Vec<char> = atom_range
            .clone()
            .map(|idx| mol.atoms[idx].alt)
            .filter(|&c| c != ' ')
            .collect();
        alt_labels.sort_unstable();
        alt_labels.dedup();

        if alt_labels.is_empty() {
            alt_labels.push(' ');
        }

        let template_id = next_template_id;
        let mut applied = false;
        for &alt in &alt_labels {
            let Some(bond_pairs) =
                collect_amino_acid_template_pairs(mol, &atom_range, alt, template)
            else {
                continue;
            };

            for &(a1, a2, order) in &bond_pairs {
                push_bond_if_absent(mol, a1, a2, order);
                template_ids[a1.as_usize()] = template_id;
                template_ids[a2.as_usize()] = template_id;
            }
            applied = true;
        }

        if applied {
            next_template_id += 1;
        }
    }
}

fn collect_amino_acid_template_pairs(
    mol: &ObjectMolecule,
    atom_range: &Range<usize>,
    alt: char,
    template: &[ProteinTemplateBond],
) -> Option<Vec<(AtomIndex, AtomIndex, BondOrder)>> {
    let mut bond_pairs = Vec::with_capacity(template.len());
    for bond in template {
        let a1 = compatible_atom_index(mol, atom_range, alt, bond.atom1)?;
        let a2 = compatible_atom_index(mol, atom_range, alt, bond.atom2)?;
        bond_pairs.push((AtomIndex(a1 as u32), AtomIndex(a2 as u32), bond.order));
    }
    Some(bond_pairs)
}

fn compatible_atom_index(
    mol: &ObjectMolecule,
    atom_range: &Range<usize>,
    alt: char,
    name: &str,
) -> Option<usize> {
    atom_range.clone().find(|&idx| {
        let atom = &mol.atoms[idx];
        &*atom.name == name && (atom.alt == ' ' || atom.alt == alt)
    })
}

fn needs_distance_probe(template_id: u32, atom_name: &str) -> bool {
    template_id == 0 || matches!(atom_name, "N" | "C" | "SG")
}

/// Apply CCD template bonds for known HETATM residues.
///
/// For each residue with a CCD template whose heavy atoms all match,
/// creates bonds with correct bond orders. Marks bonded atoms in `ccd_bonded`
/// so the distance-based pass can skip intra-residue pairs.
fn apply_ccd_template_bonds(mol: &mut ObjectMolecule, ccd_bonded: &mut [bool]) {
    // Collect residue info first (can't iterate residues while mutating mol)
    let residue_info: Vec<(String, std::ops::Range<usize>)> = mol
        .residues()
        .filter(|r| {
            // Only HETATM residues that have a CCD template
            r.atoms.first().is_some_and(|a| a.state.hetatm)
                && ccd::get_template(&r.key.resn).is_some()
        })
        .map(|r| (r.key.resn.clone(), r.atom_range.clone()))
        .collect();

    for (resn, atom_range) in residue_info {
        let template = match ccd::get_template(&resn) {
            Some(t) => t,
            None => continue,
        };

        // Collect distinct non-blank altloc labels in this residue.
        // When altlocs are present, we must apply template bonds separately
        // per group (blank-alt atoms + one altloc label) to avoid collisions
        // in the name→index map (e.g., two "CA" atoms with alt A and B).
        let mut alt_labels: Vec<char> = atom_range
            .clone()
            .map(|idx| mol.atoms[idx].alt)
            .filter(|&c| c != ' ')
            .collect();
        alt_labels.sort_unstable();
        alt_labels.dedup();

        // If no altlocs, use a single pass with blank sentinel
        if alt_labels.is_empty() {
            alt_labels.push(' ');
        }

        let mut residue_bonded = false;

        for &alt in &alt_labels {
            // Build name→index map for atoms compatible with this altloc group:
            // atoms with alt==' ' (shared) plus atoms with this specific label.
            let name_to_idx: HashMap<String, usize> = atom_range
                .clone()
                .filter(|&idx| {
                    let a = mol.atoms[idx].alt;
                    a == ' ' || a == alt
                })
                .map(|idx| (mol.atoms[idx].name.to_string(), idx))
                .collect();

            if !template.validate_heavy_atoms(|n| name_to_idx.contains_key(n)) {
                continue;
            }

            let bond_pairs: Vec<(AtomIndex, AtomIndex, BondOrder)> = template
                .iter_bonds()
                .filter_map(|(a1_name, a2_name, order)| {
                    let &a1 = name_to_idx.get(a1_name)?;
                    let &a2 = name_to_idx.get(a2_name)?;
                    Some((AtomIndex(a1 as u32), AtomIndex(a2 as u32), order))
                })
                .collect();

            for (a1, a2, order) in bond_pairs {
                push_bond(mol, a1, a2, order);
            }

            residue_bonded = true;
        }

        // Mark all atoms in this residue as CCD-bonded
        if residue_bonded {
            for idx in atom_range {
                ccd_bonded[idx] = true;
            }
        }
    }
}

/// Generate bonds based on inter-atomic distances for a specific state.
///
/// Uses covalent radii with an additive tolerance:
///   bond if distance < cov_radius(a) + cov_radius(b) + tolerance
///
/// If CCD cache is loaded, known HETATM residues get template-based bonds
/// instead of distance-based (CCD pre-pass). Distance-based bonding is used
/// for all other atoms and for inter-residue bonds.
///
/// Spatial hash grid gives O(n) performance instead of O(n²) all-pairs.
pub(crate) fn generate_bonds_for_state(mol: &mut ObjectMolecule, tolerance: f32, state: usize) {
    if mol.coord_sets.get(state).is_none() {
        return;
    }

    let n_atoms = mol.atoms.len();
    if n_atoms < 2 {
        return;
    }

    let mut amino_acid_template_ids = vec![0u32; n_atoms];
    apply_amino_acid_template_bonds(mol, &mut amino_acid_template_ids);

    // CCD template bonds for known HETATM residues. This must run before
    // borrowing coord_set because it mutates the molecule.
    let ccd_bonded = if ccd::is_loaded() {
        let mut ccd_bonded = vec![false; n_atoms];
        apply_ccd_template_bonds(mol, &mut ccd_bonded);
        Some(ccd_bonded)
    } else {
        None
    };

    // Distance-based bonds for remaining atoms.
    let coord_set = &mol.coord_sets[state];

    struct AtomData {
        coord: Vec3,
        cov: f32,
        template_id: u32,
        distance_probe: bool,
    }

    let atom_data: Vec<AtomData> = (0..n_atoms)
        .map(|i| {
            let coord = coord_set
                .get_atom_coord(AtomIndex(i as u32))
                .unwrap_or(Vec3::new(0.0, 0.0, 0.0));
            let atom = &mol.atoms[i];
            let template_id = amino_acid_template_ids[i];
            AtomData {
                coord,
                cov: atom.cov_radius(),
                template_id,
                distance_probe: needs_distance_probe(template_id, &atom.name),
            }
        })
        .collect();

    let max_cov = atom_data.iter().map(|a| a.cov).fold(0.0f32, f32::max);
    let max_cutoff = 2.0 * max_cov + tolerance;

    let mut grid = SpatialGrid::with_capacity(max_cutoff, n_atoms);
    for (i, ad) in atom_data.iter().enumerate() {
        grid.insert(ad.coord, i);
    }

    let mut new_bonds: Vec<(AtomIndex, AtomIndex)> = Vec::with_capacity(n_atoms);
    let mut neighbors = Vec::new();

    for i in 0..n_atoms {
        let a1 = &atom_data[i];
        if !a1.distance_probe {
            continue;
        }

        grid.query_neighbors(a1.coord, &mut neighbors);

        for &j in &neighbors {
            if j == i {
                continue;
            }

            let a2 = &atom_data[j];
            if a2.distance_probe && j < i {
                continue;
            }

            if a1.template_id != 0 && a1.template_id == a2.template_id {
                continue;
            }

            let dx = a1.coord.x - a2.coord.x;
            let dy = a1.coord.y - a2.coord.y;
            let dz = a1.coord.z - a2.coord.z;
            let dist_sq = dx * dx + dy * dy + dz * dz;

            let threshold = a1.cov + a2.cov + tolerance;
            let threshold_sq = threshold * threshold;

            if dist_sq < threshold_sq {
                let atom_i = &mol.atoms[i];
                let atom_j = &mol.atoms[j];

                // Skip bonds between atoms with incompatible alternate locations.
                // Atoms with alt=' ' (no altloc) are compatible with everything.
                let alt_i = atom_i.alt;
                let alt_j = atom_j.alt;
                if alt_i != ' ' && alt_j != ' ' && alt_i != alt_j {
                    continue;
                }

                // For solvent/ion residues, only allow intra-residue bonds
                if (is_water(&atom_i.residue.resn)
                    || is_ion(&atom_i.residue.resn)
                    || is_water(&atom_j.residue.resn)
                    || is_ion(&atom_j.residue.resn))
                    && !atoms_same_residue(atom_i, atom_j)
                {
                    continue;
                }
                // Skip intra-residue bonds for CCD-handled residues
                let ccd_pair_bonded = ccd_bonded
                    .as_ref()
                    .is_some_and(|bonded| bonded[i] && bonded[j]);
                if ccd_pair_bonded && atoms_same_residue(atom_i, atom_j) {
                    continue;
                }
                new_bonds.push((AtomIndex(i as u32), AtomIndex(j as u32)));
            }
        }
    }

    for (a1, a2) in new_bonds {
        push_bond(mol, a1, a2, BondOrder::Single);
    }
}

/// Assign bond orders for known protein residues based on atom names.
pub(crate) fn assign_known_residue_bond_orders(mol: &mut ObjectMolecule) {
    use crate::residue::{atoms_same_residue, is_amino_acid};

    for bond in &mut mol.bonds {
        let atom1 = match mol.atoms.get(bond.atom1.as_usize()) {
            Some(a) => a,
            None => continue,
        };
        let atom2 = match mol.atoms.get(bond.atom2.as_usize()) {
            Some(a) => a,
            None => continue,
        };

        if !atoms_same_residue(atom1, atom2) {
            continue;
        }

        if !is_amino_acid(&atom1.residue.resn) {
            continue;
        }

        let name1 = &*atom1.name;
        let name2 = &*atom2.name;
        let resn = atom1.residue.resn.as_str();

        let order = get_protein_bond_order(resn, name1, name2);
        if order != BondOrder::Single {
            bond.order = order;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coordset::CoordSet;
    use crate::element::Element;
    use crate::molecule::ObjectMolecule;
    use crate::{Atom, AtomResidue};
    use lin_alg::f32::Vec3;
    use std::sync::Arc;

    fn residue_atom(name: &str, element: Element, resn: &str, alt: char) -> Atom {
        let mut atom = Atom::new(name, element);
        atom.alt = alt;
        atom.residue = Arc::new(AtomResidue::from_parts("A", resn, 1, ' ', ""));
        atom
    }

    fn add_residue_atoms(
        mol: &mut ObjectMolecule,
        resn: &str,
        atoms: &[(&str, Element, char, Vec3)],
    ) {
        let mut positions = Vec::with_capacity(atoms.len());
        for &(name, element, alt, pos) in atoms {
            mol.add_atom(residue_atom(name, element, resn, alt));
            positions.push(pos);
        }
        mol.add_coord_set(CoordSet::from_vec3(&positions));
    }

    fn atom_index(mol: &ObjectMolecule, name: &str, alt: char) -> AtomIndex {
        mol.atoms_indexed()
            .find(|(_, atom)| &*atom.name == name && atom.alt == alt)
            .map(|(idx, _)| idx)
            .unwrap_or_else(|| panic!("missing atom {name}/{alt}"))
    }

    fn bond_order(mol: &ObjectMolecule, atom1: &str, atom2: &str) -> BondOrder {
        let a1 = atom_index(mol, atom1, ' ');
        let a2 = atom_index(mol, atom2, ' ');
        let bond_idx = mol
            .find_bond(a1, a2)
            .unwrap_or_else(|| panic!("missing bond {atom1}-{atom2}"));
        mol.get_bond(bond_idx).expect("bond index is valid").order
    }

    fn count_bonds_between(mol: &ObjectMolecule, a1: AtomIndex, a2: AtomIndex) -> usize {
        mol.bonds().filter(|bond| bond.connects(a1, a2)).count()
    }

    #[test]
    fn test_amino_acid_template_creates_exact_heavy_topology() {
        let mut mol = ObjectMolecule::new("ala");
        let origin = Vec3::new(0.0, 0.0, 0.0);
        add_residue_atoms(
            &mut mol,
            "ALA",
            &[
                ("N", Element::Nitrogen, ' ', origin),
                ("CA", Element::Carbon, ' ', origin),
                ("C", Element::Carbon, ' ', origin),
                ("O", Element::Oxygen, ' ', origin),
                ("CB", Element::Carbon, ' ', origin),
            ],
        );

        mol.generate_bonds(DEFAULT_BOND_TOLERANCE);

        assert_eq!(mol.bond_count(), 4);
        assert_eq!(bond_order(&mol, "C", "O"), BondOrder::Double);
        assert!(
            mol.find_bond(atom_index(&mol, "N", ' '), atom_index(&mol, "C", ' '))
                .is_none(),
            "template-handled heavy atoms should not get extra distance bonds"
        );
    }

    #[test]
    fn test_incomplete_amino_acid_residue_uses_distance_fallback() {
        let mut mol = ObjectMolecule::new("ala_incomplete");
        add_residue_atoms(
            &mut mol,
            "ALA",
            &[
                ("N", Element::Nitrogen, ' ', Vec3::new(0.0, 0.0, 0.0)),
                ("CA", Element::Carbon, ' ', Vec3::new(1.45, 0.0, 0.0)),
                ("C", Element::Carbon, ' ', Vec3::new(2.90, 0.0, 0.0)),
                ("O", Element::Oxygen, ' ', Vec3::new(4.13, 0.0, 0.0)),
            ],
        );

        mol.generate_bonds(DEFAULT_BOND_TOLERANCE);

        assert_eq!(mol.bond_count(), 3);
        assert_eq!(
            bond_order(&mol, "C", "O"),
            BondOrder::Single,
            "incomplete residues should not receive template bond orders"
        );
    }

    #[test]
    fn test_hydrogen_and_oxt_stay_distance_based() {
        let mut mol = ObjectMolecule::new("ala_terminal");
        add_residue_atoms(
            &mut mol,
            "ALA",
            &[
                ("N", Element::Nitrogen, ' ', Vec3::new(0.0, 0.0, 0.0)),
                ("CA", Element::Carbon, ' ', Vec3::new(1.45, 0.0, 0.0)),
                ("C", Element::Carbon, ' ', Vec3::new(2.90, 0.0, 0.0)),
                ("O", Element::Oxygen, ' ', Vec3::new(3.50, 1.10, 0.0)),
                ("CB", Element::Carbon, ' ', Vec3::new(1.45, 1.54, 0.0)),
                ("H", Element::Hydrogen, ' ', Vec3::new(-0.90, 0.0, 0.0)),
                ("OXT", Element::Oxygen, ' ', Vec3::new(3.50, -1.10, 0.0)),
            ],
        );

        mol.generate_bonds(DEFAULT_BOND_TOLERANCE);

        assert_eq!(mol.bond_count(), 6);
        assert!(mol
            .find_bond(atom_index(&mol, "N", ' '), atom_index(&mol, "H", ' '))
            .is_some());
        assert!(mol
            .find_bond(atom_index(&mol, "C", ' '), atom_index(&mol, "OXT", ' '))
            .is_some());
    }

    #[test]
    fn test_amino_acid_template_altlocs_do_not_duplicate_blank_bonds() {
        let mut mol = ObjectMolecule::new("ala_altloc");
        let origin = Vec3::new(0.0, 0.0, 0.0);
        add_residue_atoms(
            &mut mol,
            "ALA",
            &[
                ("N", Element::Nitrogen, ' ', origin),
                ("CA", Element::Carbon, 'A', origin),
                ("CA", Element::Carbon, 'B', origin),
                ("C", Element::Carbon, ' ', origin),
                ("O", Element::Oxygen, ' ', origin),
                ("CB", Element::Carbon, 'A', origin),
                ("CB", Element::Carbon, 'B', origin),
            ],
        );

        mol.generate_bonds(DEFAULT_BOND_TOLERANCE);

        let c_idx = atom_index(&mol, "C", ' ');
        let o_idx = atom_index(&mol, "O", ' ');
        assert_eq!(mol.bond_count(), 7);
        assert_eq!(count_bonds_between(&mol, c_idx, o_idx), 1);
    }

    /// Fe–N coordination bond (~2.0 Å) must be detected with covalent radii.
    /// Fe(cov=1.32) + N(cov=0.71) + 0.45 = 2.48 Å threshold → 2.0 Å bond detected.
    #[test]
    fn test_fe_n_coordination_bond() {
        let mut mol = ObjectMolecule::new("heme");
        mol.add_atom(Atom::new("FE", Element::Iron));
        mol.add_atom(Atom::new("NA", Element::Nitrogen));

        // Fe–N distance ~2.0 Å (typical heme coordination)
        let cs = CoordSet::from_vec3(&[Vec3::new(0.0, 0.0, 0.0), Vec3::new(2.0, 0.0, 0.0)]);
        mol.add_coord_set(cs);

        mol.generate_bonds(DEFAULT_BOND_TOLERANCE);
        assert_eq!(
            mol.bond_count(),
            1,
            "Fe-N coordination bond should be detected"
        );
    }

    /// Two distant atoms should not be bonded.
    #[test]
    fn test_no_bond_at_large_distance() {
        let mut mol = ObjectMolecule::new("far");
        mol.add_atom(Atom::new("C1", Element::Carbon));
        mol.add_atom(Atom::new("C2", Element::Carbon));

        // C–C at 5.0 Å — well beyond cov threshold (0.76+0.76+0.45 = 1.97)
        let cs = CoordSet::from_vec3(&[Vec3::new(0.0, 0.0, 0.0), Vec3::new(5.0, 0.0, 0.0)]);
        mol.add_coord_set(cs);

        mol.generate_bonds(DEFAULT_BOND_TOLERANCE);
        assert_eq!(mol.bond_count(), 0);
    }

    /// Standard C–C bond at 1.54 Å should be detected.
    #[test]
    fn test_cc_bond() {
        let mut mol = ObjectMolecule::new("ethane");
        mol.add_atom(Atom::new("C1", Element::Carbon));
        mol.add_atom(Atom::new("C2", Element::Carbon));

        let cs = CoordSet::from_vec3(&[Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.54, 0.0, 0.0)]);
        mol.add_coord_set(cs);

        mol.generate_bonds(DEFAULT_BOND_TOLERANCE);
        assert_eq!(mol.bond_count(), 1);
    }

    /// Atoms at bonding distance with different altlocs must NOT be bonded.
    #[test]
    fn test_no_bond_between_different_altlocs() {
        let mut mol = ObjectMolecule::new("altloc");
        let mut a1 = Atom::new("CA", Element::Carbon);
        a1.alt = 'A';
        let mut a2 = Atom::new("CB", Element::Carbon);
        a2.alt = 'B';
        mol.add_atom(a1);
        mol.add_atom(a2);

        let cs = CoordSet::from_vec3(&[Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.54, 0.0, 0.0)]);
        mol.add_coord_set(cs);

        mol.generate_bonds(DEFAULT_BOND_TOLERANCE);
        assert_eq!(mol.bond_count(), 0, "no bond between incompatible altlocs");
    }

    /// Atoms with the same altloc should bond normally.
    #[test]
    fn test_bond_between_same_altlocs() {
        let mut mol = ObjectMolecule::new("altloc");
        let mut a1 = Atom::new("CA", Element::Carbon);
        a1.alt = 'A';
        let mut a2 = Atom::new("CB", Element::Carbon);
        a2.alt = 'A';
        mol.add_atom(a1);
        mol.add_atom(a2);

        let cs = CoordSet::from_vec3(&[Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.54, 0.0, 0.0)]);
        mol.add_coord_set(cs);

        mol.generate_bonds(DEFAULT_BOND_TOLERANCE);
        assert_eq!(mol.bond_count(), 1, "same altloc atoms should bond");
    }

    /// An atom with altloc should bond to a blank-altloc atom (shared position).
    #[test]
    fn test_bond_altloc_to_blank() {
        let mut mol = ObjectMolecule::new("altloc");
        let mut a1 = Atom::new("CA", Element::Carbon);
        a1.alt = 'A';
        let a2 = Atom::new("N", Element::Nitrogen); // alt=' ' (default)
        mol.add_atom(a1);
        mol.add_atom(a2);

        // C-N bond ~1.47 Å
        let cs = CoordSet::from_vec3(&[Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.47, 0.0, 0.0)]);
        mol.add_coord_set(cs);

        mol.generate_bonds(DEFAULT_BOND_TOLERANCE);
        assert_eq!(
            mol.bond_count(),
            1,
            "altloc atom should bond to blank-alt atom"
        );
    }
}
