//! Bond generation and bond order assignment
//!
//! Extracted from molecule.rs — provides distance-based bond generation
//! using a spatial hash grid and protein-specific bond order assignment.

use lin_alg::f32::Vec3;

use crate::bond::{Bond, BondOrder};
use crate::index::{AtomIndex, BondIndex};
use crate::molecule::ObjectMolecule;
use crate::spatial::SpatialGrid;

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
    let is_bond = |a: &str, b: &str| -> bool {
        (name1 == a && name2 == b) || (name1 == b && name2 == a)
    };

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

/// Generate bonds based on inter-atomic distances for a specific state.
///
/// Uses a spatial hash grid for O(n) performance instead of O(n²) all-pairs.
pub(crate) fn generate_bonds_for_state(mol: &mut ObjectMolecule, tolerance: f32, state: usize) {
    let coord_set = match mol.coord_sets.get(state) {
        Some(cs) => cs,
        None => return,
    };

    let n_atoms = mol.atoms.len();
    if n_atoms < 2 {
        return;
    }

    struct AtomData {
        coord: Vec3,
        vdw: f32,
        is_hydrogen: bool,
    }

    let atom_data: Vec<AtomData> = (0..n_atoms)
        .map(|i| {
            let coord = coord_set
                .get_atom_coord(AtomIndex(i as u32))
                .unwrap_or(Vec3::new(0.0, 0.0, 0.0));
            let atom = &mol.atoms[i];
            AtomData {
                coord,
                vdw: atom.effective_vdw(),
                is_hydrogen: atom.is_hydrogen(),
            }
        })
        .collect();

    let max_vdw = atom_data
        .iter()
        .map(|a| a.vdw)
        .fold(0.0f32, f32::max);
    let max_cutoff = 2.0 * max_vdw * tolerance;

    let mut grid = SpatialGrid::with_capacity(max_cutoff, n_atoms);
    for (i, ad) in atom_data.iter().enumerate() {
        grid.insert(ad.coord, i);
    }

    let mut new_bonds: Vec<(AtomIndex, AtomIndex)> = Vec::new();
    let mut neighbors = Vec::new();

    for i in 0..n_atoms {
        let a1 = &atom_data[i];
        grid.query_neighbors(a1.coord, &mut neighbors);

        for &j in &neighbors {
            if j <= i {
                continue;
            }

            let a2 = &atom_data[j];

            let dx = a1.coord.x - a2.coord.x;
            let dy = a1.coord.y - a2.coord.y;
            let dz = a1.coord.z - a2.coord.z;
            let dist_sq = dx * dx + dy * dy + dz * dz;

            let effective_tolerance = if a1.is_hydrogen && a2.is_hydrogen {
                tolerance * 0.5
            } else {
                tolerance
            };

            let threshold = (a1.vdw + a2.vdw) * effective_tolerance;
            let threshold_sq = threshold * threshold;

            if dist_sq < threshold_sq {
                new_bonds.push((AtomIndex(i as u32), AtomIndex(j as u32)));
            }
        }
    }

    for (a1, a2) in new_bonds {
        let bond = Bond::new(a1, a2, BondOrder::Single);
        let bond_index = BondIndex(mol.bonds.len() as u32);
        mol.bonds.push(bond);
        mol.atom_bonds[a1.as_usize()].push(bond_index);
        mol.atom_bonds[a2.as_usize()].push(bond_index);
        mol.atoms[a1.as_usize()].state.bonded = true;
        mol.atoms[a2.as_usize()].state.bonded = true;
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
