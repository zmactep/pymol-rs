//! Side-chain completion by self-templating.
//!
//! Missing side-chain atoms are rebuilt by copying geometry from the most
//! complete residue of the same type in the same object: its backbone frame
//! (N/CA/C) is rigidly mapped onto the target residue and the missing atoms are
//! transferred. This is correct by construction for every residue type
//! (including rings) and needs no bundled rotamer library, at the cost of
//! requiring a complete donor of each type to be present.

use std::collections::{HashMap, HashSet};

use lin_alg::f32::Vec3;
use patinae_mol::{Element, ObjectMolecule};

type V = [f64; 3];

fn v(p: Vec3) -> V {
    [f64::from(p.x), f64::from(p.y), f64::from(p.z)]
}

fn sub(a: V, b: V) -> V {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: V, b: V) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: V, b: V) -> V {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn normalize(a: V) -> V {
    let n = dot(a, a).sqrt().max(1e-9);
    [a[0] / n, a[1] / n, a[2] / n]
}

/// Orthonormal backbone frame: e1 along N→CA, e3 ⟂ the N-CA-C plane, e2 = e3×e1.
struct Frame {
    origin: V,
    e1: V,
    e2: V,
    e3: V,
}

impl Frame {
    fn build(n: V, ca: V, c: V) -> Self {
        let e1 = normalize(sub(n, ca));
        let e3 = normalize(cross(e1, sub(c, ca)));
        let e2 = cross(e3, e1);
        Frame {
            origin: ca,
            e1,
            e2,
            e3,
        }
    }

    fn to_local(&self, p: V) -> V {
        let d = sub(p, self.origin);
        [dot(d, self.e1), dot(d, self.e2), dot(d, self.e3)]
    }

    fn to_world(&self, l: V) -> V {
        [
            self.origin[0] + self.e1[0] * l[0] + self.e2[0] * l[1] + self.e3[0] * l[2],
            self.origin[1] + self.e1[1] * l[0] + self.e2[1] * l[1] + self.e3[1] * l[2],
            self.origin[2] + self.e1[2] * l[0] + self.e2[2] * l[1] + self.e3[2] * l[2],
        ]
    }
}

#[derive(Debug, Clone)]
struct ResidueAtoms {
    resn: String,
    chain: String,
    resv: i32,
    inscode: char,
    /// name → (atom index, coord, element)
    atoms: HashMap<String, (usize, Vec3, Element)>,
}

impl ResidueAtoms {
    fn backbone(&self) -> Option<(V, V, V)> {
        let n = self.atoms.get("N")?;
        let ca = self.atoms.get("CA")?;
        let c = self.atoms.get("C")?;
        Some((v(n.1), v(ca.1), v(c.1)))
    }
}

/// One atom to add to a target residue.
#[derive(Debug, Clone)]
pub struct BuiltAtom {
    pub name: String,
    pub element: Element,
    pub coord: Vec3,
}

/// A residue whose side chain is completed from a donor.
#[derive(Debug, Clone)]
pub struct ResidueBuild {
    pub chain: String,
    pub resv: i32,
    pub inscode: char,
    pub resn: String,
    pub added: Vec<BuiltAtom>,
    /// Donor intra-residue connectivity (atom name pairs) to recreate.
    pub bonds: Vec<(String, String)>,
}

/// Side-chain build plan plus a list of residues that could not be completed.
#[derive(Debug, Clone, Default)]
pub struct SideChainPlan {
    pub residues: Vec<ResidueBuild>,
    pub skipped: Vec<String>,
}

fn residue_key(chain: &str, resv: i32, inscode: char) -> String {
    format!("{chain}/{resv}{inscode}")
}

/// Plans the atoms to add so every selected residue matches the most complete
/// instance of its type. `selected` is the set of selected atom indices.
pub fn plan_side_chains(molecule: &ObjectMolecule, selected: &HashSet<usize>) -> SideChainPlan {
    // Index every residue and which residue each atom belongs to.
    let mut residues: HashMap<String, ResidueAtoms> = HashMap::new();
    let mut atom_residue: HashMap<usize, String> = HashMap::new();
    for (idx, atom) in molecule.atoms_indexed() {
        let key = residue_key(&atom.residue.chain, atom.residue.resv, atom.residue.inscode);
        let coord = molecule
            .get_coord(idx, 0)
            .unwrap_or_else(|| Vec3::new(0.0, 0.0, 0.0));
        atom_residue.insert(idx.as_usize(), key.clone());
        residues
            .entry(key)
            .or_insert_with(|| ResidueAtoms {
                resn: atom.residue.resn.clone(),
                chain: atom.residue.chain.clone(),
                resv: atom.residue.resv,
                inscode: atom.residue.inscode,
                atoms: HashMap::new(),
            })
            .atoms
            .insert(atom.name.to_string(), (idx.as_usize(), coord, atom.element));
    }

    // Best donor per residue type: most atoms, must have a backbone frame.
    let mut donors: HashMap<String, String> = HashMap::new();
    for (key, res) in &residues {
        if res.backbone().is_none() {
            continue;
        }
        donors
            .entry(res.resn.clone())
            .and_modify(|best| {
                if residues[best].atoms.len() < res.atoms.len() {
                    *best = key.clone();
                }
            })
            .or_insert_with(|| key.clone());
    }

    // Donor intra-residue bonds (name pairs), computed lazily per donor key.
    let mut donor_bonds: HashMap<String, Vec<(String, String)>> = HashMap::new();
    let mut bonds_for = |donor_key: &str| -> Vec<(String, String)> {
        if let Some(found) = donor_bonds.get(donor_key) {
            return found.clone();
        }
        let mut pairs = Vec::new();
        for bond in molecule.bonds() {
            let a = bond.atom1.as_usize();
            let b = bond.atom2.as_usize();
            if atom_residue.get(&a).map(String::as_str) == Some(donor_key)
                && atom_residue.get(&b).map(String::as_str) == Some(donor_key)
            {
                if let (Some(an), Some(bn)) =
                    (molecule.get_atom(bond.atom1), molecule.get_atom(bond.atom2))
                {
                    pairs.push((an.name.to_string(), bn.name.to_string()));
                }
            }
        }
        donor_bonds.insert(donor_key.to_string(), pairs.clone());
        pairs
    };

    // Residues containing a selected atom are candidates.
    let mut candidate_keys: Vec<String> = residues
        .iter()
        .filter(|(_, res)| res.atoms.values().any(|(idx, _, _)| selected.contains(idx)))
        .map(|(key, _)| key.clone())
        .collect();
    candidate_keys.sort();

    let mut plan = SideChainPlan::default();
    for key in candidate_keys {
        let res = &residues[&key];
        let Some(donor_key) = donors.get(&res.resn) else {
            continue;
        };
        if donor_key == &key {
            continue; // the donor itself is complete
        }
        let donor = &residues[donor_key];
        let missing: Vec<&String> = donor
            .atoms
            .keys()
            .filter(|name| !res.atoms.contains_key(*name))
            .collect();
        if missing.is_empty() {
            continue;
        }
        let (Some(target_bb), Some(donor_bb)) = (res.backbone(), donor.backbone()) else {
            plan.skipped.push(format!("{} {}", res.resn, key));
            continue;
        };
        let donor_frame = Frame::build(donor_bb.0, donor_bb.1, donor_bb.2);
        let target_frame = Frame::build(target_bb.0, target_bb.1, target_bb.2);

        let mut added = Vec::new();
        for name in missing {
            let (_, coord, element) = donor.atoms[name];
            let local = donor_frame.to_local(v(coord));
            let world = target_frame.to_world(local);
            added.push(BuiltAtom {
                name: name.clone(),
                element,
                coord: Vec3::new(world[0] as f32, world[1] as f32, world[2] as f32),
            });
        }
        plan.residues.push(ResidueBuild {
            chain: res.chain.clone(),
            resv: res.resv,
            inscode: res.inscode,
            resn: res.resn.clone(),
            added,
            bonds: bonds_for(donor_key),
        });
    }

    plan
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn frame_round_trips_a_point() {
        let frame = Frame::build([0.0, 0.0, 0.0], [1.46, 0.0, 0.0], [2.0, 1.4, 0.0]);
        let p = [3.0, -1.0, 2.0];
        let local = frame.to_local(p);
        let world = frame.to_world(local);
        for i in 0..3 {
            assert!((world[i] - p[i]).abs() < 1e-6);
        }
    }
}
