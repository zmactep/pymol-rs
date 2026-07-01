//! Idealized residue side-chain geometry.
//!
//! Builds a canonical side chain (heavy atoms, CB outward) for an amino acid
//! from standard internal coordinates via NeRF placement, on a fixed backbone
//! seed (N/CA/C). Used by mutagenesis to graft a target residue's side chain
//! regardless of whether that type is present in the structure.
//!
//! Coverage: the 15 acyclic residues. Ring residues (HIS/PHE/TYR/TRP) and PRO
//! return `None` for now (a dedicated planar-ring builder is pending).

use lin_alg::f32::Vec3;
use patinae_mol::Element;

type V = [f64; 3];

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

/// NeRF: place D bonded to C, angle B-C-D = `theta_deg`, dihedral A-B-C-D = `phi_deg`.
fn place(a: V, b: V, c: V, r: f64, theta_deg: f64, phi_deg: f64) -> V {
    let theta = theta_deg.to_radians();
    let phi = phi_deg.to_radians();
    let bc = normalize(sub(c, b));
    let n = normalize(cross(sub(b, a), bc));
    let nbc = cross(n, bc);
    let d = [
        -r * theta.cos(),
        r * theta.sin() * phi.cos(),
        r * theta.sin() * phi.sin(),
    ];
    [
        c[0] + bc[0] * d[0] + nbc[0] * d[1] + n[0] * d[2],
        c[1] + bc[1] * d[0] + nbc[1] * d[1] + n[1] * d[2],
        c[2] + bc[2] * d[0] + nbc[2] * d[1] + n[2] * d[2],
    ]
}

/// One internal-coordinate entry: place `name` bonded to `c`, with bond `r`,
/// angle (`b`-`c`-`name`) and dihedral (`a`-`b`-`c`-`name`).
struct Ic {
    name: &'static str,
    element: Element,
    a: &'static str,
    b: &'static str,
    c: &'static str,
    r: f64,
    theta: f64,
    phi: f64,
}

/// A built residue side chain: atoms (name, element, coord) and bonds (names).
#[derive(Debug, Clone)]
pub struct IdealSideChain {
    pub atoms: Vec<(String, Element, Vec3)>,
    pub bonds: Vec<(String, String)>,
}

// Canonical backbone seed (Å). CA at origin, C on +x, N in the xy-plane at the
// standard N-CA-C angle (111°).
const CA: V = [0.0, 0.0, 0.0];
const C: V = [1.525, 0.0, 0.0];
const N: V = [-0.5225, 1.3612, 0.0];

#[allow(clippy::too_many_arguments)]
fn ic(
    name: &'static str,
    element: Element,
    a: &'static str,
    b: &'static str,
    c: &'static str,
    r: f64,
    theta: f64,
    phi: f64,
) -> Ic {
    Ic {
        name,
        element,
        a,
        b,
        c,
        r,
        theta,
        phi,
    }
}

/// CB placed off the backbone with L-amino-acid chirality.
fn cb_ic() -> Ic {
    ic("CB", Element::Carbon, "C", "N", "CA", 1.530, 110.4, 122.6)
}

/// Internal-coordinate table for a residue's side chain (CB first), or `None`
/// for unsupported types (rings / PRO / GLY).
fn side_chain_ics(resn: &str) -> Option<Vec<Ic>> {
    use Element::{Carbon as Cc, Nitrogen as Nn, Oxygen as Oo, Sulfur as Ss};
    let table: Vec<Ic> = match resn {
        "ALA" => vec![cb_ic()],
        "SER" => vec![cb_ic(), ic("OG", Oo, "N", "CA", "CB", 1.417, 110.5, 180.0)],
        "CYS" => vec![cb_ic(), ic("SG", Ss, "N", "CA", "CB", 1.808, 114.0, 180.0)],
        "THR" => vec![
            cb_ic(),
            ic("OG1", Oo, "N", "CA", "CB", 1.433, 109.6, 180.0),
            ic("CG2", Cc, "N", "CA", "CB", 1.521, 111.0, -58.0),
        ],
        "VAL" => vec![
            cb_ic(),
            ic("CG1", Cc, "N", "CA", "CB", 1.521, 110.5, 180.0),
            ic("CG2", Cc, "N", "CA", "CB", 1.521, 110.5, -122.0),
        ],
        "LEU" => vec![
            cb_ic(),
            ic("CG", Cc, "N", "CA", "CB", 1.530, 116.3, 180.0),
            ic("CD1", Cc, "CA", "CB", "CG", 1.521, 110.7, 180.0),
            ic("CD2", Cc, "CA", "CB", "CG", 1.521, 110.7, -120.0),
        ],
        "ILE" => vec![
            cb_ic(),
            ic("CG1", Cc, "N", "CA", "CB", 1.530, 110.4, 180.0),
            ic("CG2", Cc, "N", "CA", "CB", 1.521, 110.5, -122.0),
            ic("CD1", Cc, "CA", "CB", "CG1", 1.513, 113.8, 180.0),
        ],
        "MET" => vec![
            cb_ic(),
            ic("CG", Cc, "N", "CA", "CB", 1.520, 114.1, 180.0),
            ic("SD", Ss, "CA", "CB", "CG", 1.803, 112.7, 180.0),
            ic("CE", Cc, "CB", "CG", "SD", 1.791, 100.9, 180.0),
        ],
        "ASP" => vec![
            cb_ic(),
            ic("CG", Cc, "N", "CA", "CB", 1.516, 112.6, 180.0),
            ic("OD1", Oo, "CA", "CB", "CG", 1.249, 118.4, 0.0),
            ic("OD2", Oo, "CA", "CB", "CG", 1.249, 118.4, 180.0),
        ],
        "ASN" => vec![
            cb_ic(),
            ic("CG", Cc, "N", "CA", "CB", 1.516, 112.6, 180.0),
            ic("OD1", Oo, "CA", "CB", "CG", 1.231, 120.8, 0.0),
            ic("ND2", Nn, "CA", "CB", "CG", 1.328, 116.5, 180.0),
        ],
        "GLU" => vec![
            cb_ic(),
            ic("CG", Cc, "N", "CA", "CB", 1.520, 113.8, 180.0),
            ic("CD", Cc, "CA", "CB", "CG", 1.516, 112.6, 180.0),
            ic("OE1", Oo, "CB", "CG", "CD", 1.249, 118.4, 0.0),
            ic("OE2", Oo, "CB", "CG", "CD", 1.249, 118.4, 180.0),
        ],
        "GLN" => vec![
            cb_ic(),
            ic("CG", Cc, "N", "CA", "CB", 1.520, 113.8, 180.0),
            ic("CD", Cc, "CA", "CB", "CG", 1.516, 112.6, 180.0),
            ic("OE1", Oo, "CB", "CG", "CD", 1.231, 120.8, 0.0),
            ic("NE2", Nn, "CB", "CG", "CD", 1.328, 116.5, 180.0),
        ],
        "LYS" => vec![
            cb_ic(),
            ic("CG", Cc, "N", "CA", "CB", 1.520, 113.8, 180.0),
            ic("CD", Cc, "CA", "CB", "CG", 1.520, 111.3, 180.0),
            ic("CE", Cc, "CB", "CG", "CD", 1.520, 111.3, 180.0),
            ic("NZ", Nn, "CG", "CD", "CE", 1.489, 111.9, 180.0),
        ],
        "ARG" => vec![
            cb_ic(),
            ic("CG", Cc, "N", "CA", "CB", 1.520, 113.8, 180.0),
            ic("CD", Cc, "CA", "CB", "CG", 1.520, 111.3, 180.0),
            ic("NE", Nn, "CB", "CG", "CD", 1.460, 112.0, 180.0),
            ic("CZ", Cc, "CG", "CD", "NE", 1.329, 124.2, 180.0),
            ic("NH1", Nn, "CD", "NE", "CZ", 1.326, 120.0, 0.0),
            ic("NH2", Nn, "CD", "NE", "CZ", 1.326, 120.0, 180.0),
        ],
        // Phenyl ring (closure CE2–CZ).
        "PHE" => vec![
            cb_ic(),
            ic("CG", Cc, "N", "CA", "CB", 1.502, 113.8, 180.0),
            ic("CD1", Cc, "CA", "CB", "CG", 1.384, 120.0, 90.0),
            ic("CD2", Cc, "CA", "CB", "CG", 1.384, 120.0, 270.0),
            ic("CE1", Cc, "CB", "CG", "CD1", 1.382, 120.4, 180.0),
            ic("CE2", Cc, "CB", "CG", "CD2", 1.382, 120.4, 180.0),
            ic("CZ", Cc, "CG", "CD1", "CE1", 1.382, 120.0, 0.0),
        ],
        // Tyrosine: phenol (closure CE2–CZ) + hydroxyl.
        "TYR" => vec![
            cb_ic(),
            ic("CG", Cc, "N", "CA", "CB", 1.512, 113.8, 180.0),
            ic("CD1", Cc, "CA", "CB", "CG", 1.387, 120.0, 90.0),
            ic("CD2", Cc, "CA", "CB", "CG", 1.387, 120.0, 270.0),
            ic("CE1", Cc, "CB", "CG", "CD1", 1.382, 120.4, 180.0),
            ic("CE2", Cc, "CB", "CG", "CD2", 1.382, 120.4, 180.0),
            ic("CZ", Cc, "CG", "CD1", "CE1", 1.378, 120.0, 0.0),
            ic("OH", Oo, "CD1", "CE1", "CZ", 1.375, 119.9, 180.0),
        ],
        // Imidazole, built around the ring (closure CD2–CG).
        "HIS" => vec![
            cb_ic(),
            ic("CG", Cc, "N", "CA", "CB", 1.497, 113.8, 180.0),
            ic("ND1", Nn, "CA", "CB", "CG", 1.378, 121.6, 90.0),
            ic("CE1", Cc, "CB", "CG", "ND1", 1.321, 105.6, 180.0),
            ic("NE2", Nn, "CG", "ND1", "CE1", 1.319, 111.7, 0.0),
            ic("CD2", Cc, "ND1", "CE1", "NE2", 1.374, 107.2, 0.0),
        ],
        // Indole: 5-ring chain (closure CD2–CG), then the fused 6-ring chain
        // off CD2 (closure CZ2–CE2). All dihedrals planar.
        "TRP" => vec![
            cb_ic(),
            ic("CG", Cc, "N", "CA", "CB", 1.498, 113.6, 180.0),
            ic("CD1", Cc, "CA", "CB", "CG", 1.365, 126.8, 90.0),
            ic("NE1", Nn, "CB", "CG", "CD1", 1.374, 110.1, 180.0),
            ic("CE2", Cc, "CG", "CD1", "NE1", 1.370, 108.9, 0.0),
            ic("CD2", Cc, "CD1", "NE1", "CE2", 1.409, 107.4, 0.0),
            ic("CE3", Cc, "NE1", "CE2", "CD2", 1.400, 122.4, 180.0),
            ic("CZ3", Cc, "CE2", "CD2", "CE3", 1.385, 118.7, 0.0),
            ic("CH2", Cc, "CD2", "CE3", "CZ3", 1.400, 121.1, 0.0),
            ic("CZ2", Cc, "CE3", "CZ3", "CH2", 1.368, 121.5, 0.0),
        ],
        // Pyrrolidine ring back to backbone N (closure CD–N).
        "PRO" => vec![
            cb_ic(),
            ic("CG", Cc, "N", "CA", "CB", 1.492, 104.5, 30.0),
            ic("CD", Cc, "CA", "CB", "CG", 1.503, 106.1, -35.0),
        ],
        _ => return None,
    };
    Some(table)
}

/// Extra bonds that close rings (beyond the internal-coordinate tree). For PRO
/// the closure attaches to the backbone N of the target residue.
fn ring_closures(resn: &str) -> &'static [(&'static str, &'static str)] {
    match resn {
        "PHE" | "TYR" => &[("CE2", "CZ")],
        "HIS" => &[("CD2", "CG")],
        "TRP" => &[("CD2", "CG"), ("CZ2", "CE2")],
        "PRO" => &[("CD", "N")],
        _ => &[],
    }
}

/// Builds the idealized side chain for `resn`, or `None` if unsupported.
/// GLY has no side chain and returns an empty build.
pub fn ideal_side_chain(resn: &str) -> Option<IdealSideChain> {
    if resn == "GLY" {
        return Some(IdealSideChain {
            atoms: Vec::new(),
            bonds: Vec::new(),
        });
    }
    let ics = side_chain_ics(resn)?;
    let mut coords: Vec<(&'static str, V)> = vec![("N", N), ("CA", CA), ("C", C)];
    let lookup = |coords: &[(&'static str, V)], name: &str| -> Option<V> {
        coords.iter().find(|(n, _)| *n == name).map(|(_, p)| *p)
    };

    let mut atoms = Vec::new();
    let mut bonds = Vec::new();
    for entry in &ics {
        let (a, b, c) = (
            lookup(&coords, entry.a)?,
            lookup(&coords, entry.b)?,
            lookup(&coords, entry.c)?,
        );
        let p = place(a, b, c, entry.r, entry.theta, entry.phi);
        coords.push((entry.name, p));
        atoms.push((
            entry.name.to_string(),
            entry.element,
            Vec3::new(p[0] as f32, p[1] as f32, p[2] as f32),
        ));
        bonds.push((entry.c.to_string(), entry.name.to_string()));
    }
    for (a, b) in ring_closures(resn) {
        bonds.push(((*a).to_string(), (*b).to_string()));
    }
    Some(IdealSideChain { atoms, bonds })
}

/// Returns whether mutation to `resn` is supported by the ideal-geometry table.
pub fn is_supported(resn: &str) -> bool {
    resn == "GLY" || side_chain_ics(resn).is_some()
}

fn vd(p: Vec3) -> V {
    [f64::from(p.x), f64::from(p.y), f64::from(p.z)]
}

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
        Frame { origin: ca, e1, e2, e3 }
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

/// Like [`ideal_side_chain`], but with coordinates mapped onto a target
/// residue's backbone frame (actual N/CA/C positions), ready to graft.
pub fn ideal_side_chain_on(resn: &str, n: Vec3, ca: Vec3, c: Vec3) -> Option<IdealSideChain> {
    let base = ideal_side_chain(resn)?;
    let canonical = Frame::build(N, CA, C);
    let target = Frame::build(vd(n), vd(ca), vd(c));
    let atoms = base
        .atoms
        .into_iter()
        .map(|(name, element, coord)| {
            let world = target.to_world(canonical.to_local(vd(coord)));
            (
                name,
                element,
                Vec3::new(world[0] as f32, world[1] as f32, world[2] as f32),
            )
        })
        .collect();
    Some(IdealSideChain {
        atoms,
        bonds: base.bonds,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn d(a: Vec3, b: Vec3) -> f64 {
        f64::from((a.x - b.x).powi(2) + (a.y - b.y).powi(2) + (a.z - b.z).powi(2)).sqrt()
    }

    fn coord(sc: &IdealSideChain, name: &str) -> Vec3 {
        // Side-chain bonds may reference the backbone seed atom CA (CB–CA).
        match name {
            "CA" => Vec3::new(CA[0] as f32, CA[1] as f32, CA[2] as f32),
            "N" => Vec3::new(N[0] as f32, N[1] as f32, N[2] as f32),
            "C" => Vec3::new(C[0] as f32, C[1] as f32, C[2] as f32),
            _ => sc.atoms.iter().find(|(n, _, _)| n == name).unwrap().2,
        }
    }

    const ALL: &[&str] = &[
        "ALA", "SER", "CYS", "THR", "VAL", "LEU", "ILE", "MET", "ASP", "ASN", "GLU", "GLN", "LYS",
        "ARG", "PHE", "TYR", "HIS", "TRP", "PRO",
    ];

    #[test]
    fn all_residues_build() {
        for resn in ALL {
            let sc = ideal_side_chain(resn).unwrap_or_else(|| panic!("{resn} failed"));
            assert!(!sc.atoms.is_empty(), "{resn} empty");
            assert!(is_supported(resn));
        }
        assert!(ideal_side_chain("GLY").unwrap().atoms.is_empty());
    }

    #[test]
    fn bond_lengths_are_physical() {
        // Includes ring-closure bonds.
        for resn in ALL {
            let sc = ideal_side_chain(resn).unwrap();
            for (a, b) in &sc.bonds {
                let len = d(coord(&sc, a), coord(&sc, b));
                assert!((1.0..2.0).contains(&len), "{resn} bond {a}-{b} length {len}");
            }
        }
    }

    #[test]
    fn no_internal_clashes() {
        for resn in ALL {
            let sc = ideal_side_chain(resn).unwrap();
            for i in 0..sc.atoms.len() {
                for j in (i + 1)..sc.atoms.len() {
                    let dist = d(sc.atoms[i].2, sc.atoms[j].2);
                    assert!(dist > 1.0, "{resn} clash {} {} = {dist}", sc.atoms[i].0, sc.atoms[j].0);
                }
            }
        }
    }

    #[test]
    fn phenyl_ring_is_planar() {
        let sc = ideal_side_chain("PHE").unwrap();
        let ring = ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"];
        let p: Vec<V> = ring
            .iter()
            .map(|n| {
                let c = coord(&sc, n);
                [f64::from(c.x), f64::from(c.y), f64::from(c.z)]
            })
            .collect();
        // Plane normal from the first three; every atom must lie near it.
        let normal = normalize(cross(sub(p[1], p[0]), sub(p[2], p[0])));
        for q in &p {
            let dev = dot(normal, sub(*q, p[0])).abs();
            assert!(dev < 0.1, "PHE ring out of plane by {dev} Å");
        }
    }

    #[test]
    fn cb_has_l_chirality() {
        // Improper volume sign of (N-CA, C-CA, CB-CA) is negative for L-amino acids.
        let sc = ideal_side_chain("ALA").unwrap();
        let cb = coord(&sc, "CB");
        let nca = sub(N, CA);
        let cca = sub(C, CA);
        let cbca = [
            f64::from(cb.x) - CA[0],
            f64::from(cb.y) - CA[1],
            f64::from(cb.z) - CA[2],
        ];
        let vol = dot(cross(nca, cca), cbca);
        assert!(vol < 0.0, "CB chirality volume {vol} (expected L, < 0)");
    }
}
