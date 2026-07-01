//! Terminal capping geometry.
//!
//! Builds neutral caps for polypeptide termini from idealized internal
//! coordinates (NeRF placement) relative to the terminal residue's backbone:
//! an acetyl group (`ACE`) on the N-terminus and an N-methyl amide (`NME`) on
//! the C-terminus. Pure geometry — no force field required.

use lin_alg::f32::Vec3;
use patinae_mol::Element;

type V = [f64; 3];

fn v(p: Vec3) -> V {
    [f64::from(p.x), f64::from(p.y), f64::from(p.z)]
}

fn to_vec3(p: V) -> Vec3 {
    Vec3::new(p[0] as f32, p[1] as f32, p[2] as f32)
}

fn sub(a: V, b: V) -> V {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn cross(a: V, b: V) -> V {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn norm(a: V) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

fn normalize(a: V) -> V {
    let n = norm(a).max(1e-9);
    [a[0] / n, a[1] / n, a[2] / n]
}

/// Places atom D bonded to C, with bond length `r` (C–D), bond angle `theta_deg`
/// (B–C–D) and dihedral `phi_deg` (A–B–C–D). The classic NeRF construction.
fn place(a: V, b: V, c: V, r: f64, theta_deg: f64, phi_deg: f64) -> V {
    let theta = theta_deg.to_radians();
    let phi = phi_deg.to_radians();
    let bc = normalize(sub(c, b));
    let n = normalize(cross(sub(b, a), bc));
    let nbc = cross(n, bc);
    let d2 = [
        -r * theta.cos(),
        r * theta.sin() * phi.cos(),
        r * theta.sin() * phi.sin(),
    ];
    [
        c[0] + bc[0] * d2[0] + nbc[0] * d2[1] + n[0] * d2[2],
        c[1] + bc[1] * d2[0] + nbc[1] * d2[1] + n[1] * d2[2],
        c[2] + bc[2] * d2[0] + nbc[2] * d2[1] + n[2] * d2[2],
    ]
}

/// One atom of a cap residue.
#[derive(Debug, Clone)]
pub struct CapAtom {
    pub name: &'static str,
    pub element: Element,
    pub coord: Vec3,
}

/// A cap residue to graft onto a chain terminus.
#[derive(Debug, Clone)]
pub struct CapBuild {
    pub resn: &'static str,
    pub chain: String,
    pub resv: i32,
    pub atoms: Vec<CapAtom>,
    /// Bonds within the cap, as indices into `atoms`.
    pub internal_bonds: Vec<(usize, usize)>,
    /// Indices into `atoms` of the cap atoms bonded to the terminal atom.
    pub parent_bonds: Vec<usize>,
    /// Atom index (in the source molecule) the cap bonds to.
    pub parent_atom_index: usize,
    /// When true, the new atoms join the terminal residue (e.g. NH3⁺ hydrogens,
    /// a carboxylate OXT) instead of forming their own cap residue.
    pub extend_parent_residue: bool,
}

/// Builds an acetyl (`ACE`) cap on the N-terminus, from the first residue's
/// backbone N/CA/C positions. `n_atom_index` is the N atom in the molecule.
pub fn build_ace(n: Vec3, ca: Vec3, c: Vec3, chain: String, resv: i32, n_atom_index: usize) -> CapBuild {
    let (n, ca, c) = (v(n), v(ca), v(c));
    let c_ace = place(c, ca, n, 1.335, 121.0, 180.0);
    let o = place(ca, n, c_ace, 1.229, 122.0, 0.0);
    let ct = place(ca, n, c_ace, 1.500, 116.0, 180.0);
    let h1 = place(n, c_ace, ct, 1.090, 109.5, 60.0);
    let h2 = place(n, c_ace, ct, 1.090, 109.5, 180.0);
    let h3 = place(n, c_ace, ct, 1.090, 109.5, 300.0);
    CapBuild {
        resn: "ACE",
        chain,
        resv,
        atoms: vec![
            CapAtom { name: "C", element: Element::Carbon, coord: to_vec3(c_ace) },
            CapAtom { name: "O", element: Element::Oxygen, coord: to_vec3(o) },
            CapAtom { name: "CH3", element: Element::Carbon, coord: to_vec3(ct) },
            CapAtom { name: "HH31", element: Element::Hydrogen, coord: to_vec3(h1) },
            CapAtom { name: "HH32", element: Element::Hydrogen, coord: to_vec3(h2) },
            CapAtom { name: "HH33", element: Element::Hydrogen, coord: to_vec3(h3) },
        ],
        internal_bonds: vec![(0, 1), (0, 2), (2, 3), (2, 4), (2, 5)],
        parent_bonds: vec![0],
        parent_atom_index: n_atom_index,
        extend_parent_residue: false,
    }
}

/// Builds a charged free N-terminus (NH3⁺): three hydrogens on the terminal N,
/// staggered relative to the backbone carbonyl. The hydrogens join the terminal
/// residue rather than forming a new cap.
pub fn build_nterm_nh3(
    n: Vec3,
    ca: Vec3,
    c: Vec3,
    chain: String,
    resv: i32,
    n_atom_index: usize,
) -> CapBuild {
    let (n, ca, c) = (v(n), v(ca), v(c));
    let h1 = place(c, ca, n, 1.010, 109.5, 60.0);
    let h2 = place(c, ca, n, 1.010, 109.5, 180.0);
    let h3 = place(c, ca, n, 1.010, 109.5, 300.0);
    CapBuild {
        resn: "",
        chain,
        resv,
        atoms: vec![
            CapAtom { name: "H1", element: Element::Hydrogen, coord: to_vec3(h1) },
            CapAtom { name: "H2", element: Element::Hydrogen, coord: to_vec3(h2) },
            CapAtom { name: "H3", element: Element::Hydrogen, coord: to_vec3(h3) },
        ],
        internal_bonds: vec![],
        parent_bonds: vec![0, 1, 2],
        parent_atom_index: n_atom_index,
        extend_parent_residue: true,
    }
}

/// Builds an N-methyl amide (`NME`) cap on the C-terminus, from the last
/// residue's backbone C/CA/O positions. `c_atom_index` is the C atom.
pub fn build_nme(c: Vec3, ca: Vec3, o: Vec3, chain: String, resv: i32, c_atom_index: usize) -> CapBuild {
    let (c, ca, o) = (v(c), v(ca), v(o));
    let n = place(o, ca, c, 1.335, 116.0, 180.0);
    let h = place(ca, c, n, 1.010, 119.0, 0.0);
    let ct = place(ca, c, n, 1.448, 121.0, 180.0);
    let h1 = place(c, n, ct, 1.090, 109.5, 60.0);
    let h2 = place(c, n, ct, 1.090, 109.5, 180.0);
    let h3 = place(c, n, ct, 1.090, 109.5, 300.0);
    CapBuild {
        resn: "NME",
        chain,
        resv,
        atoms: vec![
            CapAtom { name: "N", element: Element::Nitrogen, coord: to_vec3(n) },
            CapAtom { name: "H", element: Element::Hydrogen, coord: to_vec3(h) },
            CapAtom { name: "CH3", element: Element::Carbon, coord: to_vec3(ct) },
            CapAtom { name: "HH31", element: Element::Hydrogen, coord: to_vec3(h1) },
            CapAtom { name: "HH32", element: Element::Hydrogen, coord: to_vec3(h2) },
            CapAtom { name: "HH33", element: Element::Hydrogen, coord: to_vec3(h3) },
        ],
        internal_bonds: vec![(0, 1), (0, 2), (2, 3), (2, 4), (2, 5)],
        parent_bonds: vec![0],
        parent_atom_index: c_atom_index,
        extend_parent_residue: false,
    }
}

/// Builds an amidated C-terminus (`NH2`): a carboxamide nitrogen with two
/// hydrogens grafted onto the terminal carbonyl carbon, as its own residue.
pub fn build_cterm_amide(
    c: Vec3,
    ca: Vec3,
    o: Vec3,
    chain: String,
    resv: i32,
    c_atom_index: usize,
) -> CapBuild {
    let (c, ca, o) = (v(c), v(ca), v(o));
    let n = place(o, ca, c, 1.335, 116.0, 180.0);
    let h1 = place(ca, c, n, 1.010, 119.0, 0.0);
    let h2 = place(ca, c, n, 1.010, 121.0, 180.0);
    CapBuild {
        resn: "NH2",
        chain,
        resv,
        atoms: vec![
            CapAtom { name: "N", element: Element::Nitrogen, coord: to_vec3(n) },
            CapAtom { name: "H1", element: Element::Hydrogen, coord: to_vec3(h1) },
            CapAtom { name: "H2", element: Element::Hydrogen, coord: to_vec3(h2) },
        ],
        internal_bonds: vec![(0, 1), (0, 2)],
        parent_bonds: vec![0],
        parent_atom_index: c_atom_index,
        extend_parent_residue: false,
    }
}

/// Builds a charged free C-terminus (COO⁻): a second carboxyl oxygen (`OXT`)
/// on the terminal carbon, joining the terminal residue.
pub fn build_cterm_carboxylate(
    c: Vec3,
    ca: Vec3,
    o: Vec3,
    chain: String,
    resv: i32,
    c_atom_index: usize,
) -> CapBuild {
    let (c, ca, o) = (v(c), v(ca), v(o));
    let oxt = place(o, ca, c, 1.251, 118.0, 180.0);
    CapBuild {
        resn: "",
        chain,
        resv,
        atoms: vec![CapAtom { name: "OXT", element: Element::Oxygen, coord: to_vec3(oxt) }],
        internal_bonds: vec![],
        parent_bonds: vec![0],
        parent_atom_index: c_atom_index,
        extend_parent_residue: true,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dist(a: Vec3, b: Vec3) -> f64 {
        f64::from((a.x - b.x).powi(2) + (a.y - b.y).powi(2) + (a.z - b.z).powi(2)).sqrt()
    }

    #[test]
    fn ace_bond_lengths_are_physical() {
        let n = Vec3::new(0.0, 0.0, 0.0);
        let ca = Vec3::new(1.46, 0.0, 0.0);
        let c = Vec3::new(2.0, 1.4, 0.0);
        let cap = build_ace(n, ca, c, "A".into(), 0, 0);
        let carbonyl = cap.atoms[0].coord;
        // C(ace)–N bond ~1.335 Å
        assert!((dist(carbonyl, n) - 1.335).abs() < 0.05, "C-N = {}", dist(carbonyl, n));
        // C=O bond ~1.229 Å
        assert!((dist(carbonyl, cap.atoms[1].coord) - 1.229).abs() < 0.05);
        // C-CH3 bond ~1.5 Å
        assert!((dist(carbonyl, cap.atoms[2].coord) - 1.5).abs() < 0.05);
        assert_eq!(cap.resn, "ACE");
        assert_eq!(cap.parent_atom_index, 0);
    }

    #[test]
    fn nme_bond_lengths_are_physical() {
        let c = Vec3::new(0.0, 0.0, 0.0);
        let ca = Vec3::new(-1.52, 0.0, 0.0);
        let o = Vec3::new(0.5, 1.2, 0.0);
        let cap = build_nme(c, ca, o, "A".into(), 100, 5);
        let amide_n = cap.atoms[0].coord;
        assert!((dist(amide_n, c) - 1.335).abs() < 0.05, "N-C = {}", dist(amide_n, c));
        assert!((dist(amide_n, cap.atoms[1].coord) - 1.01).abs() < 0.05);
        assert_eq!(cap.resn, "NME");
    }

    #[test]
    fn nh3_terminus_adds_three_hydrogens_to_n() {
        let n = Vec3::new(0.0, 0.0, 0.0);
        let ca = Vec3::new(1.46, 0.0, 0.0);
        let c = Vec3::new(2.0, 1.4, 0.0);
        let cap = build_nterm_nh3(n, ca, c, "A".into(), 5, 0);
        assert!(cap.extend_parent_residue, "NH3⁺ joins the terminal residue");
        assert_eq!(cap.atoms.len(), 3);
        assert_eq!(cap.parent_bonds, vec![0, 1, 2], "all three H bond to N");
        for h in &cap.atoms {
            assert_eq!(h.element, Element::Hydrogen);
            assert!((dist(h.coord, n) - 1.01).abs() < 0.05, "N-H length off");
        }
    }

    #[test]
    fn carboxylate_terminus_adds_oxt() {
        let c = Vec3::new(0.0, 0.0, 0.0);
        let ca = Vec3::new(-1.52, 0.0, 0.0);
        let o = Vec3::new(0.5, 1.2, 0.0);
        let cap = build_cterm_carboxylate(c, ca, o, "A".into(), 100, 5);
        assert!(cap.extend_parent_residue);
        assert_eq!(cap.atoms.len(), 1);
        assert_eq!(cap.atoms[0].name, "OXT");
        assert!((dist(cap.atoms[0].coord, c) - 1.251).abs() < 0.05, "C-OXT length off");
    }

    #[test]
    fn amide_terminus_is_its_own_residue() {
        let c = Vec3::new(0.0, 0.0, 0.0);
        let ca = Vec3::new(-1.52, 0.0, 0.0);
        let o = Vec3::new(0.5, 1.2, 0.0);
        let cap = build_cterm_amide(c, ca, o, "A".into(), 101, 5);
        assert!(!cap.extend_parent_residue);
        assert_eq!(cap.resn, "NH2");
        assert_eq!(cap.atoms.len(), 3, "amide N + 2 H");
        assert!((dist(cap.atoms[0].coord, c) - 1.335).abs() < 0.05);
    }
}
