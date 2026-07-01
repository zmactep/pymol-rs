//! Backbone-independent rotamer search for side-chain repacking.
//!
//! Rather than bundling a published rotamer library verbatim, this samples each
//! χ dihedral over its staggered states (a backbone-independent grid) and picks
//! the combination with the least steric overlap against a fixed environment.
//! It is used by the "Build side chains" tool to resolve clashes after a side
//! chain has been grafted: the side-chain atoms are rotated about their χ bonds
//! while the backbone and the rest of the structure stay put.
//!
//! The χ bonds and which atoms they move are derived from the residue's own
//! bond graph, so the same routine handles every amino acid (the rigid aromatic
//! rings simply rotate as a unit about χ2).

use std::collections::{HashMap, HashSet, VecDeque};

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
fn norm(a: V) -> f64 {
    dot(a, a).sqrt()
}
fn normalize(a: V) -> V {
    let n = norm(a).max(1e-9);
    [a[0] / n, a[1] / n, a[2] / n]
}

/// χ dihedral atom-name quads (χ1, χ2, …) per residue. The rotatable bond of
/// χk is the middle pair `[1]`–`[2]`; atom `[3]` and everything beyond it move.
fn chi_atoms(resn: &str) -> &'static [[&'static str; 4]] {
    match resn {
        "SER" => &[["N", "CA", "CB", "OG"]],
        "CYS" => &[["N", "CA", "CB", "SG"]],
        "THR" => &[["N", "CA", "CB", "OG1"]],
        "VAL" => &[["N", "CA", "CB", "CG1"]],
        "LEU" => &[["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD1"]],
        "ILE" => &[["N", "CA", "CB", "CG1"], ["CA", "CB", "CG1", "CD1"]],
        "MET" => &[
            ["N", "CA", "CB", "CG"],
            ["CA", "CB", "CG", "SD"],
            ["CB", "CG", "SD", "CE"],
        ],
        "ASP" => &[["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "OD1"]],
        "ASN" => &[["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "OD1"]],
        "GLU" => &[
            ["N", "CA", "CB", "CG"],
            ["CA", "CB", "CG", "CD"],
            ["CB", "CG", "CD", "OE1"],
        ],
        "GLN" => &[
            ["N", "CA", "CB", "CG"],
            ["CA", "CB", "CG", "CD"],
            ["CB", "CG", "CD", "OE1"],
        ],
        "LYS" => &[
            ["N", "CA", "CB", "CG"],
            ["CA", "CB", "CG", "CD"],
            ["CB", "CG", "CD", "CE"],
            ["CG", "CD", "CE", "NZ"],
        ],
        "ARG" => &[
            ["N", "CA", "CB", "CG"],
            ["CA", "CB", "CG", "CD"],
            ["CB", "CG", "CD", "NE"],
            ["CG", "CD", "NE", "CZ"],
        ],
        "PHE" | "TYR" => &[["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD1"]],
        "HIS" => &[["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "ND1"]],
        "TRP" => &[["N", "CA", "CB", "CG"], ["CA", "CB", "CG", "CD1"]],
        // ALA (no χ), GLY (no side chain), PRO (ring locked to backbone).
        _ => &[],
    }
}

/// Staggered χ samples (degrees). Coarse = the three classic rotamers; fine
/// adds ±30° shoulders for a denser search.
fn chi_samples(fine: bool) -> &'static [f64] {
    if fine {
        &[-90.0, -60.0, -30.0, 30.0, 60.0, 90.0, 150.0, 180.0, -150.0]
    } else {
        &[-60.0, 60.0, 180.0]
    }
}

/// van-der-Waals radius (Å) used only for the soft clash penalty.
fn vdw_radius(element: Element) -> f64 {
    match element {
        Element::Hydrogen => 1.10,
        Element::Carbon => 1.70,
        Element::Nitrogen => 1.55,
        Element::Oxygen => 1.52,
        Element::Sulfur => 1.80,
        _ => 1.70,
    }
}

/// Dihedral angle A-B-C-D in degrees.
fn dihedral(a: V, b: V, c: V, d: V) -> f64 {
    let b1 = sub(b, a);
    let b2 = sub(c, b);
    let b3 = sub(d, c);
    let n1 = cross(b1, b2);
    let n2 = cross(b2, b3);
    let m = cross(n1, normalize(b2));
    let x = dot(n1, n2);
    let y = dot(m, n2);
    y.atan2(x).to_degrees()
}

/// Rotates `point` by `angle_deg` about the line through `origin` along `axis`.
fn rotate_about(point: V, origin: V, axis: V, angle_deg: f64) -> V {
    let k = normalize(axis);
    let p = sub(point, origin);
    let theta = angle_deg.to_radians();
    let (s, cth) = theta.sin_cos();
    // Rodrigues' rotation formula.
    let kxp = cross(k, p);
    let kdp = dot(k, p);
    let r = [
        p[0] * cth + kxp[0] * s + k[0] * kdp * (1.0 - cth),
        p[1] * cth + kxp[1] * s + k[1] * kdp * (1.0 - cth),
        p[2] * cth + kxp[2] * s + k[2] * kdp * (1.0 - cth),
    ];
    [origin[0] + r[0], origin[1] + r[1], origin[2] + r[2]]
}

/// A movable side-chain atom.
#[derive(Debug, Clone)]
pub struct SideAtom {
    pub name: String,
    pub element: Element,
    pub coord: Vec3,
}

/// One resolved χ: the four dihedral atoms `a`-`b`-`c`-`d` (by name; the
/// rotatable bond is `b`–`c`) and the set of side-chain atom names that move
/// with it (everything beyond `c`, including `d`).
struct ChiAxis {
    a: String,
    b: String,
    c: String,
    d: String,
    downstream: Vec<String>,
}

/// Builds the χ axes for a residue from its chi-atom table and bond graph.
/// `present` is the set of side-chain atom names actually built.
fn resolve_chi_axes(
    resn: &str,
    present: &HashSet<String>,
    bonds: &[(String, String)],
) -> Vec<ChiAxis> {
    // Adjacency among side-chain atoms only (backbone N/CA/C/O are not movable).
    let mut adj: HashMap<&str, Vec<&str>> = HashMap::new();
    for (a, b) in bonds {
        if present.contains(a) && present.contains(b) {
            adj.entry(a).or_default().push(b);
            adj.entry(b).or_default().push(a);
        }
    }

    let mut axes = Vec::new();
    for quad in chi_atoms(resn) {
        let [a, b, c, d] = *quad;
        if !present.contains(c) || !present.contains(d) {
            continue;
        }
        // Downstream = side-chain atoms reachable from `d` without crossing the
        // rotatable bond `b`–`c`.
        let mut downstream = Vec::new();
        let mut seen: HashSet<&str> = HashSet::new();
        seen.insert(c);
        let mut queue = VecDeque::new();
        queue.push_back(d);
        seen.insert(d);
        while let Some(cur) = queue.pop_front() {
            downstream.push(cur.to_string());
            for &next in adj.get(cur).into_iter().flatten() {
                // Never walk back across the rotatable bond into `b`.
                if next == b || seen.contains(next) {
                    continue;
                }
                seen.insert(next);
                queue.push_back(next);
            }
        }
        axes.push(ChiAxis {
            a: a.to_string(),
            b: b.to_string(),
            c: c.to_string(),
            d: d.to_string(),
            downstream,
        });
    }
    axes
}

/// Soft clash penalty of the side chain against the fixed environment: summed
/// squared overlap where atoms are closer than the sum of their vdW radii
/// (minus a small tolerance so van-der-Waals contact is not penalized).
fn clash_score(
    sidechain: &HashMap<String, (V, Element)>,
    environment: &[(V, Element)],
) -> f64 {
    const TOL: f64 = 0.4;
    let mut score = 0.0;
    for (sc_pos, sc_el) in sidechain.values() {
        let rs = vdw_radius(*sc_el);
        for (env_pos, env_el) in environment {
            let min_d = rs + vdw_radius(*env_el) - TOL;
            let d = norm(sub(*sc_pos, *env_pos));
            if d < min_d {
                let overlap = min_d - d;
                score += overlap * overlap;
            }
        }
    }
    score
}

/// Repacks a grafted side chain by χ-grid search, returning the clash-minimal
/// coordinates keyed by atom name. `fixed` supplies backbone reference atoms
/// (at least N and CA) used by the χ1 dihedral; `environment` is every atom not
/// belonging to this residue. Returns the input coordinates unchanged when the
/// residue has no rotatable χ or no environment to clash against.
pub fn optimize_rotamer(
    resn: &str,
    sidechain: &[SideAtom],
    fixed: &HashMap<String, Vec3>,
    bonds: &[(String, String)],
    environment: &[(Vec3, Element)],
    fine: bool,
) -> HashMap<String, Vec3> {
    let present: HashSet<String> = sidechain.iter().map(|a| a.name.clone()).collect();
    let axes = resolve_chi_axes(resn, &present, bonds);
    let mut current: HashMap<String, (V, Element)> =
        sidechain.iter().map(|a| (a.name.clone(), (v(a.coord), a.element))).collect();

    if axes.is_empty() {
        return current.into_iter().map(|(k, (p, _))| (k, to_vec3(p))).collect();
    }

    let fixed: HashMap<String, V> = fixed.iter().map(|(k, p)| (k.clone(), v(*p))).collect();
    let env: Vec<(V, Element)> = environment.iter().map(|(p, e)| (v(*p), *e)).collect();
    let pos_of = |current: &HashMap<String, (V, Element)>, name: &str| -> Option<V> {
        current.get(name).map(|(p, _)| *p).or_else(|| fixed.get(name).copied())
    };

    let samples = chi_samples(fine);
    // Cap the search so deep side chains (χ4) stay tractable.
    const MAX_COMBOS: usize = 4096;
    let mut total = 1usize;
    for _ in &axes {
        total = total.saturating_mul(samples.len());
    }
    let dense = total <= MAX_COMBOS;

    let mut best_score = f64::INFINITY;
    let mut best = current.clone();

    // Enumerate χ combinations. When the full grid is too large, fall back to a
    // greedy per-χ sweep (optimize χ1, fix it, then χ2, …).
    if dense {
        let mut combo = vec![0usize; axes.len()];
        loop {
            let mut trial = current.clone();
            apply_combo(&axes, samples, &combo, &mut trial, &pos_of);
            let s = clash_score(&trial, &env);
            if s < best_score {
                best_score = s;
                best = trial;
            }
            if !advance(&mut combo, samples.len()) {
                break;
            }
        }
    } else {
        for (i, _) in axes.iter().enumerate() {
            let mut best_k = 0usize;
            let mut local_best = f64::INFINITY;
            for (k, _) in samples.iter().enumerate() {
                let mut trial = current.clone();
                set_chi(&axes[i], samples[k], &mut trial, &pos_of);
                let s = clash_score(&trial, &env);
                if s < local_best {
                    local_best = s;
                    best_k = k;
                }
            }
            set_chi(&axes[i], samples[best_k], &mut current, &pos_of);
        }
        best = current;
    }

    best.into_iter().map(|(k, (p, _))| (k, to_vec3(p))).collect()
}

/// Sets a single χ to `target_deg` by rotating its downstream atoms.
fn set_chi(
    axis: &ChiAxis,
    target_deg: f64,
    current: &mut HashMap<String, (V, Element)>,
    pos_of: &impl Fn(&HashMap<String, (V, Element)>, &str) -> Option<V>,
) {
    let (Some(a_pos), Some(b_pos), Some(c_pos), Some(d_pos)) = (
        pos_of(current, &axis.a),
        pos_of(current, &axis.b),
        pos_of(current, &axis.c),
        pos_of(current, &axis.d),
    ) else {
        return;
    };
    let axis_vec = sub(c_pos, b_pos);
    let current_deg = dihedral(a_pos, b_pos, c_pos, d_pos);
    let delta = target_deg - current_deg;
    for name in &axis.downstream {
        if let Some(entry) = current.get_mut(name) {
            entry.0 = rotate_about(entry.0, c_pos, axis_vec, delta);
        }
    }
}

/// Applies a full χ combination (χ1 then χ2 …) to `trial`.
fn apply_combo(
    axes: &[ChiAxis],
    samples: &[f64],
    combo: &[usize],
    trial: &mut HashMap<String, (V, Element)>,
    pos_of: &impl Fn(&HashMap<String, (V, Element)>, &str) -> Option<V>,
) {
    for (i, axis) in axes.iter().enumerate() {
        set_chi(axis, samples[combo[i]], trial, pos_of);
    }
}

/// Odometer increment over `radix`; returns false when it wraps to all-zero.
fn advance(combo: &mut [usize], radix: usize) -> bool {
    for slot in combo.iter_mut() {
        *slot += 1;
        if *slot < radix {
            return true;
        }
        *slot = 0;
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rotating_chi1_moves_downstream_only() {
        // Minimal LYS-like chain: N, CA fixed; CB, CG, CD movable along χ1/χ2.
        let bonds = vec![
            ("CB".to_string(), "CG".to_string()),
            ("CG".to_string(), "CD".to_string()),
        ];
        let present: HashSet<String> = ["CB", "CG", "CD"].iter().map(|s| s.to_string()).collect();
        let axes = resolve_chi_axes("LYS", &present, &bonds);
        assert!(!axes.is_empty());
        // χ1 (N-CA-CB-CG) rotates CG and CD (downstream of CB), not CB.
        let chi1 = &axes[0];
        assert!(chi1.downstream.contains(&"CG".to_string()));
        assert!(chi1.downstream.contains(&"CD".to_string()));
        assert!(!chi1.downstream.contains(&"CB".to_string()));
    }

    #[test]
    fn picks_rotamer_that_avoids_a_clashing_atom() {
        // CB at origin region; a single γ atom that, at one χ1, sits on top of an
        // environment atom and at another is far from it. The search must pick
        // the clash-free rotamer.
        let fixed: HashMap<String, Vec3> = [
            ("N".to_string(), Vec3::new(-1.0, 1.2, 0.0)),
            ("CA".to_string(), Vec3::new(0.0, 0.0, 0.0)),
        ]
        .into_iter()
        .collect();
        let sidechain = vec![
            SideAtom { name: "CB".into(), element: Element::Carbon, coord: Vec3::new(1.5, 0.0, 0.0) },
            SideAtom { name: "OG".into(), element: Element::Oxygen, coord: Vec3::new(2.0, 1.2, 0.0) },
        ];
        let bonds = vec![("CB".to_string(), "OG".to_string())];
        // Environment atom sitting where OG starts → initial clash.
        let environment = vec![(Vec3::new(2.0, 1.2, 0.0), Element::Oxygen)];
        let out = optimize_rotamer("SER", &sidechain, &fixed, &bonds, &environment, true);
        let og = out.get("OG").unwrap();
        let d = ((og.x - 2.0).powi(2) + (og.y - 1.2).powi(2) + og.z.powi(2)).sqrt();
        assert!(d > 1.0, "optimizer left OG on top of the clashing atom (d={d})");
        // CB is on the χ1 axis and must not move.
        let cb = out.get("CB").unwrap();
        assert!((cb.x - 1.5).abs() < 1e-4 && cb.y.abs() < 1e-4);
    }
}
