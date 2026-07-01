//! Hydrogen-bond network optimization for rotatable polar hydrogens.
//!
//! Hydroxyl (Ser/Thr/Tyr) and thiol (Cys) hydrogens, and any single H on a
//! polar heavy atom with one heavy neighbour, are free to rotate about the
//! donor–pivot bond. After protonation they are placed at an arbitrary
//! dihedral; this module reorients each one to point at a nearby acceptor,
//! forming hydrogen bonds while avoiding clashes.
//!
//! "Local" optimization orients each H independently against the fixed
//! acceptors; "global" runs several sweeps so hydrogens that can donate to one
//! another settle into a consistent network and do not pile onto the same spot.

use lin_alg::f32::Vec3;

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

fn rotate_about(point: V, origin: V, axis: V, angle_deg: f64) -> V {
    let k = normalize(axis);
    let p = sub(point, origin);
    let theta = angle_deg.to_radians();
    let (s, cth) = theta.sin_cos();
    let kxp = cross(k, p);
    let kdp = dot(k, p);
    let r = [
        p[0] * cth + kxp[0] * s + k[0] * kdp * (1.0 - cth),
        p[1] * cth + kxp[1] * s + k[1] * kdp * (1.0 - cth),
        p[2] * cth + kxp[2] * s + k[2] * kdp * (1.0 - cth),
    ];
    [origin[0] + r[0], origin[1] + r[1], origin[2] + r[2]]
}

/// One rotatable polar hydrogen: its current position `h`, the polar heavy atom
/// `donor` it is bonded to, and the `pivot` heavy atom across the rotatable
/// `donor`–`pivot` bond (e.g. CB for a Ser OG-HG hydroxyl).
#[derive(Debug, Clone, Copy)]
pub struct PolarH {
    pub h: Vec3,
    pub donor: Vec3,
    pub pivot: Vec3,
}

/// Ideal H⋯acceptor distance (Å) and the half-width of the reward well.
const HB_DIST: f64 = 1.9;
const HB_WIDTH: f64 = 0.5;
/// Below this H⋯atom distance a steric clash is penalized.
const CLASH_DIST: f64 = 1.5;

/// Scores a candidate H position against acceptors (heavy O/N) and other
/// hydrogens. Rewards near-linear H-bonds; penalizes clashes and H/H overlap.
fn score_h(h: V, donor: V, acceptors: &[V], others: &[V]) -> f64 {
    let dir_hd = normalize(sub(donor, h));
    let mut score = 0.0;
    for &a in acceptors {
        let d = norm(sub(a, h));
        // Skip the donor's own heavy atoms (acceptor list may include them).
        if !(0.3..=3.2).contains(&d) {
            continue;
        }
        let dir_ha = normalize(sub(a, h));
        // Linear D-H⋯A: donor and acceptor on opposite sides of H.
        let linear = (-dot(dir_hd, dir_ha)).max(0.0);
        let well = (-((d - HB_DIST) / HB_WIDTH).powi(2)).exp();
        score += linear * well;
        if d < CLASH_DIST {
            score -= (CLASH_DIST - d) * 2.0;
        }
    }
    for &o in others {
        let d = norm(sub(o, h));
        if d < CLASH_DIST {
            score -= (CLASH_DIST - d) * 1.0;
        }
    }
    score
}

/// Reorients each rotatable polar H about its donor–pivot bond to maximize the
/// hydrogen-bond network score against `acceptors`. `passes` ≥ 2 enables global
/// mode, where each sweep sees the other hydrogens' current positions so the
/// network settles. Returns the optimized H positions in input order.
pub fn optimize_polar_hydrogens(hs: &[PolarH], acceptors: &[Vec3], passes: usize) -> Vec<Vec3> {
    let acc: Vec<V> = acceptors.iter().map(|p| v(*p)).collect();
    let mut pos: Vec<V> = hs.iter().map(|p| v(p.h)).collect();
    let passes = passes.max(1);

    for _ in 0..passes {
        for (i, ph) in hs.iter().enumerate() {
            let donor = v(ph.donor);
            let axis = sub(v(ph.pivot), donor);
            if norm(axis) < 1e-6 {
                continue;
            }
            // Other optimized hydrogens act as soft repulsors during global sweeps.
            let others: Vec<V> = pos
                .iter()
                .enumerate()
                .filter(|(j, _)| *j != i)
                .map(|(_, p)| *p)
                .collect();

            let mut best = pos[i];
            let mut best_score = f64::NEG_INFINITY;
            // Sweep the H around the donor–pivot bond.
            for step in 0..36 {
                let angle = step as f64 * 10.0;
                let cand = rotate_about(pos[i], donor, axis, angle);
                let s = score_h(cand, donor, &acc, &others);
                if s > best_score {
                    best_score = s;
                    best = cand;
                }
            }
            pos[i] = best;
        }
    }

    pos.into_iter().map(to_vec3).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hydroxyl_h_swings_toward_a_lone_acceptor() {
        // Donor OG at origin, pivot CB along -x: the H rotates on a cone about
        // the +x axis (its x stays 0.3). It starts pointing -y; an off-axis
        // acceptor sits at +y, so optimization should swing the H toward it.
        let donor = Vec3::new(0.0, 0.0, 0.0);
        let pivot = Vec3::new(-1.4, 0.0, 0.0);
        let h = Vec3::new(0.3, -0.95, 0.0); // ~1.0 Å O-H, pointing down
        let acceptor = Vec3::new(0.6, 2.7, 0.0); // off-axis, reachable at +y
        let out = optimize_polar_hydrogens(&[PolarH { h, donor, pivot }], &[acceptor], 3);
        let p = out[0];
        let d_start = ((h.x - acceptor.x).powi(2) + (h.y - acceptor.y).powi(2)).sqrt();
        let d_end = ((p.x - acceptor.x).powi(2) + (p.y - acceptor.y).powi(2)).sqrt();
        assert!(
            d_end < d_start,
            "H did not move toward the acceptor (start {d_start}, end {d_end})"
        );
        assert!(
            p.y > 0.0,
            "H should have swung to the +y side toward the acceptor"
        );
        // Bond length to the donor is preserved by the rigid rotation.
        let r = ((p.x - donor.x).powi(2) + (p.y - donor.y).powi(2) + p.z.powi(2)).sqrt();
        assert!((r - 1.0).abs() < 0.05, "O-H length not preserved: {r}");
    }
}
