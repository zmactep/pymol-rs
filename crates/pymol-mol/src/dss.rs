//! Define Secondary Structure (DSS) algorithm
//!
//! Implements PyMOL's DSS algorithm for calculating secondary structure
//! using backbone hydrogen bonds as the primary mechanism with phi/psi
//! dihedral angles as a secondary filter.
//!
//! # References
//!
//! - PyMOL's layer3/Selector.cpp - SelectorAssignSS function
//! - PyMOL's layer2/ObjectMolecule2.cpp - ObjectMoleculeGetCheckHBond function
//! - PyMOL's layer1/SettingInfo.h - Default angle thresholds

use bitflags::bitflags;
use lin_alg::f32::Vec3;

use crate::spatial::SpatialGrid;


use crate::coordset::CoordSet;
use crate::index::AtomIndex;
use crate::molecule::ObjectMolecule;
use crate::residue::ResidueView;
use crate::secondary::SecondaryStructure;

// ============================================================================
// DSS Settings (from PyMOL's SettingInfo.h)
// ============================================================================

/// Settings for DSS secondary structure assignment
#[derive(Debug, Clone)]
pub struct DssSettings {
    /// Target phi angle for helix (degrees)
    pub helix_phi_target: f32,
    /// Target psi angle for helix (degrees)
    pub helix_psi_target: f32,
    /// Phi angle tolerance for helix inclusion (degrees)
    pub helix_phi_include: f32,
    /// Psi angle tolerance for helix inclusion (degrees)
    pub helix_psi_include: f32,
    /// Phi angle tolerance for helix exclusion (degrees)
    pub helix_phi_exclude: f32,
    /// Psi angle tolerance for helix exclusion (degrees)
    pub helix_psi_exclude: f32,
    /// Target phi angle for strand (degrees)
    pub strand_phi_target: f32,
    /// Target psi angle for strand (degrees)
    pub strand_psi_target: f32,
    /// Phi angle tolerance for strand inclusion (degrees)
    pub strand_phi_include: f32,
    /// Psi angle tolerance for strand inclusion (degrees)
    pub strand_psi_include: f32,
    /// Phi angle tolerance for strand exclusion (degrees)
    pub strand_phi_exclude: f32,
    /// Psi angle tolerance for strand exclusion (degrees)
    pub strand_psi_exclude: f32,
}

impl Default for DssSettings {
    fn default() -> Self {
        Self {
            helix_phi_target: -57.0,
            helix_psi_target: -48.0,
            helix_phi_include: 55.0,
            helix_psi_include: 55.0,
            helix_phi_exclude: 85.0,
            helix_psi_exclude: 85.0,
            strand_phi_target: -129.0,
            strand_psi_target: 124.0,
            strand_phi_include: 40.0,
            strand_psi_include: 40.0,
            strand_phi_exclude: 100.0,
            strand_psi_exclude: 90.0,
        }
    }
}

// ============================================================================
// Dihedral Angle Calculation
// ============================================================================

/// Calculate the dihedral angle between four points
///
/// Based on PyMOL's get_dihedral3f function from layer0/Vector.cpp
///
/// # Returns
/// Dihedral angle in radians, range [-PI, PI]
fn dihedral_angle(v0: Vec3, v1: Vec3, v2: Vec3, v3: Vec3) -> f32 {
    let d21 = v2 - v1;
    let d01 = v0 - v1;
    let d32 = v3 - v2;

    let d21_len = d21.magnitude();
    if d21_len < 1e-6 {
        return angle_between(d01, d32);
    }

    let dd1 = d21.cross(d01);
    let dd3 = d21.cross(d32);

    let dd1_len = dd1.magnitude();
    let dd3_len = dd3.magnitude();

    if dd1_len < 1e-6 || dd3_len < 1e-6 {
        return angle_between(d01, d32);
    }

    let mut result = angle_between(dd1, dd3);

    let pos_d = d21.cross(dd1);
    if dd3.dot(pos_d) < 0.0 {
        result = -result;
    }

    result
}

fn angle_between(v1: Vec3, v2: Vec3) -> f32 {
    let len1 = v1.magnitude();
    let len2 = v2.magnitude();

    if len1 < 1e-6 || len2 < 1e-6 {
        return 0.0;
    }

    let cos_angle = v1.dot(v2) / (len1 * len2);
    cos_angle.clamp(-1.0, 1.0).acos()
}

// ============================================================================
// Phi/Psi Calculation
// ============================================================================

/// Phi and psi angles for a residue
#[derive(Debug, Clone, Copy)]
pub struct PhiPsi {
    pub phi: f32,
    pub psi: f32,
}

// ============================================================================
// Classification Flags (matching PyMOL's cSS* flags from Selector.cpp)
// ============================================================================

bitflags! {
    #[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
    struct SsFlags: u32 {
        const HELIX_3_HBOND         = 0x0001;
        const HELIX_4_HBOND         = 0x0002;
        const HELIX_5_HBOND         = 0x0004;
        const GOT_PHI_PSI           = 0x0008;
        const PHI_PSI_HELIX         = 0x0010;
        const PHI_PSI_NOT_HELIX     = 0x0020;
        const PHI_PSI_STRAND        = 0x0040;
        const PHI_PSI_NOT_STRAND    = 0x0080;
        const ANTI_STRAND_SINGLE_HB = 0x0100;
        const ANTI_STRAND_DOUBLE_HB = 0x0200;
        const ANTI_STRAND_BULGE_HB  = 0x0400;
        const ANTI_STRAND_SKIP      = 0x0800;
        const PARA_STRAND_SINGLE_HB = 0x1000;
        const PARA_STRAND_DOUBLE_HB = 0x2000;
        const PARA_STRAND_SKIP      = 0x4000;
    }
}

const HELIX_HBOND_FLAGS: SsFlags = SsFlags::HELIX_3_HBOND
    .union(SsFlags::HELIX_4_HBOND)
    .union(SsFlags::HELIX_5_HBOND);

const SS_MAX_HBOND: usize = 6;
const SS_BREAK_SIZE: usize = 5;

// H-bond criteria (Kabsch & Sander inspired, from PyMOL's SelectorAssignSS)
const HBOND_MAX_ANGLE: f32 = 63.0;
const HBOND_MAX_DIST_AT_MAX_ANGLE: f32 = 3.2;
const HBOND_MAX_DIST_AT_ZERO: f32 = 4.0;
const HBOND_POWER_A: f32 = 1.6;
const HBOND_POWER_B: f32 = 5.0;

// ============================================================================
// Spatial Grid for H-Bond Detection — see crate::spatial::SpatialGrid
// ============================================================================

/// Internal structure for residue data during DSS calculation
#[derive(Debug, Clone)]
struct ResidueData {
    ca_idx: Option<AtomIndex>,
    n_idx: Option<AtomIndex>,
    c_idx: Option<AtomIndex>,
    o_idx: Option<AtomIndex>,
    chain: String,
    resv: i32,
    real: bool,
    phi_psi: Option<PhiPsi>,
    flags: SsFlags,
    ss: u8,
    acc: Vec<usize>,
    don: Vec<usize>,
}

impl Default for ResidueData {
    fn default() -> Self {
        Self {
            ca_idx: None,
            n_idx: None,
            c_idx: None,
            o_idx: None,
            chain: String::new(),
            resv: 0,
            real: false,
            phi_psi: None,
            flags: SsFlags::empty(),
            ss: b'L',
            acc: Vec::new(),
            don: Vec::new(),
        }
    }
}

// ============================================================================
// DSS Algorithm
// ============================================================================

/// Assign secondary structure to a molecule using the DSS algorithm
///
/// This implements PyMOL's DSS algorithm which uses hydrogen bonds as the
/// primary mechanism for secondary structure detection, with phi/psi dihedral
/// angles as a secondary filter.
pub fn assign_secondary_structure(
    molecule: &mut ObjectMolecule,
    state: usize,
    settings: &DssSettings,
) -> usize {
    let coord_set = match molecule.get_coord_set(state) {
        Some(cs) => cs,
        None => return 0,
    };

    // Step 1: Collect residue data and insert chain breaks
    let mut res = collect_residue_data_with_breaks(molecule);
    let n_res = res.len();

    if n_res < 2 * SS_BREAK_SIZE + 1 {
        return 0;
    }

    // Step 2: Verify coordinates are present
    verify_coordinates(&mut res, coord_set);

    // Step 3: Detect backbone hydrogen bonds
    detect_hydrogen_bonds(&mut res, coord_set, molecule);

    // Step 4: Compute phi/psi angles and classify
    calculate_all_phi_psi(&mut res, coord_set);
    classify_residues(&mut res, settings);

    // Step 5: Detect H-bond patterns (helix, sheet)
    detect_hbond_patterns(&mut res);

    // Step 6: Assign secondary structure
    assign_helices(&mut res);
    assign_sheets(&mut res);

    // Step 7: Validate and clean up
    validate_assignments(&mut res);

    // Step 8: Apply to molecule
    apply_ss_assignments(molecule, &res)
}

// ============================================================================
// Residue Collection with Chain Breaks
// ============================================================================

fn collect_residue_data_with_breaks(molecule: &ObjectMolecule) -> Vec<ResidueData> {
    let raw = collect_raw_residues(molecule);
    if raw.is_empty() {
        return Vec::new();
    }

    let mut result = Vec::with_capacity(raw.len() + 2 * SS_BREAK_SIZE + 10);

    // Add initial break
    for _ in 0..SS_BREAK_SIZE {
        result.push(ResidueData::default());
    }

    for (i, r) in raw.iter().enumerate() {
        if i > 0 {
            let need_break = if raw[i - 1].chain != r.chain {
                true
            } else {
                // Check for peptide bond: prev C bonded to curr N
                match (raw[i - 1].c_idx, r.n_idx) {
                    (Some(c), Some(n)) => molecule.find_bond(c, n).is_none(),
                    _ => true,
                }
            };

            if need_break {
                for _ in 0..SS_BREAK_SIZE {
                    result.push(ResidueData::default());
                }
            }
        }

        result.push(r.clone());
    }

    // Add final break
    for _ in 0..SS_BREAK_SIZE {
        result.push(ResidueData::default());
    }

    result
}

fn collect_raw_residues(molecule: &ObjectMolecule) -> Vec<ResidueData> {
    let mut residue_data = Vec::new();

    for residue in molecule.residues() {
        if !crate::residue::is_amino_acid(&residue.key.resn) {
            continue;
        }

        let ca_idx = find_atom_by_name(&residue, "CA");
        let n_idx = find_atom_by_name(&residue, "N");
        let c_idx = find_atom_by_name(&residue, "C");
        let o_idx = find_atom_by_name(&residue, "O");

        // Require CA, N, C, O for full backbone
        if ca_idx.is_some() && n_idx.is_some() && c_idx.is_some() && o_idx.is_some() {
            residue_data.push(ResidueData {
                ca_idx,
                n_idx,
                c_idx,
                o_idx,
                chain: residue.key.chain.clone(),
                resv: residue.key.resv,
                real: true,
                phi_psi: None,
                flags: SsFlags::empty(),
                ss: b'L',
                acc: Vec::new(),
                don: Vec::new(),
            });
        }
    }

    residue_data
}

fn find_atom_by_name(residue: &ResidueView, name: &str) -> Option<AtomIndex> {
    for (idx, atom) in residue.iter_indexed() {
        if &*atom.name == name {
            return Some(idx);
        }
    }
    None
}

// ============================================================================
// Coordinate Verification
// ============================================================================

fn verify_coordinates(res: &mut [ResidueData], coord_set: &CoordSet) {
    for r in res.iter_mut() {
        if !r.real {
            continue;
        }
        // Verify all backbone atoms have coordinates
        let has_coords = [r.n_idx, r.o_idx, r.c_idx, r.ca_idx]
            .iter()
            .all(|idx| idx.and_then(|i| coord_set.get_atom_coord(i)).is_some());
        if !has_coords {
            r.real = false;
        }
    }
}

// ============================================================================
// Hydrogen Bond Detection
// ============================================================================

fn detect_hydrogen_bonds(
    res: &mut [ResidueData],
    coord_set: &CoordSet,
    molecule: &ObjectMolecule,
) {
    let n_res = res.len();
    let cutoff = HBOND_MAX_DIST_AT_ZERO;

    // Precompute factor_a and factor_b for angle-based cutoff interpolation
    let factor_a = 0.5 / (HBOND_MAX_ANGLE as f64).powf(HBOND_POWER_A as f64);
    let factor_b = 0.5 / (HBOND_MAX_ANGLE as f64).powf(HBOND_POWER_B as f64);

    // Precompute O and N positions for all residues
    let o_positions: Vec<Option<Vec3>> = (0..n_res)
        .map(|i| {
            if !res[i].real {
                return None;
            }
            res[i].o_idx.and_then(|idx| coord_set.get_atom_coord(idx))
        })
        .collect();

    let n_positions: Vec<Option<Vec3>> = (0..n_res)
        .map(|i| {
            if !res[i].real {
                return None;
            }
            res[i].n_idx.and_then(|idx| coord_set.get_atom_coord(idx))
        })
        .collect();

    // Precompute N-H directions for all residues (avoids repeated bonded_atoms() calls)
    // If an actual H atom is bonded to N, use N→H directly.
    // Otherwise, approximate H direction as opposite to average of heavy atom bonds.
    let nh_directions: Vec<Option<Vec3>> = (0..n_res)
        .map(|i| {
            if !res[i].real {
                return None;
            }
            let n_idx = res[i].n_idx?;
            let n_pos = coord_set.get_atom_coord(n_idx)?;
            let bonded: smallvec::SmallVec<[AtomIndex; 4]> =
                molecule.bonded_atoms_iter(n_idx).collect();
            if bonded.is_empty() {
                return None;
            }

            // Check if there's a bonded hydrogen — use its position directly
            for &b_idx in &bonded {
                if let Some(b_atom) = molecule.get_atom(b_idx) {
                    if b_atom.is_hydrogen() {
                        if let Some(h_pos) = coord_set.get_atom_coord(b_idx) {
                            let d = h_pos - n_pos;
                            let len = d.magnitude();
                            if len > 1e-6 {
                                return Some(d * (1.0 / len));
                            }
                        }
                    }
                }
            }

            // No hydrogen found — approximate from heavy atom geometry
            let mut avg_dir = Vec3::new(0.0, 0.0, 0.0);
            let mut count = 0;
            for &b_idx in &bonded {
                if let Some(b_pos) = coord_set.get_atom_coord(b_idx) {
                    let d = b_pos - n_pos;
                    let len = d.magnitude();
                    if len > 1e-6 {
                        avg_dir = avg_dir + d * (1.0 / len);
                        count += 1;
                    }
                }
            }
            if count == 0 {
                return None;
            }
            let h_dir = avg_dir * -1.0;
            let len = h_dir.magnitude();
            if len < 1e-6 {
                return None;
            }
            Some(h_dir * (1.0 / len))
        })
        .collect();

    // Build spatial grid indexed by O atom positions (acceptors)
    let mut grid = SpatialGrid::with_capacity(cutoff, o_positions.len());
    for (a0, o_pos) in o_positions.iter().enumerate() {
        if let Some(pos) = o_pos {
            grid.insert(*pos, a0);
        }
    }

    // For each donor N, query grid for nearby acceptor O atoms
    let mut neighbors = Vec::new();
    for a1 in 0..n_res {
        let n_pos = match n_positions[a1] {
            Some(p) => p,
            None => continue,
        };

        grid.query_neighbors(n_pos, &mut neighbors);

        for &a0 in &neighbors {
            // Exclude adjacent residues (within 2 positions in the padded array)
            let dist_idx = if a0 > a1 { a0 - a1 } else { a1 - a0 };
            if dist_idx <= 2 {
                continue;
            }

            // Both acc and don lists are already at capacity
            if res[a0].acc.len() >= SS_MAX_HBOND && res[a1].don.len() >= SS_MAX_HBOND {
                continue;
            }

            let o_pos = o_positions[a0].unwrap(); // safe: grid only contains valid positions

            let don_to_acc = o_pos - n_pos;
            let dist = don_to_acc.magnitude();
            if dist > cutoff || dist < 1e-6 {
                continue;
            }

            if let Some(h_dir) = nh_directions[a1] {
                // Compute angle between N-H direction and N→O direction
                let don_to_acc_norm = don_to_acc * (1.0 / dist);
                let cos_angle = h_dir.dot(don_to_acc_norm);
                let angle = if cos_angle >= 1.0 {
                    0.0
                } else if cos_angle <= 0.0 {
                    90.0
                } else {
                    cos_angle.acos().to_degrees() as f64
                };

                if angle > HBOND_MAX_ANGLE as f64 {
                    continue;
                }

                // Interpolated distance cutoff
                let curve = (angle.powf(HBOND_POWER_A as f64) * factor_a)
                    + (angle.powf(HBOND_POWER_B as f64) * factor_b);
                let dist_cutoff = (HBOND_MAX_DIST_AT_MAX_ANGLE as f64 * curve)
                    + (HBOND_MAX_DIST_AT_ZERO as f64 * (1.0 - curve));

                if (dist as f64) > dist_cutoff {
                    continue;
                }
            } else {
                // No H direction available — fall back to distance-only check
                if dist > 3.5 {
                    continue;
                }
            }

            // Store acceptor relationship (a0 is acceptor via O, a1 is donor via N)
            if res[a0].acc.len() < SS_MAX_HBOND {
                res[a0].acc.push(a1);
            }
            if res[a1].don.len() < SS_MAX_HBOND {
                res[a1].don.push(a0);
            }
        }
    }
}

// ============================================================================
// Phi/Psi Calculation and Classification
// ============================================================================

fn calculate_all_phi_psi(res: &mut [ResidueData], coord_set: &CoordSet) {
    let n_res = res.len();

    for i in 0..n_res {
        if !res[i].real {
            continue;
        }
        // Need previous residue to be real for phi/psi
        if i == 0 || !res[i - 1].real {
            continue;
        }

        // Get backbone positions for phi/psi
        let c_prev = res[i - 1]
            .c_idx
            .and_then(|idx| coord_set.get_atom_coord(idx));
        let n_curr = res[i]
            .n_idx
            .and_then(|idx| coord_set.get_atom_coord(idx));
        let ca_curr = res[i]
            .ca_idx
            .and_then(|idx| coord_set.get_atom_coord(idx));
        let c_curr = res[i]
            .c_idx
            .and_then(|idx| coord_set.get_atom_coord(idx));

        // For psi, we need N of next real residue
        let n_next = if i + 1 < n_res && res[i + 1].real {
            res[i + 1]
                .n_idx
                .and_then(|idx| coord_set.get_atom_coord(idx))
        } else {
            None
        };

        if let (Some(c_prev), Some(n_curr), Some(ca_curr), Some(c_curr)) =
            (c_prev, n_curr, ca_curr, c_curr)
        {
            // We can compute phi even without n_next
            let phi = dihedral_angle(c_prev, n_curr, ca_curr, c_curr).to_degrees();

            if let Some(n_next) = n_next {
                let psi = dihedral_angle(n_curr, ca_curr, c_curr, n_next).to_degrees();
                res[i].phi_psi = Some(PhiPsi { phi, psi });
                res[i].flags |= SsFlags::GOT_PHI_PSI;
            }
        }
    }
}

fn classify_residues(res: &mut [ResidueData], settings: &DssSettings) {
    for r in res.iter_mut() {
        if let Some(phi_psi) = r.phi_psi {
            let helix_phi_delta = angle_delta(phi_psi.phi, settings.helix_phi_target);
            let helix_psi_delta = angle_delta(phi_psi.psi, settings.helix_psi_target);
            let strand_phi_delta = angle_delta(phi_psi.phi, settings.strand_phi_target);
            let strand_psi_delta = angle_delta(phi_psi.psi, settings.strand_psi_target);

            if helix_phi_delta > settings.helix_phi_exclude
                || helix_psi_delta > settings.helix_psi_exclude
            {
                r.flags |= SsFlags::PHI_PSI_NOT_HELIX;
            } else if helix_phi_delta < settings.helix_phi_include
                && helix_psi_delta < settings.helix_psi_include
            {
                r.flags |= SsFlags::PHI_PSI_HELIX;
            }

            if strand_phi_delta > settings.strand_phi_exclude
                || strand_psi_delta > settings.strand_psi_exclude
            {
                r.flags |= SsFlags::PHI_PSI_NOT_STRAND;
            } else if strand_phi_delta < settings.strand_phi_include
                && strand_psi_delta < settings.strand_psi_include
            {
                r.flags |= SsFlags::PHI_PSI_STRAND;
            }
        }
    }
}

fn angle_delta(angle: f32, target: f32) -> f32 {
    let mut delta = (angle - target).abs();
    if delta > 180.0 {
        delta = 360.0 - delta;
    }
    delta
}

// ============================================================================
// H-Bond Pattern Detection
// ============================================================================

fn detect_hbond_patterns(res: &mut [ResidueData]) {
    let n_res = res.len();

    for a in SS_BREAK_SIZE..n_res.saturating_sub(SS_BREAK_SIZE) {
        if !res[a].real {
            continue;
        }

        // Helix patterns: check if residue is acceptor for i+3,4,5 or donor for i-3,4,5
        let acc_list: Vec<usize> = res[a].acc.clone();
        for &acc in &acc_list {
            if acc == a + 3 {
                res[a].flags |= SsFlags::HELIX_3_HBOND;
            }
            if acc == a + 4 {
                res[a].flags |= SsFlags::HELIX_4_HBOND;
            }
            if acc == a + 5 {
                res[a].flags |= SsFlags::HELIX_5_HBOND;
            }
        }

        let don_list: Vec<usize> = res[a].don.clone();
        for &don in &don_list {
            if don == a.wrapping_sub(3) {
                res[a].flags |= SsFlags::HELIX_3_HBOND;
            }
            if don == a.wrapping_sub(4) {
                res[a].flags |= SsFlags::HELIX_4_HBOND;
            }
            if don == a.wrapping_sub(5) {
                res[a].flags |= SsFlags::HELIX_5_HBOND;
            }
        }

        // Antiparallel double H-bond: i accepts j, j accepts i
        for &acc_j in &acc_list {
            if acc_j >= n_res || !res[acc_j].real {
                continue;
            }
            let acc_j_list: Vec<usize> = res[acc_j].acc.clone();
            for &acc_k in &acc_j_list {
                if acc_k == a {
                    res[a].flags |= SsFlags::ANTI_STRAND_DOUBLE_HB;
                    res[acc_j].flags |= SsFlags::ANTI_STRAND_DOUBLE_HB;
                }
            }
        }

        // Antiparallel bulge: i accepts j, j+1 accepts i
        for &acc_j in &acc_list {
            if acc_j >= n_res || !res[acc_j].real {
                continue;
            }
            let j_plus_1 = acc_j + 1;
            if j_plus_1 < n_res && res[j_plus_1].real {
                let acc_jp1_list: Vec<usize> = res[j_plus_1].acc.clone();
                for &acc_k in &acc_jp1_list {
                    if acc_k == a {
                        res[a].flags |= SsFlags::ANTI_STRAND_DOUBLE_HB;
                        res[j_plus_1].flags |= SsFlags::ANTI_STRAND_BULGE_HB;
                        res[acc_j].flags |= SsFlags::ANTI_STRAND_BULGE_HB;
                    }
                }
            }
        }

        // Antiparallel ladder: i accepts j, (j-2) accepts (i+2)
        if a + 2 < n_res && res[a + 1].real && res[a + 2].real {
            for &acc_j in &acc_list {
                if acc_j < 2 || !res[acc_j].real {
                    continue;
                }
                let j_minus_2 = acc_j - 2;
                if !res[j_minus_2].real {
                    continue;
                }
                let acc_jm2_list: Vec<usize> = res[j_minus_2].acc.clone();
                for &acc_k in &acc_jm2_list {
                    if acc_k == a + 2 {
                        res[a].flags |= SsFlags::ANTI_STRAND_SINGLE_HB;
                        res[a + 1].flags |= SsFlags::ANTI_STRAND_SKIP;
                        res[a + 2].flags |= SsFlags::ANTI_STRAND_SINGLE_HB;
                        res[j_minus_2].flags |= SsFlags::ANTI_STRAND_SINGLE_HB;
                        if acc_j >= j_minus_2 + 2 {
                            res[j_minus_2 + 1].flags |= SsFlags::ANTI_STRAND_SKIP;
                        }
                        res[acc_j].flags |= SsFlags::ANTI_STRAND_SINGLE_HB;
                    }
                }
            }
        }

        // Parallel ladder: i accepts j, j accepts (i+2)
        if a + 2 < n_res && res[a + 1].real && res[a + 2].real {
            for &acc_j in &acc_list {
                if acc_j >= n_res || !res[acc_j].real {
                    continue;
                }
                let acc_j_list2: Vec<usize> = res[acc_j].acc.clone();
                for &acc_k in &acc_j_list2 {
                    if acc_k == a + 2 {
                        res[a].flags |= SsFlags::PARA_STRAND_SINGLE_HB;
                        res[a + 1].flags |= SsFlags::PARA_STRAND_SKIP;
                        res[a + 2].flags |= SsFlags::PARA_STRAND_SINGLE_HB;
                        res[acc_j].flags |= SsFlags::PARA_STRAND_DOUBLE_HB;
                    }
                }
            }
        }
    }
}

// ============================================================================
// Secondary Structure Assignment
// ============================================================================

fn assign_helices(res: &mut [ResidueData]) {
    let n_res = res.len();

    // Pass 1: Clean internal helical residues using H-bonds
    // Three consecutive residues with helix H-bonds and acceptable phi/psi.
    // Also require that neighbors don't have explicitly non-helical geometry,
    // unless they also have helical geometry (matching PyMOL's neighbor check).
    for a in SS_BREAK_SIZE..n_res.saturating_sub(SS_BREAK_SIZE) {
        if !res[a].real {
            continue;
        }

        if res[a - 1].flags.intersects(HELIX_HBOND_FLAGS)
            && res[a].flags.intersects(HELIX_HBOND_FLAGS)
            && res[a + 1].flags.intersects(HELIX_HBOND_FLAGS)
            && !res[a].flags.contains(SsFlags::PHI_PSI_NOT_HELIX)
            && !(res[a - 1].flags.contains(SsFlags::PHI_PSI_NOT_HELIX)
                && !res[a - 1].flags.contains(SsFlags::PHI_PSI_HELIX))
            && !(res[a + 1].flags.contains(SsFlags::PHI_PSI_NOT_HELIX)
                && !res[a + 1].flags.contains(SsFlags::PHI_PSI_HELIX))
        {
            res[a].ss = b'H';
        }
    }

    // Pass 2: Fill gaps — one missing H-bond residue surrounded by helix H-bonds + good geometry
    for a in SS_BREAK_SIZE..n_res.saturating_sub(SS_BREAK_SIZE) {
        if !res[a].real {
            continue;
        }

        if res[a - 2].flags.intersects(HELIX_HBOND_FLAGS)
            && res[a - 1].flags.intersects(HELIX_HBOND_FLAGS)
            && res[a - 1].flags.contains(SsFlags::PHI_PSI_HELIX)
            && res[a].flags.contains(SsFlags::PHI_PSI_HELIX)
            && res[a + 1].flags.intersects(HELIX_HBOND_FLAGS)
            && res[a + 1].flags.contains(SsFlags::PHI_PSI_HELIX)
            && res[a + 2].flags.intersects(HELIX_HBOND_FLAGS)
        {
            res[a].ss = b'h';
        }
    }

    // Promote 'h' to 'H' and add H-bond flags so cap extension can work
    for a in SS_BREAK_SIZE..n_res.saturating_sub(SS_BREAK_SIZE) {
        if res[a].real && res[a].ss == b'h' {
            res[a].flags |= HELIX_HBOND_FLAGS;
            res[a].ss = b'H';
        }
    }

    // Pass 3: Cap extension — extend helix boundaries using geometry.
    // Only extends if the candidate residue is not already excluded by phi/psi
    // AND its outer neighbor (the one not yet in the helix) is not anti-helical.
    for a in SS_BREAK_SIZE..n_res.saturating_sub(SS_BREAK_SIZE) {
        if !res[a].real || res[a].ss == b'H' {
            continue;
        }

        // Extend forward: a is N-terminal cap, a+1 is already H
        if res[a].flags.intersects(HELIX_HBOND_FLAGS)
            && res[a].flags.contains(SsFlags::PHI_PSI_HELIX)
            && !res[a].flags.contains(SsFlags::PHI_PSI_NOT_HELIX)
            && res[a + 1].flags.intersects(HELIX_HBOND_FLAGS)
            && res[a + 1].flags.contains(SsFlags::PHI_PSI_HELIX)
            && res[a + 2].flags.intersects(HELIX_HBOND_FLAGS)
            && res[a + 2].flags.contains(SsFlags::PHI_PSI_HELIX)
            && res[a + 1].ss == b'H'
            && !(res[a - 1].flags.contains(SsFlags::PHI_PSI_NOT_HELIX)
                && !res[a - 1].flags.contains(SsFlags::PHI_PSI_HELIX))
        {
            res[a].ss = b'H';
        }

        // Extend backward: a is C-terminal cap, a-1 is already H
        if res[a].flags.intersects(HELIX_HBOND_FLAGS)
            && res[a].flags.contains(SsFlags::PHI_PSI_HELIX)
            && !res[a].flags.contains(SsFlags::PHI_PSI_NOT_HELIX)
            && res[a - 1].flags.intersects(HELIX_HBOND_FLAGS)
            && res[a - 1].flags.contains(SsFlags::PHI_PSI_HELIX)
            && res[a - 2].flags.intersects(HELIX_HBOND_FLAGS)
            && res[a - 2].flags.contains(SsFlags::PHI_PSI_HELIX)
            && res[a - 1].ss == b'H'
            && !(res[a + 1].flags.contains(SsFlags::PHI_PSI_NOT_HELIX)
                && !res[a + 1].flags.contains(SsFlags::PHI_PSI_HELIX))
        {
            res[a].ss = b'H';
        }
    }
}

fn assign_sheets(res: &mut [ResidueData]) {
    let n_res = res.len();

    for a in SS_BREAK_SIZE..n_res.saturating_sub(SS_BREAK_SIZE) {
        if !res[a].real {
            continue;
        }

        // Antiparallel: double H-bond + phi/psi not excluded
        if res[a].flags.contains(SsFlags::ANTI_STRAND_DOUBLE_HB)
            && !res[a].flags.contains(SsFlags::PHI_PSI_NOT_STRAND)
        {
            res[a].ss = b'S';
        }

        // Antiparallel bulge: consecutive bulge pair
        if res[a].flags.contains(SsFlags::ANTI_STRAND_BULGE_HB)
            && a + 1 < n_res
            && res[a + 1].flags.contains(SsFlags::ANTI_STRAND_BULGE_HB)
        {
            res[a].ss = b'S';
            res[a + 1].ss = b'S';
        }

        // Antiparallel skip: between double HB anchors + phi/psi not excluded
        if res[a].flags.contains(SsFlags::ANTI_STRAND_SKIP)
            && !res[a].flags.contains(SsFlags::PHI_PSI_NOT_STRAND)
        {
            if res[a - 1]
                .flags
                .contains(SsFlags::ANTI_STRAND_DOUBLE_HB)
                && res[a + 1]
                    .flags
                    .intersects(SsFlags::ANTI_STRAND_SINGLE_HB | SsFlags::ANTI_STRAND_DOUBLE_HB)
            {
                res[a].ss = b'S';
            }
            if res[a - 1]
                .flags
                .intersects(SsFlags::ANTI_STRAND_SINGLE_HB | SsFlags::ANTI_STRAND_DOUBLE_HB)
                && res[a + 1]
                    .flags
                    .contains(SsFlags::ANTI_STRAND_DOUBLE_HB)
            {
                res[a].ss = b'S';
            }
        }

        // Antiparallel: open ladders with PHI_PSI_STRAND support
        if res[a - 1]
            .flags
            .intersects(SsFlags::ANTI_STRAND_SINGLE_HB | SsFlags::ANTI_STRAND_DOUBLE_HB)
            && res[a - 1].flags.contains(SsFlags::PHI_PSI_STRAND)
            && !res[a - 1].flags.contains(SsFlags::PHI_PSI_NOT_STRAND)
            && res[a].flags.contains(SsFlags::PHI_PSI_STRAND)
            && !res[a].flags.contains(SsFlags::PHI_PSI_NOT_STRAND)
            && res[a + 1]
                .flags
                .intersects(SsFlags::ANTI_STRAND_SINGLE_HB | SsFlags::ANTI_STRAND_DOUBLE_HB)
            && res[a + 1].flags.contains(SsFlags::PHI_PSI_STRAND)
        {
            res[a - 1].ss = b'S';
            res[a].ss = b'S';
            res[a + 1].ss = b'S';
        }

        // Parallel: double H-bond + phi/psi not excluded
        if res[a].flags.contains(SsFlags::PARA_STRAND_DOUBLE_HB)
            && !res[a].flags.contains(SsFlags::PHI_PSI_NOT_STRAND)
        {
            res[a].ss = b'S';
        }

        // Parallel skip: between HB anchors + phi/psi not excluded
        if res[a].flags.contains(SsFlags::PARA_STRAND_SKIP)
            && !res[a].flags.contains(SsFlags::PHI_PSI_NOT_STRAND)
        {
            if res[a - 1].flags.contains(SsFlags::PARA_STRAND_DOUBLE_HB)
                && res[a + 1]
                    .flags
                    .intersects(SsFlags::PARA_STRAND_SINGLE_HB | SsFlags::PARA_STRAND_DOUBLE_HB)
            {
                res[a].ss = b'S';
            }
            if res[a - 1]
                .flags
                .intersects(SsFlags::PARA_STRAND_SINGLE_HB | SsFlags::PARA_STRAND_DOUBLE_HB)
                && res[a + 1]
                    .flags
                    .contains(SsFlags::PARA_STRAND_DOUBLE_HB)
            {
                res[a].ss = b'S';
            }
        }

        // Parallel: open ladders with PHI_PSI_STRAND support
        if res[a - 1]
            .flags
            .intersects(SsFlags::PARA_STRAND_SINGLE_HB | SsFlags::PARA_STRAND_DOUBLE_HB)
            && res[a - 1].flags.contains(SsFlags::PHI_PSI_STRAND)
            && res[a].flags.contains(SsFlags::PARA_STRAND_SKIP)
            && res[a].flags.contains(SsFlags::PHI_PSI_STRAND)
            && res[a + 1]
                .flags
                .intersects(SsFlags::PARA_STRAND_SINGLE_HB | SsFlags::PARA_STRAND_DOUBLE_HB)
            && res[a + 1].flags.contains(SsFlags::PHI_PSI_STRAND)
        {
            res[a - 1].ss = b'S';
            res[a].ss = b'S';
            res[a + 1].ss = b'S';
        }
    }
}

// ============================================================================
// Validation / Cleanup
// ============================================================================

fn validate_assignments(res: &mut [ResidueData]) {
    let n_res = res.len();
    let mut repeat = true;

    while repeat {
        repeat = false;

        for a in SS_BREAK_SIZE..n_res.saturating_sub(SS_BREAK_SIZE) {
            if !res[a].real {
                continue;
            }

            // Remove 2-residue segments
            if res[a].ss == b'S'
                && res[a + 1].ss == b'S'
                && res[a - 1].ss != b'S'
                && res[a + 2].ss != b'S'
            {
                res[a].ss = b'L';
                res[a + 1].ss = b'L';
                repeat = true;
            }
            if res[a].ss == b'H'
                && res[a + 1].ss == b'H'
                && res[a - 1].ss != b'H'
                && res[a + 2].ss != b'H'
            {
                res[a].ss = b'L';
                res[a + 1].ss = b'L';
                repeat = true;
            }

            // Remove 1-residue segments
            if res[a].ss == b'S' && res[a - 1].ss != b'S' && res[a + 1].ss != b'S' {
                res[a].ss = b'L';
                repeat = true;
            }
            if res[a].ss == b'H' && res[a - 1].ss != b'H' && res[a + 1].ss != b'H' {
                res[a].ss = b'L';
                repeat = true;
            }

            // Terminal strand residues must have H-bond partner that's also 'S'
            if res[a].ss == b'S' && (res[a - 1].ss != b'S' || res[a + 1].ss != b'S') {
                let mut found = false;

                // Check acceptors
                for &acc in &res[a].acc {
                    if acc < n_res && res[acc].ss == b'S' {
                        found = true;
                        break;
                    }
                }

                // Check donors
                if !found {
                    for &don in &res[a].don {
                        if don < n_res && res[don].ss == b'S' {
                            found = true;
                            break;
                        }
                    }
                }

                // Allow skip residues if neighbor has partner
                if !found
                    && res[a]
                        .flags
                        .intersects(SsFlags::ANTI_STRAND_SKIP | SsFlags::PARA_STRAND_SKIP)
                {
                    if res[a + 1].ss == res[a].ss {
                        for &acc in &res[a + 1].acc {
                            if acc < n_res && res[acc].ss == b'S' {
                                found = true;
                                break;
                            }
                        }
                    }

                    if !found && res[a - 1].ss == res[a].ss {
                        for &don in &res[a - 1].don {
                            if don < n_res && res[don].ss == b'S' {
                                found = true;
                                break;
                            }
                        }
                    }
                }

                if !found {
                    res[a].ss = b'L';
                    repeat = true;
                }
            }
        }
    }
}

// ============================================================================
// Apply Assignments
// ============================================================================

fn apply_ss_assignments(molecule: &mut ObjectMolecule, res: &[ResidueData]) -> usize {
    let mut updated_count = 0;

    let mut ss_map = ahash::AHashMap::new();
    for r in res {
        if !r.real {
            continue;
        }
        let ss = match r.ss {
            b'H' => SecondaryStructure::Helix,
            b'S' => SecondaryStructure::Sheet,
            _ => SecondaryStructure::Loop,
        };
        ss_map.insert((r.chain.clone(), r.resv), ss);
    }

    for atom_idx in 0..molecule.atom_count() {
        let idx = AtomIndex(atom_idx as u32);
        if let Some(atom) = molecule.get_atom(idx) {
            let key = (atom.residue.chain.clone(), atom.residue.resv);
            if let Some(&ss) = ss_map.get(&key) {
                if let Some(atom_mut) = molecule.get_atom_mut(idx) {
                    if atom_mut.ss_type != ss {
                        atom_mut.ss_type = ss;
                        updated_count += 1;
                    }
                }
            }
        }
    }

    updated_count
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dihedral_angle_planar() {
        let v0 = Vec3::new(0.0, 0.0, 0.0);
        let v1 = Vec3::new(1.0, 0.0, 0.0);
        let v2 = Vec3::new(2.0, 0.0, 0.0);
        let v3 = Vec3::new(3.0, 0.0, 0.0);

        let angle = dihedral_angle(v0, v1, v2, v3).to_degrees();
        assert!(angle.abs() < 10.0 || (angle.abs() - 180.0).abs() < 10.0);
    }

    #[test]
    fn test_dihedral_angle_90() {
        let v0 = Vec3::new(0.0, 1.0, 0.0);
        let v1 = Vec3::new(0.0, 0.0, 0.0);
        let v2 = Vec3::new(1.0, 0.0, 0.0);
        let v3 = Vec3::new(1.0, 0.0, 1.0);

        let angle = dihedral_angle(v0, v1, v2, v3).to_degrees();
        assert!((angle.abs() - 90.0).abs() < 1.0);
    }

    #[test]
    fn test_angle_delta_wraparound() {
        assert!((angle_delta(170.0, -170.0) - 20.0).abs() < 0.01);
        assert!((angle_delta(-170.0, 170.0) - 20.0).abs() < 0.01);
        assert!((angle_delta(10.0, -10.0) - 20.0).abs() < 0.01);
    }

    #[test]
    fn test_dss_settings_default() {
        let settings = DssSettings::default();
        assert!((settings.helix_phi_target - (-57.0)).abs() < 0.01);
        assert!((settings.helix_psi_target - (-48.0)).abs() < 0.01);
        assert!((settings.strand_phi_target - (-129.0)).abs() < 0.01);
        assert!((settings.strand_psi_target - 124.0).abs() < 0.01);
    }

    #[test]
    fn test_ss_flags() {
        let mut flags = SsFlags::empty();
        assert!(!flags.intersects(HELIX_HBOND_FLAGS));

        flags |= SsFlags::HELIX_4_HBOND;
        assert!(flags.intersects(HELIX_HBOND_FLAGS));

        flags |= SsFlags::PHI_PSI_HELIX;
        assert!(flags.contains(SsFlags::PHI_PSI_HELIX));
    }

    #[test]
    fn test_hbond_flags() {
        let mut flags = SsFlags::empty();
        flags |= SsFlags::ANTI_STRAND_DOUBLE_HB;
        assert!(flags.contains(SsFlags::ANTI_STRAND_DOUBLE_HB));
        assert!(!flags.contains(SsFlags::PARA_STRAND_DOUBLE_HB));

        flags |= SsFlags::PARA_STRAND_SINGLE_HB;
        assert!(flags.intersects(SsFlags::PARA_STRAND_SINGLE_HB | SsFlags::PARA_STRAND_DOUBLE_HB));
    }
}
