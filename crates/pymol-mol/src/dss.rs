//! Define Secondary Structure (DSS) algorithm
//!
//! Implements PyMOL's DSS algorithm for calculating secondary structure
//! based on backbone phi/psi dihedral angles.
//!
//! This algorithm analyzes the backbone geometry of protein residues and assigns
//! secondary structure types (helix, sheet, loop) based on the dihedral angles.
//!
//! # References
//!
//! - PyMOL's layer3/Selector.cpp - SelectorAssignSS function
//! - PyMOL's layer2/ObjectMolecule.cpp - ObjectMoleculeGetPhiPsi function

use lin_alg::f32::Vec3;
use std::f32::consts::PI;

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
        // PyMOL default values from layer1/SettingInfo.h
        Self {
            helix_phi_target: -57.0,
            helix_psi_target: -48.0,
            helix_phi_include: 55.0,
            helix_psi_include: 55.0,
            helix_phi_exclude: 85.0,
            helix_psi_exclude: 85.0,
            strand_phi_target: -124.0,
            strand_psi_target: 124.0,
            strand_phi_include: 40.0,
            strand_psi_include: 40.0,
            strand_phi_exclude: 90.0,
            strand_psi_exclude: 90.0,
        }
    }
}

// ============================================================================
// Dihedral Angle Calculation
// ============================================================================

/// Calculate the dihedral angle between four points
///
/// The dihedral angle is the angle between the planes defined by
/// (v0, v1, v2) and (v1, v2, v3).
///
/// Based on PyMOL's get_dihedral3f function from layer0/Vector.cpp
///
/// # Returns
/// Dihedral angle in radians, range [-PI, PI]
pub fn dihedral_angle(v0: Vec3, v1: Vec3, v2: Vec3, v3: Vec3) -> f32 {
    // Vectors along the bonds
    let d21 = v2 - v1;
    let d01 = v0 - v1;
    let d32 = v3 - v2;
    
    let d21_len = d21.magnitude();
    if d21_len < 1e-6 {
        // Degenerate case: fall back to angle between vectors
        return angle_between(d01, d32);
    }
    
    // Normal vectors to the two planes
    let dd1 = d21.cross(d01);
    let dd3 = d21.cross(d32);
    
    let dd1_len = dd1.magnitude();
    let dd3_len = dd3.magnitude();
    
    if dd1_len < 1e-6 || dd3_len < 1e-6 {
        // Degenerate case: fall back to angle between vectors
        return angle_between(d01, d32);
    }
    
    // Angle between the normal vectors
    let mut result = angle_between(dd1, dd3);
    
    // Determine sign of the dihedral angle
    let pos_d = d21.cross(dd1);
    if dd3.dot(pos_d) < 0.0 {
        result = -result;
    }
    
    result
}

/// Calculate the angle between two vectors in radians
fn angle_between(v1: Vec3, v2: Vec3) -> f32 {
    let len1 = v1.magnitude();
    let len2 = v2.magnitude();
    
    if len1 < 1e-6 || len2 < 1e-6 {
        return 0.0;
    }
    
    let cos_angle = v1.dot(v2) / (len1 * len2);
    // Clamp to [-1, 1] to avoid NaN from acos
    cos_angle.clamp(-1.0, 1.0).acos()
}

/// Convert radians to degrees
#[inline]
pub fn rad_to_deg(rad: f32) -> f32 {
    rad * 180.0 / PI
}

/// Convert degrees to radians
#[inline]
pub fn deg_to_rad(deg: f32) -> f32 {
    deg * PI / 180.0
}

// ============================================================================
// Phi/Psi Calculation
// ============================================================================

/// Phi and psi angles for a residue
#[derive(Debug, Clone, Copy)]
pub struct PhiPsi {
    /// Phi angle in degrees (C-1, N, CA, C)
    pub phi: f32,
    /// Psi angle in degrees (N, CA, C, N+1)
    pub psi: f32,
}

/// Internal structure for residue data during DSS calculation
#[derive(Debug, Clone)]
struct ResidueData {
    /// Atom indices for backbone atoms
    ca_idx: Option<AtomIndex>,
    n_idx: Option<AtomIndex>,
    c_idx: Option<AtomIndex>,
    /// Chain ID
    chain: String,
    /// Residue number
    resv: i32,
    /// Phi/psi angles (if calculable)
    phi_psi: Option<PhiPsi>,
    /// Classification flags
    flags: u32,
}

// Classification flags (matching PyMOL's cSS* flags)
const SS_GOT_PHI_PSI: u32 = 0x0001;
const SS_PHI_PSI_HELIX: u32 = 0x0002;
const SS_PHI_PSI_NOT_HELIX: u32 = 0x0004;
const SS_PHI_PSI_STRAND: u32 = 0x0008;
const SS_PHI_PSI_NOT_STRAND: u32 = 0x0010;
const SS_HELIX: u32 = 0x0020;
const SS_STRAND: u32 = 0x0040;

// Break size for helix/strand detection (PyMOL's cSSBreakSize = 3)
const SS_BREAK_SIZE: usize = 3;

/// Calculate phi/psi angles for a residue given the backbone atom positions
///
/// # Arguments
/// - `c_prev`: Position of C from previous residue
/// - `n`: Position of N from current residue
/// - `ca`: Position of CA from current residue
/// - `c`: Position of C from current residue
/// - `n_next`: Position of N from next residue
///
/// # Returns
/// Phi and psi angles in degrees, or None if calculation failed
pub fn calculate_phi_psi(
    c_prev: Vec3,
    n: Vec3,
    ca: Vec3,
    c: Vec3,
    n_next: Vec3,
) -> PhiPsi {
    // Phi: dihedral(C-1, N, CA, C)
    let phi = rad_to_deg(dihedral_angle(c_prev, n, ca, c));
    
    // Psi: dihedral(N, CA, C, N+1)
    let psi = rad_to_deg(dihedral_angle(n, ca, c, n_next));
    
    PhiPsi { phi, psi }
}

// ============================================================================
// DSS Algorithm
// ============================================================================

/// Assign secondary structure to a molecule using the DSS algorithm
///
/// This implements PyMOL's DSS (Define Secondary Structure) algorithm which:
/// 1. Calculates phi/psi angles for each residue
/// 2. Classifies residues based on angle deviation from ideal values
/// 3. Assigns secondary structure based on consecutive residue patterns
///
/// # Arguments
/// - `molecule`: The molecule to process
/// - `state`: The state (coordinate set) to use (0 = first state)
/// - `settings`: DSS settings for angle thresholds
///
/// # Returns
/// Number of atoms that had their secondary structure updated
pub fn assign_secondary_structure(
    molecule: &mut ObjectMolecule,
    state: usize,
    settings: &DssSettings,
) -> usize {
    // Get the coordinate set
    let coord_set = match molecule.get_coord_set(state) {
        Some(cs) => cs,
        None => return 0,
    };
    
    // Collect residue data
    let mut residue_data = collect_residue_data(molecule, coord_set);
    
    if residue_data.len() < 3 {
        return 0; // Need at least 3 residues
    }
    
    // Calculate phi/psi angles for each residue
    calculate_all_phi_psi(&mut residue_data, molecule, coord_set);
    
    // Classify residues based on phi/psi angles
    classify_residues(&mut residue_data, settings);
    
    // Assign secondary structure based on consecutive patterns
    assign_ss_from_classification(&mut residue_data);
    
    // Apply the assignments to the molecule
    apply_ss_assignments(molecule, &residue_data)
}

/// Collect backbone atom information for all protein residues
fn collect_residue_data(molecule: &ObjectMolecule, _coord_set: &CoordSet) -> Vec<ResidueData> {
    let mut residue_data = Vec::new();
    
    for residue in molecule.residues() {
        // Skip non-amino acids
        if !crate::residue::is_amino_acid(&residue.key.resn) {
            continue;
        }
        
        // Find backbone atoms
        let ca_idx = find_atom_by_name(&residue, "CA");
        let n_idx = find_atom_by_name(&residue, "N");
        let c_idx = find_atom_by_name(&residue, "C");
        
        // Only include if we have at least CA
        if ca_idx.is_some() {
            residue_data.push(ResidueData {
                ca_idx,
                n_idx,
                c_idx,
                chain: residue.key.chain.clone(),
                resv: residue.key.resv,
                phi_psi: None,
                flags: 0,
            });
        }
    }
    
    residue_data
}

/// Find an atom by name within a residue
fn find_atom_by_name(residue: &ResidueView, name: &str) -> Option<AtomIndex> {
    for (idx, atom) in residue.iter_indexed() {
        if atom.name == name {
            return Some(idx);
        }
    }
    None
}

/// Calculate phi/psi angles for all residues
fn calculate_all_phi_psi(
    residue_data: &mut [ResidueData],
    _molecule: &ObjectMolecule,
    coord_set: &CoordSet,
) {
    let n_res = residue_data.len();
    
    for i in 1..n_res.saturating_sub(1) {
        let prev = &residue_data[i - 1];
        let curr = &residue_data[i];
        let next = &residue_data[i + 1];
        
        // Check chain continuity
        if prev.chain != curr.chain || curr.chain != next.chain {
            continue;
        }
        
        // Check residue number continuity (allow gaps up to 1)
        if curr.resv - prev.resv > 1 || next.resv - curr.resv > 1 {
            continue;
        }
        
        // Get atom positions
        let c_prev = get_atom_pos(prev.c_idx, coord_set);
        let n_curr = get_atom_pos(curr.n_idx, coord_set);
        let ca_curr = get_atom_pos(curr.ca_idx, coord_set);
        let c_curr = get_atom_pos(curr.c_idx, coord_set);
        let n_next = get_atom_pos(next.n_idx, coord_set);
        
        // Need all positions for phi/psi calculation
        if let (Some(c_prev), Some(n_curr), Some(ca_curr), Some(c_curr), Some(n_next)) =
            (c_prev, n_curr, ca_curr, c_curr, n_next)
        {
            let phi_psi = calculate_phi_psi(c_prev, n_curr, ca_curr, c_curr, n_next);
            residue_data[i].phi_psi = Some(phi_psi);
            residue_data[i].flags |= SS_GOT_PHI_PSI;
        }
    }
}

/// Get atom position from coordinate set
fn get_atom_pos(atom_idx: Option<AtomIndex>, coord_set: &CoordSet) -> Option<Vec3> {
    atom_idx.and_then(|idx| coord_set.get_atom_coord(idx))
}

/// Classify residues based on phi/psi angles
fn classify_residues(residue_data: &mut [ResidueData], settings: &DssSettings) {
    for res in residue_data.iter_mut() {
        if let Some(phi_psi) = res.phi_psi {
            // Calculate angle deltas (handling wraparound)
            let helix_phi_delta = angle_delta(phi_psi.phi, settings.helix_phi_target);
            let helix_psi_delta = angle_delta(phi_psi.psi, settings.helix_psi_target);
            let strand_phi_delta = angle_delta(phi_psi.phi, settings.strand_phi_target);
            let strand_psi_delta = angle_delta(phi_psi.psi, settings.strand_psi_target);
            
            // Helix classification
            if helix_phi_delta > settings.helix_phi_exclude
                || helix_psi_delta > settings.helix_psi_exclude
            {
                res.flags |= SS_PHI_PSI_NOT_HELIX;
            } else if helix_phi_delta < settings.helix_phi_include
                && helix_psi_delta < settings.helix_psi_include
            {
                res.flags |= SS_PHI_PSI_HELIX;
            }
            
            // Strand classification
            if strand_phi_delta > settings.strand_phi_exclude
                || strand_psi_delta > settings.strand_psi_exclude
            {
                res.flags |= SS_PHI_PSI_NOT_STRAND;
            } else if strand_phi_delta < settings.strand_phi_include
                && strand_psi_delta < settings.strand_psi_include
            {
                res.flags |= SS_PHI_PSI_STRAND;
            }
        }
    }
}

/// Calculate the angular difference, handling wraparound at 180°
fn angle_delta(angle: f32, target: f32) -> f32 {
    let mut delta = (angle - target).abs();
    if delta > 180.0 {
        delta = 360.0 - delta;
    }
    delta
}

/// Assign secondary structure based on consecutive residue patterns
///
/// PyMOL requires at least 3 consecutive residues with helix/strand characteristics
/// to assign secondary structure
fn assign_ss_from_classification(residue_data: &mut [ResidueData]) {
    let n_res = residue_data.len();
    
    // Skip if too few residues
    if n_res < 2 * SS_BREAK_SIZE + 1 {
        return;
    }
    
    // First pass: tentatively assign helices
    // Look for runs of residues with helix-compatible phi/psi
    for i in SS_BREAK_SIZE..n_res.saturating_sub(SS_BREAK_SIZE) {
        let flags = residue_data[i].flags;
        
        // Check if this residue is a potential helix residue
        if (flags & SS_PHI_PSI_HELIX) != 0 {
            // Check if we have a run of helix-compatible residues
            let mut helix_run = true;
            for j in (i.saturating_sub(1))..=(i + 1).min(n_res - 1) {
                if j != i {
                    let neighbor_flags = residue_data[j].flags;
                    // Allow if neighbor is helix-compatible OR not explicitly non-helix
                    if (neighbor_flags & SS_PHI_PSI_NOT_HELIX) != 0 
                        && (neighbor_flags & SS_PHI_PSI_HELIX) == 0 
                    {
                        helix_run = false;
                        break;
                    }
                }
            }
            
            if helix_run {
                residue_data[i].flags |= SS_HELIX;
            }
        }
    }
    
    // Second pass: tentatively assign strands
    for i in SS_BREAK_SIZE..n_res.saturating_sub(SS_BREAK_SIZE) {
        let flags = residue_data[i].flags;
        
        // Don't overwrite helix assignments
        if (flags & SS_HELIX) != 0 {
            continue;
        }
        
        if (flags & SS_PHI_PSI_STRAND) != 0 {
            let mut strand_run = true;
            for j in (i.saturating_sub(1))..=(i + 1).min(n_res - 1) {
                if j != i {
                    let neighbor_flags = residue_data[j].flags;
                    if (neighbor_flags & SS_PHI_PSI_NOT_STRAND) != 0
                        && (neighbor_flags & SS_PHI_PSI_STRAND) == 0
                    {
                        strand_run = false;
                        break;
                    }
                }
            }
            
            if strand_run {
                residue_data[i].flags |= SS_STRAND;
            }
        }
    }
}

/// Apply secondary structure assignments to the molecule
fn apply_ss_assignments(molecule: &mut ObjectMolecule, residue_data: &[ResidueData]) -> usize {
    let mut updated_count = 0;
    
    // Build a map from residue (chain, resv) to secondary structure
    let mut ss_map = std::collections::HashMap::new();
    for res in residue_data {
        let ss = if (res.flags & SS_HELIX) != 0 {
            SecondaryStructure::Helix
        } else if (res.flags & SS_STRAND) != 0 {
            SecondaryStructure::Sheet
        } else {
            SecondaryStructure::Loop
        };
        ss_map.insert((res.chain.clone(), res.resv), ss);
    }
    
    // Apply to all atoms in matching residues
    for atom_idx in 0..molecule.atom_count() {
        let idx = AtomIndex(atom_idx as u32);
        if let Some(atom) = molecule.get_atom(idx) {
            let key = (atom.chain.clone(), atom.resv);
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
        // Four points in a plane, 180° dihedral
        let v0 = Vec3::new(0.0, 0.0, 0.0);
        let v1 = Vec3::new(1.0, 0.0, 0.0);
        let v2 = Vec3::new(2.0, 0.0, 0.0);
        let v3 = Vec3::new(3.0, 0.0, 0.0);
        
        let angle = rad_to_deg(dihedral_angle(v0, v1, v2, v3));
        // Collinear points should give 0 (degenerate case)
        assert!(angle.abs() < 10.0 || (angle.abs() - 180.0).abs() < 10.0);
    }

    #[test]
    fn test_dihedral_angle_90() {
        // Test a 90° dihedral
        let v0 = Vec3::new(0.0, 1.0, 0.0);
        let v1 = Vec3::new(0.0, 0.0, 0.0);
        let v2 = Vec3::new(1.0, 0.0, 0.0);
        let v3 = Vec3::new(1.0, 0.0, 1.0);
        
        let angle = rad_to_deg(dihedral_angle(v0, v1, v2, v3));
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
    }
}
