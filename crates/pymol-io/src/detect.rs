//! File format detection
//!
//! Detects molecular file formats from file extensions and content.

use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use crate::traits::FileFormat;

/// Detect file format from a file path (by extension)
pub fn detect_from_path(path: &Path) -> FileFormat {
    FileFormat::from_path(path)
}

/// Detect file format from file content
///
/// Reads the first few lines of the file to determine the format.
/// Returns the detected format and the content that was read (for re-parsing).
pub fn detect_from_content<R: Read>(reader: R) -> (FileFormat, String) {
    let mut buf_reader = BufReader::new(reader);
    let mut content = String::new();
    let mut lines = Vec::new();

    // Read up to 10 lines for detection
    for _ in 0..10 {
        let mut line = String::new();
        match buf_reader.read_line(&mut line) {
            Ok(0) => break,
            Ok(_) => {
                content.push_str(&line);
                lines.push(line);
            }
            Err(_) => break,
        }
    }

    // Read the rest into content
    let _ = buf_reader.read_to_string(&mut content);

    let format = detect_from_lines(&lines);
    (format, content)
}

/// Detect file format from the first few lines
fn detect_from_lines(lines: &[String]) -> FileFormat {
    if lines.is_empty() {
        return FileFormat::Unknown;
    }

    let first_line = lines[0].trim();

    // Check for PDB format
    if is_pdb_format(lines) {
        return FileFormat::Pdb;
    }

    // Check for mmCIF format
    if first_line.starts_with("data_") || first_line.starts_with("_") {
        return FileFormat::Cif;
    }

    // Check for MOL2 format
    if first_line.starts_with("@<TRIPOS>") || lines.iter().any(|l| l.starts_with("@<TRIPOS>")) {
        return FileFormat::Mol2;
    }

    // Check for SDF/MOL format (connection table)
    if is_sdf_format(lines) {
        return FileFormat::Sdf;
    }

    // Check for XYZ format
    if is_xyz_format(lines) {
        return FileFormat::Xyz;
    }

    FileFormat::Unknown
}

/// Check if content looks like PDB format
fn is_pdb_format(lines: &[String]) -> bool {
    // PDB files typically start with HEADER, REMARK, ATOM, HETATM, or other keywords
    let pdb_keywords = [
        "HEADER", "TITLE", "COMPND", "SOURCE", "KEYWDS", "EXPDTA", "AUTHOR", "REVDAT", "JRNL",
        "REMARK", "DBREF", "SEQRES", "HET", "HETNAM", "FORMUL", "HELIX", "SHEET", "SITE", "CRYST1",
        "ORIGX", "SCALE", "MTRIX", "MODEL", "ATOM", "HETATM", "TER", "ENDMDL", "CONECT", "MASTER",
        "END",
    ];

    for line in lines {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        for keyword in &pdb_keywords {
            if trimmed.starts_with(keyword) {
                return true;
            }
        }
    }
    false
}

/// Check if content looks like SDF/MOL format
fn is_sdf_format(lines: &[String]) -> bool {
    // SDF/MOL files have a specific structure:
    // Line 1: Molecule name
    // Line 2: Program/timestamp line
    // Line 3: Comment
    // Line 4: Counts line (aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv)
    if lines.len() < 4 {
        return false;
    }

    // Check the counts line (4th line, 0-indexed as 3)
    let counts_line = lines[3].trim();
    if counts_line.len() < 6 {
        return false;
    }

    // The counts line should have numbers in specific positions
    // First 3 chars: atom count, next 3 chars: bond count
    let atom_count = counts_line.get(0..3).and_then(|s| s.trim().parse::<u32>().ok());
    let bond_count = counts_line.get(3..6).and_then(|s| s.trim().parse::<u32>().ok());

    // Also check if line contains "V2000" or "V3000"
    let has_version = counts_line.contains("V2000") || counts_line.contains("V3000");

    atom_count.is_some() && bond_count.is_some() && (has_version || counts_line.len() >= 33)
}

/// Check if content looks like XYZ format
fn is_xyz_format(lines: &[String]) -> bool {
    if lines.len() < 3 {
        return false;
    }

    // First line should be an integer (atom count)
    let first_line = lines[0].trim();
    if first_line.parse::<u32>().is_err() {
        return false;
    }

    // Third line (first atom line) should be: element x y z
    let atom_line = lines[2].trim();
    let parts: Vec<&str> = atom_line.split_whitespace().collect();

    if parts.len() < 4 {
        return false;
    }

    // First part should be an element symbol (1-2 letters)
    let element = parts[0];
    if element.is_empty() || element.len() > 2 {
        return false;
    }

    // Remaining parts should be numbers
    parts[1..4].iter().all(|p| p.parse::<f64>().is_ok())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_detect_pdb() {
        let lines = vec![
            "HEADER    test".to_string(),
            "ATOM      1  N   ALA A   1".to_string(),
        ];
        assert_eq!(detect_from_lines(&lines), FileFormat::Pdb);
    }

    #[test]
    fn test_detect_sdf() {
        let lines = vec![
            "Molecule Name".to_string(),
            "  Program".to_string(),
            "Comment".to_string(),
            "  3  2  0  0  0  0  0  0  0  0999 V2000".to_string(),
        ];
        assert_eq!(detect_from_lines(&lines), FileFormat::Sdf);
    }

    #[test]
    fn test_detect_mol2() {
        let lines = vec!["@<TRIPOS>MOLECULE".to_string(), "test".to_string()];
        assert_eq!(detect_from_lines(&lines), FileFormat::Mol2);
    }

    #[test]
    fn test_detect_xyz() {
        let lines = vec![
            "3".to_string(),
            "Water molecule".to_string(),
            "O  0.0 0.0 0.0".to_string(),
        ];
        assert_eq!(detect_from_lines(&lines), FileFormat::Xyz);
    }

    #[test]
    fn test_detect_cif() {
        let lines = vec!["data_1ABC".to_string(), "_cell.length_a 10.0".to_string()];
        assert_eq!(detect_from_lines(&lines), FileFormat::Cif);
    }
}
