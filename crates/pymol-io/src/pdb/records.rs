//! PDB record types
//!
//! Defines the data structures for PDB records.

use pymol_mol::Element;

/// Parsed ATOM or HETATM record
#[derive(Debug, Clone)]
pub struct AtomRecord {
    /// Record type: true for HETATM, false for ATOM
    pub hetatm: bool,
    /// Atom serial number (1-99999)
    pub serial: i32,
    /// Atom name (4 characters, may include leading space)
    pub name: String,
    /// Alternate location indicator
    pub alt_loc: char,
    /// Residue name (3 characters)
    pub resn: String,
    /// Chain identifier
    pub chain: String,
    /// Residue sequence number
    pub resv: i32,
    /// Insertion code
    pub icode: char,
    /// X coordinate (Angstroms)
    pub x: f32,
    /// Y coordinate (Angstroms)
    pub y: f32,
    /// Z coordinate (Angstroms)
    pub z: f32,
    /// Occupancy
    pub occupancy: f32,
    /// Temperature factor (B-factor)
    pub b_factor: f32,
    /// Segment identifier
    pub segi: String,
    /// Element symbol
    pub element: String,
    /// Formal charge
    pub charge: String,
}

impl Default for AtomRecord {
    fn default() -> Self {
        AtomRecord {
            hetatm: false,
            serial: 0,
            name: String::new(),
            alt_loc: ' ',
            resn: String::new(),
            chain: String::new(),
            resv: 0,
            icode: ' ',
            x: 0.0,
            y: 0.0,
            z: 0.0,
            occupancy: 1.0,
            b_factor: 0.0,
            segi: String::new(),
            element: String::new(),
            charge: String::new(),
        }
    }
}

impl AtomRecord {
    /// Get the element from the record
    pub fn get_element(&self) -> Element {
        // Try to parse element from the element column first
        if !self.element.is_empty() {
            if let Some(elem) = Element::from_symbol(self.element.trim()) {
                return elem;
            }
        }

        // Fall back to inferring from atom name
        infer_element_from_name(&self.name)
    }

    /// Parse the formal charge from the charge column
    pub fn get_formal_charge(&self) -> i8 {
        parse_pdb_charge(&self.charge)
    }
}

/// Infer element from atom name
pub fn infer_element_from_name(name: &str) -> Element {
    let name = name.trim();
    if name.is_empty() {
        return Element::Unknown;
    }

    // Common atom names that map directly to elements
    // First, try the first character if the name starts with a space
    let chars: Vec<char> = name.chars().collect();

    // PDB atom names are justified based on element:
    // - 1-letter elements: right-justified in columns 13-14, name starts in col 13
    // - 2-letter elements: left-justified in columns 13-14

    // If first char is space or digit, look at second char
    let start_idx = if chars[0] == ' ' || chars[0].is_ascii_digit() {
        1
    } else {
        0
    };

    if start_idx >= chars.len() {
        return Element::Unknown;
    }

    // Try two-letter element first
    if chars.len() > start_idx + 1 {
        let two_letter: String = chars[start_idx..=start_idx + 1].iter().collect();
        // Only try if second char is lowercase (true 2-letter element)
        // or both are uppercase letters that form a valid element
        if let Some(elem) = Element::from_symbol(&two_letter) {
            // For common atoms like CA (calcium vs C-alpha), prefer single letter
            // if the atom name looks like a protein atom name
            let is_common_protein_name = matches!(
                name.trim(),
                "CA" | "CB" | "CG" | "CD" | "CE" | "CZ" | "CH" | "NE" | "NH" | "NZ" | "OG" | "OH"
                    | "OE" | "OD" | "SD" | "SG"
            );
            if !is_common_protein_name {
                return elem;
            }
        }
    }

    // Try single-letter element
    let one_letter = chars[start_idx].to_ascii_uppercase().to_string();
    Element::from_symbol(&one_letter).unwrap_or(Element::Unknown)
}

/// Parse PDB-style charge string (e.g., "2+", "1-", "+", "-")
pub fn parse_pdb_charge(charge_str: &str) -> i8 {
    let s = charge_str.trim();
    if s.is_empty() {
        return 0;
    }

    // Format can be "2+", "2-", "+2", "-2", "+", "-"
    let chars: Vec<char> = s.chars().collect();

    match chars.as_slice() {
        ['+'] => 1,
        ['-'] => -1,
        [d, '+'] if d.is_ascii_digit() => d.to_digit(10).unwrap_or(0) as i8,
        [d, '-'] if d.is_ascii_digit() => -(d.to_digit(10).unwrap_or(0) as i8),
        ['+', d] if d.is_ascii_digit() => d.to_digit(10).unwrap_or(0) as i8,
        ['-', d] if d.is_ascii_digit() => -(d.to_digit(10).unwrap_or(0) as i8),
        _ => 0,
    }
}

/// CONECT record for bond connectivity
#[derive(Debug, Clone, Default)]
pub struct ConectRecord {
    /// Serial number of the central atom
    pub atom: i32,
    /// Serial numbers of bonded atoms (up to 4 per record)
    pub bonded: Vec<i32>,
}

/// CRYST1 record for unit cell parameters
#[derive(Debug, Clone)]
pub struct Cryst1Record {
    /// Unit cell a dimension (Angstroms)
    pub a: f32,
    /// Unit cell b dimension (Angstroms)
    pub b: f32,
    /// Unit cell c dimension (Angstroms)
    pub c: f32,
    /// Unit cell alpha angle (degrees)
    pub alpha: f32,
    /// Unit cell beta angle (degrees)
    pub beta: f32,
    /// Unit cell gamma angle (degrees)
    pub gamma: f32,
    /// Space group
    pub space_group: String,
    /// Z value
    pub z: i32,
}

impl Default for Cryst1Record {
    fn default() -> Self {
        Cryst1Record {
            a: 1.0,
            b: 1.0,
            c: 1.0,
            alpha: 90.0,
            beta: 90.0,
            gamma: 90.0,
            space_group: "P 1".to_string(),
            z: 1,
        }
    }
}

/// HELIX record for helix secondary structure
#[derive(Debug, Clone)]
pub struct HelixRecord {
    /// Helix serial number
    pub serial: i32,
    /// Helix identifier
    pub helix_id: String,
    /// Initial residue name
    pub init_resn: String,
    /// Initial chain ID
    pub init_chain: String,
    /// Initial residue sequence number
    pub init_seq: i32,
    /// Initial insertion code
    pub init_icode: char,
    /// Terminal residue name
    pub end_resn: String,
    /// Terminal chain ID
    pub end_chain: String,
    /// Terminal residue sequence number
    pub end_seq: i32,
    /// Terminal insertion code
    pub end_icode: char,
    /// Helix class (1-10)
    pub helix_class: i32,
}

/// SHEET record for beta sheet secondary structure
#[derive(Debug, Clone)]
pub struct SheetRecord {
    /// Strand number
    pub strand: i32,
    /// Sheet identifier
    pub sheet_id: String,
    /// Number of strands in sheet
    pub num_strands: i32,
    /// Initial residue name
    pub init_resn: String,
    /// Initial chain ID
    pub init_chain: String,
    /// Initial residue sequence number
    pub init_seq: i32,
    /// Initial insertion code
    pub init_icode: char,
    /// Terminal residue name
    pub end_resn: String,
    /// Terminal chain ID
    pub end_chain: String,
    /// Terminal residue sequence number
    pub end_seq: i32,
    /// Terminal insertion code
    pub end_icode: char,
    /// Sense of strand with respect to previous
    pub sense: i32,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_infer_element() {
        assert_eq!(infer_element_from_name(" CA "), Element::Carbon);
        assert_eq!(infer_element_from_name("FE  "), Element::Iron);
        assert_eq!(infer_element_from_name(" N  "), Element::Nitrogen);
        assert_eq!(infer_element_from_name(" O  "), Element::Oxygen);
        assert_eq!(infer_element_from_name("1H  "), Element::Hydrogen);
    }

    #[test]
    fn test_parse_charge() {
        assert_eq!(parse_pdb_charge("2+"), 2);
        assert_eq!(parse_pdb_charge("2-"), -2);
        assert_eq!(parse_pdb_charge("+"), 1);
        assert_eq!(parse_pdb_charge("-"), -1);
        assert_eq!(parse_pdb_charge(""), 0);
    }
}
