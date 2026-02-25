//! Space group symmetry operation lookup and parsing
//!
//! Maps space group names to lists of 4×4 transformation matrices (in fractional
//! coordinates). The symmetry operation data is ported from PyMOL's `xray.py`.

use super::space_groups_data::{ALIASES, SYM_OPS};

/// Get the symmetry operation matrices for a space group.
///
/// Each matrix is a 4×4 row-major array representing a symmetry operation
/// in fractional coordinates. Returns `None` if the space group is unknown.
pub fn get_symops(space_group: &str) -> Option<Vec<[f32; 16]>> {
    let canonical = canonicalize(space_group);
    let ops = SYM_OPS.get(canonical.as_str())?;
    Some(ops.iter().map(|op| parse_symop(op)).collect())
}

/// Canonicalize a space group name: normalize whitespace, uppercase, resolve aliases
pub fn canonicalize(sg: &str) -> String {
    // Normalize whitespace and uppercase
    let normalized: String = sg
        .split_whitespace()
        .collect::<Vec<_>>()
        .join(" ")
        .to_uppercase();

    // Check alias table
    if let Some(canonical) = ALIASES.get(normalized.as_str()) {
        return canonical.to_string();
    }

    // Try with spaces stripped (e.g., "P212121" -> check alias)
    let no_spaces: String = normalized.chars().filter(|c| *c != ' ').collect();
    if no_spaces != normalized {
        if let Some(canonical) = ALIASES.get(no_spaces.as_str()) {
            return canonical.to_string();
        }
    }

    normalized
}

/// Parse a symmetry operation string like `"-x+1/2,y,-z"` into a 4×4 row-major matrix
///
/// Each comma-separated component defines one row of the 3×3 rotation part
/// plus the translation. For example, `-x+1/2` means coefficient -1 for x,
/// 0 for y, 0 for z, and translation 1/2.
pub fn parse_symop(op: &str) -> [f32; 16] {
    let mut mat = [0.0f32; 16];
    mat[15] = 1.0;

    for (row, expr) in op.split(',').enumerate() {
        if row >= 3 {
            break;
        }
        let (coeffs, trans) = parse_expr(expr.trim());
        mat[row * 4] = coeffs[0]; // x coefficient
        mat[row * 4 + 1] = coeffs[1]; // y coefficient
        mat[row * 4 + 2] = coeffs[2]; // z coefficient
        mat[row * 4 + 3] = trans; // translation
    }

    mat
}

/// Parse a single symmetry expression component like `-x+1/2` or `y-z`
///
/// Returns ([coeff_x, coeff_y, coeff_z], translation)
fn parse_expr(expr: &str) -> ([f32; 3], f32) {
    let mut coeffs = [0.0f32; 3];
    let mut trans = 0.0f32;

    let chars: Vec<char> = expr.chars().collect();
    let len = chars.len();
    let mut i = 0;

    while i < len {
        let c = chars[i];

        // Determine sign
        let sign = if c == '-' {
            i += 1;
            -1.0
        } else if c == '+' {
            i += 1;
            1.0
        } else {
            1.0
        };

        if i >= len {
            break;
        }

        let c = chars[i];

        match c {
            'x' | 'X' => {
                coeffs[0] = sign;
                i += 1;
            }
            'y' | 'Y' => {
                coeffs[1] = sign;
                i += 1;
            }
            'z' | 'Z' => {
                coeffs[2] = sign;
                i += 1;
            }
            '0'..='9' => {
                // Parse a number, possibly a fraction like "1/2" or "1/3"
                let (val, consumed) = parse_number(&chars[i..]);
                trans += sign * val;
                i += consumed;
            }
            ' ' => {
                i += 1;
            }
            _ => {
                i += 1;
            }
        }
    }

    (coeffs, trans)
}

/// Parse a number from a character slice. Handles integers and fractions like "1/2".
/// Returns (value, chars_consumed).
fn parse_number(chars: &[char]) -> (f32, usize) {
    let mut i = 0;
    let mut numerator = 0.0f32;

    // Parse integer part
    while i < chars.len() && chars[i].is_ascii_digit() {
        numerator = numerator * 10.0 + (chars[i] as u32 - '0' as u32) as f32;
        i += 1;
    }

    // Check for fraction
    if i < chars.len() && chars[i] == '/' {
        i += 1;
        let mut denominator = 0.0f32;
        while i < chars.len() && chars[i].is_ascii_digit() {
            denominator = denominator * 10.0 + (chars[i] as u32 - '0' as u32) as f32;
            i += 1;
        }
        if denominator != 0.0 {
            return (numerator / denominator, i);
        }
    }

    (numerator, i)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_identity() {
        let mat = parse_symop("x,y,z");
        // Should be identity
        assert_eq!(mat[0], 1.0); // x -> x
        assert_eq!(mat[5], 1.0); // y -> y
        assert_eq!(mat[10], 1.0); // z -> z
        // No translation
        assert_eq!(mat[3], 0.0);
        assert_eq!(mat[7], 0.0);
        assert_eq!(mat[11], 0.0);
    }

    #[test]
    fn test_parse_negative_with_translation() {
        let mat = parse_symop("-x+1/2,-y+1/2,z");
        assert_eq!(mat[0], -1.0); // -x
        assert_eq!(mat[3], 0.5); // +1/2
        assert_eq!(mat[5], -1.0); // -y
        assert_eq!(mat[7], 0.5); // +1/2
        assert_eq!(mat[10], 1.0); // z
        assert_eq!(mat[11], 0.0);
    }

    #[test]
    fn test_parse_mixed() {
        let mat = parse_symop("-x,y+1/2,-z+1/2");
        assert_eq!(mat[0], -1.0);
        assert_eq!(mat[3], 0.0);
        assert_eq!(mat[5], 1.0);
        assert_eq!(mat[7], 0.5);
        assert_eq!(mat[10], -1.0);
        assert_eq!(mat[11], 0.5);
    }

    #[test]
    fn test_parse_third_fraction() {
        let mat = parse_symop("x+1/3,y+2/3,z+2/3");
        assert_eq!(mat[0], 1.0);
        assert!((mat[3] - 1.0 / 3.0).abs() < 1e-6);
        assert_eq!(mat[5], 1.0);
        assert!((mat[7] - 2.0 / 3.0).abs() < 1e-6);
        assert_eq!(mat[10], 1.0);
        assert!((mat[11] - 2.0 / 3.0).abs() < 1e-6);
    }

    #[test]
    fn test_get_symops_p1() {
        let ops = get_symops("P 1").unwrap();
        assert_eq!(ops.len(), 1);
        // Identity
        assert_eq!(ops[0][0], 1.0);
        assert_eq!(ops[0][5], 1.0);
        assert_eq!(ops[0][10], 1.0);
    }

    #[test]
    fn test_get_symops_p212121() {
        let ops = get_symops("P 21 21 21").unwrap();
        assert_eq!(ops.len(), 4);
    }

    #[test]
    fn test_get_symops_alias() {
        // "P212121" should resolve to "P 21 21 21"
        let ops = get_symops("P212121").unwrap();
        assert_eq!(ops.len(), 4);
    }

    #[test]
    fn test_get_symops_case_insensitive() {
        let ops = get_symops("p 21 21 21").unwrap();
        assert_eq!(ops.len(), 4);
    }

    #[test]
    fn test_get_symops_unknown() {
        assert!(get_symops("NONEXISTENT").is_none());
    }

    #[test]
    fn test_canonicalize() {
        assert_eq!(canonicalize("P212121"), "P 21 21 21");
        assert_eq!(canonicalize("p 1"), "P 1");
        assert_eq!(canonicalize("  P  21  21  21  "), "P 21 21 21");
    }

    #[test]
    fn test_p_minus_1() {
        let ops = get_symops("P -1").unwrap();
        assert_eq!(ops.len(), 2);
        // Second op should be -x,-y,-z (inversion)
        assert_eq!(ops[1][0], -1.0);
        assert_eq!(ops[1][5], -1.0);
        assert_eq!(ops[1][10], -1.0);
    }
}
