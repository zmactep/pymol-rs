//! Hybrid-36 encoding/decoding for PDB serial numbers and residue sequence numbers.
//!
//! The hybrid-36 system extends the PDB format to handle very large structures
//! whose atom count or residue sequence numbers exceed the decimal field limit:
//!
//! | Field           | Width | Decimal range | Uppercase range  | Lowercase range  | Total        |
//! |-----------------|-------|---------------|------------------|------------------|--------------|
//! | Atom serial     |   5   | 1 – 99 999    | A0000 – ZZZZZ    | a0000 – zzzzz    | 87 440 031   |
//! | Residue seq     |   4   | 1 – 9 999     | A000 – ZZZZ      | a000 – zzzz      |  2 436 111   |
//!
//! Programs that support hybrid-36 remain interoperable with programs that do
//! not, as long as fewer than 100 000 atoms (or 10 000 residues) are present.
//!
//! Reference: <https://cci.lbl.gov/hybrid_36/>

/// Base-36 digit set — uppercase variant ("0-9 A-Z").
const DIGITS_UPPER: &[u8] = b"0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
/// Base-36 digit set — lowercase variant ("0-9 a-z").
const DIGITS_LOWER: &[u8] = b"0123456789abcdefghijklmnopqrstuvwxyz";

/// Decode a hybrid-36 field of `width` characters into an integer.
///
/// * Whitespace-only fields decode to `0`.
/// * Fields whose first non-space character is a digit, `+`, or `-` are
///   parsed as plain decimal integers (standard PDB behaviour).
/// * Fields starting with an uppercase letter are decoded as upper-case
///   base-36 (values ≥ 10^width).
/// * Fields starting with a lowercase letter are decoded as lower-case
///   base-36 (values ≥ 10^width + 26 × 36^(width-1)).
///
/// Returns `None` if any character is outside the expected digit set.
pub fn hy36decode(width: usize, s: &str) -> Option<i32> {
    let s = s.trim();
    if s.is_empty() {
        return Some(0);
    }

    let first = s.as_bytes()[0];

    if first.is_ascii_digit() || first == b'+' || first == b'-' {
        // Plain decimal (standard PDB range)
        s.parse::<i32>().ok()
    } else if first.is_ascii_uppercase() {
        // Upper-case base-36 range
        let mut value: i64 = 0;
        for &b in s.as_bytes() {
            let digit = DIGITS_UPPER.iter().position(|&d| d == b)?;
            value = value * 36 + digit as i64;
        }
        let w = width as u32;
        value -= 10 * 36i64.pow(w - 1);
        value += 10i64.pow(w);
        i32::try_from(value).ok()
    } else if first.is_ascii_lowercase() {
        // Lower-case base-36 range
        let mut value: i64 = 0;
        for &b in s.as_bytes() {
            let digit = DIGITS_LOWER.iter().position(|&d| d == b)?;
            value = value * 36 + digit as i64;
        }
        let w = width as u32;
        let width_pow = 36i64.pow(w - 1);
        value -= 10 * width_pow;
        value += 10i64.pow(w);
        value += 26 * width_pow;
        i32::try_from(value).ok()
    } else {
        None
    }
}

/// Encode `value` as a hybrid-36 string of exactly `width` characters.
///
/// * Negative values and values in the range 0 – 10^width − 1 are returned
///   as right-justified decimal strings (space-padded on the left), matching
///   standard PDB convention.
/// * Values from 10^width up to 10^width + 26 × 36^(width-1) − 1 are encoded
///   as upper-case base-36.
/// * Values from there up to 10^width + 2 × 26 × 36^(width-1) − 1 are encoded
///   as lower-case base-36.
///
/// Returns `None` if `value` exceeds the representable range for `width`.
pub fn hy36encode(width: usize, value: i32) -> Option<String> {
    let w = width as u32;
    let dec_limit = 10i32.pow(w);
    let base36_range = 26 * 36i32.pow(w - 1);

    if value < dec_limit {
        // Plain decimal — space-padded right-justified (standard PDB style).
        // Negative sequence numbers are valid PDB and handled here too.
        return Some(format!("{:>width$}", value, width = width));
    }

    let mut v = value - dec_limit;

    if v < base36_range {
        // Upper-case base-36
        let mut result = vec![0u8; width];
        let mut temp = v;
        for i in (1..width).rev() {
            result[i] = DIGITS_UPPER[(temp % 36) as usize];
            temp /= 36;
        }
        result[0] = DIGITS_UPPER[(10 + temp) as usize];
        return String::from_utf8(result).ok();
    }

    v -= base36_range;

    if v < base36_range {
        // Lower-case base-36
        let mut result = vec![0u8; width];
        let mut temp = v;
        for i in (1..width).rev() {
            result[i] = DIGITS_LOWER[(temp % 36) as usize];
            temp /= 36;
        }
        result[0] = DIGITS_LOWER[(10 + temp) as usize];
        return String::from_utf8(result).ok();
    }

    None // Value exceeds hybrid-36 capacity for this width
}

#[cfg(test)]
mod tests {
    use super::*;

    // --- decode ---

    #[test]
    fn test_decode_decimal() {
        assert_eq!(hy36decode(5, "    1"), Some(1));
        assert_eq!(hy36decode(5, "99999"), Some(99999));
        assert_eq!(hy36decode(4, "   1"), Some(1));
        assert_eq!(hy36decode(4, "9999"), Some(9999));
        assert_eq!(hy36decode(5, "     "), Some(0));
        assert_eq!(hy36decode(4, "    "), Some(0));
    }

    #[test]
    fn test_decode_negative_resv() {
        // Negative residue sequence numbers are valid PDB
        assert_eq!(hy36decode(4, "  -1"), Some(-1));
        assert_eq!(hy36decode(4, "-999"), Some(-999));
    }

    #[test]
    fn test_decode_uppercase() {
        // A0000 is the first uppercase value for width=5: 100000
        assert_eq!(hy36decode(5, "A0000"), Some(100000));
        // A0001 → 100001
        assert_eq!(hy36decode(5, "A0001"), Some(100001));
        // A000 is the first uppercase value for width=4: 10000
        assert_eq!(hy36decode(4, "A000"), Some(10000));
        assert_eq!(hy36decode(4, "A001"), Some(10001));
    }

    #[test]
    fn test_decode_lowercase() {
        // a0000 is the first lowercase value for width=5: 100000 + 26*36^4 = 43770016
        let first_lower_5 = 100000 + 26 * 36i32.pow(4);
        assert_eq!(hy36decode(5, "a0000"), Some(first_lower_5));
        // a000 for width=4: 10000 + 26*36^3 = 1223056
        let first_lower_4 = 10000 + 26 * 36i32.pow(3);
        assert_eq!(hy36decode(4, "a000"), Some(first_lower_4));
    }

    #[test]
    fn test_decode_invalid() {
        assert_eq!(hy36decode(5, "!0000"), None);
        assert_eq!(hy36decode(5, "A00!0"), None);
    }

    // --- encode ---

    #[test]
    fn test_encode_decimal() {
        assert_eq!(hy36encode(5, 1), Some("    1".to_string()));
        assert_eq!(hy36encode(5, 99999), Some("99999".to_string()));
        assert_eq!(hy36encode(4, 1), Some("   1".to_string()));
        assert_eq!(hy36encode(4, 9999), Some("9999".to_string()));
        assert_eq!(hy36encode(4, 0), Some("   0".to_string()));
    }

    #[test]
    fn test_encode_negative() {
        assert_eq!(hy36encode(4, -1), Some("  -1".to_string()));
        assert_eq!(hy36encode(4, -999), Some("-999".to_string()));
    }

    #[test]
    fn test_encode_uppercase() {
        assert_eq!(hy36encode(5, 100000), Some("A0000".to_string()));
        assert_eq!(hy36encode(5, 100001), Some("A0001".to_string()));
        assert_eq!(hy36encode(4, 10000), Some("A000".to_string()));
        assert_eq!(hy36encode(4, 10001), Some("A001".to_string()));
    }

    #[test]
    fn test_encode_lowercase() {
        let first_lower_5 = 100000 + 26 * 36i32.pow(4);
        assert_eq!(hy36encode(5, first_lower_5), Some("a0000".to_string()));
        let first_lower_4 = 10000 + 26 * 36i32.pow(3);
        assert_eq!(hy36encode(4, first_lower_4), Some("a000".to_string()));
    }

    #[test]
    fn test_encode_out_of_range() {
        let max_5 = 100000 + 2 * 26 * 36i32.pow(4) - 1;
        assert!(hy36encode(5, max_5).is_some());
        assert!(hy36encode(5, max_5 + 1).is_none());
    }

    // --- round-trip ---

    #[test]
    fn test_roundtrip_serial() {
        for v in [1, 99999, 100000, 100001, 43770015, 43770016, 87440031] {
            let encoded = hy36encode(5, v).unwrap_or_else(|| panic!("encode failed for {}", v));
            assert_eq!(encoded.len(), 5, "encoded '{}' is not 5 chars", encoded);
            let decoded = hy36decode(5, &encoded).unwrap_or_else(|| panic!("decode failed for '{}'", encoded));
            assert_eq!(decoded, v, "round-trip failed for {}", v);
        }
    }

    #[test]
    fn test_roundtrip_resv() {
        for v in [1, 9999, 10000, 10001, 1223055, 1223056, 2436111] {
            let encoded = hy36encode(4, v).unwrap_or_else(|| panic!("encode failed for {}", v));
            assert_eq!(encoded.len(), 4, "encoded '{}' is not 4 chars", encoded);
            let decoded = hy36decode(4, &encoded).unwrap_or_else(|| panic!("decode failed for '{}'", encoded));
            assert_eq!(decoded, v, "round-trip failed for {}", v);
        }
    }
}
