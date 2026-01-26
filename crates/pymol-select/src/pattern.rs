//! Pattern matching for selection values
//!
//! Provides pattern types for matching atom properties like names, residues,
//! chains, etc. Supports exact matches, wildcards, lists, and ranges.

/// A pattern for matching string values
///
/// Used for atom names, residue names, chain IDs, element symbols, etc.
#[derive(Debug, Clone, PartialEq)]
pub enum Pattern {
    /// Exact match (case-insensitive for most properties)
    Exact(String),

    /// Wildcard pattern using * for any characters
    /// - `C*` matches CA, CB, CD, etc.
    /// - `*H` matches any name ending in H
    /// - `*H*` matches any name containing H
    Wildcard(String),

    /// List of patterns (OR semantics)
    /// - `CA+CB+CD` matches CA or CB or CD
    List(Vec<Pattern>),

    /// Range of values (for chains: A:C matches A, B, C)
    Range(String, String),
}

impl Pattern {
    /// Create an exact match pattern
    pub fn exact(s: impl Into<String>) -> Self {
        Pattern::Exact(s.into())
    }

    /// Create a wildcard pattern
    pub fn wildcard(s: impl Into<String>) -> Self {
        Pattern::Wildcard(s.into())
    }

    /// Create a list pattern from multiple patterns
    pub fn list(patterns: Vec<Pattern>) -> Self {
        Pattern::List(patterns)
    }

    /// Check if a value matches this pattern
    ///
    /// # Arguments
    /// * `value` - The value to match against
    /// * `case_sensitive` - Whether to perform case-sensitive matching
    pub fn matches(&self, value: &str, case_sensitive: bool) -> bool {
        match self {
            Pattern::Exact(pattern) => {
                if case_sensitive {
                    value == pattern
                } else {
                    value.eq_ignore_ascii_case(pattern)
                }
            }
            Pattern::Wildcard(pattern) => {
                Self::match_wildcard(pattern, value, case_sensitive)
            }
            Pattern::List(patterns) => {
                patterns.iter().any(|p| p.matches(value, case_sensitive))
            }
            Pattern::Range(start, end) => {
                let value_cmp = if case_sensitive {
                    value.to_string()
                } else {
                    value.to_uppercase()
                };
                let start_cmp = if case_sensitive {
                    start.clone()
                } else {
                    start.to_uppercase()
                };
                let end_cmp = if case_sensitive {
                    end.clone()
                } else {
                    end.to_uppercase()
                };
                value_cmp >= start_cmp && value_cmp <= end_cmp
            }
        }
    }

    /// Match a wildcard pattern against a value
    fn match_wildcard(pattern: &str, value: &str, case_sensitive: bool) -> bool {
        let pattern = if case_sensitive {
            pattern.to_string()
        } else {
            pattern.to_uppercase()
        };
        let value = if case_sensitive {
            value.to_string()
        } else {
            value.to_uppercase()
        };

        Self::match_wildcard_impl(&pattern, &value)
    }

    /// Recursive wildcard matching implementation
    fn match_wildcard_impl(pattern: &str, value: &str) -> bool {
        let mut p_chars = pattern.chars().peekable();
        let mut v_chars = value.chars().peekable();

        while let Some(pc) = p_chars.next() {
            match pc {
                '*' => {
                    // * matches zero or more characters
                    // Check if the rest of the pattern matches any suffix
                    let remaining_pattern: String = p_chars.collect();
                    if remaining_pattern.is_empty() {
                        return true; // * at end matches everything
                    }

                    // Try matching remaining pattern at each position
                    let remaining_value: String = v_chars.collect();
                    for i in 0..=remaining_value.len() {
                        if Self::match_wildcard_impl(&remaining_pattern, &remaining_value[i..]) {
                            return true;
                        }
                    }
                    return false;
                }
                '?' => {
                    // ? matches exactly one character
                    if v_chars.next().is_none() {
                        return false;
                    }
                }
                c => {
                    // Regular character must match exactly
                    if v_chars.next() != Some(c) {
                        return false;
                    }
                }
            }
        }

        // Pattern exhausted, value must also be exhausted
        v_chars.next().is_none()
    }

    /// Check if this pattern contains wildcards
    pub fn has_wildcards(&self) -> bool {
        match self {
            Pattern::Exact(_) => false,
            Pattern::Wildcard(_) => true,
            Pattern::List(patterns) => patterns.iter().any(|p| p.has_wildcards()),
            Pattern::Range(_, _) => false,
        }
    }
}

impl std::fmt::Display for Pattern {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Pattern::Exact(s) => write!(f, "{}", s),
            Pattern::Wildcard(s) => write!(f, "{}", s),
            Pattern::List(patterns) => {
                let strings: Vec<String> = patterns.iter().map(|p| p.to_string()).collect();
                write!(f, "{}", strings.join("+"))
            }
            Pattern::Range(start, end) => write!(f, "{}:{}", start, end),
        }
    }
}

/// Specification for integer values (like resi, index, id)
#[derive(Debug, Clone, PartialEq)]
pub struct IntSpec {
    /// List of items to match
    pub items: Vec<IntItem>,
}

impl IntSpec {
    /// Create from a single value
    pub fn single(value: i32) -> Self {
        IntSpec {
            items: vec![IntItem::Single(value)],
        }
    }

    /// Create from a range
    pub fn range(start: i32, end: i32) -> Self {
        IntSpec {
            items: vec![IntItem::Range(start, end)],
        }
    }

    /// Check if a value matches this specification
    pub fn matches(&self, value: i32) -> bool {
        self.items.iter().any(|item| item.matches(value))
    }
}

/// A single item in an integer specification
#[derive(Debug, Clone, PartialEq)]
pub enum IntItem {
    /// Single value (e.g., `100`)
    Single(i32),
    /// Range of values (e.g., `100-200`)
    Range(i32, i32),
}

impl IntItem {
    /// Check if a value matches this item
    pub fn matches(&self, value: i32) -> bool {
        match self {
            IntItem::Single(v) => value == *v,
            IntItem::Range(start, end) => value >= *start && value <= *end,
        }
    }
}

/// Specification for residue identifiers
///
/// Supports residue numbers, ranges, and insertion codes.
#[derive(Debug, Clone, PartialEq)]
pub struct ResiSpec {
    /// List of residue items to match
    pub items: Vec<ResiItem>,
}

impl ResiSpec {
    /// Create from a single residue number
    pub fn single(resv: i32) -> Self {
        ResiSpec {
            items: vec![ResiItem::Single(resv)],
        }
    }

    /// Create from a range
    pub fn range(start: i32, end: i32) -> Self {
        ResiSpec {
            items: vec![ResiItem::Range(start, end)],
        }
    }

    /// Check if a residue matches this specification
    ///
    /// # Arguments
    /// * `resv` - Residue sequence number
    /// * `inscode` - Insertion code (or space if none)
    pub fn matches(&self, resv: i32, inscode: char) -> bool {
        self.items.iter().any(|item| item.matches(resv, inscode))
    }
}

/// A single item in a residue specification
#[derive(Debug, Clone, PartialEq)]
pub enum ResiItem {
    /// Single residue number (e.g., `100`)
    Single(i32),

    /// Range of residue numbers (e.g., `100-200`)
    Range(i32, i32),

    /// Residue with insertion code (e.g., `100A`)
    InsCode(i32, char),

    /// Range with insertion codes (e.g., `100A-100D`)
    InsCodeRange(i32, char, i32, char),
}

impl ResiItem {
    /// Check if a residue matches this item
    pub fn matches(&self, resv: i32, inscode: char) -> bool {
        match self {
            ResiItem::Single(v) => resv == *v,
            ResiItem::Range(start, end) => resv >= *start && resv <= *end,
            ResiItem::InsCode(v, code) => {
                resv == *v && (inscode == *code || (*code == ' ' && inscode == ' '))
            }
            ResiItem::InsCodeRange(start_v, start_c, end_v, end_c) => {
                // This is a simplified implementation
                // In PyMOL, this compares lexicographically
                if resv < *start_v || resv > *end_v {
                    return false;
                }
                if resv == *start_v && inscode < *start_c {
                    return false;
                }
                if resv == *end_v && inscode > *end_c {
                    return false;
                }
                true
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pattern_exact() {
        let pattern = Pattern::exact("CA");
        assert!(pattern.matches("CA", false));
        assert!(pattern.matches("ca", false));
        assert!(pattern.matches("Ca", false));
        assert!(!pattern.matches("CB", false));

        // Case sensitive
        assert!(pattern.matches("CA", true));
        assert!(!pattern.matches("ca", true));
    }

    #[test]
    fn test_pattern_wildcard() {
        let pattern = Pattern::wildcard("C*");
        assert!(pattern.matches("CA", false));
        assert!(pattern.matches("CB", false));
        assert!(pattern.matches("C", false));
        assert!(!pattern.matches("N", false));

        let pattern = Pattern::wildcard("*H");
        assert!(pattern.matches("H", false));
        assert!(pattern.matches("1H", false));
        assert!(pattern.matches("OH", false));
        assert!(!pattern.matches("HA", false));

        let pattern = Pattern::wildcard("*H*");
        assert!(pattern.matches("H", false));
        assert!(pattern.matches("HA", false));
        assert!(pattern.matches("OH", false));
        assert!(pattern.matches("1H2", false));
    }

    #[test]
    fn test_pattern_list() {
        let pattern = Pattern::list(vec![
            Pattern::exact("CA"),
            Pattern::exact("CB"),
            Pattern::exact("CD"),
        ]);
        assert!(pattern.matches("CA", false));
        assert!(pattern.matches("CB", false));
        assert!(pattern.matches("CD", false));
        assert!(!pattern.matches("CE", false));
    }

    #[test]
    fn test_pattern_range() {
        let pattern = Pattern::Range("A".to_string(), "C".to_string());
        assert!(pattern.matches("A", false));
        assert!(pattern.matches("B", false));
        assert!(pattern.matches("C", false));
        assert!(!pattern.matches("D", false));
    }

    #[test]
    fn test_int_spec() {
        let spec = IntSpec {
            items: vec![IntItem::Single(100), IntItem::Range(200, 210)],
        };
        assert!(spec.matches(100));
        assert!(spec.matches(200));
        assert!(spec.matches(205));
        assert!(spec.matches(210));
        assert!(!spec.matches(101));
        assert!(!spec.matches(199));
        assert!(!spec.matches(211));
    }

    #[test]
    fn test_resi_spec() {
        let spec = ResiSpec {
            items: vec![ResiItem::Single(100), ResiItem::Range(200, 210)],
        };
        assert!(spec.matches(100, ' '));
        assert!(spec.matches(205, ' '));
        assert!(!spec.matches(150, ' '));
    }

    #[test]
    fn test_resi_inscode() {
        let spec = ResiSpec {
            items: vec![ResiItem::InsCode(100, 'A')],
        };
        assert!(spec.matches(100, 'A'));
        assert!(!spec.matches(100, 'B'));
        assert!(!spec.matches(100, ' '));
    }

    #[test]
    fn test_pattern_display() {
        assert_eq!(Pattern::exact("CA").to_string(), "CA");
        assert_eq!(Pattern::wildcard("C*").to_string(), "C*");
        assert_eq!(
            Pattern::list(vec![Pattern::exact("CA"), Pattern::exact("CB")]).to_string(),
            "CA+CB"
        );
        assert_eq!(
            Pattern::Range("A".to_string(), "C".to_string()).to_string(),
            "A:C"
        );
    }
}
