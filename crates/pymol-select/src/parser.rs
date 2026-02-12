//! Parser for the selection language
//!
//! Converts a token stream into an AST using recursive descent parsing
//! with operator precedence handling.

use crate::ast::{CompareOp, MacroSpec, PropertyValue, SelectionExpr};
use crate::error::ParseError;
use crate::keywords::{lookup, Keyword, KeywordType};
use crate::lexer::{Token, TokenStream};
use crate::pattern::{IntItem, IntSpec, Pattern, ResiItem, ResiSpec};

/// Parse a selection string into an AST
pub fn parse_selection(input: &str) -> Result<SelectionExpr, ParseError> {
    let mut stream = TokenStream::from_str(input)?;
    let expr = parse_expr(&mut stream, 0)?;

    // Ensure we consumed all input
    if !stream.is_eof() {
        return Err(ParseError::UnexpectedToken(format!(
            "{:?}",
            stream.peek()
        )));
    }

    Ok(expr)
}

/// Parse an expression with minimum precedence
fn parse_expr(stream: &mut TokenStream, min_prec: u8) -> Result<SelectionExpr, ParseError> {
    let mut left = parse_unary(stream)?;

    loop {
        let Some(tok) = stream.peek().cloned() else {
            break;
        };

        // Check for binary operators
        let (op_kw, prec) = match &tok {
            Token::Ident(s) => match lookup(s) {
                Some(kw @ Keyword::And) => (Some(kw), kw.precedence()),
                Some(kw @ Keyword::Or) => (Some(kw), kw.precedence()),
                Some(kw @ Keyword::In) => (Some(kw), kw.precedence()),
                Some(kw @ Keyword::Like) => (Some(kw), kw.precedence()),
                _ => break, // Not a binary operator
            },
            Token::Ampersand => (Some(Keyword::And), Keyword::And.precedence()),
            Token::Pipe => (Some(Keyword::Or), Keyword::Or.precedence()),
            Token::Plus => (Some(Keyword::Or), Keyword::Or.precedence()),
            Token::Minus => {
                // Minus is AND NOT
                (Some(Keyword::And), Keyword::And.precedence())
            }
            _ => break, // Not a binary operator
        };

        let Some(op_kw) = op_kw else {
            break;
        };

        if prec < min_prec {
            break;
        }

        // Consume the operator
        stream.next();

        // Handle minus as AND NOT
        let is_and_not = matches!(tok, Token::Minus);

        // Parse the right side with higher precedence
        let right = parse_expr(stream, prec + 1)?;

        // Build the expression
        left = match op_kw {
            Keyword::And => {
                if is_and_not {
                    SelectionExpr::And(Box::new(left), Box::new(SelectionExpr::Not(Box::new(right))))
                } else {
                    SelectionExpr::And(Box::new(left), Box::new(right))
                }
            }
            Keyword::Or => SelectionExpr::Or(Box::new(left), Box::new(right)),
            Keyword::In => SelectionExpr::In(Box::new(left), Box::new(right)),
            Keyword::Like => SelectionExpr::Like(Box::new(left), Box::new(right)),
            _ => unreachable!(),
        };
    }

    Ok(left)
}

/// Parse a unary expression (NOT, prefix operators)
fn parse_unary(stream: &mut TokenStream) -> Result<SelectionExpr, ParseError> {
    let Some(tok) = stream.peek().cloned() else {
        return Err(ParseError::UnexpectedEof);
    };

    match &tok {
        Token::Ident(s) => {
            if let Some(kw) = lookup(s) {
                match kw.keyword_type() {
                    KeywordType::Opr1 => {
                        return parse_prefix_operator(stream, kw);
                    }
                    KeywordType::Op22 => {
                        return parse_distance_binary(stream, kw);
                    }
                    _ => {}
                }
            }
        }
        Token::Exclamation => {
            stream.next();
            let inner = parse_unary(stream)?;
            return Ok(SelectionExpr::Not(Box::new(inner)));
        }
        _ => {}
    }

    parse_primary(stream)
}

/// Parse prefix operators (not, byres, bychain, etc.)
fn parse_prefix_operator(
    stream: &mut TokenStream,
    kw: Keyword,
) -> Result<SelectionExpr, ParseError> {
    stream.next(); // Consume the keyword

    match kw {
        Keyword::Not => {
            let inner = parse_unary(stream)?;
            Ok(SelectionExpr::Not(Box::new(inner)))
        }
        Keyword::ByRes => {
            let inner = parse_unary(stream)?;
            Ok(SelectionExpr::ByRes(Box::new(inner)))
        }
        Keyword::ByChain => {
            let inner = parse_unary(stream)?;
            Ok(SelectionExpr::ByChain(Box::new(inner)))
        }
        Keyword::ByObject => {
            let inner = parse_unary(stream)?;
            Ok(SelectionExpr::ByObject(Box::new(inner)))
        }
        Keyword::ByMolecule => {
            let inner = parse_unary(stream)?;
            Ok(SelectionExpr::ByMolecule(Box::new(inner)))
        }
        Keyword::BySegment => {
            let inner = parse_unary(stream)?;
            Ok(SelectionExpr::BySegment(Box::new(inner)))
        }
        Keyword::ByFragment => {
            let inner = parse_unary(stream)?;
            Ok(SelectionExpr::ByFragment(Box::new(inner)))
        }
        Keyword::ByCAlpha => {
            let inner = parse_unary(stream)?;
            Ok(SelectionExpr::ByCAlpha(Box::new(inner)))
        }
        Keyword::ByRing => {
            let inner = parse_unary(stream)?;
            Ok(SelectionExpr::ByRing(Box::new(inner)))
        }
        Keyword::ByCell => {
            let inner = parse_unary(stream)?;
            Ok(SelectionExpr::ByCell(Box::new(inner)))
        }
        Keyword::Neighbor => {
            let inner = parse_unary(stream)?;
            Ok(SelectionExpr::Neighbor(Box::new(inner)))
        }
        Keyword::BoundTo => {
            let inner = parse_unary(stream)?;
            Ok(SelectionExpr::BoundTo(Box::new(inner)))
        }
        Keyword::First => {
            let inner = parse_unary(stream)?;
            Ok(SelectionExpr::First(Box::new(inner)))
        }
        Keyword::Last => {
            let inner = parse_unary(stream)?;
            Ok(SelectionExpr::Last(Box::new(inner)))
        }
        Keyword::Around => {
            let dist = parse_number(stream)?;
            let inner = parse_unary(stream)?;
            Ok(SelectionExpr::Around(dist, Box::new(inner)))
        }
        Keyword::Expand => {
            let dist = parse_number(stream)?;
            let inner = parse_unary(stream)?;
            Ok(SelectionExpr::Expand(dist, Box::new(inner)))
        }
        Keyword::Extend => {
            let n = parse_int(stream)? as u32;
            let inner = parse_unary(stream)?;
            Ok(SelectionExpr::Extend(n, Box::new(inner)))
        }
        Keyword::Gap => {
            let dist = parse_number(stream)?;
            let inner = parse_unary(stream)?;
            Ok(SelectionExpr::Gap(dist, Box::new(inner)))
        }
        _ => Err(ParseError::UnknownKeyword(format!("{:?}", kw))),
    }
}

/// Parse distance binary operators (within, beyond, near_to)
fn parse_distance_binary(
    stream: &mut TokenStream,
    kw: Keyword,
) -> Result<SelectionExpr, ParseError> {
    stream.next(); // Consume the keyword

    let dist = parse_number(stream)?;

    // Expect "of"
    match stream.peek() {
        Some(Token::Of) => {
            stream.next();
        }
        _ => {
            return Err(ParseError::Expected {
                expected: "of".to_string(),
                found: format!("{:?}", stream.peek()),
            });
        }
    }

    let inner = parse_unary(stream)?;

    match kw {
        Keyword::Within => Ok(SelectionExpr::Within(
            dist,
            Box::new(SelectionExpr::All),
            Box::new(inner),
        )),
        Keyword::Beyond => Ok(SelectionExpr::Beyond(
            dist,
            Box::new(SelectionExpr::All),
            Box::new(inner),
        )),
        Keyword::NearTo => Ok(SelectionExpr::NearTo(
            dist,
            Box::new(SelectionExpr::All),
            Box::new(inner),
        )),
        _ => Err(ParseError::UnknownKeyword(format!("{:?}", kw))),
    }
}

/// Parse a primary expression (atom, parenthesized, property selector)
fn parse_primary(stream: &mut TokenStream) -> Result<SelectionExpr, ParseError> {
    let Some(tok) = stream.peek().cloned() else {
        return Err(ParseError::UnexpectedEof);
    };

    match tok {
        Token::LParen => {
            stream.next();
            let inner = parse_expr(stream, 0)?;
            match stream.next() {
                Some(Token::RParen) => Ok(inner),
                _ => Err(ParseError::UnmatchedParen),
            }
        }
        Token::Slash => parse_macro(stream),
        Token::Percent => {
            stream.next();
            // Selection reference
            match stream.next() {
                Some(Token::Ident(name)) => Ok(SelectionExpr::Selection(name)),
                _ => Err(ParseError::MissingArgument("selection name".to_string())),
            }
        }
        Token::Ident(ref s) => {
            if let Some(kw) = lookup(s) {
                match kw.keyword_type() {
                    KeywordType::Sel0 => parse_sel0(stream, kw),
                    KeywordType::Sel1 => parse_sel1(stream, kw),
                    KeywordType::Sel2 => parse_sel2(stream, kw),
                    KeywordType::Sel3 => parse_sel3(stream, kw),
                    KeywordType::Prp1 => parse_prefix_operator(stream, kw),
                    KeywordType::Opr1 => parse_prefix_operator(stream, kw),
                    KeywordType::Op22 => parse_distance_binary(stream, kw),
                    KeywordType::Opr2 => {
                        // Binary operators shouldn't appear here
                        Err(ParseError::UnexpectedToken(s.clone()))
                    }
                }
            } else {
                // Treat as a selection name or model name
                stream.next();
                Ok(SelectionExpr::Selection(s.clone()))
            }
        }
        Token::Integer(n) => {
            // Could be a residue number in shorthand
            stream.next();
            // Try to parse as shorthand selection
            if let Some(Token::Slash) = stream.peek() {
                // This is macro notation
                stream.set_position(stream.position() - 1);
                parse_macro(stream)
            } else {
                // Treat as resi shorthand
                Ok(SelectionExpr::Resi(ResiSpec::single(n)))
            }
        }
        Token::Asterisk => {
            stream.next();
            Ok(SelectionExpr::All)
        }
        _ => Err(ParseError::UnexpectedToken(format!("{:?}", tok))),
    }
}

/// Parse zero-argument selections
fn parse_sel0(stream: &mut TokenStream, kw: Keyword) -> Result<SelectionExpr, ParseError> {
    stream.next(); // Consume the keyword

    match kw {
        Keyword::All => Ok(SelectionExpr::All),
        Keyword::None => Ok(SelectionExpr::None),
        Keyword::Hetatm => Ok(SelectionExpr::Hetatm),
        Keyword::Hydrogens => Ok(SelectionExpr::Hydrogens),
        Keyword::Visible => Ok(SelectionExpr::Visible),
        Keyword::Enabled => Ok(SelectionExpr::Enabled),
        Keyword::Bonded => Ok(SelectionExpr::Bonded),
        Keyword::Polymer => Ok(SelectionExpr::Polymer),
        Keyword::PolymerProtein => Ok(SelectionExpr::PolymerProtein),
        Keyword::PolymerNucleic => Ok(SelectionExpr::PolymerNucleic),
        Keyword::Organic => Ok(SelectionExpr::Organic),
        Keyword::Inorganic => Ok(SelectionExpr::Inorganic),
        Keyword::Solvent => Ok(SelectionExpr::Solvent),
        Keyword::Metals => Ok(SelectionExpr::Metals),
        Keyword::Backbone => Ok(SelectionExpr::Backbone),
        Keyword::Sidechain => Ok(SelectionExpr::Sidechain),
        Keyword::Donors => Ok(SelectionExpr::Donors),
        Keyword::Acceptors => Ok(SelectionExpr::Acceptors),
        Keyword::HBondAcceptors => Ok(SelectionExpr::Acceptors),
        Keyword::HBondDonors => Ok(SelectionExpr::Donors),
        Keyword::Delocalized => Ok(SelectionExpr::Delocalized),
        Keyword::Fixed => Ok(SelectionExpr::Fixed),
        Keyword::Restrained => Ok(SelectionExpr::Restrained),
        Keyword::Masked => Ok(SelectionExpr::Masked),
        Keyword::Protected => Ok(SelectionExpr::Protected),
        Keyword::Present => Ok(SelectionExpr::Present),
        Keyword::Guide => Ok(SelectionExpr::Guide),
        Keyword::Origin => Ok(SelectionExpr::Origin),
        Keyword::Center => Ok(SelectionExpr::Center),
        _ => Err(ParseError::UnknownKeyword(format!("{:?}", kw))),
    }
}

/// Parse one-argument property selections
fn parse_sel1(stream: &mut TokenStream, kw: Keyword) -> Result<SelectionExpr, ParseError> {
    stream.next(); // Consume the keyword

    match kw {
        Keyword::Name => {
            let pattern = parse_pattern(stream)?;
            Ok(SelectionExpr::Name(pattern))
        }
        Keyword::Element => {
            let pattern = parse_pattern(stream)?;
            Ok(SelectionExpr::Elem(pattern))
        }
        Keyword::Resi => {
            let spec = parse_resi_spec(stream)?;
            Ok(SelectionExpr::Resi(spec))
        }
        Keyword::Resn => {
            let pattern = parse_pattern(stream)?;
            Ok(SelectionExpr::Resn(pattern))
        }
        Keyword::Chain => {
            let pattern = parse_pattern(stream)?;
            Ok(SelectionExpr::Chain(pattern))
        }
        Keyword::Segi => {
            let pattern = parse_pattern(stream)?;
            Ok(SelectionExpr::Segi(pattern))
        }
        Keyword::Alt => {
            let pattern = parse_pattern(stream)?;
            Ok(SelectionExpr::Alt(pattern))
        }
        Keyword::Flag => {
            let spec = parse_int_spec(stream)?;
            Ok(SelectionExpr::Flag(spec))
        }
        Keyword::TextType => {
            let pattern = parse_pattern(stream)?;
            Ok(SelectionExpr::TextType(pattern))
        }
        Keyword::NumericType => {
            let pattern = parse_pattern(stream)?;
            Ok(SelectionExpr::NumericType(pattern))
        }
        Keyword::Id => {
            let spec = parse_int_spec(stream)?;
            Ok(SelectionExpr::Id(spec))
        }
        Keyword::Index => {
            let spec = parse_int_spec(stream)?;
            Ok(SelectionExpr::Index(spec))
        }
        Keyword::Rank => {
            let spec = parse_int_spec(stream)?;
            Ok(SelectionExpr::Rank(spec))
        }
        Keyword::Model => {
            let pattern = parse_pattern(stream)?;
            Ok(SelectionExpr::Model(pattern))
        }
        Keyword::State => {
            let spec = parse_int_spec(stream)?;
            Ok(SelectionExpr::State(spec))
        }
        Keyword::SecondaryStructure => {
            let pattern = parse_pattern(stream)?;
            Ok(SelectionExpr::SecondaryStructure(pattern))
        }
        Keyword::Rep => {
            let pattern = parse_pattern(stream)?;
            Ok(SelectionExpr::Rep(pattern))
        }
        Keyword::Color => {
            let pattern = parse_pattern(stream)?;
            Ok(SelectionExpr::Color(pattern))
        }
        Keyword::CartoonColor => {
            let pattern = parse_pattern(stream)?;
            Ok(SelectionExpr::CartoonColor(pattern))
        }
        Keyword::RibbonColor => {
            let pattern = parse_pattern(stream)?;
            Ok(SelectionExpr::RibbonColor(pattern))
        }
        Keyword::PepSeq => {
            let pattern = parse_pattern(stream)?;
            Ok(SelectionExpr::PepSeq(pattern))
        }
        Keyword::Custom => {
            let pattern = parse_pattern(stream)?;
            Ok(SelectionExpr::Custom(pattern))
        }
        Keyword::Label => {
            let pattern = parse_pattern(stream)?;
            Ok(SelectionExpr::Label(pattern))
        }
        Keyword::Stereo => {
            let pattern = parse_pattern(stream)?;
            Ok(SelectionExpr::Stereo(pattern))
        }
        Keyword::Selection => {
            match stream.next() {
                Some(Token::Ident(name)) => Ok(SelectionExpr::Selection(name)),
                _ => Err(ParseError::MissingArgument("selection name".to_string())),
            }
        }
        _ => Err(ParseError::UnknownKeyword(format!("{:?}", kw))),
    }
}

/// Parse two-argument numeric comparison selections
fn parse_sel2(stream: &mut TokenStream, kw: Keyword) -> Result<SelectionExpr, ParseError> {
    stream.next(); // Consume the keyword

    let op = parse_compare_op(stream)?;
    let value = parse_number(stream)?;

    match kw {
        Keyword::BFactor => Ok(SelectionExpr::BFactor(op, value)),
        Keyword::Occupancy => Ok(SelectionExpr::Occupancy(op, value)),
        Keyword::PartialCharge => Ok(SelectionExpr::PartialCharge(op, value)),
        Keyword::FormalCharge => Ok(SelectionExpr::FormalCharge(op, value as i32)),
        Keyword::X => Ok(SelectionExpr::X(op, value)),
        Keyword::Y => Ok(SelectionExpr::Y(op, value)),
        Keyword::Z => Ok(SelectionExpr::Z(op, value)),
        _ => Err(ParseError::UnknownKeyword(format!("{:?}", kw))),
    }
}

/// Parse three-argument property access (p.property op value)
fn parse_sel3(stream: &mut TokenStream, _kw: Keyword) -> Result<SelectionExpr, ParseError> {
    stream.next(); // Consume "p."

    // Get the property name
    let prop_name = match stream.next() {
        Some(Token::Ident(name)) => name,
        _ => return Err(ParseError::MissingArgument("property name".to_string())),
    };

    let op = parse_compare_op(stream)?;

    // Get the value
    let value = match stream.next() {
        Some(Token::Integer(n)) => PropertyValue::Int(n),
        Some(Token::Float(f)) => PropertyValue::Float(f),
        Some(Token::Ident(s)) => PropertyValue::String(s),
        Some(Token::QuotedString(s)) => PropertyValue::String(s),
        _ => return Err(ParseError::MissingArgument("property value".to_string())),
    };

    Ok(SelectionExpr::Property(prop_name, op, value))
}

/// Parse a comparison operator
fn parse_compare_op(stream: &mut TokenStream) -> Result<CompareOp, ParseError> {
    match stream.next() {
        Some(Token::Equals) => Ok(CompareOp::Eq),
        Some(Token::NotEquals) => Ok(CompareOp::Ne),
        Some(Token::LessThan) => Ok(CompareOp::Lt),
        Some(Token::LessOrEqual) => Ok(CompareOp::Le),
        Some(Token::GreaterThan) => Ok(CompareOp::Gt),
        Some(Token::GreaterOrEqual) => Ok(CompareOp::Ge),
        other => Err(ParseError::Expected {
            expected: "comparison operator".to_string(),
            found: format!("{:?}", other),
        }),
    }
}

/// Parse a number (int or float as f32)
fn parse_number(stream: &mut TokenStream) -> Result<f32, ParseError> {
    match stream.next() {
        Some(Token::Integer(n)) => Ok(n as f32),
        Some(Token::Float(f)) => Ok(f),
        other => Err(ParseError::Expected {
            expected: "number".to_string(),
            found: format!("{:?}", other),
        }),
    }
}

/// Parse an integer
fn parse_int(stream: &mut TokenStream) -> Result<i32, ParseError> {
    match stream.next() {
        Some(Token::Integer(n)) => Ok(n),
        other => Err(ParseError::Expected {
            expected: "integer".to_string(),
            found: format!("{:?}", other),
        }),
    }
}

/// Parse a pattern (identifier, wildcard, list, or range)
fn parse_pattern(stream: &mut TokenStream) -> Result<Pattern, ParseError> {
    let first = parse_single_pattern(stream)?;

    // Check for list (+ separator) or range (: separator)
    match stream.peek() {
        Some(Token::Plus) => {
            let mut patterns = vec![first];
            while matches!(stream.peek(), Some(Token::Plus)) {
                stream.next(); // Consume +
                patterns.push(parse_single_pattern(stream)?);
            }
            Ok(Pattern::List(patterns))
        }
        Some(Token::Colon) => {
            stream.next(); // Consume :
            let second = parse_single_pattern(stream)?;
            match (first, second) {
                (Pattern::Exact(start), Pattern::Exact(end)) => Ok(Pattern::Range(start, end)),
                _ => Err(ParseError::InvalidPattern(
                    "range requires exact patterns".to_string(),
                )),
            }
        }
        _ => Ok(first),
    }
}

/// Parse a single pattern element
fn parse_single_pattern(stream: &mut TokenStream) -> Result<Pattern, ParseError> {
    match stream.next() {
        Some(Token::Ident(mut s)) => {
            // Check for trailing wildcards
            loop {
                match stream.peek() {
                    Some(Token::Asterisk) => {
                        stream.next();
                        s.push('*');
                    }
                    Some(Token::Question) => {
                        stream.next();
                        s.push('?');
                    }
                    _ => break,
                }
            }
            if s.contains('*') || s.contains('?') {
                Ok(Pattern::Wildcard(s))
            } else {
                Ok(Pattern::Exact(s))
            }
        }
        Some(Token::QuotedString(s)) => Ok(Pattern::Exact(s)),
        Some(Token::Integer(n)) => {
            // Check for wildcard after number
            let mut s = n.to_string();
            loop {
                match stream.peek() {
                    Some(Token::Asterisk) => {
                        stream.next();
                        s.push('*');
                    }
                    Some(Token::Question) => {
                        stream.next();
                        s.push('?');
                    }
                    _ => break,
                }
            }
            if s.contains('*') || s.contains('?') {
                Ok(Pattern::Wildcard(s))
            } else {
                Ok(Pattern::Exact(s))
            }
        }
        Some(Token::Asterisk) => {
            // Could be just "*" or "*" followed by more characters
            let mut s = "*".to_string();
            // Check for identifier following (like *H for ending in H)
            if let Some(Token::Ident(id)) = stream.peek().cloned() {
                stream.next();
                s.push_str(&id);
                // Check for more wildcards
                loop {
                    match stream.peek() {
                        Some(Token::Asterisk) => {
                            stream.next();
                            s.push('*');
                        }
                        Some(Token::Question) => {
                            stream.next();
                            s.push('?');
                        }
                        _ => break,
                    }
                }
            }
            Ok(Pattern::Wildcard(s))
        }
        other => Err(ParseError::Expected {
            expected: "pattern".to_string(),
            found: format!("{:?}", other),
        }),
    }
}

/// Parse an integer specification (for index, id, etc.)
fn parse_int_spec(stream: &mut TokenStream) -> Result<IntSpec, ParseError> {
    let mut items = vec![parse_int_item(stream)?];

    // Check for list (+ separator)
    while matches!(stream.peek(), Some(Token::Plus)) {
        stream.next(); // Consume +
        items.push(parse_int_item(stream)?);
    }

    Ok(IntSpec { items })
}

/// Parse a single integer item (value or range)
fn parse_int_item(stream: &mut TokenStream) -> Result<IntItem, ParseError> {
    let first = parse_int(stream)?;

    // Check for range (- or :)
    match stream.peek() {
        Some(Token::Minus) | Some(Token::Colon) => {
            stream.next(); // Consume - or :
            let second = parse_int(stream)?;
            Ok(IntItem::Range(first, second))
        }
        _ => Ok(IntItem::Single(first)),
    }
}

/// Parse a residue specification
fn parse_resi_spec(stream: &mut TokenStream) -> Result<ResiSpec, ParseError> {
    let mut items = vec![parse_resi_item(stream)?];

    // Check for list (+ separator)
    while matches!(stream.peek(), Some(Token::Plus)) {
        stream.next(); // Consume +
        items.push(parse_resi_item(stream)?);
    }

    Ok(ResiSpec { items })
}

/// Parse a single residue item
fn parse_resi_item(stream: &mut TokenStream) -> Result<ResiItem, ParseError> {
    // First token could be an integer or identifier (for insertion codes)
    let first_tok = stream.next();
    let (first_val, first_ins) = match first_tok {
        Some(Token::Integer(n)) => (n, None),
        Some(Token::Ident(s)) => {
            // Try to parse as number with insertion code (e.g., "100A")
            parse_resi_with_inscode(&s)?
        }
        other => {
            return Err(ParseError::Expected {
                expected: "residue number".to_string(),
                found: format!("{:?}", other),
            });
        }
    };

    // Check for insertion code after integer
    let first_ins = if first_ins.is_none() {
        let maybe_inscode = if let Some(Token::Ident(s)) = stream.peek() {
            if s.len() == 1 && s.chars().next().unwrap().is_alphabetic() {
                Some(s.chars().next().unwrap())
            } else {
                None
            }
        } else {
            None
        };
        if maybe_inscode.is_some() {
            stream.next();
        }
        maybe_inscode
    } else {
        first_ins
    };

    // Check for range (- or :)
    match stream.peek() {
        Some(Token::Minus) | Some(Token::Colon) => {
            stream.next(); // Consume - or :

            // Parse second part
            let second_tok = stream.next();
            let (second_val, second_ins) = match second_tok {
                Some(Token::Integer(n)) => (n, None),
                Some(Token::Ident(s)) => parse_resi_with_inscode(&s)?,
                other => {
                    return Err(ParseError::Expected {
                        expected: "residue number".to_string(),
                        found: format!("{:?}", other),
                    });
                }
            };

            // Check for insertion code after second integer
            let second_ins = if second_ins.is_none() {
                let maybe_inscode = if let Some(Token::Ident(s)) = stream.peek() {
                    if s.len() == 1 && s.chars().next().unwrap().is_alphabetic() {
                        Some(s.chars().next().unwrap())
                    } else {
                        None
                    }
                } else {
                    None
                };
                if maybe_inscode.is_some() {
                    stream.next();
                }
                maybe_inscode
            } else {
                second_ins
            };

            match (first_ins, second_ins) {
                (Some(fi), Some(si)) => Ok(ResiItem::InsCodeRange(first_val, fi, second_val, si)),
                (None, None) => Ok(ResiItem::Range(first_val, second_val)),
                _ => {
                    // Mixed - just use range
                    Ok(ResiItem::Range(first_val, second_val))
                }
            }
        }
        _ => {
            // Single value
            match first_ins {
                Some(c) => Ok(ResiItem::InsCode(first_val, c)),
                None => Ok(ResiItem::Single(first_val)),
            }
        }
    }
}

/// Parse a residue number with optional insertion code (e.g., "100A")
fn parse_resi_with_inscode(s: &str) -> Result<(i32, Option<char>), ParseError> {
    let mut num_part = String::new();
    let mut ins_code = None;

    for c in s.chars() {
        if c.is_ascii_digit() || (num_part.is_empty() && c == '-') {
            num_part.push(c);
        } else if c.is_alphabetic() && ins_code.is_none() {
            ins_code = Some(c);
        } else {
            return Err(ParseError::InvalidResidue(s.to_string()));
        }
    }

    let num: i32 = num_part
        .parse()
        .map_err(|_| ParseError::InvalidResidue(s.to_string()))?;

    Ok((num, ins_code))
}

/// Parse slash macro notation
fn parse_macro(stream: &mut TokenStream) -> Result<SelectionExpr, ParseError> {
    let mut spec = MacroSpec::default();

    // Consume leading slash if present
    if matches!(stream.peek(), Some(Token::Slash)) {
        stream.next();
    }

    // Parse components: model/segi/chain/resn/name
    // Empty fields are wildcards

    let mut components: Vec<Option<String>> = Vec::new();
    let mut current = String::new();
    let mut has_content = false;

    loop {
        // Clone the peeked token to avoid borrow issues
        let peeked = stream.peek().cloned();
        match peeked {
            Some(Token::Slash) => {
                stream.next();
                if has_content {
                    components.push(Some(current.clone()));
                } else {
                    components.push(None);
                }
                current.clear();
                has_content = false;
            }
            Some(Token::Ident(s)) => {
                stream.next();
                current.push_str(&s);
                has_content = true;
            }
            Some(Token::Integer(n)) => {
                stream.next();
                current.push_str(&n.to_string());
                has_content = true;
            }
            Some(Token::Backtick) => {
                stream.next();
                current.push('`');
                has_content = true;
            }
            Some(Token::Asterisk) => {
                stream.next();
                current.push('*');
                has_content = true;
            }
            _ => {
                // End of macro
                if has_content {
                    components.push(Some(current.clone()));
                }
                break;
            }
        }
    }

    // Assign components to spec fields based on position
    // Format: model/segi/chain/resn/name
    // But if started with /, first component is model
    let mut idx = 0;
    for component in components {
        let pattern = component.map(|s| {
            if s.contains('*') || s.contains('?') {
                Pattern::Wildcard(s)
            } else {
                Pattern::Exact(s)
            }
        });

        match idx {
            0 => spec.model = pattern,
            1 => spec.segi = pattern,
            2 => spec.chain = pattern,
            3 => {
                // This could be resn or resn`resi
                if let Some(ref p) = pattern {
                    if let Pattern::Exact(ref s) = p {
                        if let Some(pos) = s.find('`') {
                            let (resn_part, resi_part) = s.split_at(pos);
                            let resi_part = &resi_part[1..]; // Skip backtick
                            spec.resn = if resn_part.is_empty() {
                                None
                            } else {
                                Some(Pattern::Exact(resn_part.to_string()))
                            };
                            if let Ok((resi_num, ins_code)) = parse_resi_with_inscode(resi_part) {
                                spec.resi = Some(ResiSpec {
                                    items: vec![match ins_code {
                                        Some(c) => ResiItem::InsCode(resi_num, c),
                                        None => ResiItem::Single(resi_num),
                                    }],
                                });
                            }
                        } else {
                            // Could be resn or resi depending on content
                            if s.chars().next().map(|c| c.is_ascii_digit()).unwrap_or(false) {
                                // Starts with digit - treat as resi
                                if let Ok((resi_num, ins_code)) = parse_resi_with_inscode(s) {
                                    spec.resi = Some(ResiSpec {
                                        items: vec![match ins_code {
                                            Some(c) => ResiItem::InsCode(resi_num, c),
                                            None => ResiItem::Single(resi_num),
                                        }],
                                    });
                                }
                            } else {
                                spec.resn = pattern;
                            }
                        }
                    } else {
                        spec.resn = pattern;
                    }
                }
            }
            4 => {
                // name or name`alt
                if let Some(ref p) = pattern {
                    if let Pattern::Exact(ref s) = p {
                        if let Some(pos) = s.find('`') {
                            let (name_part, alt_part) = s.split_at(pos);
                            let alt_part = &alt_part[1..]; // Skip backtick
                            spec.name = if name_part.is_empty() {
                                None
                            } else {
                                Some(Pattern::Exact(name_part.to_string()))
                            };
                            spec.alt = if alt_part.is_empty() {
                                None
                            } else {
                                Some(Pattern::Exact(alt_part.to_string()))
                            };
                        } else {
                            spec.name = pattern;
                        }
                    } else {
                        spec.name = pattern;
                    }
                }
            }
            _ => {} // Ignore extra components
        }
        idx += 1;
    }

    Ok(SelectionExpr::Macro(spec))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_all() {
        let expr = parse_selection("all").unwrap();
        assert!(matches!(expr, SelectionExpr::All));
    }

    #[test]
    fn test_parse_name() {
        let expr = parse_selection("name CA").unwrap();
        assert!(matches!(expr, SelectionExpr::Name(Pattern::Exact(ref s)) if s == "CA"));
    }

    #[test]
    fn test_parse_name_wildcard() {
        let expr = parse_selection("name C*").unwrap();
        assert!(matches!(expr, SelectionExpr::Name(Pattern::Wildcard(ref s)) if s == "C*"));
    }

    #[test]
    fn test_parse_and() {
        let expr = parse_selection("name CA and chain A").unwrap();
        assert!(matches!(expr, SelectionExpr::And(_, _)));
    }

    #[test]
    fn test_parse_or() {
        let expr = parse_selection("name CA or name CB").unwrap();
        assert!(matches!(expr, SelectionExpr::Or(_, _)));
    }

    #[test]
    fn test_parse_not() {
        let expr = parse_selection("not hydrogens").unwrap();
        assert!(matches!(expr, SelectionExpr::Not(_)));
    }

    #[test]
    fn test_parse_parens() {
        let expr = parse_selection("(name CA)").unwrap();
        assert!(matches!(expr, SelectionExpr::Name(_)));
    }

    #[test]
    fn test_parse_complex() {
        let expr = parse_selection("(name CA and chain A) or name CB").unwrap();
        assert!(matches!(expr, SelectionExpr::Or(_, _)));
    }

    #[test]
    fn test_parse_parens_with_or() {
        let expr = parse_selection("(chain A and resi 29) or (chain A and resi 31)").unwrap();
        assert!(matches!(expr, SelectionExpr::Or(_, _)));
    }

    #[test]
    fn test_parse_resi_single() {
        let expr = parse_selection("resi 100").unwrap();
        if let SelectionExpr::Resi(spec) = expr {
            assert_eq!(spec.items.len(), 1);
            assert!(matches!(spec.items[0], ResiItem::Single(100)));
        } else {
            panic!("Expected Resi");
        }
    }

    #[test]
    fn test_parse_resi_range() {
        let expr = parse_selection("resi 100-200").unwrap();
        if let SelectionExpr::Resi(spec) = expr {
            assert_eq!(spec.items.len(), 1);
            assert!(matches!(spec.items[0], ResiItem::Range(100, 200)));
        } else {
            panic!("Expected Resi");
        }
    }

    #[test]
    fn test_parse_byres() {
        let expr = parse_selection("byres name CA").unwrap();
        assert!(matches!(expr, SelectionExpr::ByRes(_)));
    }

    #[test]
    fn test_parse_within() {
        let expr = parse_selection("within 5 of name CA").unwrap();
        if let SelectionExpr::Within(dist, _, _) = expr {
            assert!((dist - 5.0).abs() < 0.01);
        } else {
            panic!("Expected Within");
        }
    }

    #[test]
    fn test_parse_b_factor() {
        let expr = parse_selection("b > 50").unwrap();
        if let SelectionExpr::BFactor(op, val) = expr {
            assert_eq!(op, CompareOp::Gt);
            assert!((val - 50.0).abs() < 0.01);
        } else {
            panic!("Expected BFactor");
        }
    }

    #[test]
    fn test_parse_name_list() {
        let expr = parse_selection("name CA+CB+CD").unwrap();
        if let SelectionExpr::Name(Pattern::List(patterns)) = expr {
            assert_eq!(patterns.len(), 3);
        } else {
            panic!("Expected Name with List pattern");
        }
    }

    #[test]
    fn test_parse_selection_reference() {
        let expr = parse_selection("%sele1").unwrap();
        if let SelectionExpr::Selection(name) = expr {
            assert_eq!(name, "sele1");
        } else {
            panic!("Expected Selection");
        }
    }

    #[test]
    fn test_parse_macro() {
        let expr = parse_selection("/protein//A/100/CA").unwrap();
        if let SelectionExpr::Macro(spec) = expr {
            assert!(matches!(spec.model, Some(Pattern::Exact(ref s)) if s == "protein"));
            assert!(matches!(spec.chain, Some(Pattern::Exact(ref s)) if s == "A"));
        } else {
            panic!("Expected Macro");
        }
    }
}
