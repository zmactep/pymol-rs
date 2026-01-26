//! CIF/STAR format lexer
//!
//! Tokenizes CIF format files.

use nom::{
    branch::alt,
    bytes::complete::{tag, take_while, take_while1},
    character::complete::{char, multispace1, not_line_ending},
    combinator::value,
    multi::many0,
    sequence::preceded,
    IResult,
};

/// CIF token types
#[derive(Debug, Clone, PartialEq)]
pub enum Token {
    /// data_ block header
    DataBlock(String),
    /// loop_ keyword
    Loop,
    /// _category.item data name
    DataName(String),
    /// Unquoted value
    Value(String),
    /// Single-quoted string 'value'
    SingleQuoted(String),
    /// Double-quoted string "value"
    DoubleQuoted(String),
    /// Semicolon-delimited text block
    TextField(String),
    /// Missing value (.)
    Missing,
    /// Unknown value (?)
    Unknown,
    /// End of input
    Eof,
}

impl Token {
    /// Get the string value of a token (for value tokens)
    #[allow(dead_code)]
    pub fn as_str(&self) -> Option<&str> {
        match self {
            Token::Value(s) | Token::SingleQuoted(s) | Token::DoubleQuoted(s) | Token::TextField(s) => {
                Some(s)
            }
            _ => None,
        }
    }

    /// Check if this token represents a missing or unknown value
    #[allow(dead_code)]
    pub fn is_missing(&self) -> bool {
        matches!(self, Token::Missing | Token::Unknown)
    }
}

/// Skip whitespace and comments
fn skip_ws_comments(input: &str) -> IResult<&str, ()> {
    let (input, _) = many0(alt((
        value((), multispace1),
        value((), preceded(char('#'), not_line_ending)),
    )))(input)?;
    Ok((input, ()))
}

/// Parse a data_ block header
fn parse_data_block(input: &str) -> IResult<&str, Token> {
    let (input, _) = tag("data_")(input)?;
    let (input, name) = take_while(|c: char| !c.is_whitespace())(input)?;
    Ok((input, Token::DataBlock(name.to_string())))
}

/// Parse loop_ keyword
fn parse_loop(input: &str) -> IResult<&str, Token> {
    let (input, _) = tag("loop_")(input)?;
    Ok((input, Token::Loop))
}

/// Parse a data name (_category.item)
fn parse_data_name(input: &str) -> IResult<&str, Token> {
    let (input, _) = char('_')(input)?;
    let (input, name) = take_while1(|c: char| !c.is_whitespace())(input)?;
    Ok((input, Token::DataName(format!("_{}", name))))
}

/// Parse a single-quoted string
fn parse_single_quoted(input: &str) -> IResult<&str, Token> {
    let (input, _) = char('\'')(input)?;
    let (input, content) = take_while(|c| c != '\'')(input)?;
    let (input, _) = char('\'')(input)?;
    Ok((input, Token::SingleQuoted(content.to_string())))
}

/// Parse a double-quoted string
fn parse_double_quoted(input: &str) -> IResult<&str, Token> {
    let (input, _) = char('"')(input)?;
    let (input, content) = take_while(|c| c != '"')(input)?;
    let (input, _) = char('"')(input)?;
    Ok((input, Token::DoubleQuoted(content.to_string())))
}

/// Parse a semicolon-delimited text field
fn parse_text_field(input: &str) -> IResult<&str, Token> {
    // Text fields start with ; at beginning of line and end with ; at beginning of line
    let (input, _) = char(';')(input)?;
    
    // Find the closing semicolon (at start of a line)
    let mut end_pos = 0;
    let mut found = false;
    for (i, c) in input.char_indices() {
        if c == '\n' {
            if input.get(i + 1..i + 2) == Some(";") {
                end_pos = i;
                found = true;
                break;
            }
        }
    }
    
    if !found {
        // If no closing ; found, take until end
        return Ok(("", Token::TextField(input.to_string())));
    }
    
    let content = &input[..end_pos];
    let rest = &input[end_pos + 2..]; // Skip \n;
    Ok((rest, Token::TextField(content.to_string())))
}

/// Parse missing value (.)
#[allow(dead_code)]
fn parse_missing(input: &str) -> IResult<&str, Token> {
    let (input, _) = char('.')(input)?;
    // Make sure it's not followed by more characters
    if input.is_empty() || input.starts_with(char::is_whitespace) {
        Ok((input, Token::Missing))
    } else {
        // It's part of a value
        Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Char,
        )))
    }
}

/// Parse unknown value (?)
#[allow(dead_code)]
fn parse_unknown(input: &str) -> IResult<&str, Token> {
    let (input, _) = char('?')(input)?;
    if input.is_empty() || input.starts_with(char::is_whitespace) {
        Ok((input, Token::Unknown))
    } else {
        Err(nom::Err::Error(nom::error::Error::new(
            input,
            nom::error::ErrorKind::Char,
        )))
    }
}

/// Parse an unquoted value
fn parse_unquoted_value(input: &str) -> IResult<&str, Token> {
    let (input, value) = take_while1(|c: char| !c.is_whitespace() && c != '#')(input)?;
    
    // Check for special values
    if value == "." {
        return Ok((input, Token::Missing));
    }
    if value == "?" {
        return Ok((input, Token::Unknown));
    }
    
    Ok((input, Token::Value(value.to_string())))
}

/// Parse a single token
pub fn parse_token(input: &str) -> IResult<&str, Token> {
    // Skip whitespace and comments first
    let (input, _) = skip_ws_comments(input)?;
    
    if input.is_empty() {
        return Ok((input, Token::Eof));
    }
    
    alt((
        parse_data_block,
        parse_loop,
        parse_data_name,
        parse_single_quoted,
        parse_double_quoted,
        parse_text_field,
        parse_unquoted_value,
    ))(input)
}

/// Tokenize a CIF file
pub fn tokenize(input: &str) -> Vec<Token> {
    let mut tokens = Vec::new();
    let mut remaining = input;
    
    loop {
        match parse_token(remaining) {
            Ok((_rest, Token::Eof)) => {
                tokens.push(Token::Eof);
                break;
            }
            Ok((rest, token)) => {
                tokens.push(token);
                remaining = rest;
            }
            Err(_) => {
                // Skip problematic character and continue
                if !remaining.is_empty() {
                    remaining = &remaining[1..];
                } else {
                    tokens.push(Token::Eof);
                    break;
                }
            }
        }
    }
    
    tokens
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_data_block() {
        let (rest, token) = parse_data_block("data_1ABC").unwrap();
        assert_eq!(token, Token::DataBlock("1ABC".to_string()));
        assert_eq!(rest, "");
    }

    #[test]
    fn test_parse_loop() {
        let (_rest, token) = parse_loop("loop_").unwrap();
        assert_eq!(token, Token::Loop);
    }

    #[test]
    fn test_parse_data_name() {
        let (_rest, token) = parse_data_name("_atom_site.id").unwrap();
        assert_eq!(token, Token::DataName("_atom_site.id".to_string()));
    }

    #[test]
    fn test_parse_quoted() {
        let (_, token) = parse_single_quoted("'hello world'").unwrap();
        assert_eq!(token, Token::SingleQuoted("hello world".to_string()));

        let (_, token) = parse_double_quoted("\"hello world\"").unwrap();
        assert_eq!(token, Token::DoubleQuoted("hello world".to_string()));
    }

    #[test]
    fn test_tokenize() {
        let cif = "data_TEST\n_cell.length_a 10.0\nloop_\n_atom_site.id\n1\n2\n";
        let tokens = tokenize(cif);
        
        assert!(matches!(tokens[0], Token::DataBlock(_)));
        assert!(matches!(tokens[1], Token::DataName(_)));
        assert!(matches!(tokens[2], Token::Value(_)));
        assert!(matches!(tokens[3], Token::Loop));
    }
}
