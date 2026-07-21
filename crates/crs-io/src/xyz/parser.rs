//! XYZ parser from the specification at: https://github.com/libAtoms/extxyz/blob/master/README.md
//!
//! Supports both standard XYZ (species + 3 coordinates) and extended XYZ
//! (arbitrary `Properties` columns, `Lattice`, and extra `info` key=value pairs).

use super::{XYZDataType, XYZFrame, XYZPropery};
use std::collections::HashMap;
use winnow::ascii::{digit1, line_ending, space0, space1};
use winnow::combinator::{alt, opt};
use winnow::error::{ContextError, ErrMode, ParserError};
use winnow::prelude::*;
use winnow::token::{none_of, take_while};

/// Convenience type alias matching winnow 1.0's `ModalResult`.
type PResult<T> = ModalResult<T>;

// ---------------------------------------------------------------------------
// Atom count (first line)
// ---------------------------------------------------------------------------

/// Parse the first line of a frame: optional leading whitespace, an integer,
/// optional trailing whitespace, then a newline.
fn parse_atom_count(input: &mut &str) -> PResult<usize> {
    let _ = space0.parse_next(input)?;
    let digits: &str = digit1.parse_next(input)?;
    let _ = space0.parse_next(input)?;
    let _ = line_ending.parse_next(input)?;
    digits
        .parse::<usize>()
        .map_err(|_| ErrMode::Backtrack(ContextError::from_input(input)))
}

// ---------------------------------------------------------------------------
// Comment-line value parsers
// ---------------------------------------------------------------------------

/// Parse a double-quoted value: `"..."`.
/// Handles escaped quotes (`\"`) inside the value.
fn parse_quoted_value(input: &mut &str) -> PResult<String> {
    let _ = '"'.parse_next(input)?;
    let mut result = String::new();
    loop {
        // Take characters that are neither `"` nor `\`
        let chunk: &str = take_while(0.., |c: char| c != '"' && c != '\\').parse_next(input)?;
        result.push_str(chunk);
        if input.starts_with("\\\"") {
            // Escaped quote
            let _ = "\\\"".parse_next(input)?;
            result.push('"');
        } else if input.starts_with('\\') {
            // Other escape — keep the backslash as-is
            let _ = '\\'.parse_next(input)?;
            result.push('\\');
        } else {
            break;
        }
    }
    let _ = '"'.parse_next(input)?;
    Ok(result)
}

/// Parse a brace-delimited value: `{...}` (no nesting).
fn parse_braced_value(input: &mut &str) -> PResult<String> {
    let _ = '{'.parse_next(input)?;
    let inner: &str = take_while(0.., |c: char| c != '}').parse_next(input)?;
    let _ = '}'.parse_next(input)?;
    Ok(inner.to_string())
}

/// Parse a bare (unquoted) value: a run of non-whitespace characters that
/// does not start with `"` or `{`.
fn parse_bare_value(input: &mut &str) -> PResult<String> {
    let v: &str = take_while(1.., |c: char| {
        !c.is_ascii_whitespace() && c != '"' && c != '{'
    })
    .parse_next(input)?;
    Ok(v.to_string())
}

/// Parse any value: quoted, braced, or bare.
fn parse_value(input: &mut &str) -> PResult<String> {
    alt((parse_quoted_value, parse_braced_value, parse_bare_value)).parse_next(input)
}

// ---------------------------------------------------------------------------
// Comment-line key parsers
// ---------------------------------------------------------------------------

/// Parse a key: one or more alphanumeric / `_` / `-` characters.
fn parse_key(input: &mut &str) -> PResult<String> {
    let k: &str =
        take_while(1.., |c: char| c.is_alphanumeric() || c == '_' || c == '-').parse_next(input)?;
    Ok(k.to_string())
}

/// Parse a single `key=value` pair **or** a standalone key (treated as boolean true).
fn parse_key_value(input: &mut &str) -> PResult<(String, String)> {
    let key = parse_key(input)?;
    // Optional spaces around `=`
    let _ = space0.parse_next(input)?;
    if input.starts_with('=') {
        let _ = '='.parse_next(input)?;
        let _ = space0.parse_next(input)?;
        let value = parse_value(input)?;
        Ok((key, value))
    } else {
        // Standalone key → boolean true
        Ok((key, "T".to_string()))
    }
}

/// Parse the entire comment line into a `HashMap<String, String>`.
/// Consumes everything up to (and including) the line ending.
fn parse_comment_line(input: &mut &str) -> PResult<HashMap<String, String>> {
    let _ = space0.parse_next(input)?;
    let mut map = HashMap::new();

    // Try to parse key=value pairs separated by whitespace until we hit a newline or EOF.
    loop {
        // Skip any leading spaces
        let _ = space0.parse_next(input)?;
        // Check for end of line
        if input.is_empty() || input.starts_with('\n') || input.starts_with("\r\n") {
            break;
        }
        let (key, value) = parse_key_value(input)?;
        map.insert(key, value);
    }

    // Consume the newline if present
    let _ = opt(line_ending).parse_next(input)?;
    Ok(map)
}

// ---------------------------------------------------------------------------
// Properties parsing
// ---------------------------------------------------------------------------

/// Default properties for standard XYZ: species (string, 1 col) + pos (real, 3 cols).
fn default_properties() -> Vec<XYZPropery> {
    vec![
        XYZPropery {
            name: "species".to_string(),
            data_type: XYZDataType::String,
            n_consecutive_cols: 1,
        },
        XYZPropery {
            name: "pos".to_string(),
            data_type: XYZDataType::Real,
            n_consecutive_cols: 3,
        },
    ]
}

/// Parse the `Properties` value string, e.g. `"species:S:1:pos:R:3:forces:R:3"`.
fn parse_properties(input: &mut &str) -> PResult<Vec<XYZPropery>> {
    let mut props = Vec::new();

    loop {
        // name
        let name: &str =
            take_while(1.., |c: char| c != ':' && !c.is_ascii_whitespace()).parse_next(input)?;
        let _ = ':'.parse_next(input)?;

        // type character
        let type_char: char = none_of([':']).parse_next(input)?;
        let data_type = match type_char {
            'S' | 's' => XYZDataType::String,
            'I' | 'i' => XYZDataType::Integer,
            'R' | 'r' => XYZDataType::Real,
            'L' | 'l' => XYZDataType::Logical,
            _ => return Err(ErrMode::Backtrack(ContextError::from_input(input))),
        };
        let _ = ':'.parse_next(input)?;

        // number of columns
        let n_str: &str = digit1.parse_next(input)?;
        let n_consecutive_cols: usize = n_str
            .parse()
            .map_err(|_| ErrMode::Backtrack(ContextError::from_input(input)))?;

        props.push(XYZPropery {
            name: name.to_string(),
            data_type,
            n_consecutive_cols,
        });

        // If next char is `:`, consume it and loop; otherwise done
        if input.starts_with(':') {
            let _ = ':'.parse_next(input)?;
        } else {
            break;
        }
    }

    Ok(props)
}

// ---------------------------------------------------------------------------
// Per-atom line parsing
// ---------------------------------------------------------------------------

/// Parse a single per-atom data line: whitespace-separated tokens.
/// Returns the tokens as a `Vec<String>`.
fn parse_atom_line(input: &mut &str, expected_cols: usize) -> PResult<Vec<String>> {
    let _ = space0.parse_next(input)?;
    let mut tokens = Vec::with_capacity(expected_cols);

    for i in 0..expected_cols {
        if i > 0 {
            let _ = space1.parse_next(input)?;
        }
        let tok: &str = take_while(1.., |c: char| !c.is_ascii_whitespace()).parse_next(input)?;
        tokens.push(tok.to_string());
    }

    let _ = space0.parse_next(input)?;
    let _ = opt(line_ending).parse_next(input)?;
    Ok(tokens)
}

// ---------------------------------------------------------------------------
// Full frame parsing
// ---------------------------------------------------------------------------

/// Parse one complete XYZ frame and return an `XYZFrame`.
///
/// Returns `Ok(None)` when the input is empty (EOF).
pub(super) fn parse_frame(input: &mut &str) -> PResult<Option<XYZFrame>> {
    // Skip any leading blank lines / whitespace
    let _ = take_while(0.., |c: char| c.is_ascii_whitespace()).parse_next(input)?;
    if input.is_empty() {
        return Ok(None);
    }

    // 1. Atom count
    let n_atoms = parse_atom_count(input)?;

    // 2. Comment line → info map
    let info = parse_comment_line(input)?;

    // 3. Determine row properties
    let row_properties = if let Some(props_str) = info.get("Properties") {
        let mut s: &str = props_str.as_str();
        parse_properties(&mut s).map_err(|_| ErrMode::Backtrack(ContextError::from_input(input)))?
    } else {
        default_properties()
    };

    let expected_cols: usize = row_properties.iter().map(|p| p.n_consecutive_cols).sum();

    // 4. Per-atom lines
    let mut rows = Vec::with_capacity(n_atoms);
    for _ in 0..n_atoms {
        let tokens = parse_atom_line(input, expected_cols)?;
        rows.push(tokens);
    }

    Ok(Some(XYZFrame::new(row_properties, info, rows)))
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // ---- atom count -------------------------------------------------------

    #[test]
    fn test_parse_atom_count() {
        let mut input = "3\n";
        assert_eq!(parse_atom_count(&mut input).unwrap(), 3);
    }

    #[test]
    fn test_parse_atom_count_with_spaces() {
        let mut input = "  42  \n";
        assert_eq!(parse_atom_count(&mut input).unwrap(), 42);
    }

    // ---- quoted value -----------------------------------------------------

    #[test]
    fn test_parse_quoted_value_simple() {
        let mut input = "\"hello world\"";
        assert_eq!(parse_quoted_value(&mut input).unwrap(), "hello world");
    }

    #[test]
    fn test_parse_quoted_value_with_escaped_quote() {
        let mut input = r#""he said \"hi\"""#;
        assert_eq!(parse_quoted_value(&mut input).unwrap(), r#"he said "hi""#);
    }

    // ---- braced value -----------------------------------------------------

    #[test]
    fn test_parse_braced_value() {
        let mut input = "{1.0 2.0 3.0}";
        assert_eq!(parse_braced_value(&mut input).unwrap(), "1.0 2.0 3.0");
    }

    // ---- bare value -------------------------------------------------------

    #[test]
    fn test_parse_bare_value() {
        let mut input = "42 rest";
        assert_eq!(parse_bare_value(&mut input).unwrap(), "42");
        assert_eq!(input, " rest");
    }

    // ---- key-value --------------------------------------------------------

    #[test]
    fn test_parse_key_value_pair() {
        let mut input = "energy=-3.14 next";
        let (k, v) = parse_key_value(&mut input).unwrap();
        assert_eq!(k, "energy");
        assert_eq!(v, "-3.14");
    }

    #[test]
    fn test_parse_key_value_standalone() {
        let mut input = "pbc\n";
        let (k, v) = parse_key_value(&mut input).unwrap();
        assert_eq!(k, "pbc");
        assert_eq!(v, "T");
    }

    #[test]
    fn test_parse_key_value_quoted() {
        let mut input = r#"Lattice="1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0""#;
        let (k, v) = parse_key_value(&mut input).unwrap();
        assert_eq!(k, "Lattice");
        assert_eq!(v, "1.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 1.0");
    }

    #[test]
    fn test_parse_key_value_with_spaces_around_equals() {
        let mut input = "energy = -3.14 next";
        let (k, v) = parse_key_value(&mut input).unwrap();
        assert_eq!(k, "energy");
        assert_eq!(v, "-3.14");
    }

    // ---- comment line -----------------------------------------------------

    #[test]
    fn test_parse_comment_line_empty() {
        let mut input = "\n";
        let map = parse_comment_line(&mut input).unwrap();
        assert!(map.is_empty());
    }

    #[test]
    fn test_parse_comment_line_whitespace_only() {
        let mut input = "   \n";
        let map = parse_comment_line(&mut input).unwrap();
        assert!(map.is_empty());
    }

    #[test]
    fn test_parse_comment_line_multiple_keys() {
        let mut input = "Properties=species:S:1:pos:R:3 Lattice=\"5.0 0.0 0.0 0.0 5.0 0.0 0.0 0.0 5.0\" energy=-3.14 pbc\n";
        let map = parse_comment_line(&mut input).unwrap();
        assert_eq!(map["Properties"], "species:S:1:pos:R:3");
        assert_eq!(map["Lattice"], "5.0 0.0 0.0 0.0 5.0 0.0 0.0 0.0 5.0");
        assert_eq!(map["energy"], "-3.14");
        assert_eq!(map["pbc"], "T");
    }

    // ---- properties -------------------------------------------------------

    #[test]
    fn test_parse_properties_basic() {
        let mut input = "species:S:1:pos:R:3";
        let props = parse_properties(&mut input).unwrap();
        assert_eq!(props.len(), 2);
        assert_eq!(props[0].name, "species");
        assert_eq!(props[0].data_type, XYZDataType::String);
        assert_eq!(props[0].n_consecutive_cols, 1);
        assert_eq!(props[1].name, "pos");
        assert_eq!(props[1].data_type, XYZDataType::Real);
        assert_eq!(props[1].n_consecutive_cols, 3);
    }

    #[test]
    fn test_parse_properties_with_forces() {
        let mut input = "species:S:1:pos:R:3:forces:R:3";
        let props = parse_properties(&mut input).unwrap();
        assert_eq!(props.len(), 3);
        assert_eq!(props[2].name, "forces");
        assert_eq!(props[2].data_type, XYZDataType::Real);
        assert_eq!(props[2].n_consecutive_cols, 3);
    }

    #[test]
    fn test_parse_properties_all_types() {
        let mut input = "label:S:1:index:I:1:pos:R:3:fixed:L:1";
        let props = parse_properties(&mut input).unwrap();
        assert_eq!(props.len(), 4);
        assert_eq!(props[0].data_type, XYZDataType::String);
        assert_eq!(props[1].data_type, XYZDataType::Integer);
        assert_eq!(props[2].data_type, XYZDataType::Real);
        assert_eq!(props[3].data_type, XYZDataType::Logical);
    }

    // ---- atom line --------------------------------------------------------

    #[test]
    fn test_parse_atom_line_basic() {
        let mut input = "Si 0.0 0.0 0.0\n";
        let tokens = parse_atom_line(&mut input, 4).unwrap();
        assert_eq!(tokens, vec!["Si", "0.0", "0.0", "0.0"]);
    }

    #[test]
    fn test_parse_atom_line_with_forces() {
        let mut input = "O 1.2 3.4 5.6 -0.1 -0.2 -0.3\n";
        let tokens = parse_atom_line(&mut input, 7).unwrap();
        assert_eq!(tokens.len(), 7);
        assert_eq!(tokens[0], "O");
        assert_eq!(tokens[6], "-0.3");
    }

    // ---- full frame -------------------------------------------------------

    #[test]
    fn test_parse_frame_plain_xyz() {
        let input_str = "\
3
water molecule
O  0.0000  0.0000  0.1173
H  0.0000  0.7572 -0.4692
H  0.0000 -0.7572 -0.4692
";
        let mut input = input_str;
        let frame = parse_frame(&mut input).unwrap().unwrap();

        // Default properties: species:S:1:pos:R:3
        assert_eq!(frame.row_properties().len(), 2);
        assert_eq!(frame.rows().len(), 3);
        assert_eq!(frame.rows()[0][0], "O");
        assert_eq!(frame.rows()[1][0], "H");
    }

    #[test]
    fn test_parse_frame_extxyz() {
        let input_str = "\
2
Properties=species:S:1:pos:R:3 Lattice=\"5.0 0.0 0.0 0.0 5.0 0.0 0.0 0.0 5.0\" energy=-3.14
Si 0.0 0.0 0.0
Si 2.5 2.5 2.5
";
        let mut input = input_str;
        let frame = parse_frame(&mut input).unwrap().unwrap();

        assert_eq!(frame.rows().len(), 2);
        assert_eq!(
            frame.get_info("Lattice").unwrap(),
            "5.0 0.0 0.0 0.0 5.0 0.0 0.0 0.0 5.0"
        );
        assert_eq!(frame.get_info("energy").unwrap(), "-3.14");
        assert_eq!(frame.row_properties().len(), 2);
    }

    #[test]
    fn test_parse_frame_extxyz_with_forces() {
        let input_str = "\
2
Properties=species:S:1:pos:R:3:forces:R:3 energy=0.5
Si 0.0 0.0 0.0 0.1 0.2 0.3
Si 2.5 2.5 2.5 -0.1 -0.2 -0.3
";
        let mut input = input_str;
        let frame = parse_frame(&mut input).unwrap().unwrap();

        assert_eq!(frame.row_properties().len(), 3);
        let expected_cols: usize = frame
            .row_properties()
            .iter()
            .map(|p| p.n_consecutive_cols)
            .sum();
        assert_eq!(expected_cols, 7);
        assert_eq!(frame.rows()[0].len(), 7);
        assert_eq!(frame.rows()[0][4], "0.1"); // first force component
    }

    #[test]
    fn test_parse_frame_eof_returns_none() {
        let mut input = "";
        assert!(parse_frame(&mut input).unwrap().is_none());
    }

    #[test]
    fn test_parse_frame_whitespace_only_eof() {
        let mut input = "  \n  \n  ";
        assert!(parse_frame(&mut input).unwrap().is_none());
    }

    #[test]
    fn test_parse_multiple_frames() {
        let input_str = "\
2
first frame
H 0.0 0.0 0.0
H 1.0 0.0 0.0
3
second frame
O 0.0 0.0 0.0
H 0.5 0.5 0.0
H -0.5 0.5 0.0
";
        let mut input = input_str;

        let frame1 = parse_frame(&mut input).unwrap().unwrap();
        assert_eq!(frame1.rows().len(), 2);

        let frame2 = parse_frame(&mut input).unwrap().unwrap();
        assert_eq!(frame2.rows().len(), 3);
        assert_eq!(frame2.rows()[0][0], "O");

        // No more frames
        assert!(parse_frame(&mut input).unwrap().is_none());
    }

    #[test]
    fn test_parse_frame_braced_lattice() {
        let input_str = "\
1
Lattice={5.0 0.0 0.0 0.0 5.0 0.0 0.0 0.0 5.0} Properties=species:S:1:pos:R:3
Si 0.0 0.0 0.0
";
        let mut input = input_str;
        let frame = parse_frame(&mut input).unwrap().unwrap();
        assert_eq!(
            frame.get_info("Lattice").unwrap(),
            "5.0 0.0 0.0 0.0 5.0 0.0 0.0 0.0 5.0"
        );
    }

    #[test]
    fn test_parse_comment_standalone_keys_mixed() {
        let mut input = "pbc energy=-1.5 periodic\n";
        let map = parse_comment_line(&mut input).unwrap();
        assert_eq!(map["pbc"], "T");
        assert_eq!(map["energy"], "-1.5");
        assert_eq!(map["periodic"], "T");
    }

    #[test]
    fn test_parse_properties_lowercase_types() {
        let mut input = "species:s:1:pos:r:3";
        let props = parse_properties(&mut input).unwrap();
        assert_eq!(props[0].data_type, XYZDataType::String);
        assert_eq!(props[1].data_type, XYZDataType::Real);
    }

    #[test]
    fn test_parse_atom_line_leading_spaces() {
        let mut input = "  C 1.0 2.0 3.0\n";
        let tokens = parse_atom_line(&mut input, 4).unwrap();
        assert_eq!(tokens, vec!["C", "1.0", "2.0", "3.0"]);
    }

    #[test]
    fn test_parse_key_with_hyphens() {
        let mut input = "config-type=bulk\n";
        let (k, v) = parse_key_value(&mut input).unwrap();
        assert_eq!(k, "config-type");
        assert_eq!(v, "bulk");
    }
}
