use super::{CifBlock, CifItem, CifValue};
use winnow::ascii::{multispace0, multispace1, space0};
use winnow::combinator::{alt, delimited, opt, preceded, repeat, separated, terminated};
use winnow::error::{ContextError, ErrMode};
use winnow::prelude::*;
use winnow::token::{any, none_of, take, take_until, take_while};

type PResult<T> = ModalResult<T>;

// ---------------------------------------------------------------------------
// Whitespace and Comments
// ---------------------------------------------------------------------------

/// Parses a comment: `#` until the end of the line.
fn comment(input: &mut &str) -> PResult<()> {
    let _ = '#'.parse_next(input)?;
    let _ = take_while(0.., |c: char| c != '\n' && c != '\r').parse_next(input)?;
    Ok(())
}

/// Parses any combination of whitespace and comments.
fn wspace0(input: &mut &str) -> PResult<()> {
    let _: () = repeat(0.., alt((multispace1.void(), comment.void()))).parse_next(input)?;
    Ok(())
}

fn wspace1(input: &mut &str) -> PResult<()> {
    let _: () = repeat(1.., alt((multispace1.void(), comment.void()))).parse_next(input)?;
    Ok(())
}

// ---------------------------------------------------------------------------
// Primitive Values
// ---------------------------------------------------------------------------

/// Parse a single-quoted string `'...'`.
fn single_quoted(input: &mut &str) -> PResult<String> {
    let _ = '\''.parse_next(input)?;
    let inner: &str = take_while(0.., |c: char| c != '\'').parse_next(input)?;
    let _ = '\''.parse_next(input)?;
    Ok(inner.to_string())
}

/// Parse a double-quoted string `"..."`.
fn double_quoted(input: &mut &str) -> PResult<String> {
    let _ = '"'.parse_next(input)?;
    let inner: &str = take_while(0.., |c: char| c != '"').parse_next(input)?;
    let _ = '"'.parse_next(input)?;
    Ok(inner.to_string())
}

/// Parse a triple-single-quoted string `'''...'''`.
fn triple_single_quoted(input: &mut &str) -> PResult<String> {
    let _ = "'''".parse_next(input)?;
    let inner: &str = take_until(0.., "'''").parse_next(input)?;
    let _ = "'''".parse_next(input)?;
    Ok(inner.to_string())
}

/// Parse a triple-double-quoted string `"""..."""`.
fn triple_double_quoted(input: &mut &str) -> PResult<String> {
    let _ = "\"\"\"".parse_next(input)?;
    let inner: &str = take_until(0.., "\"\"\"").parse_next(input)?;
    let _ = "\"\"\"".parse_next(input)?;
    Ok(inner.to_string())
}

/// Parse a CIF text block.
/// It starts with `;` as the first character of a line and ends with `;` as the first character of a line.
/// Because whitespace is consumed before values, the `;` is guaranteed to be at the start of our input.
fn text_block(input: &mut &str) -> PResult<String> {
    let _ = ';'.parse_next(input)?;

    let mut end_idx = input.find("\n;");
    if let Some(rn_idx) = input.find("\r\n;") {
        if end_idx.is_none() || rn_idx < end_idx.unwrap() {
            end_idx = Some(rn_idx);
        }
    }

    if let Some(idx) = end_idx {
        let content = &input[..idx];
        let term_len = if input[idx..].starts_with("\r\n;") {
            3
        } else {
            2
        };
        let _ = take(idx + term_len).parse_next(input)?;
        Ok(content.to_string())
    } else {
        let content = *input;
        let _ = take(input.len()).parse_next(input)?;
        Ok(content.to_string())
    }
}

/// Parse an unquoted string.
fn unquoted_string(input: &mut &str) -> PResult<String> {
    let first_char: char =
        none_of(['_', '\'', '"', '#', '$', ';', '[', ']', '{', '}']).parse_next(input)?;
    let rest: &str = take_while(0.., |c: char| {
        !c.is_ascii_whitespace() && c != '[' && c != ']' && c != '{' && c != '}'
    })
    .parse_next(input)?;

    let mut s = String::new();
    s.push(first_char);
    s.push_str(rest);

    // Check for reserved keywords
    let lower = s.to_lowercase();
    if lower.starts_with("data_")
        || lower.starts_with("save_")
        || lower == "loop_"
        || lower == "stop_"
        || lower == "global_"
    {
        return Err(ErrMode::Backtrack(ContextError::from_input(input)));
    }

    Ok(s)
}

fn primitive_value(input: &mut &str) -> PResult<CifValue> {
    alt((
        triple_double_quoted.map(CifValue::Primitive),
        triple_single_quoted.map(CifValue::Primitive),
        double_quoted.map(CifValue::Primitive),
        single_quoted.map(CifValue::Primitive),
        text_block.map(CifValue::Primitive),
        unquoted_string.map(CifValue::Primitive),
    ))
    .parse_next(input)
}

// ---------------------------------------------------------------------------
// Complex Values (Lists and Tables)
// ---------------------------------------------------------------------------

fn cif_list(input: &mut &str) -> PResult<CifValue> {
    let _ = '['.parse_next(input)?;
    let _ = wspace0.parse_next(input)?;
    let elements: Vec<CifValue> = separated(0.., cif_value, wspace1).parse_next(input)?;
    let _ = wspace0.parse_next(input)?;
    let _ = ']'.parse_next(input)?;
    Ok(CifValue::List(elements))
}

fn table_key(input: &mut &str) -> PResult<String> {
    alt((single_quoted, double_quoted, unquoted_string)).parse_next(input)
}

fn table_entry(input: &mut &str) -> PResult<(String, CifValue)> {
    let key = table_key(input)?;
    let _ = space0.parse_next(input)?;
    let _ = ':'.parse_next(input)?;
    let _ = wspace0.parse_next(input)?;
    let value = cif_value(input)?;
    Ok((key, value))
}

fn cif_table(input: &mut &str) -> PResult<CifValue> {
    let _ = '{'.parse_next(input)?;
    let _ = wspace0.parse_next(input)?;
    let entries: Vec<(String, CifValue)> =
        separated(0.., table_entry, wspace1).parse_next(input)?;
    let _ = wspace0.parse_next(input)?;
    let _ = '}'.parse_next(input)?;
    Ok(CifValue::Table(entries))
}

fn cif_value(input: &mut &str) -> PResult<CifValue> {
    alt((cif_list, cif_table, primitive_value)).parse_next(input)
}

// ---------------------------------------------------------------------------
// CIF Items and Loops
// ---------------------------------------------------------------------------

fn data_name(input: &mut &str) -> PResult<String> {
    let _ = '_'.parse_next(input)?;
    let rest: &str = take_while(1.., |c: char| !c.is_ascii_whitespace()).parse_next(input)?;
    Ok(format!("_{}", rest))
}

fn single_data_item(input: &mut &str) -> PResult<CifItem> {
    let name = data_name(input)?;
    let _ = wspace1.parse_next(input)?;
    let value = cif_value(input)?;
    Ok(CifItem::Data(name, value))
}

fn loop_item(input: &mut &str) -> PResult<CifItem> {
    let _ = alt(("loop_", "LOOP_", "Loop_")).parse_next(input)?;
    let _ = wspace1.parse_next(input)?;

    // Parse the data names
    let names: Vec<String> = separated(1.., data_name, wspace1).parse_next(input)?;
    let _ = wspace1.parse_next(input)?;

    // Parse the values
    let values: Vec<CifValue> = separated(names.len().., cif_value, wspace1).parse_next(input)?;

    // Chunk values by the number of names
    let n_cols = names.len();
    if values.len() % n_cols != 0 {
        return Err(ErrMode::Backtrack(ContextError::from_input(input)));
    }

    let rows: Vec<Vec<CifValue>> = values.chunks(n_cols).map(|chunk| chunk.to_vec()).collect();

    Ok(CifItem::Loop(names, rows))
}

// ---------------------------------------------------------------------------
// Save Frames (Ignored)
// ---------------------------------------------------------------------------

fn save_frame(input: &mut &str) -> PResult<()> {
    // Starts with `save_something`
    let _ = alt(("save_", "SAVE_", "Save_")).parse_next(input)?;
    let _ = take_while(1.., |c: char| !c.is_ascii_whitespace()).parse_next(input)?;
    let _ = wspace1.parse_next(input)?;

    // Parse internal items but discard them
    loop {
        if input.to_lowercase().starts_with("save_") {
            let _ = alt(("save_", "SAVE_", "Save_")).parse_next(input)?;
            let _ = take_while(0.., |c: char| !c.is_ascii_whitespace()).parse_next(input)?;
            break;
        }

        let _ = alt((single_data_item.void(), loop_item.void())).parse_next(input)?;
        let _ = wspace0.parse_next(input)?;
    }

    Ok(())
}

// ---------------------------------------------------------------------------
// Blocks
// ---------------------------------------------------------------------------

/// Parse an entire data block: `data_name` then items.
pub(super) fn parse_cif_block(input: &mut &str) -> PResult<Option<CifBlock>> {
    let _ = wspace0.parse_next(input)?;

    if input.is_empty() {
        return Ok(None);
    }

    // Read `data_` prefix
    let _ = alt(("data_", "DATA_", "Data_")).parse_next(input)?;
    let name: &str = take_while(1.., |c: char| !c.is_ascii_whitespace()).parse_next(input)?;

    let mut items = Vec::new();

    loop {
        let _ = wspace0.parse_next(input)?;

        if input.is_empty() || input.to_lowercase().starts_with("data_") {
            break;
        }

        if input.to_lowercase().starts_with("save_") {
            let _ = save_frame(input)?;
        } else if input.starts_with('_') {
            let item = single_data_item(input)?;
            items.push(item);
        } else if input.to_lowercase().starts_with("loop_") {
            let item = loop_item(input)?;
            items.push(item);
        } else {
            // Not recognizable as a data item, loop, or save frame.
            return Err(ErrMode::Backtrack(ContextError::from_input(input)));
        }
    }

    Ok(Some(CifBlock::new(name.to_string(), items)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_unquoted_string() {
        let mut input = "hello_world 123";
        let res = unquoted_string(&mut input).unwrap();
        assert_eq!(res, "hello_world");
        assert_eq!(input, " 123");
    }

    #[test]
    fn test_quoted_string() {
        let mut input = "'hello world' 123";
        let res = single_quoted(&mut input).unwrap();
        assert_eq!(res, "hello world");
        assert_eq!(input, " 123");
    }

    #[test]
    fn test_text_block() {
        let mut input = ";This is a
text block
; 123";
        let res = text_block(&mut input).unwrap();
        assert_eq!(res, "This is a\ntext block");
        assert_eq!(input, " 123");
    }

    #[test]
    fn test_single_data_item() {
        let mut input = "_cell_length_a 10.5";
        let item = single_data_item(&mut input).unwrap();
        if let CifItem::Data(name, CifValue::Primitive(val)) = item {
            assert_eq!(name, "_cell_length_a");
            assert_eq!(val, "10.5");
        } else {
            panic!("Expected Data item");
        }
    }

    #[test]
    fn test_loop_item() {
        let mut input = "loop_
            _atom_site_label
            _atom_site_fract_x
            _atom_site_fract_y
            _atom_site_fract_z
            O1 0.0 0.0 0.0
            O2 0.5 0.5 0.5
        ";
        let item = loop_item(&mut input).unwrap();
        if let CifItem::Loop(names, rows) = item {
            assert_eq!(names.len(), 4);
            assert_eq!(names[0], "_atom_site_label");
            assert_eq!(rows.len(), 2);
            assert_eq!(rows[0].len(), 4);
            assert_eq!(rows[0][0].as_str(), Some("O1"));
        } else {
            panic!("Expected Loop item");
        }
    }

    #[test]
    fn test_cif_block() {
        let mut input = "data_global
            _cell_length_a 10.5
            _cell_length_b 10.5
            _cell_length_c 10.5
            loop_
            _atom_site_label
            _atom_site_fract_x
            _atom_site_fract_y
            _atom_site_fract_z
            O1 0.0 0.0 0.0
            O2 0.5 0.5 0.5
        ";
        let block = parse_cif_block(&mut input).unwrap().unwrap();
        assert_eq!(block.name, "global");
        assert_eq!(block.items.len(), 4);
    }
}
