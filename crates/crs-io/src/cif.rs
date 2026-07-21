use crate::traits::{Parser, Writer};
use std::collections::HashMap;
use std::io::{BufRead, Error as IoError};
use thiserror::Error as ThisError;

mod parser;

#[derive(Debug, ThisError)]
pub enum CifError {
    #[error("I/O error: {0}")]
    IoError(#[from] IoError),
    #[error("Parse error: {0}")]
    ParseError(String),
    #[error("Empty file or string")]
    EmptyFile,
}

#[derive(Debug, Clone, PartialEq)]
pub enum CifValue {
    Primitive(String),
    List(Vec<CifValue>),
    Table(Vec<(String, CifValue)>),
}

impl CifValue {
    pub fn as_str(&self) -> Option<&str> {
        if let CifValue::Primitive(s) = self {
            Some(s.as_str())
        } else {
            None
        }
    }
}

#[derive(Debug, Clone)]
pub enum CifItem {
    Data(String, CifValue),
    Loop(Vec<String>, Vec<Vec<CifValue>>),
}

#[derive(Debug, Clone)]
pub struct CifBlock {
    pub name: String,
    pub items: Vec<CifItem>,
}

impl CifBlock {
    pub fn new(name: String, items: Vec<CifItem>) -> Self {
        Self { name, items }
    }

    /// Retrieve the value associated with a single data key (if it exists).
    pub fn get_data(&self, key: &str) -> Option<&CifValue> {
        let key = if !key.starts_with('_') {
            format!("_{}", key)
        } else {
            key.to_string()
        };
        let key_lower = key.to_lowercase();

        self.items.iter().find_map(|item| match item {
            CifItem::Data(k, v) if k.to_lowercase() == key_lower => Some(v),
            _ => None,
        })
    }
}

impl Parser<CifBlock> for CifBlock {
    type E = CifError;

    fn parse_from_reader<R: BufRead>(&self, reader: &mut R) -> Result<Option<CifBlock>, Self::E> {
        let mut buffer = String::new();
        // Since a CIF block can be very large, and `parse_from_reader` is supposed to
        // read until it completes a block. The easiest approach for now is to read the
        // entire stream into a string and parse it, but that doesn't work well if we want
        // to parse multiple blocks incrementally from a stream.
        // For a true stream parser, we'd read chunk by chunk. For simplicity (like XYZ),
        // let's read the whole file if this is the first call, or we can use winnow's
        // streaming capabilities, but `xyz` reads all into memory.

        reader.read_to_string(&mut buffer)?;
        if buffer.trim().is_empty() {
            return Ok(None);
        }

        let mut input: &str = buffer.as_str();

        match parser::parse_cif_block(&mut input) {
            Ok(Some(block)) => Ok(Some(block)),
            Ok(None) => Ok(None),
            Err(e) => Err(CifError::ParseError(e.to_string())),
        }
    }
}

impl TryFrom<&str> for CifBlock {
    type Error = CifError;

    fn try_from(s: &str) -> Result<Self, Self::Error> {
        let mut input: &str = s;
        match parser::parse_cif_block(&mut input) {
            Ok(Some(result)) => Ok(result),
            Ok(None) => Err(CifError::EmptyFile),
            Err(e) => Err(CifError::ParseError(e.to_string())),
        }
    }
}
