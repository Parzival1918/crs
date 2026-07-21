use crate::traits::{Parser, Writer};
use crs_core::au::AsymmetricUnit;
use crs_core::crystal::Crystal;
use crs_core::data::SPECIES_NAMES;
use crs_core::molecule::Molecule;
use crs_core::spacegroup::SpaceGroup;
use crs_core::traits::CellData;
use crs_core::unitcell::UnitCell;
use nalgebra::{Matrix3, MatrixXx3};
use std::collections::HashMap;
use std::{
    fmt::Write,
    hash::Hash,
    io::{BufRead, Error, ErrorKind, Result as IoResult},
};
use thiserror::Error as ThisError;

mod parser;

#[derive(Debug, ThisError)]
pub enum XYZError {
    #[error("I/O error: {0}")]
    IoError(#[from] std::io::Error),
    #[error("Invalid data")]
    InvalidData,
    #[error("Missing comment line")]
    MissingComment,
    #[error("Missing Lattice in comment")]
    MissingLattice,
    #[error("Unterminated Lattice string")]
    UnterminatedLattice,
    #[error("Lattice must have 9 elements")]
    InvalidLattice,
    #[error("Parse error: {0}")]
    ParseError(String),
}

#[derive(Debug, Clone, PartialEq)]
pub enum XYZDataType {
    String,
    Integer,
    Real,
    Logical,
}

#[derive(Debug, Clone)]
pub struct XYZPropery {
    name: String,
    data_type: XYZDataType,
    n_consecutive_cols: usize,
}

/// Structure for holding the data of a single XYZ frame.
/// It can parse both standard XYZ and extended XYZ formats.
#[derive(Debug, Clone)]
pub struct XYZFrame {
    row_properties: Vec<XYZPropery>,
    info: HashMap<String, String>,
    rows: Vec<Vec<String>>,
}

impl XYZFrame {
    pub fn new(
        row_properties: Vec<XYZPropery>,
        info: HashMap<String, String>,
        rows: Vec<Vec<String>>,
    ) -> Self {
        Self {
            row_properties,
            info,
            rows,
        }
    }

    pub fn row_properties(&self) -> &Vec<XYZPropery> {
        &self.row_properties
    }

    pub fn info(&self) -> &HashMap<String, String> {
        &self.info
    }

    pub fn rows(&self) -> &Vec<Vec<String>> {
        &self.rows
    }

    pub fn get_info(&self, key: &str) -> Option<&String> {
        self.info.get(key)
    }

    pub fn get_row_property(&self, name: &str) -> Option<&XYZPropery> {
        self.row_properties.iter().find(|prop| prop.name == name)
    }
}

impl Parser<XYZFrame> for XYZFrame {
    type E = XYZError;

    fn parse_from_reader<R: BufRead>(&self, reader: &mut R) -> Result<Option<XYZFrame>, Self::E> {
        // 1. Read the atom-count line (skip blank lines)
        let mut count_line = String::new();
        loop {
            count_line.clear();
            let bytes_read = reader.read_line(&mut count_line)?;
            if bytes_read == 0 {
                return Ok(None); // EOF
            }
            if !count_line.trim().is_empty() {
                break;
            }
        }

        // 2. Read the comment line
        let mut comment_line = String::new();
        reader.read_line(&mut comment_line)?;

        // 3. Parse atom count to know how many atom lines to read
        let n_atoms: usize = count_line.trim().parse().map_err(|_| {
            XYZError::ParseError(format!("invalid atom count: {:?}", count_line.trim()))
        })?;

        // 4. Read N atom lines
        let mut atom_lines = String::new();
        for _ in 0..n_atoms {
            reader.read_line(&mut atom_lines)?;
        }

        // 5. Assemble the full frame text and parse with winnow
        let mut frame_text = count_line;
        frame_text.push_str(&comment_line);
        frame_text.push_str(&atom_lines);

        let mut input: &str = frame_text.as_str();
        parser::parse_frame(&mut input).map_err(|e| XYZError::ParseError(e.to_string()))
    }
}

impl Writer<Molecule> for XYZFrame {
    fn write_to_writer<W: std::io::prelude::Write>(
        &self,
        data: &Molecule,
        writer: &mut W,
    ) -> IoResult<()> {
        unimplemented!("Writing XYZ frames is not yet implemented.");
    }
}

impl Writer<Crystal> for XYZFrame {
    fn write_to_writer<W: std::io::prelude::Write>(
        &self,
        data: &Crystal,
        writer: &mut W,
    ) -> IoResult<()> {
        unimplemented!("Writing XYZ frames is not yet implemented.");
    }
}
