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

use crate::traits::{Parser, Writer};

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
        unimplemented!("Parsing XYZ frames is not yet implemented.");
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
