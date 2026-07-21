use crate::traits::{Parser, Writer};
use crs_core::au::AsymmetricUnit;
use crs_core::crystal::Crystal;
use crs_core::data::SPECIES_NAMES;
use crs_core::molecule::Molecule;
use crs_core::spacegroup::SpaceGroup;
use crs_core::traits::{AtomicData, CartAtomicData, CellData};
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

impl Writer for XYZFrame {
    fn write_to_writer<W: std::io::prelude::Write>(&self, writer: &mut W) -> IoResult<()> {
        writeln!(writer, "{}", self.rows.len())?;

        let mut comment_parts = Vec::new();
        for (k, v) in &self.info {
            if v.contains(' ') || v.contains('{') || v.contains('}') || v.contains('"') {
                let escaped_v = v.replace("\"", "\\\"");
                comment_parts.push(format!("{}=\"{}\"", k, escaped_v));
            } else if v == "T" && k != "Properties" {
                comment_parts.push(k.clone());
            } else {
                comment_parts.push(format!("{}={}", k, v));
            }
        }
        writeln!(writer, "{}", comment_parts.join(" "))?;

        for row in &self.rows {
            writeln!(writer, "{}", row.join(" "))?;
        }

        Ok(())
    }
}

impl TryFrom<&str> for XYZFrame {
    type Error = XYZError;

    fn try_from(s: &str) -> Result<Self, Self::Error> {
        let mut input: &str = s;
        let presult =
            parser::parse_frame(&mut input).map_err(|e| XYZError::ParseError(e.to_string()))?;
        if let Some(result) = presult {
            Ok(result)
        } else {
            Err(XYZError::ParseError("Failed to parse XYZ frame".into()))
        }
    }
}

impl TryFrom<XYZFrame> for Molecule {
    type Error = XYZError;

    fn try_from(frame: XYZFrame) -> Result<Self, Self::Error> {
        let mut col = 0;
        let mut species_col = None;
        let mut pos_col = None;

        for prop in frame.row_properties() {
            if prop.name == "species" {
                species_col = Some(col);
            } else if prop.name == "pos" {
                pos_col = Some(col);
            }
            col += prop.n_consecutive_cols;
        }

        let species_col =
            species_col.ok_or_else(|| XYZError::ParseError("Missing 'species' property".into()))?;
        let pos_col =
            pos_col.ok_or_else(|| XYZError::ParseError("Missing 'pos' property".into()))?;

        let mut atomic_nums = Vec::with_capacity(frame.rows().len());
        let mut coords = Vec::with_capacity(frame.rows().len() * 3);

        for row in frame.rows() {
            if row.len() <= pos_col + 2 {
                return Err(XYZError::ParseError(
                    "Row has fewer columns than expected".into(),
                ));
            }

            let species_str = &row[species_col];
            let atomic_num = if let Ok(num) = species_str.parse::<u8>() {
                num
            } else {
                SPECIES_NAMES
                    .iter()
                    .position(|&s| s.eq_ignore_ascii_case(species_str))
                    .map(|i| i as u8)
                    .ok_or_else(|| {
                        XYZError::ParseError(format!("Unknown species '{}'", species_str))
                    })?
            };
            atomic_nums.push(atomic_num);

            let x: f64 = row[pos_col]
                .parse()
                .map_err(|_| XYZError::ParseError("Invalid x coordinate".into()))?;
            let y: f64 = row[pos_col + 1]
                .parse()
                .map_err(|_| XYZError::ParseError("Invalid y coordinate".into()))?;
            let z: f64 = row[pos_col + 2]
                .parse()
                .map_err(|_| XYZError::ParseError("Invalid z coordinate".into()))?;
            coords.extend_from_slice(&[x, y, z]);
        }

        let cartesian_coords = MatrixXx3::from_row_slice(&coords);
        Ok(Molecule::new(atomic_nums, cartesian_coords))
    }
}

impl TryFrom<XYZFrame> for Crystal {
    type Error = XYZError;

    fn try_from(frame: XYZFrame) -> Result<Self, Self::Error> {
        let lattice_str = frame.get_info("Lattice").ok_or(XYZError::MissingLattice)?;
        let lattice_parts: Vec<&str> = lattice_str.split_whitespace().collect();
        if lattice_parts.len() != 9 {
            return Err(XYZError::InvalidLattice);
        }

        let mut lattice_matrix = [0.0; 9];
        for (i, part) in lattice_parts.iter().enumerate() {
            lattice_matrix[i] = part
                .parse()
                .map_err(|_| XYZError::ParseError("Invalid Lattice value".into()))?;
        }

        let cell_matrix = Matrix3::from_row_slice(&lattice_matrix);
        let lengths = cell_matrix.lengths();
        let angles = cell_matrix.angles();
        let unit_cell = UnitCell::new(lengths, angles);

        let mol = Molecule::try_from(frame)?;
        let frac_coords = unit_cell.to_fractional(mol.cartesian_coords());

        let space_group = SpaceGroup::default_setting(1);
        let asymmetric_unit = AsymmetricUnit::new(mol.atomic_nums().to_vec(), frac_coords);

        Ok(Crystal::new(unit_cell, space_group, asymmetric_unit))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_try_from_molecule() {
        let rows = vec![
            vec![
                "O".to_string(),
                "0.0".to_string(),
                "0.0".to_string(),
                "0.1173".to_string(),
            ],
            vec![
                "H".to_string(),
                "0.0".to_string(),
                "0.7572".to_string(),
                "-0.4692".to_string(),
            ],
            vec![
                "H".to_string(),
                "0.0".to_string(),
                "-0.7572".to_string(),
                "-0.4692".to_string(),
            ],
        ];
        let properties = vec![
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
        ];

        let frame = XYZFrame::new(properties, HashMap::new(), rows);

        let mol = Molecule::try_from(frame).unwrap();
        assert_eq!(mol.atomic_nums(), &[8, 1, 1]); // O is 8, H is 1
        assert_eq!(mol.cartesian_coords().nrows(), 3);
    }

    #[test]
    fn test_try_from_crystal() {
        let mut info = HashMap::new();
        info.insert(
            "Lattice".to_string(),
            "5.0 0.0 0.0 0.0 5.0 0.0 0.0 0.0 5.0".to_string(),
        );

        let rows = vec![
            vec![
                "C".to_string(),
                "0.0".to_string(),
                "0.0".to_string(),
                "0.0".to_string(),
            ],
            vec![
                "C".to_string(),
                "2.5".to_string(),
                "2.5".to_string(),
                "2.5".to_string(),
            ],
        ];
        let properties = vec![
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
        ];
        let frame = XYZFrame::new(properties, info, rows);

        let crystal = Crystal::try_from(frame).unwrap();
        assert_eq!(crystal.atomic_nums(), &[6, 6]); // C is 6
        assert_eq!(crystal.lengths(), [5.0, 5.0, 5.0]);
    }

    #[test]
    fn test_xyz_writer() {
        let mut info = HashMap::new();
        info.insert(
            "Lattice".to_string(),
            "5.0 0.0 0.0 0.0 5.0 0.0 0.0 0.0 5.0".to_string(),
        );
        info.insert("energy".to_string(), "-3.14".to_string());

        let rows = vec![
            vec![
                "Si".to_string(),
                "0.0".to_string(),
                "0.0".to_string(),
                "0.0".to_string(),
            ],
            vec![
                "Si".to_string(),
                "2.5".to_string(),
                "2.5".to_string(),
                "2.5".to_string(),
            ],
        ];
        let frame = XYZFrame::new(vec![], info, rows);

        let mut out = Vec::new();
        frame.write_to_writer(&mut out).unwrap();
        let s = String::from_utf8(out).unwrap();

        assert!(s.starts_with("2\n"));
        let lines: Vec<&str> = s.lines().collect();
        assert_eq!(lines.len(), 4);
        assert!(lines[1].contains("Lattice=\"5.0 0.0 0.0 0.0 5.0 0.0 0.0 0.0 5.0\""));
        assert!(lines[1].contains("energy=-3.14"));
        assert_eq!(lines[2], "Si 0.0 0.0 0.0");
        assert_eq!(lines[3], "Si 2.5 2.5 2.5");
    }
}
