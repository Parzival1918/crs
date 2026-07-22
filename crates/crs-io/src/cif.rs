use crate::traits::{Parser, Writer};
use crs_core::au::AsymmetricUnit;
use crs_core::crystal::Crystal;
use crs_core::data::SPECIES_NAMES;
use crs_core::molecule::Molecule;
use crs_core::spacegroup::SpaceGroup;
use crs_core::traits::CellData;
use crs_core::unitcell::UnitCell;
use nalgebra::MatrixXx3;
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

    fn parse_from_reader<R: BufRead>(reader: &mut R) -> Result<Option<CifBlock>, Self::E> {
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

fn parse_cif_float(s: &str) -> Result<f64, CifError> {
    let s = s.split('(').next().unwrap_or(s);
    s.parse::<f64>()
        .map_err(|_| CifError::ParseError(format!("Invalid float value: {}", s)))
}

fn extract_atomic_num(symbol: &str) -> Result<u8, CifError> {
    let alpha_part: String = symbol.chars().take_while(|c| c.is_alphabetic()).collect();
    SPECIES_NAMES
        .iter()
        .position(|&s| s.eq_ignore_ascii_case(&alpha_part))
        .map(|i| i as u8)
        .ok_or_else(|| CifError::ParseError(format!("Unknown species '{}'", symbol)))
}

fn parse_unit_cell(block: &CifBlock) -> Result<UnitCell, CifError> {
    let a = parse_cif_float(
        block
            .get_data("_cell_length_a")
            .and_then(|v| v.as_str())
            .unwrap_or("0"),
    )?;
    let b = parse_cif_float(
        block
            .get_data("_cell_length_b")
            .and_then(|v| v.as_str())
            .unwrap_or("0"),
    )?;
    let c = parse_cif_float(
        block
            .get_data("_cell_length_c")
            .and_then(|v| v.as_str())
            .unwrap_or("0"),
    )?;
    let alpha = parse_cif_float(
        block
            .get_data("_cell_angle_alpha")
            .and_then(|v| v.as_str())
            .unwrap_or("90"),
    )?;
    let beta = parse_cif_float(
        block
            .get_data("_cell_angle_beta")
            .and_then(|v| v.as_str())
            .unwrap_or("90"),
    )?;
    let gamma = parse_cif_float(
        block
            .get_data("_cell_angle_gamma")
            .and_then(|v| v.as_str())
            .unwrap_or("90"),
    )?;
    Ok(UnitCell::new([a, b, c], [alpha, beta, gamma]))
}

fn find_atom_loop<'a>(block: &'a CifBlock) -> Option<(&'a Vec<String>, &'a Vec<Vec<CifValue>>)> {
    block.items.iter().find_map(|item| {
        if let CifItem::Loop(headers, rows) = item {
            if headers.iter().any(|h| {
                h.eq_ignore_ascii_case("_atom_site_label")
                    || h.eq_ignore_ascii_case("_atom_site_type_symbol")
            }) {
                return Some((headers, rows));
            }
        }
        None
    })
}

fn parse_space_group(block: &CifBlock) -> Result<SpaceGroup, CifError> {
    if let Some(val) = block
        .get_data("_space_group_IT_number")
        .or_else(|| block.get_data("_symmetry_Int_Tables_number"))
    {
        let s = val
            .as_str()
            .ok_or_else(|| CifError::ParseError("Invalid space group number".into()))?;
        let num = s
            .parse::<u16>()
            .map_err(|_| CifError::ParseError("Invalid space group number".into()))?;
        Ok(SpaceGroup::default_setting(num))
    } else {
        Ok(SpaceGroup::default_setting(1))
    }
}

impl TryFrom<CifBlock> for Crystal {
    type Error = CifError;

    fn try_from(block: CifBlock) -> Result<Self, Self::Error> {
        let unit_cell = parse_unit_cell(&block)?;
        let space_group = parse_space_group(&block)?;

        let (headers, rows) = find_atom_loop(&block)
            .ok_or_else(|| CifError::ParseError("No atom loop found".into()))?;

        let symbol_idx = headers
            .iter()
            .position(|h| h.eq_ignore_ascii_case("_atom_site_type_symbol"))
            .or_else(|| {
                headers
                    .iter()
                    .position(|h| h.eq_ignore_ascii_case("_atom_site_label"))
            })
            .ok_or_else(|| CifError::ParseError("Missing atom symbol column".into()))?;
        let x_idx = headers
            .iter()
            .position(|h| h.eq_ignore_ascii_case("_atom_site_fract_x"))
            .ok_or_else(|| CifError::ParseError("Missing fract x column".into()))?;
        let y_idx = headers
            .iter()
            .position(|h| h.eq_ignore_ascii_case("_atom_site_fract_y"))
            .ok_or_else(|| CifError::ParseError("Missing fract y column".into()))?;
        let z_idx = headers
            .iter()
            .position(|h| h.eq_ignore_ascii_case("_atom_site_fract_z"))
            .ok_or_else(|| CifError::ParseError("Missing fract z column".into()))?;

        let mut atomic_nums = Vec::with_capacity(rows.len());
        let mut coords = Vec::with_capacity(rows.len() * 3);

        for row in rows {
            let symbol = row[symbol_idx]
                .as_str()
                .ok_or_else(|| CifError::ParseError("Invalid atom symbol".into()))?;
            let atomic_num = extract_atomic_num(symbol)?;
            atomic_nums.push(atomic_num);

            let x = parse_cif_float(row[x_idx].as_str().unwrap_or("0"))?;
            let y = parse_cif_float(row[y_idx].as_str().unwrap_or("0"))?;
            let z = parse_cif_float(row[z_idx].as_str().unwrap_or("0"))?;
            coords.extend_from_slice(&[x, y, z]);
        }

        let frac_coords = MatrixXx3::from_row_slice(&coords);
        let asymmetric_unit = AsymmetricUnit::new(atomic_nums, frac_coords);

        Ok(Crystal::new(unit_cell, space_group, asymmetric_unit))
    }
}

impl TryFrom<CifBlock> for Molecule {
    type Error = CifError;

    fn try_from(block: CifBlock) -> Result<Self, Self::Error> {
        let (headers, rows) = find_atom_loop(&block)
            .ok_or_else(|| CifError::ParseError("No atom loop found".into()))?;

        let symbol_idx = headers
            .iter()
            .position(|h| h.eq_ignore_ascii_case("_atom_site_type_symbol"))
            .or_else(|| {
                headers
                    .iter()
                    .position(|h| h.eq_ignore_ascii_case("_atom_site_label"))
            })
            .ok_or_else(|| CifError::ParseError("Missing atom symbol column".into()))?;

        let cartn_x_idx = headers
            .iter()
            .position(|h| h.eq_ignore_ascii_case("_atom_site_Cartn_x"));
        let cartn_y_idx = headers
            .iter()
            .position(|h| h.eq_ignore_ascii_case("_atom_site_Cartn_y"));
        let cartn_z_idx = headers
            .iter()
            .position(|h| h.eq_ignore_ascii_case("_atom_site_Cartn_z"));

        let fract_x_idx = headers
            .iter()
            .position(|h| h.eq_ignore_ascii_case("_atom_site_fract_x"));
        let fract_y_idx = headers
            .iter()
            .position(|h| h.eq_ignore_ascii_case("_atom_site_fract_y"));
        let fract_z_idx = headers
            .iter()
            .position(|h| h.eq_ignore_ascii_case("_atom_site_fract_z"));

        let mut atomic_nums = Vec::with_capacity(rows.len());
        let mut coords = Vec::with_capacity(rows.len() * 3);
        let mut is_fractional = false;

        if let (Some(x), Some(y), Some(z)) = (cartn_x_idx, cartn_y_idx, cartn_z_idx) {
            for row in rows {
                let symbol = row[symbol_idx]
                    .as_str()
                    .ok_or_else(|| CifError::ParseError("Invalid atom symbol".into()))?;
                let atomic_num = extract_atomic_num(symbol)?;
                atomic_nums.push(atomic_num);

                let vx = parse_cif_float(row[x].as_str().unwrap_or("0"))?;
                let vy = parse_cif_float(row[y].as_str().unwrap_or("0"))?;
                let vz = parse_cif_float(row[z].as_str().unwrap_or("0"))?;
                coords.extend_from_slice(&[vx, vy, vz]);
            }
        } else if let (Some(x), Some(y), Some(z)) = (fract_x_idx, fract_y_idx, fract_z_idx) {
            is_fractional = true;
            for row in rows {
                let symbol = row[symbol_idx]
                    .as_str()
                    .ok_or_else(|| CifError::ParseError("Invalid atom symbol".into()))?;
                let atomic_num = extract_atomic_num(symbol)?;
                atomic_nums.push(atomic_num);

                let vx = parse_cif_float(row[x].as_str().unwrap_or("0"))?;
                let vy = parse_cif_float(row[y].as_str().unwrap_or("0"))?;
                let vz = parse_cif_float(row[z].as_str().unwrap_or("0"))?;
                coords.extend_from_slice(&[vx, vy, vz]);
            }
        } else {
            return Err(CifError::ParseError("Missing coordinate columns".into()));
        }

        let mut coords_matrix = MatrixXx3::from_row_slice(&coords);

        if is_fractional {
            let unit_cell = parse_unit_cell(&block)?;
            coords_matrix = unit_cell.to_cartesian(&coords_matrix);
        }

        Ok(Molecule::new(atomic_nums, coords_matrix))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_extract_atomic_num() {
        assert_eq!(extract_atomic_num("C").unwrap(), 6);
        assert_eq!(extract_atomic_num("C1").unwrap(), 6);
        assert_eq!(extract_atomic_num("He2+").unwrap(), 2);
        assert_eq!(extract_atomic_num("O_1").unwrap(), 8);
    }

    #[test]
    fn test_try_from_molecule() {
        let block = CifBlock::try_from(
            r#"
            data_test
            loop_
            _atom_site_label
            _atom_site_Cartn_x
            _atom_site_Cartn_y
            _atom_site_Cartn_z
            O 0.0 0.0 0.1173
            H 0.0 0.7572 -0.4692
            H 0.0 -0.7572 -0.4692
            "#,
        )
        .unwrap();

        let mol = Molecule::try_from(block).unwrap();
        use crs_core::traits::AtomicData;
        assert_eq!(mol.atomic_nums(), &[8, 1, 1]);
    }

    #[test]
    fn test_try_from_crystal() {
        let block = CifBlock::try_from(
            r#"
            data_test
            _cell_length_a 5.0
            _cell_length_b 5.0
            _cell_length_c 5.0
            _cell_angle_alpha 90
            _cell_angle_beta 90
            _cell_angle_gamma 90
            _space_group_IT_number 1
            loop_
            _atom_site_type_symbol
            _atom_site_fract_x
            _atom_site_fract_y
            _atom_site_fract_z
            C 0.0 0.0 0.0
            C 0.5 0.5 0.5
            "#,
        )
        .unwrap();

        use crs_core::traits::AtomicData;
        use crs_core::traits::CellData;
        let crystal = Crystal::try_from(block).unwrap();
        assert_eq!(crystal.atomic_nums(), &[6, 6]);
        assert_eq!(crystal.lengths(), [5.0, 5.0, 5.0]);
    }
}
