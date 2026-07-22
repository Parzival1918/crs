use crate::traits::Parser;
use crs_core::au::AsymmetricUnit;
use crs_core::crystal::Crystal;
use crs_core::data::SPECIES_NAMES;
use crs_core::spacegroup::{SpaceGroup, SymOp};
use crs_core::unitcell::UnitCell;
use nalgebra::MatrixXx3;
use std::io::{BufRead, Result as IoResult};
use std::str::FromStr;
use thiserror::Error as ThisError;

#[derive(Debug, ThisError)]
pub enum ShelxError {
    #[error("I/O error: {0}")]
    IoError(#[from] std::io::Error),
    #[error("Parse error: {0}")]
    ParseError(String),
}

#[derive(Debug, Clone, PartialEq)]
pub struct ShelxAtom {
    pub name: String,
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct ShelxFrame {
    pub title: String,
    pub cell: [f64; 6],
    pub symops: Vec<String>,
    pub latt: Option<i32>,
    pub atoms: Vec<ShelxAtom>,
}

impl Parser<ShelxFrame> for ShelxFrame {
    type E = ShelxError;

    fn parse_from_reader<R: BufRead>(reader: &mut R) -> Result<Option<ShelxFrame>, Self::E> {
        let mut frame = ShelxFrame::default();
        let mut has_data = false;

        let mut line = String::new();
        loop {
            line.clear();
            let bytes_read = reader.read_line(&mut line)?;
            if bytes_read == 0 {
                break; // EOF
            }

            let trimmed = line.trim();
            if trimmed.is_empty() {
                continue;
            }
            has_data = true;

            if trimmed.starts_with("END") || trimmed.starts_with("HKLF") {
                break;
            }

            let mut tokens = trimmed.split_whitespace();
            let kw_opt = tokens.next();
            if kw_opt.is_none() {
                continue;
            }
            let kw = kw_opt.unwrap().to_uppercase();

            match kw.as_str() {
                "TITL" => {
                    frame.title = trimmed[4..].trim().to_string();
                }
                "CELL" => {
                    // CELL lambda a b c alpha beta gamma
                    let _lambda = tokens.next();
                    let a = tokens.next().and_then(|s| parse_f64_strip_esd(s).ok());
                    let b = tokens.next().and_then(|s| parse_f64_strip_esd(s).ok());
                    let c = tokens.next().and_then(|s| parse_f64_strip_esd(s).ok());
                    let alpha = tokens.next().and_then(|s| parse_f64_strip_esd(s).ok());
                    let beta = tokens.next().and_then(|s| parse_f64_strip_esd(s).ok());
                    let gamma = tokens.next().and_then(|s| parse_f64_strip_esd(s).ok());

                    if let (Some(a), Some(b), Some(c), Some(alpha), Some(beta), Some(gamma)) =
                        (a, b, c, alpha, beta, gamma)
                    {
                        frame.cell = [a, b, c, alpha, beta, gamma];
                    }
                }
                "LATT" => {
                    if let Some(val) = tokens.next() {
                        frame.latt = val.parse::<i32>().ok();
                    }
                }
                "SYMM" => {
                    let sym_str = trimmed[4..].trim().to_string();
                    frame.symops.push(sym_str);
                }
                _ => {
                    if is_known_keyword(&kw) {
                        continue;
                    }

                    // Attempt to parse as atom
                    // Format: name sfac x y z ...
                    let tokens: Vec<&str> = trimmed.split_whitespace().collect();
                    if tokens.len() >= 5 {
                        let name = tokens[0].to_string();
                        let sfac = tokens[1].parse::<i32>();
                        let x = parse_f64_strip_esd(tokens[2]);
                        let y = parse_f64_strip_esd(tokens[3]);
                        let z = parse_f64_strip_esd(tokens[4]);

                        if let (Ok(_), Ok(x), Ok(y), Ok(z)) = (sfac, x, y, z) {
                            frame.atoms.push(ShelxAtom { name, x, y, z });
                        }
                    }
                }
            }
        }

        if has_data { Ok(Some(frame)) } else { Ok(None) }
    }
}

fn parse_f64_strip_esd(s: &str) -> Result<f64, std::num::ParseFloatError> {
    let s = if let Some(idx) = s.find('(') {
        &s[..idx]
    } else {
        s
    };
    s.parse::<f64>()
}

fn is_known_keyword(kw: &str) -> bool {
    matches!(
        kw,
        "ZERR"
            | "SFAC"
            | "UNIT"
            | "DISP"
            | "LAUE"
            | "REST"
            | "SADI"
            | "DFIX"
            | "DANG"
            | "BUMP"
            | "SAME"
            | "SINC"
            | "FLAT"
            | "DELU"
            | "SIMU"
            | "ISOR"
            | "NCSY"
            | "SUMP"
            | "CHIV"
            | "EADP"
            | "EXYZ"
            | "EXTI"
            | "SWAT"
            | "HOPE"
            | "MERG"
            | "SPEC"
            | "RESI"
            | "MOVE"
            | "ANIS"
            | "AFIX"
            | "HFIX"
            | "FRAG"
            | "FEND"
            | "EXPT"
            | "SIZE"
            | "HTAB"
            | "LIST"
            | "ACTA"
            | "BOND"
            | "CONF"
            | "MPLA"
            | "RTAB"
            | "PLAN"
            | "PRIG"
            | "WIGL"
            | "BISO"
            | "BIND"
            | "FREE"
            | "ELEM"
            | "GRID"
            | "TIME"
            | "WGHT"
            | "FVAR"
            | "OMIT"
            | "SHEL"
            | "BASF"
            | "TWIN"
            | "EQIV"
            | "CONN"
            | "PART"
            | "SPIN"
            | "LONE"
            | "REM"
            | "MOLE"
    )
}

impl TryFrom<&str> for ShelxFrame {
    type Error = ShelxError;

    fn try_from(s: &str) -> Result<Self, Self::Error> {
        let mut reader = s.as_bytes();
        let presult = Self::parse_from_reader(&mut reader)?;
        if let Some(result) = presult {
            Ok(result)
        } else {
            Err(ShelxError::ParseError("Failed to parse SHELX frame".into()))
        }
    }
}

impl TryFrom<ShelxFrame> for Crystal {
    type Error = ShelxError;

    fn try_from(frame: ShelxFrame) -> Result<Self, Self::Error> {
        let unit_cell = UnitCell::new(
            [frame.cell[0], frame.cell[1], frame.cell[2]],
            [frame.cell[3], frame.cell[4], frame.cell[5]],
        );

        let mut ops = Vec::new();
        ops.push(SymOp::from_str("x, y, z").unwrap());
        for sym_str in &frame.symops {
            if let Ok(op) = SymOp::from_str(sym_str) {
                ops.push(op);
            }
        }

        let space_group =
            SpaceGroup::try_from_symops(&ops).unwrap_or_else(|_| SpaceGroup::default_setting(1));

        let mut atomic_nums = Vec::with_capacity(frame.atoms.len());
        let mut coords = Vec::with_capacity(frame.atoms.len() * 3);

        for atom in &frame.atoms {
            let alpha_part: String = atom
                .name
                .chars()
                .take_while(|c| c.is_alphabetic())
                .collect();
            let atomic_num = SPECIES_NAMES
                .iter()
                .position(|&s| s.eq_ignore_ascii_case(&alpha_part))
                .map(|i| i as u8)
                .ok_or_else(|| {
                    ShelxError::ParseError(format!("Unknown species '{}'", atom.name))
                })?;

            atomic_nums.push(atomic_num);
            coords.extend_from_slice(&[atom.x, atom.y, atom.z]);
        }

        let frac_coords = MatrixXx3::from_row_slice(&coords);
        let asymmetric_unit = AsymmetricUnit::new(atomic_nums, frac_coords);

        Ok(Crystal::new(unit_cell, space_group, asymmetric_unit))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_shelx() {
        let input = "TITL test molecule
CELL 0.71073 10.0 10.0 10.0 90.0 90.0 90.0
ZERR 4 0.01 0.01 0.01 0.0 0.0 0.0
SFAC C H
UNIT 10 20
C1 1 0.1(1) 0.2 0.3 11 0.05
H1 2 0.4 0.5 0.6 11 0.05
END";
        let mut reader = input.as_bytes();
        let frame = ShelxFrame::parse_from_reader(&mut reader).unwrap().unwrap();
        assert_eq!(frame.title, "test molecule");
        assert_eq!(frame.cell, [10.0, 10.0, 10.0, 90.0, 90.0, 90.0]);
        assert_eq!(frame.symops.len(), 0);
        assert_eq!(frame.latt, None);
        assert_eq!(frame.atoms.len(), 2);
        assert_eq!(frame.atoms[0].name, "C1");
        assert_eq!(frame.atoms[0].x, 0.1);
        assert_eq!(frame.atoms[1].name, "H1");
    }

    #[test]
    fn test_parse_shelx_symmetry() {
        let input = "TITL test molecule
CELL 0.71073 10.0 10.0 10.0 90.0 90.0 90.0
ZERR 4 0.01 0.01 0.01 0.0 0.0 0.0
LATT -1
SYMM -x, y+1/2, -z
SYMM x, -y, z+1/2
SFAC C H
UNIT 10 20
C1 1 0.1(1) 0.2 0.3 11 0.05
H1 2 0.4 0.5 0.6 11 0.05
END";
        let mut reader = input.as_bytes();
        let frame = ShelxFrame::parse_from_reader(&mut reader).unwrap().unwrap();
        assert_eq!(frame.symops.len(), 2);
        assert_eq!(frame.symops[0], "-x, y+1/2, -z");
        assert_eq!(frame.symops[1], "x, -y, z+1/2");
        assert_eq!(frame.latt, Some(-1));
    }

    #[test]
    fn test_parse_shelx_no_esd() {
        let input = "C1 1 0.1 0.2 0.3
HKLF 4";
        let mut reader = input.as_bytes();
        let frame = ShelxFrame::parse_from_reader(&mut reader).unwrap().unwrap();
        assert_eq!(frame.atoms.len(), 1);
        assert_eq!(frame.atoms[0].x, 0.1);
    }
    
    #[test]
    fn test_try_from_crystal() {
        let input = "TITL test molecule
CELL 0.71073 10.0 10.0 10.0 90.0 90.0 90.0
SYMM -x, y+1/2, -z
SYMM x, -y, z+1/2
SFAC C H
C1 1 0.1 0.2 0.3 11 0.05
H1 2 0.4 0.5 0.6 11 0.05
END";
        let frame = ShelxFrame::try_from(input).unwrap();
        let crystal = Crystal::try_from(frame).unwrap();

        let au = crystal.asymmetric_unit();
        assert_eq!(au.atomic_nums(), &[6, 1]);
        assert_eq!(au.frac_coords().nrows(), 2);
    }
}
