use crate::errors::CoreError;
use moyo::data::{Setting, operations_from_number, hall_symbol_entry};
use nalgebra::{Matrix3, Vector3};
use std::fmt;
use std::str::FromStr;
use std::sync::OnceLock;
use std::collections::HashMap;

type SymOpCache = HashMap<Vec<String>, SpaceGroup>;
static SPACE_GROUP_CACHE: OnceLock<SymOpCache> = OnceLock::new();

fn get_space_group_cache() -> &'static SymOpCache {
    SPACE_GROUP_CACHE.get_or_init(|| {
        let mut cache = HashMap::new();
        for hall_number in 1..=530 {
            if let Some(entry) = hall_symbol_entry(hall_number) {
                let number = entry.number as u16;
                let setting = Setting::HallNumber(hall_number);
                for primitive in [true, false] {
                    if let Ok(sg) = SpaceGroup::new(number, setting, primitive) {
                        let mut ops: Vec<String> =
                            sg.operations().iter().map(|op| op.to_string()).collect();
                        ops.sort();
                        cache.insert(ops, sg);
                    }
                }
            }
        }
        cache
    })
}

#[derive(Debug, Clone)]
pub struct SymOp {
    rotation: Matrix3<f64>,
    translation: Vector3<f64>,
}

impl SymOp {
    pub fn new(rotation: Matrix3<f64>, translation: Vector3<f64>) -> Self {
        Self {
            rotation,
            translation,
        }
    }

    pub fn rotation(&self) -> &Matrix3<f64> {
        &self.rotation
    }

    pub fn translation(&self) -> &Vector3<f64> {
        &self.translation
    }
}

#[derive(Debug, Clone)]
pub struct SpaceGroup {
    number: u16,
    setting: Setting,
    primitive: bool,
    operations: Vec<SymOp>,
}

impl SpaceGroup {
    pub fn new(number: u16, setting: Setting, primitive: bool) -> Result<Self, CoreError> {
        if !(1 <= number && number <= 230) {
            return Err(CoreError::InvalidSpaceGroupNumber(number));
        }
        let operations = operations_from_number(number as i32, setting, primitive)
            .map_err(|_| CoreError::InvalidSpaceGroupSettings(setting, primitive))?
            .into_iter()
            .map(|op| SymOp::new(op.rotation.cast::<f64>(), op.translation))
            .collect();
        Ok(Self {
            number,
            setting,
            primitive,
            operations,
        })
    }

    pub fn default_setting(number: u16) -> Self {
        Self::new(number, Setting::Spglib, true).expect("Default settings should always be valid")
    }

    pub fn number(&self) -> u16 {
        self.number
    }

    // Get the symmetry operations for the space group
    pub fn operations(&self) -> &[SymOp] {
        &self.operations
    }

    pub fn n_symops(&self) -> usize {
        self.operations.len()
    }

    /// Identify a Space Group from a list of symmetry operations (as strings).
    pub fn try_from_symops(ops: &[SymOp]) -> Result<Self, CoreError> {
        let mut target_ops = Vec::with_capacity(ops.len());
        for op in ops {
            target_ops.push(op.to_string());
        }
        target_ops.sort();

        let cache = get_space_group_cache();
        if let Some(sg) = cache.get(&target_ops) {
            Ok(sg.clone())
        } else {
            Err(CoreError::SpaceGroupNotFound(target_ops))
        }
    }
}

fn parse_fraction(s: &str) -> Option<f64> {
    if let Some((num, den)) = s.split_once('/') {
        let n: f64 = num.parse().ok()?;
        let d: f64 = den.parse().ok()?;
        if d == 0.0 { None } else { Some(n / d) }
    } else {
        s.parse().ok()
    }
}

fn format_fraction(val: f64) -> String {
    const TOL: f64 = 1e-5;
    let fractions = [
        (1.0 / 2.0, "1/2"),
        (1.0 / 3.0, "1/3"),
        (2.0 / 3.0, "2/3"),
        (1.0 / 4.0, "1/4"),
        (3.0 / 4.0, "3/4"),
        (1.0 / 6.0, "1/6"),
        (5.0 / 6.0, "5/6"),
    ];
    for (f, s) in fractions {
        if (val - f).abs() < TOL {
            return s.to_string();
        }
    }
    if (val - val.round()).abs() < TOL {
        return (val.round() as i32).to_string();
    }
    format!("{val}")
}

pub fn symop_string(op: &SymOp) -> String {
    let mut op_string = String::new();
    let vars = ['x', 'y', 'z'];
    for i in 0..3 {
        let mut row_str = String::new();
        for j in 0..3 {
            let r = op.rotation[(i, j)];
            if r.abs() > 1e-5 {
                if r > 0.0 {
                    row_str.push('+');
                } else {
                    row_str.push('-');
                }

                let abs_r = r.abs();
                if (abs_r - 1.0).abs() > 1e-5 {
                    row_str.push_str(&format_fraction(abs_r));
                }
                row_str.push(vars[j]);
            }
        }

        let t = op.translation[i];
        if t.abs() > 1e-5 {
            if t > 0.0 {
                row_str.push('+');
            } else {
                row_str.push('-');
            }
            row_str.push_str(&format_fraction(t.abs()));
        }

        if row_str.is_empty() {
            row_str.push('0');
        }

        op_string.push_str(&row_str);
        if i < 2 {
            op_string.push(',');
        }
    }
    op_string
}

impl FromStr for SymOp {
    type Err = CoreError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let s_no_ws: String = s.chars().filter(|c| !c.is_whitespace()).collect();
        let parts: Vec<&str> = s_no_ws.split(',').collect();
        if parts.len() != 3 {
            return Err(CoreError::InvalidSymOpString(s.to_string()));
        }

        let mut rotation = Matrix3::zeros();
        let mut translation = Vector3::zeros();

        for (i, part) in parts.iter().enumerate() {
            let mut part_mod = String::new();
            for c in part.chars() {
                if c == '-' {
                    part_mod.push_str("+-");
                } else {
                    part_mod.push(c);
                }
            }

            let terms = part_mod.split('+').filter(|t| !t.is_empty());
            for term in terms {
                if term.contains('x') || term.contains('y') || term.contains('z') {
                    let mut sign = 1.0;
                    let mut var_idx = 0;
                    let mut coeff_str = String::new();

                    for c in term.chars() {
                        match c {
                            'x' => var_idx = 0,
                            'y' => var_idx = 1,
                            'z' => var_idx = 2,
                            '-' => sign = -1.0,
                            _ => coeff_str.push(c),
                        }
                    }

                    if coeff_str.starts_with('/') {
                        coeff_str.insert(0, '1');
                    }

                    let coeff = if coeff_str.is_empty() {
                        1.0
                    } else {
                        parse_fraction(&coeff_str)
                            .ok_or_else(|| CoreError::InvalidSymOpString(s.to_string()))?
                    };

                    rotation[(i, var_idx)] += sign * coeff;
                } else {
                    translation[i] += parse_fraction(term)
                        .ok_or_else(|| CoreError::InvalidSymOpString(s.to_string()))?;
                }
            }
        }

        Ok(SymOp::new(rotation, translation))
    }
}

impl fmt::Display for SymOp {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", symop_string(self))
    }
}

#[cfg(test)]
mod tests {
    use std::ops;
    use super::*;

    #[test]
    fn test_symop_display_and_parse() {
        let op_str = "+x,+y,+z";
        let op = SymOp::from_str(op_str).unwrap();
        assert_eq!(op.rotation(), &Matrix3::identity());
        assert_eq!(op.translation(), &Vector3::zeros());
        assert_eq!(op.to_string(), "+x,+y,+z");

        let op_str = "-x,-y,-z";
        let op = SymOp::from_str(op_str).unwrap();
        assert_eq!(
            op.rotation(),
            &Matrix3::new(-1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0)
        );
        assert_eq!(op.translation(), &Vector3::zeros());
        assert_eq!(op.to_string(), "-x,-y,-z");

        let op_str = "-y,x-y,z+1/3";
        let op = SymOp::from_str(op_str).unwrap();
        assert_eq!(
            op.rotation(),
            &Matrix3::new(0.0, -1.0, 0.0, 1.0, -1.0, 0.0, 0.0, 0.0, 1.0)
        );
        assert_eq!(op.translation(), &Vector3::new(0.0, 0.0, 1.0 / 3.0));
        // Note: the order might be different depending on exact parsing.
        // Our parsing logic maintains matrix structure so to_string will output in x, y, z order:
        assert_eq!(op.to_string(), "-y,+x-y,+z+1/3");
    }

    #[test]
    fn test_space_group_new() {
        let sg = SpaceGroup::new(1, Setting::Spglib, true).unwrap();
        assert_eq!(sg.number(), 1);
    }

    #[test]
    fn test_space_group_new_invalid_number() {
        let result = SpaceGroup::new(231, Setting::Spglib, true);
        assert!(result.is_err());
    }

    #[test]
    fn test_space_group_from_symops() {
        let ops = ["+x,+y,+z"];
        let ops: Vec<SymOp> = ops.iter().map(|s| SymOp::from_str(s).unwrap()).collect();
        let sg = SpaceGroup::try_from_symops(&ops).unwrap();
        assert_eq!(sg.number(), 1);

        let ops_p1bar = ["+x,+y,+z", "-x,-y,-z"];
        let ops_p1bar: Vec<SymOp> = ops_p1bar
            .iter()
            .map(|s| SymOp::from_str(s).unwrap())
            .collect();
        let sg = SpaceGroup::try_from_symops(&ops_p1bar).unwrap();
        assert_eq!(sg.number(), 2);

        let ops = ["-x,-y+1/2,+z"];
        let ops: Vec<SymOp> = ops.iter().map(|s| SymOp::from_str(s).unwrap()).collect();
        let sg = SpaceGroup::try_from_symops(&ops);
        assert!(sg.is_err());
    }

    #[test]
    fn test_space_group_operations() {
        let sg = SpaceGroup::new(1, Setting::Spglib, true).unwrap();
        let ops = sg.operations();
        assert_eq!(ops.len(), 1);

        let sg = SpaceGroup::new(2, Setting::Spglib, true).unwrap();
        let ops = sg.operations();
        assert_eq!(ops.len(), 2);

        let sg = SpaceGroup::new(14, Setting::Spglib, true).unwrap();
        let ops = sg.operations();
        assert_eq!(ops.len(), 4);

        let sg = SpaceGroup::new(15, Setting::Spglib, true).unwrap();
        let ops = sg.operations();
        assert_eq!(ops.len(), 4);
        let sg = SpaceGroup::new(15, Setting::Spglib, false).unwrap();
        let ops = sg.operations();
        assert_eq!(ops.len(), 8);
    }
}
