use crate::errors::CoreError;
use moyo::data::{Setting, operations_from_number};
use nalgebra::{Matrix3, Vector3};

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
}

impl SpaceGroup {
    pub fn new(number: u16, setting: Setting, primitive: bool) -> Result<Self, CoreError> {
        if !(1 <= number && number <= 230) {
            return Err(CoreError::InvalidSpaceGroupNumber(number));
        }
        Ok(Self {
            number,
            setting,
            primitive,
        })
    }

    pub fn default_setting(number: u16) -> Self {
        Self::new(number, Setting::Spglib, true).expect("Default settings should always be valid")
    }

    pub fn number(&self) -> u16 {
        self.number
    }

    // Get the symmetry operations for the space group
    // setting: The setting of the space group
    // primitive: Whether to return the operations for the primitive cell or the conventional cell
    pub fn operations(&self) -> Result<Vec<SymOp>, CoreError> {
        let ops = operations_from_number(self.number as i32, self.setting, self.primitive)
            .map_err(|e| CoreError::InvalidSpaceGroupNumber(self.number))?;
        Ok(ops
            .into_iter()
            .map(|op| SymOp::new(op.rotation.cast::<f64>(), op.translation))
            .collect())
    }

    pub fn n_symops(&self) -> Result<usize, CoreError> {
        let ops = operations_from_number(self.number as i32, self.setting, self.primitive)
            .map_err(|e| CoreError::InvalidSpaceGroupNumber(self.number))?;
        Ok(ops.len())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
    fn test_space_group_operations() {
        let sg = SpaceGroup::new(1, Setting::Spglib, true).unwrap();
        let ops = sg.operations().unwrap();
        assert_eq!(ops.len(), 1);

        let sg = SpaceGroup::new(2, Setting::Spglib, true).unwrap();
        let ops = sg.operations().unwrap();
        assert_eq!(ops.len(), 2);

        let sg = SpaceGroup::new(14, Setting::Spglib, true).unwrap();
        let ops = sg.operations().unwrap();
        assert_eq!(ops.len(), 4);

        let sg = SpaceGroup::new(15, Setting::Spglib, true).unwrap();
        let ops = sg.operations().unwrap();
        assert_eq!(ops.len(), 4);
        let sg = SpaceGroup::new(15, Setting::Spglib, false).unwrap();
        let ops = sg.operations().unwrap();
        assert_eq!(ops.len(), 8);
    }
}
