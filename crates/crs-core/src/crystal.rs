use crate::au::AsymmetricUnit;
use crate::data::{ANGSTROM3_TO_CM3, AVOGADRO_NUMBER};
use crate::spacegroup::SpaceGroup;
use crate::traits::{AtomicData, CartAtomicData, CellData, FracAtomicData};
use crate::unitcell::UnitCell;
use crate::utils::{chemical_formula, molar_mass};

#[derive(Debug, Clone)]
pub struct Crystal {
    unit_cell: UnitCell,
    space_group: SpaceGroup,
    asymmetric_unit: AsymmetricUnit,
}

impl Crystal {
    pub fn new(
        unit_cell: UnitCell,
        space_group: SpaceGroup,
        asymmetric_unit: AsymmetricUnit,
    ) -> Self {
        Self {
            unit_cell,
            space_group,
            asymmetric_unit,
        }
    }

    /// Compute the chemical formula of the full crystal contents.
    fn formula(&self) -> String {
        chemical_formula(
            self.atomic_nums(),
            self.space_group.n_symops() as usize,
        )
    }

    pub fn unit_cell(&self) -> &UnitCell {
        &self.unit_cell
    }

    pub fn space_group(&self) -> &SpaceGroup {
        &self.space_group
    }

    pub fn asymmetric_unit(&self) -> &AsymmetricUnit {
        &self.asymmetric_unit
    }

    pub fn density(&self) -> f64 {
        let au_mass = molar_mass(self.asymmetric_unit.atomic_nums()) / AVOGADRO_NUMBER; // Convert to kg
        let total_mass = au_mass * self.space_group().n_symops() as f64; // Total mass in kg
        let volume = self.unit_cell.volume(); // Volume in cubic angstroms
        total_mass / (volume * ANGSTROM3_TO_CM3) // Convert volume to cubic centimeters and calculate density in g/cm^3
    }

    /// Atomic numbers of the full cell contents.
    pub fn atomic_nums(&self) -> &[u8] {
        unimplemented!()
    }
}

impl AtomicData for Crystal {
    /// Atomic numbers of the full cell contents.
    fn atomic_nums(&self) -> &[u8] {
        self.atomic_nums()
    }
}

impl CartAtomicData for Crystal {
    /// Cartesian coordinates of the full cell contents.
    fn cartesian_coords(&self) -> &nalgebra::MatrixXx3<f64> {
        unimplemented!()
    }
}

impl FracAtomicData for Crystal {
    /// Fractional coordinates of the full cell contents.
    fn fractional_coords(&self) -> &nalgebra::MatrixXx3<f64> {
        unimplemented!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_crystal_creation() {
        let unit_cell = UnitCell::cubic(5.0);
        let space_group = SpaceGroup::default_setting(1);
        let atomic_nums = vec![6, 1, 1, 1, 1]; // CH4
        let frac_coords = nalgebra::MatrixXx3::repeat(5, 0.0);
        let asymmetric_unit = AsymmetricUnit::new(atomic_nums, frac_coords);
        let crystal = Crystal::new(unit_cell, space_group, asymmetric_unit);
        assert_eq!(
            crystal.unit_cell().cell_type(),
            crate::unitcell::CellType::Cubic
        );
        assert_eq!(crystal.space_group().number(), 1);
        assert_eq!(crystal.asymmetric_unit().atomic_nums().len(), 5);
    }

    #[test]
    fn test_crystal_density() {
        let unit_cell = UnitCell::cubic(5.0);
        let space_group = SpaceGroup::default_setting(1);
        let atomic_nums = vec![6, 1, 1, 1, 1]; // CH4
        let frac_coords = nalgebra::MatrixXx3::repeat(5, 0.0);
        let asymmetric_unit = AsymmetricUnit::new(atomic_nums, frac_coords);
        let crystal = Crystal::new(unit_cell, space_group, asymmetric_unit);
        let density = crystal.density();
        assert!(density > 0.0);
    }
}
