use crate::au::AsymmetricUnit;
use crate::data::{ANGSTROM3_TO_CM3, AVOGADRO_NUMBER};
use crate::spacegroup::SpaceGroup;
use crate::traits::{AtomicData, CartAtomicData, CellData, PeriodicAtomicData, SupercellData};
use crate::unitcell::UnitCell;
use crate::utils::{apply_symops, chemical_formula, molar_mass, wrap_coordinates_in_place, frac_to_cart, cart_to_frac};
use moyo::base::Cell;
use nalgebra::{Matrix3, MatrixXx3};

#[derive(Debug, Clone)]
pub struct Crystal {
    unit_cell: UnitCell,
    space_group: SpaceGroup,
    asymmetric_unit: AsymmetricUnit,
    atomic_nums: Vec<u8>,
    fractional_coords: MatrixXx3<f64>,
    cartesian_coords: MatrixXx3<f64>,
}

impl Crystal {
    pub fn new(
        unit_cell: UnitCell,
        space_group: SpaceGroup,
        asymmetric_unit: AsymmetricUnit,
    ) -> Self {
        let sym_ops = space_group.operations();
        let mut fractional_coords = apply_symops(asymmetric_unit.frac_coords(), sym_ops);
        wrap_coordinates_in_place(&mut fractional_coords);
        let cartesian_coords = unit_cell.to_cartesian(&fractional_coords);

        let mut atomic_nums =
            Vec::with_capacity(asymmetric_unit.atomic_nums().len() * sym_ops.len());
        for _ in 0..sym_ops.len() {
            atomic_nums.extend_from_slice(asymmetric_unit.atomic_nums());
        }

        Self {
            unit_cell,
            space_group,
            asymmetric_unit,
            atomic_nums,
            fractional_coords,
            cartesian_coords,
        }
    }

    /// Compute the chemical formula of the full crystal contents.
    fn formula(&self) -> String {
        chemical_formula(self.atomic_nums(), self.space_group.n_symops() as usize)
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
}

impl AtomicData for Crystal {
    /// Atomic numbers of the full cell contents.
    fn atomic_nums(&self) -> &[u8] {
        &self.atomic_nums
    }

    fn covalent_radii(&self) -> Vec<f64> {
        unimplemented!()
    }
}

impl CartAtomicData for Crystal {
    /// Cartesian coordinates of the full cell contents.
    fn cartesian_coords(&self) -> &MatrixXx3<f64> {
        &self.cartesian_coords
    }
}

impl CellData for Crystal {
    fn lengths(&self) -> [f64; 3] {
        self.unit_cell.lengths()
    }

    fn angles(&self) -> [f64; 3] {
        self.unit_cell.angles()
    }

    fn cell_matrix(&self) -> Matrix3<f64> {
        self.unit_cell.cell_matrix()
    }

    fn inv_cell_matrix(&self) -> Matrix3<f64> {
        self.unit_cell.inv_cell_matrix()
    }

    fn volume(&self) -> f64 {
        self.unit_cell.volume()
    }

    fn to_cartesian(&self, frac_coords: &MatrixXx3<f64>) -> MatrixXx3<f64> {
        frac_to_cart(frac_coords, self)
    }

    fn to_fractional(&self, cart_coords: &MatrixXx3<f64>) -> MatrixXx3<f64> {
        cart_to_frac(cart_coords, self)
    }
}

impl PeriodicAtomicData for Crystal {
    /// Fractional coordinates of the full cell contents.
    fn fractional_coords(&self) -> &MatrixXx3<f64> {
        &self.fractional_coords
    }
}

pub struct SupercellCrystal {
    crystal: Crystal,
    supercell_data: Box<dyn SupercellData>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_crystal_creation() {
        let unit_cell = UnitCell::cubic(5.0);
        let space_group = SpaceGroup::default_setting(1);
        let atomic_nums = vec![6, 1, 1, 1, 1]; // CH4
        let frac_coords = MatrixXx3::repeat(5, 0.0);
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
        let frac_coords = MatrixXx3::repeat(5, 0.0);
        let asymmetric_unit = AsymmetricUnit::new(atomic_nums, frac_coords);
        let crystal = Crystal::new(unit_cell, space_group, asymmetric_unit);
        let density = crystal.density();
        assert!(density > 0.0);
    }
}
