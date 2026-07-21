use crate::data::{DEFAULT_ANG_TOL, DEFAULT_DIST_TOL};
use crate::traits::CellData;
use crate::utils::{cart_to_frac, cell_matrix_and_volume, frac_to_cart};
use core::f64::consts::{FRAC_2_PI, FRAC_PI_2};
use moyo::base::Cell;
use nalgebra::{Matrix3, MatrixXx3};

#[derive(Debug, Clone, PartialEq)]
pub enum CellType {
    Cubic,
    Rhombohedral,
    Hexagonal,
    Tetragonal,
    Orthorhombic,
    Monoclinic,
    Triclinic,
}

/// Angles are in radians
#[derive(Debug, Clone)]
pub struct UnitCell {
    lengths: [f64; 3],
    angles: [f64; 3],
    matrix: Matrix3<f64>,
    inv_matrix: Matrix3<f64>,
    volume: f64,
}

impl UnitCell {
    pub fn new(lengths: [f64; 3], angles: [f64; 3]) -> Self {
        let (matrix, inv_matrix, volume) = cell_matrix_and_volume(lengths, angles);
        Self {
            lengths,
            angles,
            matrix,
            inv_matrix,
            volume,
        }
    }

    pub fn from_lengths_and_angles(
        a: f64,
        b: f64,
        c: f64,
        alpha: f64,
        beta: f64,
        gamma: f64,
    ) -> Self {
        Self::new([a, b, c], [alpha, beta, gamma])
    }

    pub fn cubic(a: f64) -> Self {
        Self::new([a, a, a], [FRAC_PI_2, FRAC_PI_2, FRAC_PI_2])
    }

    pub fn tetragonal(a: f64, c: f64) -> Self {
        Self::new([a, a, c], [FRAC_PI_2, FRAC_PI_2, FRAC_PI_2])
    }

    pub fn orthorhombic(a: f64, b: f64, c: f64) -> Self {
        Self::new([a, b, c], [FRAC_PI_2, FRAC_PI_2, FRAC_PI_2])
    }

    pub fn hexagonal(a: f64, c: f64) -> Self {
        Self::new([a, a, c], [FRAC_PI_2, FRAC_PI_2, FRAC_2_PI / 3.0])
    }

    pub fn rhombohedral(a: f64, alpha: f64) -> Self {
        Self::new([a, a, a], [alpha, alpha, alpha])
    }

    pub fn monoclinic(a: f64, b: f64, c: f64, beta: f64) -> Self {
        Self::new([a, b, c], [FRAC_PI_2, beta, FRAC_PI_2])
    }

    pub fn triclinic(a: f64, b: f64, c: f64, alpha: f64, beta: f64, gamma: f64) -> Self {
        Self::new([a, b, c], [alpha, beta, gamma])
    }
}

impl CellData for UnitCell {
    fn cell_matrix(&self) -> Matrix3<f64> {
        self.matrix.clone()
    }

    fn inv_cell_matrix(&self) -> Matrix3<f64> {
        self.inv_matrix.clone()
    }

    fn volume(&self) -> f64 {
        self.volume
    }

    fn lengths(&self) -> [f64; 3] {
        self.lengths
    }

    fn angles(&self) -> [f64; 3] {
        self.angles
    }

    fn to_cartesian(&self, frac_coords: &MatrixXx3<f64>) -> MatrixXx3<f64> {
        frac_to_cart(frac_coords, self)
    }

    fn to_fractional(&self, cart_coords: &MatrixXx3<f64>) -> MatrixXx3<f64> {
        cart_to_frac(cart_coords, self)
    }
}

impl CellData for Matrix3<f64> {
    fn cell_matrix(&self) -> Matrix3<f64> {
        self.clone()
    }

    fn inv_cell_matrix(&self) -> Matrix3<f64> {
        self.try_inverse()
            .expect("Matrix is not invertible")
            .clone()
    }

    fn volume(&self) -> f64 {
        self.determinant()
    }

    fn lengths(&self) -> [f64; 3] {
        [self.row(0).norm(), self.row(1).norm(), self.row(2).norm()]
    }

    fn angles(&self) -> [f64; 3] {
        let a = self.row(0);
        let b = self.row(1);
        let c = self.row(2);

        let alpha = (b.dot(&c) / (b.norm() * c.norm())).acos();
        let beta = (a.dot(&c) / (a.norm() * c.norm())).acos();
        let gamma = (a.dot(&b) / (a.norm() * b.norm())).acos();

        [alpha, beta, gamma]
    }

    fn to_cartesian(&self, frac_coords: &MatrixXx3<f64>) -> MatrixXx3<f64> {
        frac_to_cart(frac_coords, self)
    }

    fn to_fractional(&self, cart_coords: &MatrixXx3<f64>) -> MatrixXx3<f64> {
        cart_to_frac(cart_coords, self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::f64::consts::FRAC_PI_3;

    #[test]
    fn test_unit_cell() {
        let lengths = [1.0, 1.0, 1.0];
        let angles = [FRAC_PI_2, FRAC_PI_2, FRAC_PI_2];
        let cell = UnitCell::new(lengths, angles);
        assert_eq!(cell.cell_type(), CellType::Cubic);

        let lengths = [1.0, 1.0, 1.0];
        let angles = [FRAC_PI_3, FRAC_PI_3, FRAC_PI_3];
        let cell = UnitCell::new(lengths, angles);
        assert_eq!(cell.cell_type(), CellType::Rhombohedral);

        let cell = UnitCell::hexagonal(1.0, 2.0);
        assert_eq!(cell.cell_type(), CellType::Hexagonal);

        let cell = UnitCell::tetragonal(1.0, 2.0);
        assert_eq!(cell.cell_type(), CellType::Tetragonal);

        let cell = UnitCell::orthorhombic(1.0, 2.0, 3.0);
        assert_eq!(cell.cell_type(), CellType::Orthorhombic);
    }
}
