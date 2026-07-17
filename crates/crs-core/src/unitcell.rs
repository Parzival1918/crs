use crate::data::{DEFAULT_ANG_TOL, DEFAULT_DIST_TOL};
use crate::traits::CellData;
use crate::utils::{cart_to_frac, cell_matrix_and_volume, frac_to_cart};
use core::f64::consts::FRAC_PI_2;
use moyo::base::Cell;
use nalgebra::{Matrix3, MatrixXx3};

#[derive(Debug, Clone, PartialEq)]
pub enum CellType {
    Cubic,
    Rombohedral,
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
        Self::new(
            [a, a, c],
            [FRAC_PI_2, FRAC_PI_2, 2.0 * std::f64::consts::PI / 3.0],
        )
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

    pub fn a(&self) -> f64 {
        self.lengths[0]
    }

    pub fn b(&self) -> f64 {
        self.lengths[1]
    }

    pub fn c(&self) -> f64 {
        self.lengths[2]
    }

    pub fn alpha(&self) -> f64 {
        self.angles[0]
    }

    pub fn beta(&self) -> f64 {
        self.angles[1]
    }

    pub fn gamma(&self) -> f64 {
        self.angles[2]
    }

    pub fn matrix(&self) -> &Matrix3<f64> {
        &self.matrix
    }

    pub fn volume(&self) -> f64 {
        self.volume
    }

    fn abc_close(&self) -> bool {
        (self.a() - self.b()).abs() < DEFAULT_DIST_TOL
            && (self.a() - self.c()).abs() < DEFAULT_DIST_TOL
            && (self.b() - self.c()).abs() < DEFAULT_DIST_TOL
    }

    fn abc_different(&self) -> bool {
        (self.a() - self.b()).abs() > DEFAULT_DIST_TOL
            && (self.a() - self.c()).abs() > DEFAULT_DIST_TOL
            && (self.b() - self.c()).abs() > DEFAULT_DIST_TOL
    }

    fn ab_close_c_different(&self) -> bool {
        (self.a() - self.b()).abs() < DEFAULT_DIST_TOL
            && (self.a() - self.c()).abs() > DEFAULT_DIST_TOL
            && (self.b() - self.c()).abs() > DEFAULT_DIST_TOL
    }

    fn orthogonal_angles(&self) -> bool {
        (self.alpha() - FRAC_PI_2).abs() < DEFAULT_ANG_TOL
            && (self.beta() - FRAC_PI_2).abs() < DEFAULT_ANG_TOL
            && (self.gamma() - FRAC_PI_2).abs() < DEFAULT_ANG_TOL
    }

    fn angles_close(&self) -> bool {
        (self.alpha() - self.beta()).abs() < DEFAULT_ANG_TOL
            && (self.alpha() - self.gamma()).abs() < DEFAULT_ANG_TOL
            && (self.beta() - self.gamma()).abs() < DEFAULT_ANG_TOL
    }

    pub fn cell_type(&self) -> CellType {
        if self.abc_close() && self.orthogonal_angles() {
            CellType::Cubic
        } else if self.abc_close()
            && self.angles_close()
            && (self.alpha() - FRAC_PI_2).abs() > DEFAULT_ANG_TOL
        {
            CellType::Rombohedral
        } else if self.ab_close_c_different()
            && (self.alpha() - FRAC_PI_2).abs() < DEFAULT_ANG_TOL
            && (self.beta() - FRAC_PI_2).abs() < DEFAULT_ANG_TOL
            && (self.gamma() - FRAC_PI_2).abs() > DEFAULT_ANG_TOL
        {
            CellType::Hexagonal
        } else if self.ab_close_c_different() && self.orthogonal_angles() {
            CellType::Tetragonal
        } else if self.abc_different() && self.orthogonal_angles() {
            CellType::Orthorhombic
        } else if self.abc_different()
            && (self.alpha() - FRAC_PI_2).abs() < DEFAULT_ANG_TOL
            && (self.beta() - FRAC_PI_2).abs() > DEFAULT_ANG_TOL
            && (self.gamma() - FRAC_PI_2).abs() < DEFAULT_ANG_TOL
        {
            CellType::Monoclinic
        } else {
            CellType::Triclinic
        }
    }

    pub fn inverse_matrix(&self) -> &Matrix3<f64> {
        &self.inv_matrix
    }

    pub fn to_cartesian(&self, frac_coords: &MatrixXx3<f64>) -> MatrixXx3<f64> {
        frac_to_cart(frac_coords, self)
    }

    pub fn to_fractional(&self, cart_coords: &MatrixXx3<f64>) -> MatrixXx3<f64> {
        cart_to_frac(cart_coords, self)
    }
}

impl CellData for UnitCell {
    fn cell_matrix(&self) -> &Matrix3<f64> {
        self.matrix()
    }

    fn inv_cell_matrix(&self) -> &Matrix3<f64> {
        self.inverse_matrix()
    }

    fn volume(&self) -> f64 {
        self.volume
    }
}

impl CellData for Matrix3<f64> {
    fn cell_matrix(&self) -> &Matrix3<f64> {
        self
    }

    fn inv_cell_matrix(&self) -> &Matrix3<f64> {
        self
    }

    fn volume(&self) -> f64 {
        self.determinant()
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
        assert_eq!(cell.cell_type(), CellType::Rombohedral);

        let cell = UnitCell::hexagonal(1.0, 2.0);
        assert_eq!(cell.cell_type(), CellType::Hexagonal);

        let cell = UnitCell::tetragonal(1.0, 2.0);
        assert_eq!(cell.cell_type(), CellType::Tetragonal);

        let cell = UnitCell::orthorhombic(1.0, 2.0, 3.0);
        assert_eq!(cell.cell_type(), CellType::Orthorhombic);
    }
}
