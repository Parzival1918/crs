use crate::utils::cart_to_frac;
use nalgebra::{Matrix3, MatrixXx3};

pub trait CellData {
    fn cell_matrix(&self) -> &Matrix3<f64>;
    fn inv_cell_matrix(&self) -> &Matrix3<f64>;
    fn volume(&self) -> f64;
}

pub trait AtomicData {
    fn atomic_nums(&self) -> &[u8];
}

pub trait CartAtomicData {
    fn cartesian_coords(&self) -> &MatrixXx3<f64>;
}

pub trait FracAtomicData {
    fn fractional_coords(&self) -> &MatrixXx3<f64>;
}
