use nalgebra::{Matrix3, MatrixXx3};

pub trait CellData {
    fn cell_matrix(&self) -> &Matrix3<f64>;
    fn inv_cell_matrix(&self) -> &Matrix3<f64>;
    fn volume(&self) -> f64;
}

pub trait AtomicData {
    fn atomic_nums(&self) -> &[u8];
    fn n_atoms(&self) -> usize {
        self.atomic_nums().len()
    }
}

pub trait CartAtomicData {
    fn cartesian_coords(&self) -> &MatrixXx3<f64>;
}

pub trait FracAtomicData {
    fn fractional_coords(&self) -> &MatrixXx3<f64>;
}
