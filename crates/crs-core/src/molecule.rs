use nalgebra::MatrixXx3;
use std::collections::HashMap;

#[derive(Debug, Clone)]
pub enum BondSettings {
    /// Use the default covalent radii and a tolerance of 0.45 Å.
    Default,
    /// Use the default covalent radii and the provided tolerance.
    DefaultWithTolerance(f64),
    /// Use the provided covalent radii and tolerance.
    CustomRadiiAndTolerance(HashMap<u8, f64>, f64),
}

#[derive(Debug, Clone)]
pub struct Molecule {
    atomic_nums: Vec<u8>,
    cartesian_coords: MatrixXx3<f64>,
}

impl Molecule {
    pub fn new(atomic_nums: Vec<u8>, cartesian_coords: MatrixXx3<f64>) -> Self {
        assert_eq!(
            atomic_nums.len(),
            cartesian_coords.nrows(),
            "The number of atoms must match the number of coordinates."
        );
        Self {
            atomic_nums,
            cartesian_coords,
        }
    }
}