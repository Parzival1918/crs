use nalgebra::MatrixXx3;

#[derive(Debug, Clone)]
pub struct Molecule {
    pub atomic_nums: Vec<u8>,
    pub cartesian_coords: MatrixXx3<f64>,
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