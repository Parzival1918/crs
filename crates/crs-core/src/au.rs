use crate::data::get_covalent_radius;
use crate::utils::chemical_formula;
use nalgebra::MatrixXx3;

#[derive(Debug, Clone)]
pub struct AsymmetricUnit {
    atomic_nums: Vec<u8>,
    frac_coords: MatrixXx3<f64>,
}

impl AsymmetricUnit {
    pub fn new(atomic_nums: Vec<u8>, frac_coords: MatrixXx3<f64>) -> Self {
        assert_eq!(
            atomic_nums.len(),
            frac_coords.nrows(),
            "The number of atoms must match the number of coordinates."
        );
        Self {
            atomic_nums,
            frac_coords,
        }
    }

    pub fn atomic_nums(&self) -> &[u8] {
        &self.atomic_nums
    }

    pub fn frac_coords(&self) -> &MatrixXx3<f64> {
        &self.frac_coords
    }

    pub fn covalent_radii(&self) -> Vec<f64> {
        self.atomic_nums
            .iter()
            .map(|&num| get_covalent_radius(num))
            .collect()
    }

    pub fn chemical_formula(&self) -> String {
        chemical_formula(&self.atomic_nums, 1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::COVALENT_RADII;

    #[test]
    fn test_new() {
        let atomic_nums = vec![1, 6, 8];
        let frac_coords = MatrixXx3::repeat(3, 0.0);
        let au = AsymmetricUnit::new(atomic_nums, frac_coords);
        assert_eq!(au.atomic_nums().len(), 3);
    }

    #[test]
    fn test_covalent_radii() {
        let atomic_nums = vec![1, 6, 8];
        let frac_coords = MatrixXx3::repeat(3, 0.0);
        let au = AsymmetricUnit::new(atomic_nums, frac_coords);
        let radii = au.covalent_radii();
        assert_eq!(radii.len(), 3);
        assert_eq!(radii[0], COVALENT_RADII[1]); // H
        assert_eq!(radii[1], COVALENT_RADII[6]); // C
        assert_eq!(radii[2], COVALENT_RADII[8]); // O
    }
}
