use nalgebra::MatrixXx3;
use std::collections::HashMap;
use crate::data::get_default_covalent_radii_map;
use crate::traits::{AtomicData, CartAtomicData};

#[derive(Debug, Clone)]
pub enum BondSettings {
    /// Use the default covalent radii and a tolerance of 0.45 Å.
    Default,
    /// Use the default covalent radii and the provided tolerance.
    DefaultWithTolerance(f64),
    /// Use the provided covalent radii and tolerance.
    CustomRadiiAndTolerance(HashMap<u8, f64>, f64),
}

impl BondSettings {
    pub fn tolerance(&self) -> f64 {
        match self {
            BondSettings::Default => 0.45,
            BondSettings::DefaultWithTolerance(tol) => *tol,
            BondSettings::CustomRadiiAndTolerance(_, tol) => *tol,
        }
    }

    pub fn covalent_radii(&self) -> HashMap<u8, f64> {
        match self {
            BondSettings::Default | BondSettings::DefaultWithTolerance(_) => {
                get_default_covalent_radii_map()
            }
            BondSettings::CustomRadiiAndTolerance(radii, _) => radii.clone(),
        }
    }

    /// Returns `(radii_per_atom, tolerance, max_cutoff)` where `max_cutoff` is the
    /// largest possible bonding distance across all atom pairs.
    pub fn resolve_bond_params(&self, atomic_nums: &[u8]) -> (Vec<f64>, f64, f64) {
        let radii_map = self.covalent_radii();
        let tolerance = self.tolerance();
        let mut radii_per_atom = Vec::with_capacity(atomic_nums.len());
        let mut max_cutoff: f64 = 0.0;

        for &atomic_num in atomic_nums {
            let radius = radii_map
                .get(&atomic_num)
                .copied()
                .unwrap_or_else(|| panic!("Covalent radius not found for atomic number: {}", atomic_num));
            radii_per_atom.push(radius);
            max_cutoff = max_cutoff.max(2.0 * radius + tolerance);
        }

        (radii_per_atom, tolerance, max_cutoff)
    }
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

impl AtomicData for Molecule {
    fn atomic_nums(&self) -> &[u8] {
        &self.atomic_nums
    }
}

impl CartAtomicData for Molecule {
    fn cartesian_coords(&self) -> &MatrixXx3<f64> {
        &self.cartesian_coords
    }
}