use crate::data::{DEFAULT_ANG_TOL, DEFAULT_DIST_TOL};
use crate::unitcell::CellType;
use crate::utils::wrap_coordinates_in_place;
use core::f64::consts::{FRAC_2_PI, FRAC_PI_2};
use nalgebra::{Matrix3, MatrixXx3, RowVector3, Vector3};

pub trait CellData {
    /// Cell vectors stored in rows of a 3x3 matrix.
    fn cell_matrix(&self) -> Matrix3<f64>;
    fn inv_cell_matrix(&self) -> Matrix3<f64>;
    fn volume(&self) -> f64;
    fn vec_a(&self) -> RowVector3<f64> {
        self.cell_matrix().row(0).into()
    }
    fn vec_b(&self) -> RowVector3<f64> {
        self.cell_matrix().row(1).into()
    }
    fn vec_c(&self) -> RowVector3<f64> {
        self.cell_matrix().row(2).into()
    }
    fn lengths(&self) -> [f64; 3];
    fn angles(&self) -> [f64; 3];
    fn a(&self) -> f64 {
        self.lengths()[0]
    }
    fn b(&self) -> f64 {
        self.lengths()[1]
    }
    fn c(&self) -> f64 {
        self.lengths()[2]
    }
    fn alpha(&self) -> f64 {
        self.angles()[0]
    }
    fn beta(&self) -> f64 {
        self.angles()[1]
    }
    fn gamma(&self) -> f64 {
        self.angles()[2]
    }
    fn to_cartesian(&self, frac_coords: &MatrixXx3<f64>) -> MatrixXx3<f64>;
    fn to_fractional(&self, cart_coords: &MatrixXx3<f64>) -> MatrixXx3<f64>;
    fn cell_type(&self) -> CellType {
        let [a, b, c] = self.lengths();
        let [alpha, beta, gamma] = self.angles();
        if (alpha - FRAC_PI_2).abs() < DEFAULT_ANG_TOL
            && (beta - FRAC_PI_2).abs() < DEFAULT_ANG_TOL
            && (gamma - FRAC_PI_2).abs() < DEFAULT_ANG_TOL
        {
            if (a - b).abs() < DEFAULT_DIST_TOL && (b - c).abs() < DEFAULT_DIST_TOL {
                CellType::Cubic
            } else if (a - b).abs() < DEFAULT_DIST_TOL {
                CellType::Tetragonal
            } else {
                CellType::Orthorhombic
            }
        } else if (alpha - FRAC_PI_2).abs() > DEFAULT_ANG_TOL
            && (beta - FRAC_PI_2).abs() > DEFAULT_ANG_TOL
            && (gamma - FRAC_PI_2).abs() > DEFAULT_ANG_TOL
            && (a - b).abs() < DEFAULT_DIST_TOL
            && (b - c).abs() < DEFAULT_DIST_TOL
        {
            CellType::Rhombohedral
        } else if (alpha - FRAC_PI_2).abs() < DEFAULT_ANG_TOL
            && (beta - FRAC_PI_2).abs() < DEFAULT_ANG_TOL
            && (gamma - FRAC_2_PI / 3.0).abs() < DEFAULT_ANG_TOL
        {
            CellType::Hexagonal
        } else if (alpha - FRAC_PI_2).abs() < DEFAULT_ANG_TOL
            && (gamma - FRAC_PI_2).abs() < DEFAULT_ANG_TOL
        {
            CellType::Monoclinic
        } else {
            CellType::Triclinic
        }
    }
    fn is_cubic(&self) -> bool {
        self.cell_type() == CellType::Cubic
    }
    fn is_tetragonal(&self) -> bool {
        self.cell_type() == CellType::Tetragonal
    }
    fn is_orthorhombic(&self) -> bool {
        self.cell_type() == CellType::Orthorhombic
    }
    fn is_hexagonal(&self) -> bool {
        self.cell_type() == CellType::Hexagonal
    }
    fn is_rhombohedral(&self) -> bool {
        self.cell_type() == CellType::Rhombohedral
    }
    fn is_monoclinic(&self) -> bool {
        self.cell_type() == CellType::Monoclinic
    }
    fn is_triclinic(&self) -> bool {
        self.cell_type() == CellType::Triclinic
    }
}

pub trait SupercellData {
    /// Unit cells along a lattice vector.
    fn n(&self) -> usize;
    /// Unit cells along b lattice vector.
    fn m(&self) -> usize;
    /// Unit cells along c lattice vector.
    fn l(&self) -> usize;
}

pub trait AtomicData {
    fn atomic_nums(&self) -> &[u8];
    fn n_atoms(&self) -> usize {
        self.atomic_nums().len()
    }
    fn covalent_radii(&self) -> Vec<f64>;
}

pub trait CartAtomicData: AtomicData {
    fn cartesian_coords(&self) -> &MatrixXx3<f64>;
    fn com(&self) -> Vector3<f64> {
        let coords = self.cartesian_coords();
        let n_atoms = self.n_atoms() as f64;
        let sum_coords = coords.row_sum();
        Vector3::new(
            sum_coords[0] / n_atoms,
            sum_coords[1] / n_atoms,
            sum_coords[2] / n_atoms,
        )
    }
    /// Computes the moment of inertia tensor for the atoms.
    fn moi(&self) -> Matrix3<f64> {
        let coords = self.cartesian_coords();
        let n_atoms = self.n_atoms() as f64;
        let com = self.com();
        let mut moi = Matrix3::zeros();

        for i in 0..self.n_atoms() {
            let r = coords.row(i).transpose() - com;
            moi += r * r.transpose();
        }

        moi / n_atoms
    }
    // fn rotate(&self, rotation_matrix: &Matrix3<f64>, rotation_center: &[f64; 3]);
    // fn translate(&self, shift_vector: &[f64; 3]);
    // fn scale(&self, scale_factor: f64);
}

pub trait PeriodicAtomicData: AtomicData + CellData + CartAtomicData {
    fn fractional_coords(&self) -> &MatrixXx3<f64>;
    fn density(&self) -> f64;
    // fn wrap_fractional_coords(&mut self) {
    //     wrap_coordinates_in_place(self.fractional_coords_mut());
    // }
}
