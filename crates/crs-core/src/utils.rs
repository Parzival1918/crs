use crate::data::{get_atomic_mass, get_species_name};
use crate::spacegroup::SymOp;
use crate::traits::CellData;
use nalgebra::{Matrix3, MatrixXx3};

/// Generate a chemical formula from a list of atomic numbers
/// using the Hill system:
///
/// 1. Carbon atoms
/// 2. Hydrogen atoms
/// 3. All other elements in alphabetical order
///
/// The `mult` parameter is used to multiply the counts of each element in the formula.
pub fn chemical_formula(atomic_nums: &[u8], mult: usize) -> String {
    let mut counts = std::collections::HashMap::new();
    for &num in atomic_nums {
        *counts.entry(num).or_insert(0) += 1 * mult;
    }

    let mut formula = String::new();

    // Add carbon atoms first
    if let Some(count) = counts.remove(&6) {
        formula.push_str("C");
        if count > 1 {
            formula.push_str(&count.to_string());
        }
    }

    // Add hydrogen atoms next
    if let Some(count) = counts.remove(&1) {
        formula.push_str("H");
        if count > 1 {
            formula.push_str(&count.to_string());
        }
    }

    // Add all other elements in alphabetical order
    let mut remaining: Vec<_> = counts.iter().collect();
    remaining.sort_by_key(|&(num, _)| get_species_name(*num));
    for (num, count) in remaining {
        let species = get_species_name(*num);
        formula.push_str(species);
        if *count > 1 {
            formula.push_str(&count.to_string());
        }
    }

    formula
}

/// Calculate the molar mass of a compound given a list of atomic numbers
pub fn molar_mass(atomic_nums: &[u8]) -> f64 {
    atomic_nums.iter().map(|&num| get_atomic_mass(num)).sum()
}

/// Calculate the cell matrix, inverse matrix, and volume from lengths and angles.
///
/// # Arguments
/// * `lengths` - The lengths of the unit cell (a, b, c)
/// * `angles` - The angles of the unit cell (alpha, beta, gamma) in radians
///
/// The cell vectors are stored in the rows of the matrix.
pub fn cell_matrix_and_volume(
    lengths: [f64; 3],
    angles: [f64; 3],
) -> (Matrix3<f64>, Matrix3<f64>, f64) {
    let (a, b, c) = (lengths[0], lengths[1], lengths[2]);
    let (alpha, beta, gamma) = (angles[0], angles[1], angles[2]);

    let ca = alpha.cos();
    let cb = beta.cos();
    let cg = gamma.cos();
    let sg = gamma.sin();

    let volume = a * b * c * (1.0 - ca * ca - cb * cb - cg * cg + 2.0 * ca * cb * cg).sqrt();

    (
        Matrix3::new(
            a,
            0.0,
            0.0,
            b * cg,
            b * sg,
            0.0,
            c * cb,
            c * (ca - cb * cg) / sg,
            volume / (a * b * sg),
        ),
        Matrix3::new(
            1.0 / a,
            0.0,
            0.0,
            -cg / (a * sg),
            1.0 / (b * sg),
            0.0,
            b * c * (ca * cg - cb) / (volume * sg),
            a * c * (cb * cg - ca) / (volume * sg),
            a * b * sg / volume,
        ),
        volume,
    )
}

/// Convert fractional coordinates to Cartesian coordinates using the cell matrix.
pub fn frac_to_cart(frac_coords: &MatrixXx3<f64>, cell: &impl CellData) -> MatrixXx3<f64> {
    frac_coords * cell.cell_matrix()
}

/// Convert Cartesian coordinates to fractional coordinates using the inverse of the cell matrix.
pub fn cart_to_frac(cart_coords: &MatrixXx3<f64>, cell: &impl CellData) -> MatrixXx3<f64> {
    cart_coords * cell.inv_cell_matrix()
}

/// Apply symmetry operations to fractional coordinates.
///
/// The function applies the unit symmetry operation.
pub fn apply_symops(frac_coords: &MatrixXx3<f64>, sym_ops: &[SymOp]) -> MatrixXx3<f64> {
    let mut new_coords = MatrixXx3::zeros(frac_coords.nrows() * sym_ops.len());
    for (i, op) in sym_ops.iter().enumerate() {
        let rotated = frac_coords * op.rotation();
        let translated = rotated + op.translation().transpose();
        new_coords
            .view_mut((i * frac_coords.nrows(), 0), (frac_coords.nrows(), 3))
            .copy_from(&translated);
    }
    new_coords
}

/// Wrap fractional coordinates to the range [0, 1).
pub fn wrap_coordinates_in_place(frac_coords: &mut MatrixXx3<f64>) {
    frac_coords.apply(|x| *x -= x.floor())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::spacegroup::SpaceGroup;
    use crate::unitcell::UnitCell;
    use std::f64::consts::FRAC_PI_2;

    #[test]
    fn test_chemical_formula() {
        let atomic_nums = &[6, 1, 1, 1, 8]; // C, H, H, H, O
        let formula = chemical_formula(atomic_nums, 1);
        assert_eq!(formula, "CH3O");

        let atomic_nums = &[1, 6, 1, 8, 1]; // H, C, H, O, H
        let formula = chemical_formula(atomic_nums, 1);
        assert_eq!(formula, "CH3O");

        let atomic_nums = &[1, 6, 1, 8, 1]; // H, C, H, O, H
        let formula = chemical_formula(atomic_nums, 4);
        assert_eq!(formula, "C4H12O4");
    }

    #[test]
    fn test_molar_mass() {
        let atomic_nums = &[6, 1, 1, 1, 8]; // C, H, H, H, O
        let mass = molar_mass(atomic_nums);
        assert_eq!(mass, 12.011 + 1.008 * 3.0 + 15.999);
    }

    #[test]
    fn test_cell_matrix_and_volume() {
        let lengths = [1.0, 1.0, 1.0];
        let angles = [FRAC_PI_2, FRAC_PI_2, FRAC_PI_2]; // 90 degrees in radians
        let (matrix, inv_matrix, volume) = cell_matrix_and_volume(lengths, angles);
        assert_eq!(volume, 1.0);

        let expected_matrix = Matrix3::new(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
        // assert_eq!(matrix, expected_matrix);

        let lengths = [6.0, 10.0, 2.0];
        let angles = [FRAC_PI_2, FRAC_PI_2, FRAC_PI_2]; // 90 degrees in radians
        let (matrix, inv_matrix, volume) = cell_matrix_and_volume(lengths, angles);
        assert_eq!(volume, 6.0 * 10.0 * 2.0);

        let lengths = [5.0, 5.0, 5.0];
        let angles = [
            90.0_f64.to_radians(),
            60.0_f64.to_radians(),
            90.0_f64.to_radians(),
        ]; // 90 degrees in radians
        let (matrix, inv_matrix, volume) = cell_matrix_and_volume(lengths, angles);
        assert!((volume - 108.2532).abs() < 1e-4); // volume from https://www.gallixa.com/XRPD/
    }

    #[test]
    fn test_frac_to_cart() {
        // Check coords in https://ic50.org/fractorth/
        let lengths = [10.0, 10.0, 10.0];
        let angles = [FRAC_PI_2, FRAC_PI_2, FRAC_PI_2];
        let (cell_matrix, _, volume) = cell_matrix_and_volume(lengths, angles);

        let frac_coords = MatrixXx3::from_row_slice(&[0.5, 0.5, 0.5]);
        let cart_coords = frac_to_cart(&frac_coords, &cell_matrix);
        assert_eq!(cart_coords, MatrixXx3::from_row_slice(&[5.0, 5.0, 5.0]));

        let uc = UnitCell::monoclinic(10.0, 5.0, 8.0, 75.0_f64.to_radians());
        let frac_coords = MatrixXx3::from_row_slice(&[0.5, 0.5, 0.5]);
        let cart_coords = frac_to_cart(&frac_coords, &uc);
        let expected_cart_coords = MatrixXx3::from_row_slice(&[6.035276, 2.500000, 3.863703]);
        assert!((cart_coords - expected_cart_coords).norm() < 1e-6);
    }

    #[test]
    fn test_cart_to_frac() {
        let lengths = [10.0, 10.0, 10.0];
        let angles = [FRAC_PI_2, FRAC_PI_2, FRAC_PI_2];
        let (cell_matrix, inv_matrix, _) = cell_matrix_and_volume(lengths, angles);

        let cart_coords = MatrixXx3::from_row_slice(&[5.0, 5.0, 5.0]);
        let frac_coords = cart_to_frac(&cart_coords, &inv_matrix);
        assert!((frac_coords - MatrixXx3::from_row_slice(&[0.5, 0.5, 0.5])).norm() < 1e-6);

        let uc = UnitCell::monoclinic(10.0, 5.0, 8.0, 75.0_f64.to_radians());
        let cart_coords = MatrixXx3::from_row_slice(&[6.035276, 2.500000, 3.863703]);
        let frac_coords = cart_to_frac(&cart_coords, &uc);
        let expected_frac_coords = MatrixXx3::from_row_slice(&[0.5, 0.5, 0.5]);
        assert!((frac_coords - expected_frac_coords).norm() < 1e-6);
    }

    #[test]
    fn test_apply_symmetry_operations() {
        let lengths = [10.0, 10.0, 10.0];
        let angles = [FRAC_PI_2, FRAC_PI_2, FRAC_PI_2];
        let (cell_matrix, _, volume) = cell_matrix_and_volume(lengths, angles);

        let frac_coords = MatrixXx3::from_row_slice(&[0.5, 0.5, 0.5]);
        let sg = SpaceGroup::default_setting(1);
        let sym_ops = sg.operations();
        let new_coords = apply_symops(&frac_coords, &sym_ops);
        assert_eq!(new_coords.nrows(), 1);
        assert_eq!(new_coords, frac_coords);

        let frac_coords = MatrixXx3::from_row_slice(&[0.25, 0.25, 0.25]);
        let sg = SpaceGroup::default_setting(2);
        let sym_ops = sg.operations();
        let new_coords = apply_symops(&frac_coords, &sym_ops);
        assert_eq!(new_coords.nrows(), 2);
        println!("New coords: {}", &new_coords);
        assert_eq!(
            new_coords,
            MatrixXx3::from_row_slice(&[0.25, 0.25, 0.25, -0.25, -0.25, -0.25])
        );
    }

    #[test]
    fn test_wrap_coordinates_in_place() {
        let frac_coords = MatrixXx3::from_row_slice(&[0.5, 0.5, 0.5]);
        let sg = SpaceGroup::default_setting(2);
        let sym_ops = sg.operations();
        let mut new_coords = apply_symops(&frac_coords, &sym_ops);
        assert_eq!(new_coords.nrows(), 2);
        println!("New coords: {}", &new_coords);
        wrap_coordinates_in_place(&mut new_coords);
        assert_eq!(
            new_coords,
            MatrixXx3::from_row_slice(&[0.5, 0.5, 0.5, 0.5, 0.5, 0.5])
        );
    }
}
