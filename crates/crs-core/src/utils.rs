use crate::data::{get_atomic_mass, get_covalent_radius, get_species_name};
use crate::molecule::{BondSettings, Molecule};
use crate::spacegroup::SymOp;
use crate::traits::{AtomicData, CellData, FracAtomicData};
use kdtree::KdTree;
use kdtree::distance::squared_euclidean;
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
        let mut translated = rotated;
        let trans_row = op.translation().transpose();
        for mut row in translated.row_iter_mut() {
            row += trans_row;
        }
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

/// A simple Union-Find (disjoint set) data structure for grouping atoms into molecules.
struct UnionFind {
    parent: Vec<usize>,
    rank: Vec<usize>,
}

impl UnionFind {
    fn new(n: usize) -> Self {
        Self {
            parent: (0..n).collect(),
            rank: vec![0; n],
        }
    }

    fn find(&mut self, x: usize) -> usize {
        if self.parent[x] != x {
            self.parent[x] = self.find(self.parent[x]);
        }
        self.parent[x]
    }

    fn union(&mut self, x: usize, y: usize) {
        let rx = self.find(x);
        let ry = self.find(y);
        if rx == ry {
            return;
        }
        match self.rank[rx].cmp(&self.rank[ry]) {
            std::cmp::Ordering::Less => self.parent[rx] = ry,
            std::cmp::Ordering::Greater => self.parent[ry] = rx,
            std::cmp::Ordering::Equal => {
                self.parent[ry] = rx;
                self.rank[rx] += 1;
            }
        }
    }
}

/// Resolve the covalent radii and tolerance from `BondSettings`.
///
/// Returns `(radii_per_atom, tolerance, max_cutoff)` where `max_cutoff` is the
/// largest possible bonding distance across all atom pairs.
fn resolve_bond_params(atomic_nums: &[u8], settings: &BondSettings) -> (Vec<f64>, f64, f64) {
    let (radii, tolerance) = match settings {
        BondSettings::Default => {
            let radii: Vec<f64> = atomic_nums
                .iter()
                .map(|&z| get_covalent_radius(z))
                .collect();
            (radii, 0.45)
        }
        BondSettings::DefaultWithTolerance(tol) => {
            let radii: Vec<f64> = atomic_nums
                .iter()
                .map(|&z| get_covalent_radius(z))
                .collect();
            (radii, *tol)
        }
        BondSettings::CustomRadiiAndTolerance(custom_radii, tol) => {
            let radii: Vec<f64> = atomic_nums
                .iter()
                .map(|&z| {
                    custom_radii
                        .get(&z)
                        .copied()
                        .unwrap_or_else(|| get_covalent_radius(z))
                })
                .collect();
            (radii, *tol)
        }
    };

    let max_radius = radii.iter().cloned().fold(0.0_f64, f64::max);
    let max_cutoff = 2.0 * max_radius + tolerance;
    (radii, tolerance, max_cutoff)
}

/// Compute the fractional-space cutoff distance along each cell axis.
///
/// For each axis `i`, compute how far `max_cutoff` Å extends in fractional
/// coordinates by measuring the Cartesian length of each inverse-cell row
/// (which converts Cartesian distances to fractional ones).
fn compute_fractional_cutoffs(cell: &impl CellData, max_cutoff: f64) -> [f64; 3] {
    let inv = cell.inv_cell_matrix();
    let mut frac_cutoffs = [0.0; 3];
    for i in 0..3 {
        // The i-th column of the inverse matrix maps a unit Cartesian displacement
        // to fractional space. The norm of this column gives the conversion factor.
        let col = inv.column(i);
        frac_cutoffs[i] = max_cutoff * col.norm();
    }
    frac_cutoffs
}

/// Find molecules in a set of atoms, accept an Option<CellData> that if provided
/// will use periodicity in the search for molecules. If None is provided,
/// it will treat the atoms as a cluster and find molecules without periodicity.
pub fn find_molecules<A: AtomicData + FracAtomicData, C: CellData>(
    atoms: &A,
    cell: Option<&C>,
    settings: BondSettings,
) -> Vec<Molecule> {
    let atomic_nums = atoms.atomic_nums();
    let n_atoms = atomic_nums.len();
    if n_atoms == 0 {
        return Vec::new();
    }

    let (radii, tolerance, max_cutoff) = resolve_bond_params(atomic_nums, &settings);

    // Convert fractional coordinates to Cartesian.
    let frac_coords = atoms.fractional_coords();
    let cart_coords = match cell {
        Some(c) => frac_to_cart(frac_coords, c),
        None => frac_coords.clone(),
    };

    // Build the kd-tree with all real atoms and (optionally) periodic ghost images.
    // Each entry stores (cartesian_position, original_atom_index).
    let mut tree: KdTree<f64, usize, [f64; 3]> = KdTree::new(3);
    for i in 0..n_atoms {
        let pt = [
            cart_coords[(i, 0)],
            cart_coords[(i, 1)],
            cart_coords[(i, 2)],
        ];
        tree.add(pt, i).expect("Failed to add point to kd-tree");
    }

    // Add periodic ghost images for atoms near cell boundaries.
    if let Some(c) = cell {
        let frac_cutoffs = compute_fractional_cutoffs(c, max_cutoff);
        let cell_matrix = c.cell_matrix();

        for i in 0..n_atoms {
            let fx = frac_coords[(i, 0)];
            let fy = frac_coords[(i, 1)];
            let fz = frac_coords[(i, 2)];

            // For each axis, determine which periodic shifts are needed.
            // An atom at fractional coord `f` is near the lower boundary if `f < cutoff`,
            // and near the upper boundary if `f > 1 - cutoff`.
            let mut shifts_x: Vec<f64> = Vec::with_capacity(3);
            shifts_x.push(0.0); // always include original (no shift)
            if fx < frac_cutoffs[0] {
                shifts_x.push(1.0);
            }
            if fx > 1.0 - frac_cutoffs[0] {
                shifts_x.push(-1.0);
            }

            let mut shifts_y: Vec<f64> = Vec::with_capacity(3);
            shifts_y.push(0.0);
            if fy < frac_cutoffs[1] {
                shifts_y.push(1.0);
            }
            if fy > 1.0 - frac_cutoffs[1] {
                shifts_y.push(-1.0);
            }

            let mut shifts_z: Vec<f64> = Vec::with_capacity(3);
            shifts_z.push(0.0);
            if fz < frac_cutoffs[2] {
                shifts_z.push(1.0);
            }
            if fz > 1.0 - frac_cutoffs[2] {
                shifts_z.push(-1.0);
            }

            // Generate ghost images from the Cartesian product of shifts,
            // skipping the (0,0,0) combination which is the original atom.
            for &sx in &shifts_x {
                for &sy in &shifts_y {
                    for &sz in &shifts_z {
                        if sx == 0.0 && sy == 0.0 && sz == 0.0 {
                            continue;
                        }
                        let ghost_frac = MatrixXx3::from_row_slice(&[fx + sx, fy + sy, fz + sz]);
                        let ghost_cart = &ghost_frac * cell_matrix;
                        let pt = [ghost_cart[(0, 0)], ghost_cart[(0, 1)], ghost_cart[(0, 2)]];
                        tree.add(pt, i)
                            .expect("Failed to add ghost point to kd-tree");
                    }
                }
            }
        }
    }

    // Find bonds using kd-tree range queries + union-find.
    let mut uf = UnionFind::new(n_atoms);
    let max_cutoff_sq = max_cutoff * max_cutoff;

    for i in 0..n_atoms {
        let pt_i = [
            cart_coords[(i, 0)],
            cart_coords[(i, 1)],
            cart_coords[(i, 2)],
        ];
        let neighbors = tree
            .within(&pt_i, max_cutoff_sq, &squared_euclidean)
            .expect("kd-tree query failed");

        let r_i = radii[i];
        for &(ref dist_sq, &j) in &neighbors {
            if j <= i {
                // Avoid duplicate pair processing and self-bonds.
                continue;
            }
            let r_j = radii[j];
            let bond_cutoff = r_i + r_j + tolerance;
            if *dist_sq < bond_cutoff * bond_cutoff {
                uf.union(i, j);
            }
        }
    }

    // Group atoms by their connected-component root.
    let mut components: std::collections::HashMap<usize, Vec<usize>> =
        std::collections::HashMap::new();
    for i in 0..n_atoms {
        let root = uf.find(i);
        components.entry(root).or_default().push(i);
    }

    // Build Molecule objects from each component.
    let mut molecules: Vec<Molecule> = components
        .into_values()
        .map(|indices| {
            let mol_atomic_nums: Vec<u8> = indices.iter().map(|&i| atomic_nums[i]).collect();
            let mol_coords =
                MatrixXx3::from_fn(indices.len(), |row, col| cart_coords[(indices[row], col)]);
            Molecule::new(mol_atomic_nums, mol_coords)
        })
        .collect();

    // Sort by descending size for deterministic output.
    molecules.sort_by(|a, b| b.n_atoms().cmp(&a.n_atoms()));
    molecules
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

        let frac_coords = MatrixXx3::from_row_slice(&[0.5, 0.5, 0.5, 0.8, 0.2, 0.1]);
        let cart_coords = frac_to_cart(&frac_coords, &cell_matrix);
        assert_eq!(
            cart_coords,
            MatrixXx3::from_row_slice(&[5.0, 5.0, 5.0, 8.0, 2.0, 1.0])
        );

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

    // =========================================================================
    // find_molecules tests
    // =========================================================================

    /// A simple test struct implementing AtomicData + FracAtomicData.
    struct TestAtoms {
        atomic_nums: Vec<u8>,
        frac_coords: MatrixXx3<f64>,
    }

    impl AtomicData for TestAtoms {
        fn atomic_nums(&self) -> &[u8] {
            &self.atomic_nums
        }
    }

    impl FracAtomicData for TestAtoms {
        fn fractional_coords(&self) -> &MatrixXx3<f64> {
            &self.frac_coords
        }
    }

    #[test]
    fn test_find_molecules_empty() {
        let atoms = TestAtoms {
            atomic_nums: vec![],
            frac_coords: MatrixXx3::zeros(0),
        };
        let uc = UnitCell::cubic(10.0);
        let molecules = find_molecules(&atoms, Some(&uc), BondSettings::Default);
        assert!(molecules.is_empty());
    }

    #[test]
    fn test_find_molecules_isolated_atoms() {
        // Two atoms far apart in a large cubic cell — should be two separate molecules.
        let atoms = TestAtoms {
            atomic_nums: vec![6, 8], // C and O
            frac_coords: MatrixXx3::from_row_slice(&[
                0.1, 0.1, 0.1, // C at (1,1,1) Å in a 10 Å cell
                0.9, 0.9, 0.9, // O at (9,9,9) Å — far away
            ]),
        };
        let uc = UnitCell::cubic(10.0);
        let molecules = find_molecules(&atoms, Some(&uc), BondSettings::Default);
        assert_eq!(molecules.len(), 2);
        assert_eq!(molecules[0].n_atoms(), 1);
        assert_eq!(molecules[1].n_atoms(), 1);
    }

    #[test]
    fn test_find_molecules_simple_no_periodicity() {
        // C-O bond: covalent radii C=0.76, O=0.66 => bond if dist < 0.76+0.66+0.45 = 1.87 Å
        // Place them 1.3 Å apart (well within bonding distance).
        let uc = UnitCell::cubic(10.0);
        let atoms = TestAtoms {
            atomic_nums: vec![6, 8], // C and O
            frac_coords: MatrixXx3::from_row_slice(&[
                0.5, 0.5, 0.5, // C at center
                0.5, 0.5, 0.63, // O at 0.13*10 = 1.3 Å away along z
            ]),
        };
        let molecules = find_molecules(&atoms, Some(&uc), BondSettings::Default);
        assert_eq!(molecules.len(), 1);
        assert_eq!(molecules[0].n_atoms(), 2);
    }

    #[test]
    fn test_find_molecules_periodic_bonding() {
        // Two atoms at opposite edges of the cell that bond through the periodic boundary.
        // H-H bond: covalent radii H=0.31 => bond if dist < 0.31+0.31+0.45 = 1.07 Å
        // Place one H at frac 0.02 and another at frac 0.95 along x in a 10 Å cell.
        // Direct distance: (0.95 - 0.02) * 10 = 9.3 Å (too far)
        // Periodic distance: (1.0 - 0.95 + 0.02) * 10 = 0.7 Å (bonds!)
        let uc = UnitCell::cubic(10.0);
        let atoms = TestAtoms {
            atomic_nums: vec![1, 1], // H and H
            frac_coords: MatrixXx3::from_row_slice(&[
                0.02, 0.5, 0.5, // H near lower x boundary
                0.95, 0.5, 0.5, // H near upper x boundary
            ]),
        };
        let molecules = find_molecules(&atoms, Some(&uc), BondSettings::Default);
        assert_eq!(
            molecules.len(),
            1,
            "Atoms should bond across the periodic boundary"
        );
        assert_eq!(molecules[0].n_atoms(), 2);
    }

    #[test]
    fn test_find_molecules_periodic_no_bond_without_cell() {
        // Same atoms as above, but without a cell (no periodicity).
        // Without periodicity, fractional coords are treated as Cartesian directly.
        // The atoms at frac 0.02 and 0.95 are 0.93 units apart (treated as Å).
        // H-H bond cutoff is 1.07 Å, so they WOULD bond. Use a different setup:
        // place them far enough apart that they don't bond without periodicity.
        let uc = UnitCell::cubic(10.0);

        // With the cell, atoms bond across the boundary. Without cell, they don't
        // (because direct Cartesian distance through the cell is too large).
        let atoms = TestAtoms {
            atomic_nums: vec![1, 1],
            frac_coords: MatrixXx3::from_row_slice(&[0.02, 0.5, 0.5, 0.95, 0.5, 0.5]),
        };

        // With periodicity: should bond
        let mols_periodic = find_molecules(&atoms, Some(&uc), BondSettings::Default);
        assert_eq!(mols_periodic.len(), 1);

        // Without periodicity: coords are raw fractional values (tiny numbers),
        // the distance is 0.93 which is < 1.07 for H-H. So use a more extreme case.
        let atoms_far = TestAtoms {
            atomic_nums: vec![1, 1],
            frac_coords: MatrixXx3::from_row_slice(&[
                0.0, 0.0, 0.0, 5.0, 5.0,
                5.0, // In no-cell mode, these are Cartesian: very far apart
            ]),
        };
        let mols_no_cell = find_molecules::<_, UnitCell>(&atoms_far, None, BondSettings::Default);
        assert_eq!(
            mols_no_cell.len(),
            2,
            "Without a cell, far atoms should not bond"
        );
    }

    #[test]
    fn test_find_molecules_custom_radii_and_tolerance() {
        // Place two C atoms at distance ~1.6 Å.
        // Default C radius = 0.76, so default bond cutoff = 0.76+0.76+0.45 = 1.97 Å -> bonds.
        // With a tiny tolerance of 0.0 Å, cutoff = 0.76+0.76+0.0 = 1.52 Å -> no bond at 1.6 Å.
        let uc = UnitCell::cubic(10.0);
        let atoms = TestAtoms {
            atomic_nums: vec![6, 6],
            frac_coords: MatrixXx3::from_row_slice(&[
                0.5, 0.5, 0.5, 0.5, 0.5, 0.66, // 1.6 Å apart along z
            ]),
        };

        // Default settings: should bond
        let mols = find_molecules(&atoms, Some(&uc), BondSettings::Default);
        assert_eq!(mols.len(), 1);

        // Zero tolerance: should NOT bond at 1.6 Å (cutoff = 1.52)
        let mols = find_molecules(&atoms, Some(&uc), BondSettings::DefaultWithTolerance(0.0));
        assert_eq!(
            mols.len(),
            2,
            "With zero tolerance, 1.6 Å exceeds C-C cutoff of 1.52 Å"
        );

        // Custom radii: set C radius to 1.0 Å, tolerance 0.0 => cutoff = 2.0 Å -> bonds
        let mut custom_radii = std::collections::HashMap::new();
        custom_radii.insert(6u8, 1.0);
        let mols = find_molecules(
            &atoms,
            Some(&uc),
            BondSettings::CustomRadiiAndTolerance(custom_radii, 0.0),
        );
        assert_eq!(
            mols.len(),
            1,
            "Custom radii of 1.0 gives cutoff 2.0 Å, which bonds at 1.6 Å"
        );
    }

    #[test]
    fn test_find_molecules_multiple_molecules() {
        // Three groups of atoms in a 20 Å cubic cell:
        // Group 1: C-O pair (bonded)
        // Group 2: Isolated N
        // Group 3: H-H pair (bonded)
        let uc = UnitCell::cubic(20.0);
        let atoms = TestAtoms {
            atomic_nums: vec![6, 8, 7, 1, 1],
            frac_coords: MatrixXx3::from_row_slice(&[
                0.1, 0.1, 0.1, // C
                0.1, 0.1, 0.16, // O at 1.2 Å from C (bonds: 0.76+0.66+0.45=1.87)
                0.5, 0.5, 0.5, // N (isolated)
                0.9, 0.9, 0.1, // H
                0.9, 0.9, 0.13, // H at 0.6 Å from other H (bonds: 0.31+0.31+0.45=1.07)
            ]),
        };
        let molecules = find_molecules(&atoms, Some(&uc), BondSettings::Default);
        assert_eq!(molecules.len(), 3);
        // Sorted by descending size: two 2-atom molecules, then one 1-atom molecule
        assert_eq!(molecules[0].n_atoms(), 2);
        assert_eq!(molecules[1].n_atoms(), 2);
        assert_eq!(molecules[2].n_atoms(), 1);
    }
}
