use crate::data::{get_atomic_mass, get_covalent_radius, get_species_name};
use crate::molecule::{BondSettings, Molecule};
use crate::spacegroup::SymOp;
use crate::traits::{AtomicData, CartAtomicData, CellData};
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

/// Compute the face normals (unit vectors), perpendicular heights, and lattice
/// vectors of the unit cell. All quantities are in Cartesian space.
///
/// Returns `(normals, heights, lattice_vectors)` where:
/// - `normals[i]` is the outward unit normal to the i-th pair of opposing faces
/// - `heights[i]` is the perpendicular distance between the i-th pair of faces
/// - `lattice_vectors[i]` is the i-th lattice vector (row of cell matrix)
fn compute_cell_face_geometry(cell: &impl CellData) -> ([[f64; 3]; 3], [f64; 3], [[f64; 3]; 3]) {
    use nalgebra::Vector3;

    let m = cell.cell_matrix();
    let a = Vector3::new(m[(0, 0)], m[(0, 1)], m[(0, 2)]);
    let b = Vector3::new(m[(1, 0)], m[(1, 1)], m[(1, 2)]);
    let c = Vector3::new(m[(2, 0)], m[(2, 1)], m[(2, 2)]);

    let vol = cell.volume().abs();

    // Cross products give the (unnormalized) face normals:
    //   face 0 (bc-face, perpendicular to a): normal = b × c
    //   face 1 (ca-face, perpendicular to b): normal = c × a
    //   face 2 (ab-face, perpendicular to c): normal = a × b
    let crosses = [b.cross(&c), c.cross(&a), a.cross(&b)];

    let mut normals = [[0.0; 3]; 3];
    let mut heights = [0.0; 3];
    for i in 0..3 {
        let norm = crosses[i].norm();
        normals[i] = [
            crosses[i][0] / norm,
            crosses[i][1] / norm,
            crosses[i][2] / norm,
        ];
        heights[i] = vol / norm;
    }

    let lattice_vectors = [[a[0], a[1], a[2]], [b[0], b[1], b[2]], [c[0], c[1], c[2]]];

    (normals, heights, lattice_vectors)
}

/// Find molecules in a set of atoms, accept an Option<CellData> that if provided
/// will use periodicity in the search for molecules. If None is provided,
/// it will treat the atoms as a cluster and find molecules without periodicity.
pub fn find_molecules<A: AtomicData + CartAtomicData, C: CellData>(
    atoms: &A,
    cell: Option<&C>,
    bond_settings: BondSettings,
) -> Vec<Molecule> {
    let atomic_nums = atoms.atomic_nums();
    let n_atoms = atomic_nums.len();
    if n_atoms == 0 {
        return Vec::new();
    }

    let (radii, tolerance, max_cutoff) = bond_settings.resolve_bond_params(atomic_nums);

    let cart_coords = atoms.cartesian_coords();

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
        let (normals, heights, lattice_vecs) = compute_cell_face_geometry(c);

        for i in 0..n_atoms {
            let r = [
                cart_coords[(i, 0)],
                cart_coords[(i, 1)],
                cart_coords[(i, 2)],
            ];

            // For each pair of opposing faces, determine which periodic shifts
            // are needed by computing the perpendicular distance from the atom
            // to each face in Cartesian space.
            let mut dir_shifts: [Vec<f64>; 3] = [vec![0.0], vec![0.0], vec![0.0]];
            for d in 0..3 {
                let dist_near: f64 = r
                    .iter()
                    .zip(normals[d].iter())
                    .map(|(ri, ni)| ri * ni)
                    .sum();
                let dist_far = heights[d] - dist_near;

                if dist_near < max_cutoff {
                    dir_shifts[d].push(1.0); // ghost shifted by +lattice_vec[d]
                }
                if dist_far < max_cutoff {
                    dir_shifts[d].push(-1.0); // ghost shifted by -lattice_vec[d]
                }
            }

            // Generate ghost images from the Cartesian product of per-direction shifts,
            // skipping the (0,0,0) combination which is the original atom.
            for &s0 in &dir_shifts[0] {
                for &s1 in &dir_shifts[1] {
                    for &s2 in &dir_shifts[2] {
                        if s0 == 0.0 && s1 == 0.0 && s2 == 0.0 {
                            continue;
                        }
                        let pt = [
                            r[0] + s0 * lattice_vecs[0][0]
                                + s1 * lattice_vecs[1][0]
                                + s2 * lattice_vecs[2][0],
                            r[1] + s0 * lattice_vecs[0][1]
                                + s1 * lattice_vecs[1][1]
                                + s2 * lattice_vecs[2][1],
                            r[2] + s0 * lattice_vecs[0][2]
                                + s1 * lattice_vecs[1][2]
                                + s2 * lattice_vecs[2][2],
                        ];
                        tree.add(pt, i)
                            .expect("Failed to add ghost point to kd-tree");
                    }
                }
            }
        }
    }

    // Find bonds using kd-tree range queries + union-find.
    let mut uf = UnionFind::new(n_atoms);
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n_atoms];
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
                adj[i].push(j);
                adj[j].push(i);
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
            let mut mol_coords = MatrixXx3::zeros(indices.len());
            let mut unrolled = vec![false; n_atoms];
            let mut queue = std::collections::VecDeque::new();

            // Map original index to the row in mol_coords
            let mut idx_to_row = std::collections::HashMap::new();
            for (row, &idx) in indices.iter().enumerate() {
                idx_to_row.insert(idx, row);
            }

            let start_idx = indices[0];
            queue.push_back(start_idx);
            unrolled[start_idx] = true;

            // Set the first atom's coordinates to its original Cartesian coordinates
            let start_row = idx_to_row[&start_idx];
            mol_coords[(start_row, 0)] = cart_coords[(start_idx, 0)];
            mol_coords[(start_row, 1)] = cart_coords[(start_idx, 1)];
            mol_coords[(start_row, 2)] = cart_coords[(start_idx, 2)];

            while let Some(u) = queue.pop_front() {
                let u_row = idx_to_row[&u];
                let p_u_mat = MatrixXx3::from_row_slice(&[
                    mol_coords[(u_row, 0)],
                    mol_coords[(u_row, 1)],
                    mol_coords[(u_row, 2)],
                ]);

                for &v in &adj[u] {
                    if !unrolled[v] {
                        unrolled[v] = true;
                        queue.push_back(v);

                        let v_row = idx_to_row[&v];
                        let p_v_mat = MatrixXx3::from_row_slice(&[
                            cart_coords[(v, 0)],
                            cart_coords[(v, 1)],
                            cart_coords[(v, 2)],
                        ]);

                        if let Some(c) = cell {
                            let f_u = cart_to_frac(&p_u_mat, c);
                            let f_v = cart_to_frac(&p_v_mat, c);
                            let mut df = &f_v - &f_u;
                            df.apply(|x| *x -= x.round());
                            let f_v_unrolled = f_u + df;
                            let p_v_unrolled = frac_to_cart(&f_v_unrolled, c);

                            mol_coords[(v_row, 0)] = p_v_unrolled[(0, 0)];
                            mol_coords[(v_row, 1)] = p_v_unrolled[(0, 1)];
                            mol_coords[(v_row, 2)] = p_v_unrolled[(0, 2)];
                        } else {
                            mol_coords[(v_row, 0)] = p_v_mat[(0, 0)];
                            mol_coords[(v_row, 1)] = p_v_mat[(0, 1)];
                            mol_coords[(v_row, 2)] = p_v_mat[(0, 2)];
                        }
                    }
                }
            }

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

    /// A simple test struct implementing AtomicData + CartAtomicData.
    struct TestAtoms {
        atomic_nums: Vec<u8>,
        cart_coords: MatrixXx3<f64>,
    }

    impl AtomicData for TestAtoms {
        fn atomic_nums(&self) -> &[u8] {
            &self.atomic_nums
        }
    }

    impl CartAtomicData for TestAtoms {
        fn cartesian_coords(&self) -> &MatrixXx3<f64> {
            &self.cart_coords
        }
    }

    #[test]
    fn test_find_molecules_empty() {
        let atoms = TestAtoms {
            atomic_nums: vec![],
            cart_coords: MatrixXx3::zeros(0),
        };
        let uc = UnitCell::cubic(10.0);
        let molecules = find_molecules(&atoms, Some(&uc), BondSettings::Default);
        assert!(molecules.is_empty());
    }

    #[test]
    fn test_find_molecules_isolated_atoms() {
        // Two atoms far apart in a large cubic cell — should be two separate molecules.
        let uc = UnitCell::cubic(10.0);
        let atoms = TestAtoms {
            atomic_nums: vec![6, 8], // C and O
            cart_coords: MatrixXx3::from_row_slice(&[
                1.0, 1.0, 1.0, // C at (1,1,1) Å
                9.0, 9.0, 9.0, // O at (9,9,9) Å — far away
            ]),
        };
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
            cart_coords: MatrixXx3::from_row_slice(&[
                5.0, 5.0, 5.0, // C at center
                5.0, 5.0, 6.3, // O at 1.3 Å away along z
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
        // Place one H at cart 0.2 Å and another at 9.5 Å along x in a 10 Å cell.
        // Direct distance: 9.5 - 0.2 = 9.3 Å (too far)
        // Periodic distance: (10.0 - 9.5 + 0.2) = 0.7 Å (bonds!)
        let uc = UnitCell::cubic(10.0);
        let atoms = TestAtoms {
            atomic_nums: vec![1, 1], // H and H
            cart_coords: MatrixXx3::from_row_slice(&[
                0.2, 5.0, 5.0, // H near lower x boundary
                9.5, 5.0, 5.0, // H near upper x boundary
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
        // Two H atoms at opposite edges. With periodicity they bond (0.7 Å apart),
        // without periodicity they don't (9.3 Å apart, direct Cartesian distance).
        let uc = UnitCell::cubic(10.0);
        let atoms = TestAtoms {
            atomic_nums: vec![1, 1],
            cart_coords: MatrixXx3::from_row_slice(&[
                0.2, 5.0, 5.0, // H near lower x boundary
                9.5, 5.0, 5.0, // H near upper x boundary
            ]),
        };

        // With periodicity: should bond
        let mols_periodic = find_molecules(&atoms, Some(&uc), BondSettings::Default);
        assert_eq!(mols_periodic.len(), 1);

        // Without periodicity: direct distance is 9.3 Å, should NOT bond
        let mols_no_cell = find_molecules::<_, UnitCell>(&atoms, None, BondSettings::Default);
        assert_eq!(
            mols_no_cell.len(),
            2,
            "Without a cell, far atoms should not bond"
        );
    }

    #[test]
    fn test_find_molecules_custom_radii_and_tolerance() {
        // Place two C atoms at distance 1.6 Å.
        // Default C radius = 0.76, so default bond cutoff = 0.76+0.76+0.45 = 1.97 Å -> bonds.
        // With a tiny tolerance of 0.0 Å, cutoff = 0.76+0.76+0.0 = 1.52 Å -> no bond at 1.6 Å.
        let uc = UnitCell::cubic(10.0);
        let atoms = TestAtoms {
            atomic_nums: vec![6, 6],
            cart_coords: MatrixXx3::from_row_slice(&[
                5.0, 5.0, 5.0, // C at center
                5.0, 5.0, 6.6, // C at 1.6 Å apart along z
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
            cart_coords: MatrixXx3::from_row_slice(&[
                2.0, 2.0, 2.0, // C
                2.0, 2.0, 3.2, // O at 1.2 Å from C (bonds: 0.76+0.66+0.45=1.87)
                10.0, 10.0, 10.0, // N (isolated)
                18.0, 18.0, 2.0, // H
                18.0, 18.0, 2.6, // H at 0.6 Å from other H (bonds: 0.31+0.31+0.45=1.07)
            ]),
        };
        let molecules = find_molecules(&atoms, Some(&uc), BondSettings::Default);
        assert_eq!(molecules.len(), 3);
        // Sorted by descending size: two 2-atom molecules, then one 1-atom molecule
        assert_eq!(molecules[0].n_atoms(), 2);
        assert_eq!(molecules[1].n_atoms(), 2);
        assert_eq!(molecules[2].n_atoms(), 1);
    }
    #[test]
    fn test_find_molecules_unroll_benzene() {
        let uc = UnitCell::cubic(10.0);
        let mut atomic_nums = vec![6; 6];
        atomic_nums.extend(vec![1; 6]);

        let mut coords = vec![];

        let r_c = 1.40;
        let r_h = 2.49;

        for i in 0..6 {
            let angle = (i as f64) * std::f64::consts::PI / 3.0;
            coords.push([r_c * angle.cos(), r_c * angle.sin(), 0.0]);
        }
        for i in 0..6 {
            let angle = (i as f64) * std::f64::consts::PI / 3.0;
            coords.push([r_h * angle.cos(), r_h * angle.sin(), 0.0]);
        }
        println!("Initial atom coordinates:");
        for row in coords.iter() {
            println!("Atom coords: {:}, {:}, {:}", row[0], row[1], row[2]);
        }

        // Shift center to boundary and wrap to force splitting across periodic boundary
        for p in coords.iter_mut() {
            p[0] = (p[0] + 9.5).rem_euclid(10.0);
            p[1] = (p[1] + 9.5).rem_euclid(10.0);
            p[2] = (p[2] + 9.5).rem_euclid(10.0);
        }

        let cart_coords = MatrixXx3::from_fn(12, |r, c| coords[r][c]);
        let atoms = TestAtoms {
            atomic_nums,
            cart_coords,
        };

        let molecules = find_molecules(&atoms, Some(&uc), BondSettings::Default);

        assert_eq!(
            molecules.len(),
            1,
            "Benzene should be found as a single molecule"
        );

        let mol = &molecules[0];
        assert_eq!(mol.n_atoms(), 12);

        // Check that the unrolled molecule doesn't have artificially large internal distances
        let mol_coords = mol.cartesian_coords();
        println!("Unrolled molecule coordinates:");
        for row in mol_coords.row_iter() {
            println!("Atom coords: {:}, {:}, {:}", row[0], row[1], row[2]);
        }
        let mut max_dist = 0.0_f64;
        for i in 0..12 {
            for j in 0..i {
                let p_i = nalgebra::Vector3::new(
                    mol_coords[(i, 0)],
                    mol_coords[(i, 1)],
                    mol_coords[(i, 2)],
                );
                let p_j = nalgebra::Vector3::new(
                    mol_coords[(j, 0)],
                    mol_coords[(j, 1)],
                    mol_coords[(j, 2)],
                );
                let dist = (p_i - p_j).norm();
                if dist > max_dist {
                    max_dist = dist;
                }
            }
        }

        // Maximum distance in benzene is between opposite H atoms: 2.49 * 2 = 4.98 A.
        // Unrolling ensures atoms are grouped tightly. Without unrolling, max_dist > 8.0 A.
        assert!(
            max_dist < 5.0,
            "Molecule atoms are too far apart, unrolling failed! Max dist: {}",
            max_dist
        );
    }
}
