use std::collections::HashMap;

// Defaults
pub const DEFAULT_DIST_TOL: f64 = 1e-5;
pub const DEFAULT_ANG_TOL: f64 = 1e-3;
pub static DEFAULT_COVALENT_RADIUS: f64 = 2.00;
pub static DEFAULT_SPECIES_NAME: &str = "X";
pub static DEFAULT_ATOMIC_MASS: f64 = 0.0;

// Conversion factors
pub const AVOGADRO_NUMBER: f64 = 6.02214076e23; // mol^-1
pub const AMU_TO_KG: f64 = 1.66053906660e-27; // kg/amu
pub const ANGSTROM3_TO_CM3: f64 = 1e-24; // cm^3/Å^3

/// Data from: https://chem.libretexts.org/Ancillary_Materials/Reference/Reference_Tables/Atomic_and_Molecular_Properties/A3%3A_Covalent_Radii
/// Units: Angstroms
pub static COVALENT_RADII: [f64; 10] = [
    DEFAULT_COVALENT_RADIUS, // X
    0.31,                    // H
    0.28,                    // He
    1.28,                    // Li
    0.96,                    // Be
    0.84,                    // B
    0.76,                    // C
    0.71,                    // N
    0.66,                    // O
    0.57,                    // F
];

/// Get the covalent radius for a given atomic number.
pub fn get_covalent_radius(atomic_num: u8) -> f64 {
    if atomic_num as usize >= COVALENT_RADII.len() {
        DEFAULT_COVALENT_RADIUS
    } else {
        COVALENT_RADII[atomic_num as usize]
    }
}

pub fn get_default_covalent_radii_map() -> HashMap<u8, f64> {
    let mut radii_map = HashMap::new();
    for (atomic_num, &radius) in COVALENT_RADII.iter().enumerate() {
        radii_map.insert(atomic_num as u8, radius);
    }
    radii_map
}

/// List of species names in order of atomic numbers
pub static SPECIES_NAMES: [&str; 10] = [
    DEFAULT_SPECIES_NAME, // Unknown
    "H",                  // Hydrogen
    "He",                 // Helium
    "Li",                 // Lithium
    "Be",                 // Beryllium
    "B",                  // Boron
    "C",                  // Carbon
    "N",                  // Nitrogen
    "O",                  // Oxygen
    "F",                  // Fluorine
];

/// Get the species name for a given atomic number.
pub fn get_species_name(atomic_num: u8) -> &'static str {
    if atomic_num as usize >= SPECIES_NAMES.len() {
        DEFAULT_SPECIES_NAME
    } else {
        SPECIES_NAMES[atomic_num as usize]
    }
}

pub static ATOMIC_MASSES: [f64; 10] = [
    DEFAULT_ATOMIC_MASS, // Unknown
    1.008,               // Hydrogen
    4.0026,              // Helium
    6.94,                // Lithium
    9.0122,              // Beryllium
    10.81,               // Boron
    12.011,              // Carbon
    14.007,              // Nitrogen
    15.999,              // Oxygen
    18.998,              // Fluorine
];

/// Get the atomic mass for a given atomic number.
pub fn get_atomic_mass(atomic_num: u8) -> f64 {
    if atomic_num as usize >= ATOMIC_MASSES.len() {
        DEFAULT_ATOMIC_MASS
    } else {
        ATOMIC_MASSES[atomic_num as usize]
    }
}
