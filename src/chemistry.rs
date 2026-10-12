//! Monoisotopic masses.

pub const PROTON: f64 = 1.007_276_466_621;
pub const WATER: f64 = 18.010_564_683_7;
const ISOTOPE_SPACING: f64 = 1.003_355;

pub fn residue_mass(residue: char) -> Option<f64> {
    Some(match residue.to_ascii_uppercase() {
        'G' => 57.021_463_72,
        'A' => 71.037_113_79,
        'S' => 87.032_028_41,
        'P' => 97.052_763_85,
        'V' => 99.068_413_91,
        'T' => 101.047_678_5,
        'C' => 103.009_184_5,
        'L' | 'I' => 113.084_064_0,
        'N' => 114.042_927_4,
        'D' => 115.026_943_1,
        'Q' => 128.058_577_5,
        'K' => 128.094_963_0,
        'E' => 129.042_593_1,
        'M' => 131.040_484_6,
        'H' => 137.058_911_9,
        'F' => 147.068_413_9,
        'U' => 150.953_633_4,
        'R' => 156.101_111_0,
        'Y' => 163.063_328_6,
        'W' => 186.079_312_9,
        'O' => 237.147_726_9,
        _ => return None,
    })
}

/// Neutral mass of an unmodified peptide; `None` for residues like X, B, Z.
pub fn peptide_mass(sequence: &str) -> Option<f64> {
    sequence
        .chars()
        .map(residue_mass)
        .sum::<Option<f64>>()
        .map(|residues| residues + WATER)
}

/// m/z of a singly protonated mass ([M+H]+) at another charge.
pub fn mz_from_mh(mh: f64, charge: u8) -> f64 {
    let z = f64::from(charge.max(1));
    (mh - PROTON) / z + PROTON
}

/// Neutral mass error in daltons against an [M+H]+, allowing the isotope
/// errors engines search (−1..+3).
pub fn mass_error(mass: f64, mh: f64) -> f64 {
    let observed = mh - PROTON;
    (-1..=3)
        .map(|isotope| observed - mass - f64::from(isotope) * ISOTOPE_SPACING)
        .min_by(|a, b| a.abs().total_cmp(&b.abs()))
        .expect("five isotopes")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn peptide_masses_match_reference_values() {
        // PEPTIDE: 799.359964 Da (Expasy).
        assert!((peptide_mass("PEPTIDE").unwrap() - 799.359_964).abs() < 1e-5);
        assert_eq!(peptide_mass("PEPXIDE"), None);
        assert!((mz_from_mh(1001.0, 2) - 501.003_638).abs() < 1e-5);
        // One 13C isotope off is no error.
        assert!(mass_error(799.359_964, 799.359_964 + PROTON + ISOTOPE_SPACING).abs() < 1e-9);
    }
}
