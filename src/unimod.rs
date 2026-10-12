//! Unimod accessions for common modifications (from unimod.org's
//! `unimod.obo`). A name only maps when its mass agrees too.

/// (accession number, name, monoisotopic delta)
const TABLE: &[(u32, &str, f64)] = &[
    (1, "Acetyl", 42.010_565),
    (2, "Amidated", -0.984_016),
    (3, "Biotin", 226.077_598),
    (4, "Carbamidomethyl", 57.021_464),
    (5, "Carbamyl", 43.005_814),
    (7, "Deamidated", 0.984_016),
    (21, "Phospho", 79.966_331),
    (23, "Dehydrated", -18.010_565),
    (24, "Propionamide", 71.037_114),
    (26, "Pyro-carbamidomethyl", 39.994_915),
    (27, "Glu->pyro-Glu", -18.010_565),
    (28, "Gln->pyro-Glu", -17.026_549),
    (30, "Cation:Na", 21.981_943),
    (34, "Methyl", 14.015_650),
    (35, "Oxidation", 15.994_915),
    (36, "Dimethyl", 28.031_300),
    (37, "Trimethyl", 42.046_950),
    (39, "Methylthio", 45.987_721),
    (40, "Sulfo", 79.956_815),
    (41, "Hex", 162.052_824),
    (43, "HexNAc", 203.079_373),
    (58, "Propionyl", 56.026_215),
    (64, "Succinyl", 100.016_044),
    (105, "ICAT-C", 227.126_991),
    (108, "Nethylmaleimide", 125.047_679),
    (121, "GG", 114.042_927),
    (122, "Formyl", 27.994_915),
    (184, "Label:13C(9)", 9.030_193),
    (188, "Label:13C(6)", 6.020_129),
    (199, "Dimethyl:2H(4)", 32.056_407),
    (214, "iTRAQ4plex", 144.102_063),
    (259, "Label:13C(6)15N(2)", 8.014_199),
    (267, "Label:13C(6)15N(4)", 10.008_269),
    (275, "Nitrosyl", 28.990_164),
    (312, "Cysteinyl", 119.004_099),
    (330, "Dimethyl:2H(6)13C(2)", 36.075_670),
    (345, "Trioxidation", 47.984_744),
    (354, "Nitro", 44.985_078),
    (385, "Ammonia-loss", -17.026_549),
    (425, "Dioxidation", 31.989_829),
    (481, "Label:2H(4)", 4.025_107),
    (535, "LRGG", 383.228_103),
    (730, "iTRAQ8plex", 304.205_360),
    (737, "TMT6plex", 229.162_932),
    (738, "TMT2plex", 225.155_833),
    (739, "TMT", 224.152_478),
    (747, "Malonyl", 86.000_394),
    (765, "Met-loss", -131.040_485),
    (766, "Met-loss+Acetyl", -89.029_920),
    (1289, "Butyryl", 70.041_865),
    (1290, "Dicarbamidomethyl", 114.042_927),
    (1363, "Crotonyl", 68.026_215),
    (2016, "TMTpro", 304.207_146),
    (2017, "TMTpro_zero", 295.189_592),
];

const MASS_TOLERANCE: f64 = 0.01;

/// `("UNIMOD:35", "Oxidation")` for a modification named `Oxidation` (any
/// case) of about +15.995 Da.
pub fn lookup(name: &str, mass_delta: f64) -> Option<(String, &'static str)> {
    TABLE
        .iter()
        .find(|(_, known, mass)| {
            known.eq_ignore_ascii_case(name) && (mass - mass_delta).abs() <= MASS_TOLERANCE
        })
        .map(|(id, known, _)| (format!("UNIMOD:{id}"), *known))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn names_map_only_when_the_mass_agrees() {
        assert_eq!(
            lookup("oxidation", 15.9949),
            Some(("UNIMOD:35".into(), "Oxidation"))
        );
        assert_eq!(lookup("Oxidation", 31.99), None);
        assert_eq!(lookup("Nonsense", 1.0), None);
    }
}
