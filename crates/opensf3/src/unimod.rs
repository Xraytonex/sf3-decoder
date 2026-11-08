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

#[cfg(feature = "mzid")]
const MASS_TOLERANCE: f64 = 0.01;
/// Deltas this close to a table entry are that modification reported at
/// lower precision (3-decimal engines are within 0.0005). Small enough to
/// keep TMTpro (304.2071) and iTRAQ8plex (304.2054) apart.
pub(crate) const SNAP_TOLERANCE: f64 = 0.001;

/// Unimod's own name for a modification of about `mass_delta` whose name
/// agrees with it: the same name in any case, or the Unimod name followed
/// by `_` (Scaffold writes `Phospho_STY`).
pub(crate) fn named_by_mass(name: &str, mass_delta: f64) -> Option<&'static str> {
    let lowered = name.to_ascii_lowercase();
    TABLE
        .iter()
        .filter(|(_, _, mass)| (mass - mass_delta).abs() <= SNAP_TOLERANCE)
        .find(|(_, known, _)| {
            let known = known.to_ascii_lowercase();
            lowered == known || lowered.starts_with(&format!("{known}_"))
        })
        .map(|(_, known, _)| *known)
}

/// Unimod's name for an accession like `UNIMOD:35`, when its mass is about
/// `mass_delta`.
pub(crate) fn named_by_accession(accession: &str, mass_delta: f64) -> Option<&'static str> {
    let number: u32 = accession
        .strip_prefix("UNIMOD:")
        .or_else(|| accession.strip_prefix("unimod:"))?
        .parse()
        .ok()?;
    TABLE
        .iter()
        .find(|(id, _, mass)| *id == number && (mass - mass_delta).abs() <= SNAP_TOLERANCE)
        .map(|(_, known, _)| *known)
}

#[cfg(feature = "mzid")]
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
    #[cfg(feature = "mzid")]
    fn names_map_only_when_the_mass_agrees() {
        assert_eq!(
            lookup("oxidation", 15.9949),
            Some(("UNIMOD:35".into(), "Oxidation"))
        );
        assert_eq!(lookup("Oxidation", 31.99), None);
        assert_eq!(lookup("Nonsense", 1.0), None);
    }

    #[test]
    fn scaffold_names_find_unimod_by_mass_and_name() {
        assert_eq!(
            named_by_mass("Phospho_STY", 79.966_331_481_933_6),
            Some("Phospho")
        );
        assert_eq!(named_by_mass("acetyl", 42.010_565), Some("Acetyl"));
        // Two entries share -18.0106; the name picks one.
        assert_eq!(named_by_mass("Dehydrated", -18.010_565), Some("Dehydrated"));
        assert_eq!(named_by_mass("Acetyl", 42.02), None);
        assert_eq!(named_by_mass("Mystery", 15.994_915), None);
        assert_eq!(
            named_by_accession("UNIMOD:4", 57.021_464),
            Some("Carbamidomethyl")
        );
        assert_eq!(named_by_accession("UNIMOD:1047", 42.010_565), None);
        assert_eq!(named_by_accession("UNIMOD:35", 31.99), None);
    }
}
