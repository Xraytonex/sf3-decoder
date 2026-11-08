//! What a Scaffold project holds, decoded.

use std::borrow::Cow;
use std::collections::BTreeMap;
use std::fmt::Write as _;

use crate::chemistry::{mz_from_mh, PROTON};
use crate::sf3::peaks::PeakStore;
use crate::sfdb::peaks::SfdbPeaks;

/// Which kind of Scaffold file a project was read from.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum Format {
    /// A Scaffold 4/5 project (`.sf3`): engine candidates plus Scaffold's
    /// probabilities, protein groups and saved thresholds.
    Sf3,
    /// A Scaffold 5 or Scaffold DDA database (`.sfdb`): engine candidates
    /// and scores only; Scaffold's own verdicts are not read.
    Sfdb,
}

impl Format {
    /// `sf3` or `sfdb`.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Sf3 => "sf3",
            Self::Sfdb => "sfdb",
        }
    }
}

/// A decoded Scaffold file.
#[non_exhaustive]
pub struct Project {
    pub name: String,
    pub format: Format,
    /// As Scaffold wrote it, e.g. `Scaffold_5.3.3` or `6.5.0`.
    pub scaffold_version: Option<String>,
    pub runs: Vec<Run>,
    pub engines: Vec<Engine>,
    /// Every spectrum of every run, run by run.
    pub spectra: Vec<Spectrum>,
    pub proteins: BTreeMap<String, Protein>,
    /// Scaffold's protein groups (`.sf3` only).
    pub groups: Vec<ProteinGroup>,
    /// The thresholds saved with the project (`.sf3` only).
    pub thresholds: Option<Thresholds>,
    /// What was left out, and why.
    pub notes: Vec<String>,
    pub(crate) peaks: Option<PeakSource>,
}

/// Where a project's stored peaks are read from.
pub(crate) enum PeakSource {
    Sf3(PeakStore),
    Sfdb(SfdbPeaks),
}

/// One MS run (one raw file).
#[non_exhaustive]
pub struct Run {
    /// Unique in the project.
    pub name: String,
    pub biosample: String,
    pub category: Option<String>,
}

/// A search engine whose candidates the file holds.
#[non_exhaustive]
pub struct Engine {
    pub name: String,
    pub version: Option<String>,
    /// Score names rest on value ranges, not on Scaffold's own list.
    pub inferred_scores: bool,
}

/// One MS/MS spectrum Scaffold kept.
#[non_exhaustive]
pub struct Spectrum {
    /// Unique in the project.
    pub id: i32,
    pub run: usize,
    pub title: Option<String>,
    /// The scan the title names, when it names exactly one and no other
    /// spectrum of the run names the same.
    pub scan: Option<u32>,
    /// Precursor [M+H]+. `.sf3` stores this; for `.sfdb` it follows from
    /// the stored m/z and charge.
    pub mh: Option<f64>,
    pub charge: Option<u8>,
    pub rt_seconds: Option<f64>,
    /// Precursor m/z: stored in `.sfdb`, from [M+H]+ and charge in `.sf3`.
    pub(crate) mz: Option<f64>,
    pub(crate) peak_list: i32,
    /// Every engine's candidates, ranked per engine.
    pub matches: Vec<Match>,
}

impl Spectrum {
    /// A spectrum from an `.sf3` record, which stores [M+H]+.
    pub(crate) fn from_mh(
        id: i32,
        run: usize,
        title: Option<String>,
        mh: Option<f64>,
        charge: Option<u8>,
    ) -> Self {
        Self {
            id,
            run,
            title,
            scan: None,
            mh,
            charge,
            rt_seconds: None,
            mz: mh.zip(charge).map(|(mh, charge)| mz_from_mh(mh, charge)),
            peak_list: 0,
            matches: Vec::new(),
        }
    }

    /// A spectrum from an `.sfdb` row, which stores m/z.
    pub(crate) fn from_mz(
        id: i32,
        run: usize,
        title: Option<String>,
        mz: Option<f64>,
        charge: Option<u8>,
    ) -> Self {
        let mh = mz
            .zip(charge)
            .map(|(mz, charge)| (mz - PROTON) * f64::from(charge) + PROTON);
        Self {
            mz,
            ..Self::from_mh(id, run, title, mh, charge)
        }
    }

    #[cfg(test)]
    pub(crate) fn for_tests(run: usize, scan: Option<u32>) -> Self {
        Self {
            scan,
            ..Self::from_mh(0, run, None, None, None)
        }
    }

    pub fn precursor_mz(&self) -> Option<f64> {
        self.mz
    }
}

/// One engine's candidate peptide for a spectrum.
#[non_exhaustive]
pub struct Match {
    pub engine: usize,
    /// By the engine's own score. In `.sf3` ties share a rank; in `.sfdb`
    /// ties are broken by the file's PSM id.
    pub rank: u32,
    pub peptide: Peptide,
    /// Neutral monoisotopic mass, when every residue has one.
    pub mass: Option<f64>,
    /// The charge the engine assigned (the spectrum's, in `.sf3`).
    pub charge: Option<u8>,
    pub scores: Vec<Score>,
    pub proteins: Vec<ProteinRef>,
    /// Scaffold's peptide probability (`.sf3` only).
    pub probability: Option<f32>,
    /// The candidate Scaffold picked for this spectrum and engine (`.sf3`
    /// only; always false for `.sfdb`).
    pub chosen: bool,
}

/// One engine score of a candidate.
#[non_exhaustive]
pub struct Score {
    pub name: Cow<'static, str>,
    /// PSI-MS accession, when the score has one.
    pub accession: Option<Cow<'static, str>>,
    pub value: f64,
}

impl Score {
    pub(crate) fn new(
        name: impl Into<Cow<'static, str>>,
        accession: Option<Cow<'static, str>>,
        value: f64,
    ) -> Self {
        Self {
            name: name.into(),
            accession,
            value,
        }
    }
}

/// Where a candidate's peptide sits on one protein.
#[non_exhaustive]
pub struct ProteinRef {
    pub accession: String,
    pub pre: Option<char>,
    pub post: Option<char>,
    /// 1-based first and last residue on the protein, when known (`.sfdb`).
    pub start: Option<u32>,
    pub end: Option<u32>,
}

impl ProteinRef {
    pub(crate) fn new(accession: String, pre: Option<char>, post: Option<char>) -> Self {
        Self {
            accession,
            pre,
            post,
            start: None,
            end: None,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Site {
    NTerm,
    /// 1-based residue position.
    Residue(usize),
    CTerm,
}

#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct Modification {
    pub site: Site,
    pub mass_delta: f64,
    pub name: Option<String>,
}

#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub struct Peptide {
    /// Upper case.
    pub sequence: String,
    /// Sorted by site.
    pub modifications: Vec<Modification>,
}

impl Peptide {
    pub fn new(sequence: String, mut modifications: Vec<Modification>) -> Self {
        modifications.sort_by_key(|modification| modification.site);
        Self {
            sequence,
            modifications,
        }
    }

    /// ProForma-style text, e.g. `[Acetyl]-PEPM[Oxidation]IDE`; unnamed
    /// modifications show their mass.
    pub fn key(&self) -> String {
        let tag = |modification: &Modification| match &modification.name {
            Some(name) => format!("[{name}]"),
            None => format!("[{:+.4}]", modification.mass_delta),
        };
        let mut key = String::new();
        for modification in &self.modifications {
            if modification.site == Site::NTerm {
                let _ = write!(key, "{}-", tag(modification));
            }
        }
        for (index, residue) in self.sequence.chars().enumerate() {
            key.push(residue);
            for modification in &self.modifications {
                if modification.site == Site::Residue(index + 1) {
                    key.push_str(&tag(modification));
                }
            }
        }
        for modification in &self.modifications {
            if modification.site == Site::CTerm {
                let _ = write!(key, "-{}", tag(modification));
            }
        }
        key
    }
}

/// A protein a candidate or group names.
#[non_exhaustive]
pub struct Protein {
    pub description: Option<String>,
    pub sequence: Option<String>,
    pub decoy: bool,
}

/// One of Scaffold's protein groups.
#[non_exhaustive]
pub struct ProteinGroup {
    /// Scaffold's preferred accession first.
    pub accessions: Vec<String>,
    pub decoy: bool,
    /// Scaffold's protein probability in each run that has the group.
    pub probabilities: Vec<(usize, f64)>,
}

impl ProteinGroup {
    pub fn best_probability(&self) -> Option<f64> {
        self.probabilities
            .iter()
            .map(|&(_, probability)| probability)
            .max_by(f64::total_cmp)
    }
}

/// The thresholds saved with the project.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub struct Thresholds {
    pub peptide_probability: f64,
    pub protein_probability: f64,
    pub min_peptides: u32,
}

impl Project {
    /// A spectrum's peaks as (m/z, intensity), sorted by m/z. `Ok(None)`
    /// when the project stores none for it.
    pub fn peaks(&self, spectrum: &Spectrum) -> Result<Option<Vec<(f64, f32)>>, &'static str> {
        if spectrum.peak_list <= 0 {
            return Ok(None);
        }
        match &self.peaks {
            Some(PeakSource::Sf3(store)) => store.read(spectrum.peak_list).map(Some),
            Some(PeakSource::Sfdb(store)) => store.read(spectrum.peak_list),
            None => Ok(None),
        }
    }

    pub fn run_spectra(&self, run: usize) -> impl Iterator<Item = &Spectrum> {
        self.spectra
            .iter()
            .filter(move |spectrum| spectrum.run == run)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn peptide_keys_show_every_modification_in_place() {
        let peptide = Peptide::new(
            "PEPMK".into(),
            vec![
                Modification {
                    site: Site::Residue(4),
                    mass_delta: 15.994_915,
                    name: Some("Oxidation".into()),
                },
                Modification {
                    site: Site::NTerm,
                    mass_delta: 42.010_565,
                    name: Some("Acetyl".into()),
                },
                Modification {
                    site: Site::CTerm,
                    mass_delta: -0.984_016,
                    name: None,
                },
            ],
        );
        assert_eq!(peptide.key(), "[Acetyl]-PEPM[Oxidation]K-[-0.9840]");
    }
}
