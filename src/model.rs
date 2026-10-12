//! What a Scaffold project holds, decoded.

use std::collections::BTreeMap;
use std::fmt::Write as _;

use crate::chemistry::mz_from_mh;
use crate::sf3::peaks::PeakStore;

pub struct Project {
    pub name: String,
    /// As Scaffold wrote it, e.g. `Scaffold_5.3.3`.
    pub scaffold_version: Option<String>,
    pub runs: Vec<Run>,
    pub engines: Vec<Engine>,
    /// Every spectrum of every run, run by run.
    pub spectra: Vec<Spectrum>,
    pub proteins: BTreeMap<String, Protein>,
    pub groups: Vec<ProteinGroup>,
    pub thresholds: Option<Thresholds>,
    /// What was left out, and why.
    pub notes: Vec<String>,
    pub(crate) peaks: Option<PeakStore>,
}

pub struct Run {
    /// Unique in the project.
    pub name: String,
    pub biosample: String,
    pub category: Option<String>,
}

pub struct Engine {
    pub name: String,
    pub version: Option<String>,
    /// Score names rest on value ranges, not on Scaffold's own list.
    pub inferred_scores: bool,
}

pub struct Spectrum {
    pub id: i32,
    pub run: usize,
    pub title: Option<String>,
    /// Scaffold stores the precursor as [M+H]+, not m/z.
    pub mh: Option<f64>,
    pub charge: Option<u8>,
    pub rt_seconds: Option<f64>,
    pub(crate) peak_list: i32,
    /// Every engine's candidates, ranked per engine.
    pub matches: Vec<Match>,
}

impl Spectrum {
    pub fn precursor_mz(&self) -> Option<f64> {
        Some(mz_from_mh(self.mh?, self.charge?))
    }
}

pub struct Match {
    pub engine: usize,
    /// By the engine's own score; ties share a rank.
    pub rank: u32,
    pub peptide: Peptide,
    /// Neutral monoisotopic mass, when every residue has one.
    pub mass: Option<f64>,
    pub scores: Vec<Score>,
    pub proteins: Vec<ProteinRef>,
    /// Scaffold's peptide probability.
    pub probability: Option<f32>,
    /// The candidate Scaffold picked for this spectrum and engine.
    pub chosen: bool,
}

pub struct Score {
    pub name: &'static str,
    /// PSI-MS accession, when the score has one.
    pub accession: Option<&'static str>,
    pub value: f64,
}

pub struct ProteinRef {
    pub accession: String,
    pub pre: Option<char>,
    pub post: Option<char>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Site {
    NTerm,
    /// 1-based residue position.
    Residue(usize),
    CTerm,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Modification {
    pub site: Site,
    pub mass_delta: f64,
    pub name: Option<String>,
}

#[derive(Debug, Clone, PartialEq)]
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

pub struct Protein {
    pub description: Option<String>,
    pub sequence: Option<String>,
    pub decoy: bool,
}

/// One of Scaffold's protein groups.
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
pub struct Thresholds {
    pub peptide_probability: f64,
    pub protein_probability: f64,
    pub min_peptides: u32,
}

impl Project {
    /// A spectrum's peaks as (m/z, intensity), sorted by m/z. `Ok(None)`
    /// when the project stores none for it.
    pub fn peaks(&self, spectrum: &Spectrum) -> Result<Option<Vec<(f64, f32)>>, &'static str> {
        match &self.peaks {
            Some(store) if spectrum.peak_list > 0 => store.read(spectrum.peak_list).map(Some),
            _ => Ok(None),
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
