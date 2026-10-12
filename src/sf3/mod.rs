//! Reading a Scaffold `.sf3` project into a [`Project`].
//!
//! Every check here guards the decoding itself: a file whose structure
//! isn't the one known is refused rather than half-read, and a candidate
//! that can't be decoded faithfully is left out and counted in the notes.

mod accessions;
mod blockbased;
mod candidates;
mod java;
pub(crate) mod peaks;
mod project;
mod records;

use std::collections::{BTreeMap, HashMap, HashSet};
use std::path::Path;

use crate::chemistry::mass_error;
use crate::decoy::is_decoy_accession;
use crate::model::{Engine, Match, Peptide, Project, Protein, ProteinGroup, Run, Score, Spectrum};
use accessions::{Completions, INLINE_ACCESSION_MAX, PROTEIN_BLOCK, proteins_of};
use candidates::{Candidate, EngineScores, PEPTIDE_BLOCK, Peptides, analysis_of, engine_scores};
use peaks::PeakStore;
use records::{Container, Records, read_f32, read_f64, read_i32};

/// Daltons; far beyond any search tolerance, near any decoding error.
pub const MAX_MASS_ERROR: f64 = 1.0;
/// Share of candidates allowed to miss their spectrum's mass before the
/// file is refused as misread.
pub const MAX_MASS_FAILURES: f64 = 0.05;

const SPECTRUM_BLOCK: usize = 48;
const ANALYSIS_BLOCK: usize = 180;
const INT_LIST_BLOCK: usize = 256;

#[derive(Debug, thiserror::Error)]
pub enum Sf3Error {
    #[error("not a Scaffold .sf3 file: {0}")]
    NotSf3(String),
    #[error("cannot read the .sf3 container: {0}")]
    Container(String),
    #[error("the .sf3 file has no {0}")]
    Missing(String),
    #[error("{file}: {error}")]
    Java {
        file: String,
        error: java::JavaError,
    },
    #[error(
        "{file}: {problem}. This file's layout differs from the Scaffold versions the decoder \
         knows (4.0–5.3.4), so it was refused rather than guessed at"
    )]
    Layout { file: String, problem: String },
}

pub(crate) fn layout(file: &str, problem: &str) -> Sf3Error {
    Sf3Error::Layout {
        file: file.to_string(),
        problem: problem.to_string(),
    }
}

/// A candidate whose peptide is settled.
struct Resolved {
    peptide: Peptide,
    key: String,
    mass: Option<f64>,
    protein_record: i32,
    primary: f32,
    scores: Vec<Score>,
    probability: Option<f32>,
    chosen: bool,
    flanks: [u8; 2],
}

pub fn decode(path: &Path) -> Result<Project, Sf3Error> {
    let container = Container::open(path)?;
    // The largest stream by far; read once and dropped.
    let (runs, name, thresholds, group_entries) = {
        let experiment = container.stream("Experiment")?;
        (
            project::runs(&container, &experiment)?,
            project::experiment_name(&experiment),
            project::thresholds(&experiment),
            project::protein_groups(&experiment),
        )
    };
    let programs = project::programs(&container)?;
    let types = project::modification_types(&container)?;
    let protein_map = project::protein_map(&container)?;
    let stored = project::proteins(&container)?;
    let names = Completions::new(&stored, &protein_map);
    let spectrum_records = container.records("Spectrum", SPECTRUM_BLOCK)?;
    let titles = container.strings("SpectrumName")?;
    let analyses = container.records("SpectrumAnalysis", ANALYSIS_BLOCK)?;
    let protein_records = container.records("ProteinIdentification", PROTEIN_BLOCK)?;
    let records = container.records("PeptideIdentification", PEPTIDE_BLOCK)?;
    let tails = container.strings("PeptideSequence")?;
    let type_lists = container.records("ModificationTypeIds", INT_LIST_BLOCK)?;
    let position_lists = container.records("ModificationLocations", INT_LIST_BLOCK)?;
    let peptides = Peptides {
        tails: &tails,
        type_lists: &type_lists,
        position_lists: &position_lists,
        types: &types,
    };
    let options = peptides.residue_options(&records);

    let spectrum_run = spectrum_runs(&runs, &analyses, &spectrum_records)?;
    let mut spectra = spectra(&spectrum_run, &spectrum_records, &titles);
    let position: HashMap<i32, usize> = spectra
        .iter()
        .enumerate()
        .map(|(index, spectrum)| (spectrum.id, index))
        .collect();
    let mut notes = Notes::default();

    // (spectrum, engine) → candidates, in record order.
    let mut engines: Vec<Engine> = Vec::new();
    let mut layouts: HashMap<u8, Option<(usize, EngineScores)>> = HashMap::new();
    let mut groups: BTreeMap<(i32, usize), Vec<Candidate>> = BTreeMap::new();
    for id in records.ids() {
        let block = records
            .get(id)
            .ok_or_else(|| layout("PeptideIdentification", &format!("record {id} is missing")))?;
        let Some(analysis) = analyses.get(analysis_of(block)) else {
            return Err(layout(
                "PeptideIdentification",
                &format!("record {id} names no spectrum analysis"),
            ));
        };
        let spectrum = read_i32(analysis, 4).unwrap_or(0);
        if !spectrum_run.contains_key(&spectrum) {
            notes.without_run.insert(spectrum);
            continue;
        }
        let program = analysis[8];
        let engine = layouts.entry(program).or_insert_with(|| {
            let (name, version) = programs
                .get(&program)
                .cloned()
                .unwrap_or_else(|| (format!("program {program}"), None));
            let scores = engine_scores(&name)?;
            // Scaffold keeps one program per search; one engine is one
            // name and version.
            let known = engines
                .iter()
                .position(|engine| engine.name == name && engine.version == version);
            let index = known.unwrap_or_else(|| {
                engines.push(Engine {
                    name,
                    version,
                    inferred_scores: scores.inferred,
                });
                engines.len() - 1
            });
            Some((index, scores))
        });
        let Some((engine, scores)) = engine else {
            let name = programs
                .get(&program)
                .map_or_else(|| format!("program {program}"), |(name, _)| name.clone());
            *notes.unknown_engines.entry(name).or_default() += 1;
            continue;
        };
        let candidate = peptides.candidate(block, scores).map_err(|problem| {
            layout(
                "PeptideIdentification",
                &format!("record {id}: {}", problem.0),
            )
        })?;
        groups
            .entry((spectrum, *engine))
            .or_default()
            .push(candidate);
    }

    let mut checked = 0_usize;
    for ((spectrum_id, engine), candidates) in groups {
        let spectrum = &mut spectra[position[&spectrum_id]];
        let (Some(_), Some(mh)) = (spectrum.charge, spectrum.mh) else {
            notes.without_charge += candidates.len();
            continue;
        };
        let mut resolved = Vec::with_capacity(candidates.len());
        for candidate in candidates {
            match candidate.decoded.resolve(mh, &options, MAX_MASS_ERROR) {
                Some((peptide, mass)) => resolved.push(Resolved {
                    key: peptide.key(),
                    peptide,
                    mass,
                    protein_record: candidate.protein_record,
                    primary: candidate.primary,
                    scores: candidate.scores,
                    probability: candidate.probability,
                    chosen: candidate.chosen,
                    flanks: candidate.flanks,
                }),
                None => notes.unresolved += 1,
            }
        }
        for (rank, candidate) in rank(resolved, &mut notes) {
            if let Some(mass) = candidate.mass {
                checked += 1;
                let error = mass_error(mass, mh);
                if error.is_nan() || error.abs() > MAX_MASS_ERROR {
                    notes.mass_failures.push((candidate.key, error));
                    continue;
                }
            }
            let Some(proteins) = proteins_of(
                &candidate.peptide.sequence,
                candidate.protein_record,
                candidate.flanks,
                &protein_map,
                &protein_records,
                &names,
            ) else {
                notes.cut_accessions += 1;
                continue;
            };
            spectrum.matches.push(Match {
                engine,
                rank,
                peptide: candidate.peptide,
                mass: candidate.mass,
                scores: candidate.scores,
                proteins,
                probability: candidate.probability,
                chosen: candidate.chosen,
            });
        }
    }
    if checked > 0 && notes.mass_failures.len() as f64 > MAX_MASS_FAILURES * checked as f64 {
        return Err(layout(
            "PeptideIdentification",
            &format!(
                "{} of {checked} candidates miss their spectrum's mass by more than \
                 {MAX_MASS_ERROR} Da",
                notes.mass_failures.len()
            ),
        ));
    }

    let groups = protein_groups(
        group_entries,
        &container.protein_group_rows()?,
        &runs,
        &mut notes,
    );
    let proteins = proteins(&spectra, &groups, &stored, &mut notes);
    let peaks = match PeakStore::load(&container) {
        Ok(peaks) => peaks,
        Err(error) => {
            notes.peaks_unreadable = Some(error.to_string());
            None
        }
    };
    let name = name.unwrap_or_else(|| {
        path.file_stem().map_or_else(
            || "project".into(),
            |stem| stem.to_string_lossy().into_owned(),
        )
    });
    Ok(Project {
        name,
        scaffold_version: container.version(),
        runs: runs
            .into_iter()
            .map(|run| Run {
                name: run.name,
                biosample: run.biosample,
                category: run.category,
            })
            .collect(),
        engines,
        spectra,
        proteins,
        groups,
        thresholds,
        notes: notes.lines(),
        peaks,
    })
}

/// Every spectrum of every run, run by run, then by id.
fn spectra(
    spectrum_run: &HashMap<i32, usize>,
    records: &Records,
    titles: &Records,
) -> Vec<Spectrum> {
    let mut ids: Vec<(usize, i32)> = spectrum_run.iter().map(|(&id, &run)| (run, id)).collect();
    ids.sort_unstable();
    ids.into_iter()
        .filter_map(|(run, id)| {
            let block = records.get(id)?; // checked in spectrum_runs
            Some(Spectrum {
                id,
                run,
                title: titles
                    .string(read_i32(block, 8).unwrap_or(0))
                    .filter(|title| !title.trim().is_empty()),
                mh: read_f64(block, 12).filter(|mh| mh.is_finite() && *mh > 0.0),
                charge: read_i32(block, 20)
                    .and_then(|charge| u8::try_from(charge).ok())
                    .filter(|&charge| charge > 0),
                rt_seconds: read_f32(block, 32)
                    .filter(|rt| rt.is_finite() && *rt > 0.0)
                    .map(f64::from),
                peak_list: read_i32(block, 4).unwrap_or(0),
                matches: Vec::new(),
            })
        })
        .collect()
}

/// Rank one engine's candidates for a spectrum by its own score (ties share
/// a rank). Scaffold keeps stale copies of some candidates; a copy counts
/// once, and Scaffold's pick among copies is kept.
fn rank(candidates: Vec<Resolved>, notes: &mut Notes) -> Vec<(u32, Resolved)> {
    let mut seen: HashMap<String, usize> = HashMap::new();
    let mut unique: Vec<Resolved> = Vec::with_capacity(candidates.len());
    for candidate in candidates {
        if !candidate.primary.is_finite() {
            notes.unscored += 1;
            continue;
        }
        match seen.get(&candidate.key) {
            Some(&kept) => {
                notes.duplicates += 1;
                let kept = &mut unique[kept];
                if candidate.chosen && !kept.chosen {
                    kept.chosen = true;
                    kept.probability = candidate.probability;
                }
            }
            None => {
                seen.insert(candidate.key.clone(), unique.len());
                unique.push(candidate);
            }
        }
    }
    let scores: Vec<f32> = unique.iter().map(|candidate| candidate.primary).collect();
    let mut ranked: Vec<(u32, Resolved)> = unique
        .into_iter()
        .map(|candidate| {
            let better = scores
                .iter()
                .filter(|&&score| score > candidate.primary)
                .count();
            (u32::try_from(better + 1).unwrap_or(u32::MAX), candidate)
        })
        .collect();
    ranked.sort_by_key(|(rank, _)| *rank);
    ranked
}

/// Spectrum id → index of its run, from each run's listed analyses.
fn spectrum_runs(
    runs: &[project::Run],
    analyses: &Records,
    spectra: &Records,
) -> Result<HashMap<i32, usize>, Sf3Error> {
    let mut owner = HashMap::new();
    for (index, run) in runs.iter().enumerate() {
        for &analysis in &run.analyses {
            let block = analyses.get(analysis).ok_or_else(|| {
                layout(
                    "SpectrumCollation",
                    &format!("run {} lists missing analysis {analysis}", run.name),
                )
            })?;
            let spectrum = read_i32(block, 4)
                .filter(|&spectrum| spectra.get(spectrum).is_some())
                .ok_or_else(|| {
                    layout(
                        "SpectrumAnalysis",
                        &format!("analysis {analysis} names no spectrum"),
                    )
                })?;
            if owner
                .insert(spectrum, index)
                .is_some_and(|other| other != index)
            {
                return Err(layout(
                    "SpectrumCollation",
                    &format!("spectrum {spectrum} is in two runs"),
                ));
            }
        }
    }
    Ok(owner)
}

/// Scaffold's protein groups with their probability in each run
/// (`ProteinGroup` rows: group, MS sample, probability).
fn protein_groups(
    entries: Vec<project::GroupEntry>,
    rows: &[(i64, i64, f64)],
    runs: &[project::Run],
    notes: &mut Notes,
) -> Vec<ProteinGroup> {
    let run_of: HashMap<i64, usize> = runs
        .iter()
        .enumerate()
        .filter_map(|(index, run)| Some((run.sample_id?, index)))
        .collect();
    let mut probabilities: HashMap<i64, Vec<(usize, f64)>> = HashMap::new();
    for &(group, sample, probability) in rows {
        match run_of.get(&sample) {
            Some(&run) if probability.is_finite() => {
                probabilities
                    .entry(group)
                    .or_default()
                    .push((run, probability));
            }
            _ => notes.unplaced_group_rows += 1,
        }
    }
    let mut groups: Vec<(i64, ProteinGroup)> = entries
        .into_iter()
        .map(|entry| {
            let mut found = probabilities.remove(&entry.id).unwrap_or_default();
            found.sort_by_key(|&(run, _)| run);
            (
                entry.id,
                ProteinGroup {
                    accessions: entry.accessions,
                    decoy: entry.decoy,
                    probabilities: found,
                },
            )
        })
        .collect();
    notes.unplaced_group_rows += probabilities.values().map(Vec::len).sum::<usize>();
    groups.sort_by_key(|(id, _)| *id);
    groups.into_iter().map(|(_, group)| group).collect()
}

/// Every protein a match or group names. A stored sequence is kept only when
/// it holds all of the protein's matched peptides (I and L alike).
fn proteins(
    spectra: &[Spectrum],
    groups: &[ProteinGroup],
    stored: &HashMap<String, project::ProteinEntry>,
    notes: &mut Notes,
) -> BTreeMap<String, Protein> {
    let mut peptides: BTreeMap<&str, HashSet<&str>> = BTreeMap::new();
    for found in spectra.iter().flat_map(|spectrum| &spectrum.matches) {
        for protein in &found.proteins {
            peptides
                .entry(&protein.accession)
                .or_default()
                .insert(&found.peptide.sequence);
        }
    }
    let mut decoy_groups: HashMap<&str, bool> = HashMap::new();
    for group in groups {
        for accession in &group.accessions {
            peptides.entry(accession).or_default();
            decoy_groups.insert(accession, group.decoy);
        }
    }
    let leucine = |text: &str| text.replace('I', "L");
    peptides
        .into_iter()
        .map(|(accession, sequences)| {
            let entry = stored.get(accession);
            let sequence = entry
                .and_then(|entry| entry.sequence.clone())
                .filter(|sequence| {
                    let protein = leucine(sequence);
                    let holds = sequences
                        .iter()
                        .all(|peptide| protein.contains(&leucine(peptide)));
                    if !holds {
                        notes.dropped_sequences += 1;
                    }
                    holds
                });
            let protein = Protein {
                description: entry.and_then(|entry| entry.description.clone()),
                sequence,
                decoy: decoy_groups
                    .get(accession)
                    .copied()
                    .unwrap_or_else(|| is_decoy_accession(accession)),
            };
            (accession.to_string(), protein)
        })
        .collect()
}

#[derive(Default)]
struct Notes {
    without_run: HashSet<i32>,
    without_charge: usize,
    unknown_engines: BTreeMap<String, usize>,
    unscored: usize,
    unresolved: usize,
    duplicates: usize,
    mass_failures: Vec<(String, f64)>,
    dropped_sequences: usize,
    cut_accessions: usize,
    unplaced_group_rows: usize,
    peaks_unreadable: Option<String>,
}

impl Notes {
    fn lines(&self) -> Vec<String> {
        let mut lines = Vec::new();
        let mut note = |count: usize, text: &str| {
            if count > 0 {
                lines.push(format!("{count} {text}"));
            }
        };
        note(
            self.without_run.len(),
            "spectra belong to no MS run and their candidates were left out.",
        );
        note(
            self.without_charge,
            "candidates on spectra without a precursor charge or mass were left out.",
        );
        note(
            self.unscored,
            "candidates without their engine's score were left out.",
        );
        note(
            self.unresolved,
            "candidates with residues marked modified (lower case) but no single \
             modification that fits the precursor mass were left out.",
        );
        note(
            self.duplicates,
            "repeated copies of candidates were counted once.",
        );
        note(
            self.cut_accessions,
            &format!(
                "candidates were left out: Scaffold cut their protein's accession at \
                 {INLINE_ACCESSION_MAX} characters and nothing tells whether it names the \
                 target or its decoy."
            ),
        );
        note(
            self.dropped_sequences,
            "stored protein sequences don't contain their peptides and were left out.",
        );
        note(
            self.unplaced_group_rows,
            "protein group probabilities name no known run or group and were left out.",
        );
        for (engine, count) in &self.unknown_engines {
            lines.push(format!(
                "{engine}: {count} candidates left out; where this engine's scores are \
                 stored isn't known yet."
            ));
        }
        if !self.mass_failures.is_empty() {
            let examples: Vec<String> = self
                .mass_failures
                .iter()
                .take(3)
                .map(|(peptide, error)| format!("{peptide} ({error:+.3} Da)"))
                .collect();
            lines.push(format!(
                "{} candidates whose mass misses their spectrum's by more than \
                 {MAX_MASS_ERROR} Da were left out, e.g. {}.",
                self.mass_failures.len(),
                examples.join(", ")
            ));
        }
        if let Some(error) = &self.peaks_unreadable {
            lines.push(format!("The peaks could not be read: {error}."));
        }
        lines
    }
}
