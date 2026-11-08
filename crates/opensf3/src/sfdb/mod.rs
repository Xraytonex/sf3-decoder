//! Reading a Scaffold `.sfdb` database (Scaffold 5.1 and Scaffold DDA 6.5)
//! into a [`Project`], read-only.
//!
//! An `.sfdb` is a plain SQLite database in an mzIdentML-like layout; see
//! `docs/format/02-sfdb.md` for the tables and what was checked. The rules
//! are the `.sf3` decoder's:
//! - Only the engines' own scores are read (`ScoreIndex`,
//!   `PEPTIDE_SCORE_ARRAY`). Scaffold's probabilities and post-processor
//!   verdicts (Percolator, PeptideProphet) are left out and named in the
//!   notes, and so are Scaffold's protein groups and thresholds: what they
//!   mean in this layout isn't settled. Ranks are by the engine's main
//!   score over every PSM of a spectrum (equal scores by PSM id), before
//!   anything is left out; a spectrum whose best match can't be read is
//!   left out whole, so a lower match never stands in as rank 1.
//! - A peptide's modifications are the ones its row lists; the N-terminus
//!   each PSM carries is decided by its precursor (`peptides.rs`). Every
//!   PSM's mass must fit its precursor, or it is left out, and the file is
//!   refused when more than 5% miss.
//! - Decoys are Scaffold's `IS_DECOY` flags, else the accession convention.
//! - Runs come from the spectrum titles where they name the raw file
//!   (Scaffold 5), else from the MS sample, which is also the biosample.

pub(crate) mod peaks;
mod peptides;

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::path::Path;

use rusqlite::{Connection, OpenFlags};

use crate::chemistry::PROTON;
use crate::decoy::is_decoy_accession;
use crate::model::{
    Engine, Format, Match, PeakSource, Project, Protein, ProteinRef, Run, Score, Spectrum,
};
use crate::sf3::candidates::expectation;
use crate::titles::{drop_clashing_scans, scan_from_title};
use peaks::{Layout, SfdbPeaks};
use peptides::{peptide_neutral, Readings, Tolerance};

/// More PSMs than this missing their precursor refuses the file.
pub const MAX_MASS_MISSES: f64 = 0.05;
/// Programs whose output is a verdict on the engines, not engine output:
/// their scores are left out, and they are never an engine.
const VERDICTS: [&str; 3] = ["scaffold", "percolator", "peptideprophet"];
/// Score array slots at most.
const MAX_SCORES: i64 = 64;

/// Why an `.sfdb` couldn't be read.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum SfdbError {
    #[error("not a Scaffold .sfdb file: {0}")]
    NotSfdb(String),
    #[error("cannot read the .sfdb file: {0}")]
    Sqlite(#[from] rusqlite::Error),
    #[error(
        "{0}. The file's layout differs from the Scaffold .sfdb versions the decoder knows \
         (Scaffold 5.1, Scaffold DDA 6.5), so it was refused rather than guessed at"
    )]
    Layout(String),
}

/// Tables every `.sfdb` the decoder reads has.
pub(crate) const REQUIRED: [&str; 9] = [
    "PeptideSpectrumMatch",
    "Peptide",
    "PeptideEvidence",
    "Protein",
    "Spectrum",
    "MSSample",
    "ScoreIndex",
    "Modification",
    "AnalysisSoftware",
];

/// Read an `.sfdb` file.
pub fn decode(path: &Path) -> Result<Project, SfdbError> {
    let conn = Connection::open_with_flags(
        path,
        OpenFlags::SQLITE_OPEN_READ_ONLY | OpenFlags::SQLITE_OPEN_NO_MUTEX,
    )
    .map_err(|error| SfdbError::NotSfdb(error.to_string()))?;
    let tables = table_names(&conn).map_err(|error| SfdbError::NotSfdb(error.to_string()))?;
    if let Some(missing) = REQUIRED.iter().find(|table| !tables.contains(**table)) {
        return Err(SfdbError::NotSfdb(format!("no {missing} table")));
    }
    let version: Option<String> = conn
        .query_row("SELECT version_name FROM Version LIMIT 1", [], |row| {
            row.get(0)
        })
        .ok();
    let mut notes = Notes::default();
    let (engines, engine_of) = engines(&conn, &tables)?;
    let scores = scores(&conn, &mut notes)?;
    let modifications = peptides::modifications(&conn)?;
    let samples = samples(&conn)?;
    let stored = spectra(&conn, &samples)?;
    let peptides = peptides::peptides(&conn, &modifications, &tables)?;
    let evidence = evidence(&conn)?;
    let tolerance = if tables.contains("IdentificationProtocol") {
        Tolerance::read(&conn)
    } else {
        Tolerance::default()
    };

    let raws = rank(read_psms(&conn, &engine_of, &scores, &mut notes)?, &scores);
    let context = Context {
        spectra: &stored,
        peptides: &peptides,
        evidence: &evidence,
        scores: &scores,
        tolerance,
    };
    let mut drafts: Vec<Draft> = Vec::new();
    let mut lost: HashSet<(usize, i64, &str)> = HashSet::new();
    for raw in &raws {
        match context.draft(raw, &mut notes) {
            Ok(draft) => drafts.push(draft),
            Err(reason) => {
                if raw.rank == 1 {
                    lost.insert((raw.engine, raw.sample, raw.spectrum_id.as_str()));
                }
                notes.left_out(reason);
            }
        }
    }
    if notes.mass_misses as f64 > MAX_MASS_MISSES * notes.checked.max(1) as f64 {
        return Err(SfdbError::Layout(format!(
            "{} of {} PSMs miss their precursor mass (e.g. {}), so the modifications were \
             probably misread",
            notes.mass_misses,
            notes.checked,
            notes.mass_examples.as_deref().unwrap_or("")
        )));
    }
    let before = drafts.len();
    drafts.retain(|draft| !lost.contains(&(draft.engine, draft.sample, draft.spectrum_id)));
    notes.lost_spectra = lost.len();
    notes.lost_psms = before - drafts.len();
    if drafts.is_empty() && !raws.is_empty() {
        return Err(SfdbError::Layout("no PSMs could be read".into()));
    }

    let (runs, mut spectra, keys) = arrange(stored, &samples)?;
    notes.scan_clashes = drop_clashing_scans(&mut spectra);
    let proteins = proteins(&conn, &drafts, &mut notes)?;
    notes.no_decoys = proteins.values().all(|protein| !protein.decoy);
    let position: HashMap<(i64, &str), usize> = keys
        .iter()
        .enumerate()
        .map(|(index, (sample, id))| ((*sample, id.as_str()), index))
        .collect();
    for draft in drafts {
        let index = position[&(draft.sample, draft.spectrum_id)];
        spectra[index].matches.push(draft.found);
    }
    for spectrum in &mut spectra {
        spectrum
            .matches
            .sort_by_key(|found| (found.engine, found.rank));
    }
    let layout = if tables.contains("SpectrumData") && stored_blobs(&conn) {
        Layout::Blob
    } else {
        Layout::Text
    };
    let name = path.file_stem().map_or_else(
        || "project".into(),
        |stem| stem.to_string_lossy().into_owned(),
    );
    Ok(Project {
        name,
        format: Format::Sfdb,
        scaffold_version: version,
        runs,
        engines,
        spectra,
        proteins,
        groups: Vec::new(),
        thresholds: None,
        notes: notes.lines(&scores),
        peaks: Some(PeakSource::Sfdb(SfdbPeaks::new(conn, layout, keys))),
    })
}

fn table_names(conn: &Connection) -> rusqlite::Result<HashSet<String>> {
    let mut statement = conn.prepare("SELECT name FROM sqlite_master WHERE type = 'table'")?;
    let names = statement
        .query_map([], |row| row.get(0))?
        .collect::<rusqlite::Result<_>>()?;
    Ok(names)
}

fn is_verdict(name: &str) -> bool {
    let lowered = name.to_ascii_lowercase();
    VERDICTS.iter().any(|prefix| lowered.starts_with(prefix))
}

/// One `PeptideSpectrumMatch` row with a main score, ranked.
struct Raw {
    id: i64,
    engine: usize,
    sample: i64,
    spectrum_id: String,
    peptide_ref: String,
    charge: i64,
    values: Vec<f64>,
    main: f64,
    rank: u32,
}

/// A PSM that was read, before it joins its spectrum.
struct Draft<'a> {
    engine: usize,
    sample: i64,
    spectrum_id: &'a str,
    found: Match,
}

/// What every PSM row is read against.
struct Context<'a> {
    spectra: &'a Stored,
    peptides: &'a HashMap<String, Readings>,
    evidence: &'a HashMap<String, Vec<Place>>,
    scores: &'a Scores,
    tolerance: Tolerance,
}

impl Context<'_> {
    /// The PSM a row is, or why it is left out (a mass miss is counted here).
    fn draft<'r>(&self, raw: &'r Raw, notes: &mut Notes) -> Result<Draft<'r>, String> {
        let precursor_mz = self
            .spectra
            .by_key
            .get(&(raw.sample, raw.spectrum_id.clone()))
            .and_then(|&row| self.spectra.rows.get(row))
            .ok_or("its spectrum isn't stored")?
            .mz;
        let readings = match self.peptides.get(&raw.peptide_ref) {
            Some(Ok(readings)) => readings,
            Some(Err(reason)) => return Err(reason.clone()),
            None => return Err("its peptide isn't stored".into()),
        };
        let charge = u8::try_from(raw.charge)
            .ok()
            .filter(|charge| *charge >= 1)
            .ok_or("it has no charge")?;
        let precursor = precursor_mz.map(|mz| (mz - PROTON) * f64::from(charge));
        let masses: Vec<Option<f64>> = readings
            .iter()
            .map(|reading| peptide_neutral(&reading.peptide))
            .collect();
        let chosen = match (precursor, masses.first().copied().flatten()) {
            (Some(precursor), Some(_)) => {
                notes.checked += 1;
                let fitting = masses
                    .iter()
                    .position(|mass| mass.is_some_and(|mass| self.tolerance.fits(mass, precursor)));
                let Some(index) = fitting else {
                    notes.mass_misses += 1;
                    notes.mass_examples.get_or_insert_with(|| {
                        format!(
                            "{} {charge}+ at m/z {:.4}",
                            raw.peptide_ref,
                            precursor_mz.unwrap_or_default()
                        )
                    });
                    return Err("it misses its precursor mass".into());
                };
                index
            }
            _ => {
                notes.unchecked += 1;
                0
            }
        };
        let reading = readings.get(chosen).ok_or("its peptide has no reading")?;
        if !reading.named {
            notes.other_terminus += 1;
        }
        if reading.met_lost {
            notes.met_lost += 1;
        }
        let proteins = self
            .evidence
            .get(&raw.peptide_ref)
            .map(|places| places.iter().map(Place::to_ref).collect())
            .unwrap_or_default();
        Ok(Draft {
            engine: raw.engine,
            sample: raw.sample,
            spectrum_id: &raw.spectrum_id,
            found: Match {
                engine: raw.engine,
                rank: raw.rank,
                peptide: reading.peptide.clone(),
                mass: masses.get(chosen).copied().flatten(),
                charge: Some(charge),
                scores: engine_scores(self.scores, &raw.values),
                proteins,
                probability: None,
                chosen: false,
            },
        })
    }
}

/// The search engines, and identification protocol -> engine. A protocol
/// credited to Scaffold itself (Scaffold 5 does that) belongs to the file's
/// only other program; with several, the file can't say which ran, and is
/// refused.
fn engines(
    conn: &Connection,
    tables: &HashSet<String>,
) -> Result<(Vec<Engine>, HashMap<i64, usize>), SfdbError> {
    let mut statement =
        conn.prepare("SELECT ID, SOFTWARE_NAME, SOFTWARE_VERSION FROM AnalysisSoftware")?;
    let rows: Vec<(i64, String, Option<String>)> = statement
        .query_map([], |row| Ok((row.get(0)?, row.get(1)?, row.get(2)?)))?
        .collect::<rusqlite::Result<_>>()?;
    let versions: HashMap<String, Option<String>> = rows
        .iter()
        .map(|(_, name, version)| (name.clone(), version.clone()))
        .collect();
    let software: HashMap<i64, String> = rows.into_iter().map(|(id, name, _)| (id, name)).collect();
    let others: BTreeSet<&String> = software.values().filter(|name| !is_verdict(name)).collect();
    let sole = if others.len() == 1 {
        others.iter().next().map(|name| (*name).clone())
    } else {
        None
    };
    let mut statement = conn.prepare("SELECT DISTINCT PROTOCOL_ID FROM PeptideSpectrumMatch")?;
    let used: Vec<i64> = statement
        .query_map([], |row| row.get(0))?
        .collect::<rusqlite::Result<_>>()?;
    let mut protocols: HashMap<i64, i64> = HashMap::new();
    if tables.contains("IdentificationProtocol") {
        let mut statement =
            conn.prepare("SELECT ID, ANALYSIS_SOFTWARE_ID FROM IdentificationProtocol")?;
        protocols = statement
            .query_map([], |row| Ok((row.get(0)?, row.get(1)?)))?
            .collect::<rusqlite::Result<_>>()?;
    }
    let mut named: BTreeMap<i64, String> = BTreeMap::new();
    for protocol in used {
        let engine = protocols
            .get(&protocol)
            .and_then(|software_id| software.get(software_id))
            .filter(|name| !is_verdict(name))
            .cloned()
            .or_else(|| sole.clone());
        let Some(engine) = engine else {
            return Err(SfdbError::Layout(format!(
                "identification protocol {protocol} isn't tied to one search engine (the file \
                 lists {})",
                others
                    .iter()
                    .map(|name| name.as_str())
                    .collect::<Vec<_>>()
                    .join(", ")
            )));
        };
        named.insert(protocol, engine);
    }
    let names: BTreeSet<&String> = named.values().collect();
    let engines: Vec<Engine> = names
        .iter()
        .map(|name| Engine {
            name: (*name).clone(),
            version: versions.get(*name).cloned().flatten(),
            inferred_scores: false,
        })
        .collect();
    let index: HashMap<&String, usize> = names
        .iter()
        .enumerate()
        .map(|(index, name)| (*name, index))
        .collect();
    let engine_of = named
        .iter()
        .map(|(protocol, name)| (*protocol, index[name]))
        .collect();
    Ok((engines, engine_of))
}

/// The engine scores in `PEPTIDE_SCORE_ARRAY` order.
#[derive(Default)]
struct Scores {
    /// (index into the array, name, PSI-MS accession)
    kept: Vec<(usize, String, Option<String>)>,
    left_out: Vec<String>,
    /// Index of the score ranks are by.
    main: usize,
    lower_is_better: bool,
    /// Indices of Mascot's ion score and identity threshold, for its
    /// expectation value.
    mascot: Option<(usize, usize)>,
}

fn scores(conn: &Connection, notes: &mut Notes) -> Result<Scores, SfdbError> {
    let mut statement = conn.prepare(
        "SELECT SCORE_INDEX, SCORE_NAME, SCORE_ACCESSION, SCORE_DECREASING FROM ScoreIndex \
         WHERE SCORE_TYPE = 'Peptide Score' ORDER BY SCORE_INDEX",
    )?;
    let rows: Vec<(i64, String, Option<String>, i64)> = statement
        .query_map([], |row| {
            Ok((row.get(0)?, row.get(1)?, row.get(2)?, row.get(3)?))
        })?
        .collect::<rusqlite::Result<_>>()?;
    let mut scores = Scores::default();
    let mut main = None;
    let mut seen = HashSet::new();
    for (index, name, accession, decreasing) in rows {
        let Some(index) = usize::try_from(index)
            .ok()
            .filter(|_| (0..MAX_SCORES).contains(&index))
        else {
            return Err(SfdbError::Layout(format!("score {name} has index {index}")));
        };
        if !seen.insert(index) {
            return Err(SfdbError::Layout(format!(
                "two peptide scores share index {index}"
            )));
        }
        if is_verdict(&name) {
            scores.left_out.push(name);
            continue;
        }
        let lowered = name.to_ascii_lowercase();
        let ranks_on = !lowered.contains("threshold") && !lowered.contains("rank");
        if ranks_on && main.is_none() {
            main = Some(index);
            scores.lower_is_better = decreasing == 1;
        }
        scores.kept.push((index, name, accession));
    }
    let find = |wanted: &str| {
        scores
            .kept
            .iter()
            .find(|(_, name, _)| name.eq_ignore_ascii_case(wanted))
            .map(|(index, _, _)| *index)
    };
    scores.mascot = find("Mascot:score").zip(find("Mascot:identity threshold"));
    let Some(main) = main else {
        return Err(SfdbError::Layout(format!(
            "no engine score to rank by among {}",
            scores
                .kept
                .iter()
                .map(|(_, name, _)| name.as_str())
                .chain(scores.left_out.iter().map(String::as_str))
                .collect::<Vec<_>>()
                .join(", ")
        )));
    };
    scores.main = main;
    notes.scores_read = scores
        .kept
        .iter()
        .map(|(_, name, _)| name.clone())
        .collect();
    Ok(scores)
}

fn parse_scores(array: &str) -> Vec<f64> {
    array
        .split(',')
        .map(|value| value.trim().parse::<f64>().unwrap_or(f64::NAN))
        .collect()
}

fn value_at(values: &[f64], index: usize) -> Option<f64> {
    values.get(index).copied().filter(|value| value.is_finite())
}

fn engine_scores(scores: &Scores, values: &[f64]) -> Vec<Score> {
    let mut out: Vec<Score> = scores
        .kept
        .iter()
        .filter_map(|(index, name, accession)| {
            value_at(values, *index)
                .map(|value| Score::new(name.clone(), accession.clone().map(Into::into), value))
        })
        .collect();
    let mascot = scores
        .mascot
        .and_then(|(score, identity)| value_at(values, score).zip(value_at(values, identity)));
    if let Some((score, identity)) = mascot {
        let value = expectation(score, identity);
        if value.is_finite() {
            out.push(Score::new(
                "Mascot:expectation value",
                Some("MS:1001172".into()),
                value,
            ));
        }
    }
    out
}

/// Every PSM row. A spectrum with a match lacking the main score can't be
/// ranked, so all its matches are left out, counted.
fn read_psms(
    conn: &Connection,
    engine_of: &HashMap<i64, usize>,
    scores: &Scores,
    notes: &mut Notes,
) -> Result<Vec<Raw>, SfdbError> {
    let mut statement = conn.prepare(
        "SELECT ID, MSSAMPLE_ID, SPECTRUM_ID, PEPTIDE_REF, ASSIGNED_CHARGE, \
         PEPTIDE_SCORE_ARRAY, PROTOCOL_ID FROM PeptideSpectrumMatch",
    )?;
    let mut rows = statement.query([])?;
    let mut raws = Vec::new();
    let mut unranked: HashSet<(usize, i64, String)> = HashSet::new();
    while let Some(row) = rows.next()? {
        let id: i64 = row.get(0)?;
        let protocol: i64 = row.get(6)?;
        let Some(&engine) = engine_of.get(&protocol) else {
            return Err(SfdbError::Layout(format!(
                "PSM {id} names identification protocol {protocol}, which the file doesn't \
                 tie to a search engine"
            )));
        };
        let values = parse_scores(&row.get::<_, String>(5)?);
        let sample: i64 = row.get(1)?;
        let spectrum_id: String = row.get(2)?;
        let main = value_at(&values, scores.main);
        if main.is_none() {
            unranked.insert((engine, sample, spectrum_id.clone()));
        }
        raws.push(Raw {
            id,
            engine,
            sample,
            spectrum_id,
            peptide_ref: row.get(3)?,
            charge: row.get(4)?,
            values,
            main: main.unwrap_or(f64::NAN),
            rank: 0,
        });
    }
    let before = raws.len();
    raws.retain(|raw| !unranked.contains(&(raw.engine, raw.sample, raw.spectrum_id.clone())));
    if before > raws.len() {
        notes.left_out_n(
            "a match of their spectrum lacks the engine's main score, so ranks are unknown",
            before - raws.len(),
        );
    }
    Ok(raws)
}

/// Rank each spectrum's PSMs (per engine) by the main score; equal scores
/// keep the file's order (by PSM id).
fn rank(mut raws: Vec<Raw>, scores: &Scores) -> Vec<Raw> {
    raws.sort_by(|a, b| {
        (a.engine, a.sample, &a.spectrum_id)
            .cmp(&(b.engine, b.sample, &b.spectrum_id))
            .then_with(|| {
                let order = a.main.total_cmp(&b.main);
                if scores.lower_is_better {
                    order
                } else {
                    order.reverse()
                }
            })
            .then(a.id.cmp(&b.id))
    });
    let mut previous: Option<(usize, i64, String, u32)> = None;
    for raw in &mut raws {
        let rank = match &previous {
            Some((engine, sample, spectrum, rank))
                if (*engine, *sample, spectrum) == (raw.engine, raw.sample, &raw.spectrum_id) =>
            {
                rank + 1
            }
            _ => 1,
        };
        raw.rank = rank;
        previous = Some((raw.engine, raw.sample, raw.spectrum_id.clone(), rank));
    }
    raws
}

/// MS sample id -> name (Scaffold's biosample).
fn samples(conn: &Connection) -> Result<HashMap<i64, String>, SfdbError> {
    let mut statement = conn.prepare("SELECT ID, SAMPLE_NAME FROM MSSample")?;
    let samples: HashMap<i64, String> = statement
        .query_map([], |row| Ok((row.get(0)?, row.get(1)?)))?
        .collect::<rusqlite::Result<_>>()?;
    if samples.is_empty() {
        return Err(SfdbError::Layout("no MS samples".into()));
    }
    Ok(samples)
}

/// One stored spectrum row.
struct Row {
    sample: i64,
    spectrum_id: String,
    run: String,
    scan: Option<u32>,
    title: Option<String>,
    mz: Option<f64>,
    charge: Option<u8>,
    rt_seconds: Option<f64>,
}

/// Every `Spectrum` row in file order, and (sample, id) -> row.
struct Stored {
    rows: Vec<Row>,
    by_key: HashMap<(i64, String), usize>,
}

fn spectra(conn: &Connection, samples: &HashMap<i64, String>) -> Result<Stored, SfdbError> {
    let mut statement = conn.prepare(
        "SELECT MSSAMPLE_ID, SPECTRUM_ID, SPECTRUM_TITLE, ACQUISITION_MASS, ACQUISITION_CHARGE, \
         SCAN_START_TIME_IN_SEC FROM Spectrum ORDER BY rowid",
    )?;
    let mut rows = statement.query([])?;
    let mut stored = Stored {
        rows: Vec::new(),
        by_key: HashMap::new(),
    };
    while let Some(row) = rows.next()? {
        let sample: i64 = row.get(0)?;
        let spectrum_id: String = row.get(1)?;
        let title: Option<String> = row.get(2)?;
        let Some(sample_name) = samples.get(&sample) else {
            return Err(SfdbError::Layout(format!(
                "spectrum {spectrum_id} is in MS sample {sample}, which isn't listed"
            )));
        };
        let title = title.filter(|title| !title.trim().is_empty());
        let (run, scan) = locate(title.as_deref(), sample_name);
        let key = (sample, spectrum_id.clone());
        if stored.by_key.insert(key, stored.rows.len()).is_some() {
            return Err(SfdbError::Layout(format!(
                "spectrum {spectrum_id} of MS sample {sample} is stored twice"
            )));
        }
        stored.rows.push(Row {
            sample,
            spectrum_id,
            run,
            scan,
            title,
            mz: row
                .get::<_, Option<f64>>(3)?
                .filter(|mz| mz.is_finite() && *mz > 0.0),
            charge: row
                .get::<_, Option<i64>>(4)?
                .and_then(|charge| u8::try_from(charge).ok())
                .filter(|charge| *charge > 0),
            rt_seconds: row
                .get::<_, Option<f64>>(5)?
                .filter(|rt| rt.is_finite() && *rt >= 0.0),
        });
    }
    Ok(stored)
}

/// The run and scan of a spectrum from its title: Scaffold 5's
/// `File: "...\X.raw"; SpectrumID: "1445"; scans: "5791"` names both;
/// Scaffold DDA's `1-7_7` (id, first and last scan) names the scan and the
/// MS sample is the run; anything else as `.sf3` titles are read.
fn locate(title: Option<&str>, sample: &str) -> (String, Option<u32>) {
    let title = title.map(str::trim);
    if let Some(file) = title.and_then(|title| quoted(title, "File:")) {
        let scan = title
            .and_then(|title| quoted(title, "scans:"))
            .and_then(|scans| scans.parse::<u32>().ok());
        return (run_name(&file), scan);
    }
    let scan = title.and_then(|title| dda_scan(title).or_else(|| scan_from_title(title)));
    (sample.to_string(), scan)
}

/// The value after `label` in `label "value"`.
fn quoted(title: &str, label: &str) -> Option<String> {
    let rest = &title[title.find(label)? + label.len()..];
    let rest = rest.trim_start().strip_prefix('"')?;
    Some(rest[..rest.find('"')?].to_string())
}

/// `ID-FIRST_LAST` with one scan.
fn dda_scan(title: &str) -> Option<u32> {
    let (id, scans) = title.split_once('-')?;
    let (first, last) = scans.split_once('_')?;
    id.parse::<u32>().ok()?;
    let (first, last) = (first.parse::<u32>().ok()?, last.parse::<u32>().ok()?);
    (first == last).then_some(first)
}

/// The run a raw file path names: its file name without a spectra
/// extension (`C:\raw\HeLa_01.RAW` -> `HeLa_01`).
fn run_name(location: &str) -> String {
    const SPECTRA_EXTENSIONS: &[&str] = &[
        ".mzml", ".mzxml", ".mgf", ".raw", ".d", ".wiff", ".ms2", ".mzml.gz", ".mgf.gz",
    ];
    let name = location.rsplit(['/', '\\']).next().unwrap_or(location);
    let lowered = name.to_ascii_lowercase();
    SPECTRA_EXTENSIONS
        .iter()
        .filter(|extension| lowered.ends_with(*extension))
        .max_by_key(|extension| extension.len())
        .map_or(name, |extension| &name[..name.len() - extension.len()])
        .to_string()
}

/// Runs (by name), and the spectra run by run in file order, with each
/// spectrum's (MS sample, spectrum id). A run must lie in one MS sample.
#[allow(clippy::type_complexity)]
fn arrange(
    stored: Stored,
    samples: &HashMap<i64, String>,
) -> Result<(Vec<Run>, Vec<Spectrum>, Vec<(i64, String)>), SfdbError> {
    let mut run_sample: BTreeMap<&str, i64> = BTreeMap::new();
    for row in &stored.rows {
        match run_sample.get(row.run.as_str()) {
            Some(&other) if other != row.sample => {
                return Err(SfdbError::Layout(format!(
                    "run {} is in two MS samples ({}, {})",
                    row.run,
                    samples.get(&other).map_or("?", String::as_str),
                    samples.get(&row.sample).map_or("?", String::as_str)
                )));
            }
            Some(_) => {}
            None => {
                run_sample.insert(&row.run, row.sample);
            }
        }
    }
    let runs: Vec<Run> = run_sample
        .iter()
        .map(|(name, sample)| Run {
            name: (*name).to_string(),
            biosample: samples.get(sample).cloned().unwrap_or_default(),
            category: None,
        })
        .collect();
    let run_index: HashMap<String, usize> = runs
        .iter()
        .enumerate()
        .map(|(index, run)| (run.name.clone(), index))
        .collect();
    let mut order: Vec<(usize, usize)> = stored
        .rows
        .iter()
        .enumerate()
        .map(|(position, row)| (run_index[&row.run], position))
        .collect();
    order.sort_unstable();
    let mut rows: Vec<Option<Row>> = stored.rows.into_iter().map(Some).collect();
    let mut spectra = Vec::with_capacity(rows.len());
    let mut keys = Vec::with_capacity(rows.len());
    for (run, position) in order {
        let Some(row) = rows.get_mut(position).and_then(Option::take) else {
            continue;
        };
        let id = i32::try_from(spectra.len() + 1)
            .map_err(|_| SfdbError::Layout("more spectra than the decoder can number".into()))?;
        let mut spectrum = Spectrum::from_mz(id, run, row.title, row.mz, row.charge);
        spectrum.scan = row.scan;
        spectrum.rt_seconds = row.rt_seconds;
        spectrum.peak_list = id;
        spectra.push(spectrum);
        keys.push((row.sample, row.spectrum_id));
    }
    Ok((runs, spectra, keys))
}

/// Where a peptide sits on one protein.
struct Place {
    accession: String,
    start: Option<u32>,
    end: Option<u32>,
    pre: Option<char>,
    post: Option<char>,
}

impl Place {
    fn to_ref(&self) -> ProteinRef {
        let mut place = ProteinRef::new(self.accession.clone(), self.pre, self.post);
        place.start = self.start;
        place.end = self.end;
        place
    }
}

/// Peptide ref -> where it sits on each protein, once per protein (a second
/// place on the same protein, a repeat, is not kept). A start below 1
/// (Scaffold DDA 7.1 is reported to write -1) is rebuilt from the stop and
/// the peptide's length.
fn evidence(conn: &Connection) -> Result<HashMap<String, Vec<Place>>, SfdbError> {
    let mut statement = conn.prepare(
        "SELECT e.PEPTIDE_REF, e.PROTEIN_ACCESSION, e.PROTEIN_START_INDEX, \
         e.PROTEIN_STOP_INDEX, e.PRE_AA, e.POST_AA, length(p.PEPTIDE_SEQUENCE) \
         FROM PeptideEvidence e LEFT JOIN Peptide p ON p.PEPTIDE_REF = e.PEPTIDE_REF \
         ORDER BY e.PEPTIDE_REF, e.PROTEIN_ACCESSION, e.PROTEIN_START_INDEX",
    )?;
    let mut rows = statement.query([])?;
    let mut out: HashMap<String, Vec<Place>> = HashMap::new();
    while let Some(row) = rows.next()? {
        let reference: String = row.get(0)?;
        let accession: String = row.get(1)?;
        let position = |value: Option<i64>| {
            value
                .and_then(|value| u32::try_from(value).ok())
                .filter(|value| *value > 0)
        };
        let end = position(row.get(3)?);
        let length = position(row.get(6)?);
        let start = position(row.get(2)?).or_else(|| {
            end.zip(length)
                .filter(|(end, length)| end >= length)
                .map(|(end, length)| end - length + 1)
        });
        // `-` marks the protein's terminus.
        let residue = |value: Option<String>| {
            value
                .and_then(|value| value.chars().next())
                .map(|residue| residue.to_ascii_uppercase())
                .filter(|residue| residue.is_ascii_uppercase() || *residue == '-')
        };
        let pre = residue(row.get(4)?);
        let post = residue(row.get(5)?);
        let list = out.entry(reference).or_default();
        if list.iter().any(|place| place.accession == accession) {
            continue;
        }
        list.push(Place {
            accession,
            start,
            end,
            pre,
            post,
        });
    }
    Ok(out)
}

/// The proteins PSMs name. A protein Scaffold flags `IS_DECOY` 1 is a decoy;
/// one flagged 0 follows the accession convention, since the flag is only
/// Scaffold's own decoy pattern, not the engine's statement. A stored
/// sequence is kept only if it holds every peptide placed on it (I = L).
fn proteins(
    conn: &Connection,
    drafts: &[Draft],
    notes: &mut Notes,
) -> Result<BTreeMap<String, Protein>, SfdbError> {
    let mut peptides_of: HashMap<&str, HashSet<&str>> = HashMap::new();
    for draft in drafts {
        for place in &draft.found.proteins {
            peptides_of
                .entry(place.accession.as_str())
                .or_default()
                .insert(draft.found.peptide.sequence.as_str());
        }
    }
    let mut statement = conn.prepare(
        "SELECT ACCESSION_NUMBER, ORIGINAL_PROTEIN_NAME, PROTEIN_SEQUENCE, IS_DECOY FROM Protein",
    )?;
    let mut rows = statement.query([])?;
    let mut out = BTreeMap::new();
    while let Some(row) = rows.next()? {
        let accession: String = row.get(0)?;
        let Some(peptides) = peptides_of.get(accession.as_str()) else {
            continue;
        };
        if out.contains_key(&accession) {
            continue;
        }
        let description: Option<String> = row.get(1)?;
        let sequence: Option<String> = row.get(2)?;
        let flagged: Option<i64> = row.get(3)?;
        let sequence = sequence
            .map(|sequence| sequence.trim().to_ascii_uppercase())
            .filter(|sequence| !sequence.is_empty())
            .filter(|sequence| {
                let leucine = sequence.replace('I', "L");
                let holds = peptides
                    .iter()
                    .all(|peptide| leucine.contains(&peptide.replace('I', "L")));
                if !holds {
                    notes.dropped_sequences += 1;
                }
                holds
            });
        let decoy = flagged == Some(1) || is_decoy_accession(&accession);
        out.insert(
            accession.clone(),
            Protein {
                description: description.filter(|text| !text.is_empty() && *text != accession),
                sequence,
                decoy,
            },
        );
    }
    for accession in peptides_of.keys() {
        if !out.contains_key(*accession) {
            notes.unlisted_proteins += 1;
            out.insert(
                (*accession).to_string(),
                Protein {
                    description: None,
                    sequence: None,
                    decoy: is_decoy_accession(accession),
                },
            );
        }
    }
    Ok(out)
}

fn stored_blobs(conn: &Connection) -> bool {
    conn.query_row("SELECT 1 FROM SpectrumData LIMIT 1", [], |row| {
        row.get::<_, i64>(0)
    })
    .is_ok()
}

#[derive(Default)]
struct Notes {
    left_out: BTreeMap<String, usize>,
    checked: usize,
    unchecked: usize,
    mass_misses: usize,
    mass_examples: Option<String>,
    other_terminus: usize,
    met_lost: usize,
    lost_spectra: usize,
    lost_psms: usize,
    scan_clashes: usize,
    dropped_sequences: usize,
    unlisted_proteins: usize,
    scores_read: Vec<String>,
    no_decoys: bool,
}

impl Notes {
    fn left_out(&mut self, reason: String) {
        self.left_out_n(&reason, 1);
    }

    fn left_out_n(&mut self, reason: &str, count: usize) {
        *self.left_out.entry(reason.to_string()).or_default() += count;
    }

    fn lines(&self, scores: &Scores) -> Vec<String> {
        let mut lines = vec![format!(
            "Read the engines' own scores ({}). Scaffold's probabilities, protein groups, \
             thresholds and other verdicts{} were not read: what they mean in an .sfdb isn't \
             settled.",
            self.scores_read.join(", "),
            if scores.left_out.is_empty() {
                String::new()
            } else {
                format!(" ({})", scores.left_out.join(", "))
            }
        )];
        if self.no_decoys {
            lines.push(
                "The file holds no decoy proteins (no decoy search, or Scaffold kept none).".into(),
            );
        }
        for (reason, count) in &self.left_out {
            lines.push(format!("{count} PSMs were left out: {reason}."));
        }
        if self.mass_misses > 0 {
            lines.push(format!(
                "{} of {} PSMs checked miss their precursor mass under every reading of their \
                 peptide (e.g. {}).",
                self.mass_misses,
                self.checked,
                self.mass_examples.as_deref().unwrap_or("")
            ));
        }
        if self.lost_spectra > 0 {
            lines.push(format!(
                "{} spectra whose best match was left out lost all their matches of that \
                 engine ({} more PSMs), so a lower match never counts as rank 1.",
                self.lost_spectra, self.lost_psms
            ));
        }
        if self.unchecked > 0 {
            lines.push(format!(
                "{} PSMs couldn't be checked against their precursor mass (no precursor m/z, \
                 or a residue without a mass).",
                self.unchecked
            ));
        }
        if self.other_terminus > 0 {
            lines.push(format!(
                "{} PSMs fit their precursor only with N-terminal modifications other than \
                 their Scaffold peptide name shows ({} of them with the initial methionine \
                 lost); Scaffold files such PSMs under one peptide.",
                self.other_terminus, self.met_lost
            ));
        }
        if self.scan_clashes > 0 {
            lines.push(format!(
                "{} spectra name a scan another spectrum of their run names too; their scan \
                 number is left out.",
                self.scan_clashes
            ));
        }
        if self.dropped_sequences > 0 {
            lines.push(format!(
                "{} stored protein sequences don't contain their peptides and were left out.",
                self.dropped_sequences
            ));
        }
        if self.unlisted_proteins > 0 {
            lines.push(format!(
                "{} proteins PSMs name aren't in the file's protein table; they are kept by \
                 accession alone (decoys by the accession convention).",
                self.unlisted_proteins
            ));
        }
        lines
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn titles_name_run_and_scan() {
        let (run, scan) = locate(
            Some(r#"File: "Z:\CORE\RAW\P706_WT_7T.raw"; SpectrumID: "1445"; scans: "5791""#),
            "mascot_daemon_merge (F015729)",
        );
        assert_eq!((run.as_str(), scan), ("P706_WT_7T", Some(5791)));
        assert_eq!(
            locate(Some("1-7_7"), "Elute-1-1"),
            ("Elute-1-1".to_string(), Some(7))
        );
        // Summed scans have no single scan.
        assert_eq!(locate(Some("4-110_112"), "Elute-1-1").1, None);
        assert_eq!(locate(None, "s"), ("s".to_string(), None));
    }

    #[test]
    fn run_names_drop_folders_and_spectra_extensions() {
        assert_eq!(run_name("/data/run1.mzML"), "run1");
        assert_eq!(run_name("C:\\raw\\HeLa_01.RAW"), "HeLa_01");
        assert_eq!(run_name("/data/runA"), "runA");
        assert_eq!(run_name("x.mzML.gz"), "x");
    }

    fn raw(id: i64, spectrum: &str, main: f64) -> Raw {
        Raw {
            id,
            engine: 0,
            sample: 1,
            spectrum_id: spectrum.into(),
            peptide_ref: String::new(),
            charge: 2,
            values: Vec::new(),
            main,
            rank: 0,
        }
    }

    #[test]
    fn ranks_follow_the_main_score_per_spectrum() {
        let higher = Scores::default();
        let ranked = rank(
            vec![
                raw(1, "a", 10.0),
                raw(2, "a", 30.0),
                raw(3, "b", 5.0),
                raw(4, "a", 30.0),
            ],
            &higher,
        );
        let ranks: Vec<(i64, u32)> = ranked.iter().map(|r| (r.id, r.rank)).collect();
        // Equal scores keep the file's order.
        assert_eq!(ranks, [(2, 1), (4, 2), (1, 3), (3, 1)]);
        let lower = Scores {
            lower_is_better: true,
            ..Scores::default()
        };
        let ranked = rank(vec![raw(1, "a", 0.9), raw(2, "a", 0.01)], &lower);
        assert_eq!(ranked[0].id, 2);
    }
}
