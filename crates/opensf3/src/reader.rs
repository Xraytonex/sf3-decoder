//! Opening a Scaffold file of either kind, and its stored peak lists as
//! [`SpectrumSource`]s.
//!
//! A Scaffold file holds many MS runs, and an mzML file describes one, so
//! the spectra are offered one source per run ([`Reader::run_spectra`]).
//! Each run's spectra come in the order its MGF lists them, and spectrum
//! `N` (from 0) of a run is `index=N` in that run's MGF, mzML, and in the
//! mzIdentML's `spectrumID`, so the three files line up.

use std::io::Read;
use std::path::{Path, PathBuf};

use openmassspec_core::{
    CvTerm, PrecursorInfo, RunMetadata, ScanMode, SpectrumIter, SpectrumRecord, SpectrumSource,
};
use rusqlite::{Connection, OpenFlags};

use crate::error::Error;
use crate::model::{Format, Project, Run, Spectrum};

const SQLITE_MAGIC: &[u8; 16] = b"SQLite format 3\0";

/// An opened Scaffold `.sf3` or `.sfdb` file, decoded.
pub struct Reader {
    path: PathBuf,
    project: Project,
}

impl std::fmt::Debug for Reader {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Reader")
            .field("path", &self.path)
            .field("format", &self.project.format)
            .field("runs", &self.project.runs.len())
            .field("spectra", &self.project.spectra.len())
            .finish()
    }
}

impl Reader {
    /// Open and decode a Scaffold file. Which kind it is is told from its
    /// content, never its name: both are SQLite databases, an `.sf3` with
    /// a `Files` table of compressed inner files, an `.sfdb` with
    /// mzIdentML-like tables (`PeptideSpectrumMatch`, `Spectrum`, ...).
    pub fn open(path: impl AsRef<Path>) -> Result<Self, Error> {
        let path = path.as_ref();
        let project = match detect(path)? {
            Format::Sf3 => crate::sf3::decode(path)?,
            Format::Sfdb => crate::sfdb::decode(path)?,
        };
        Ok(Self {
            path: path.to_path_buf(),
            project,
        })
    }

    pub fn path(&self) -> &Path {
        &self.path
    }

    pub fn format(&self) -> Format {
        self.project.format
    }

    /// The identifications and everything else the file holds.
    pub fn project(&self) -> &Project {
        &self.project
    }

    pub fn into_project(self) -> Project {
        self.project
    }

    pub fn runs(&self) -> &[Run] {
        &self.project.runs
    }

    /// One run's stored peak lists as a [`SpectrumSource`].
    pub fn run_spectra(&self, run: usize) -> Result<RunSpectra<'_>, Error> {
        if run >= self.project.runs.len() {
            return Err(Error::NoSuchRun {
                run,
                runs: self.project.runs.len(),
            });
        }
        Ok(RunSpectra {
            project: &self.project,
            run,
            file_name: self
                .path
                .file_name()
                .map_or_else(String::new, |name| name.to_string_lossy().into_owned()),
        })
    }
}

/// Which kind of Scaffold file `path` is, from its content.
pub fn detect(path: &Path) -> Result<Format, Error> {
    let mut magic = [0u8; 16];
    std::fs::File::open(path)?
        .read_exact(&mut magic)
        .map_err(|_| Error::Unrecognized("too short to be a SQLite database".into()))?;
    if &magic != SQLITE_MAGIC {
        return Err(Error::Unrecognized("not a SQLite database".into()));
    }
    let conn = Connection::open_with_flags(
        path,
        OpenFlags::SQLITE_OPEN_READ_ONLY | OpenFlags::SQLITE_OPEN_NO_MUTEX,
    )
    .map_err(|error| Error::Unrecognized(error.to_string()))?;
    let tables: Vec<String> = conn
        .prepare("SELECT name FROM sqlite_master WHERE type = 'table'")
        .and_then(|mut statement| {
            statement
                .query_map([], |row| row.get(0))?
                .collect::<rusqlite::Result<_>>()
        })
        .map_err(|error| Error::Unrecognized(error.to_string()))?;
    let has = |table: &str| tables.iter().any(|name| name == table);
    if has("Files") {
        Ok(Format::Sf3)
    } else if crate::sfdb::REQUIRED.iter().all(|table| has(table)) {
        Ok(Format::Sfdb)
    } else {
        Err(Error::Unrecognized(
            "a SQLite database without the tables of an .sf3 or .sfdb".into(),
        ))
    }
}

/// The stored MS/MS peak lists of one run.
///
/// Every spectrum becomes an MS2, centroided [`SpectrumRecord`]:
/// - `native_id` is `index=N`, its position in the run (the PSI-MS
///   "multiple peak list nativeID format"), as in the run's MGF;
/// - `scan_number` is the scan its title names, 0 when it names none;
/// - the precursor has the stored m/z and charge (none when the file
///   stores no m/z, as for an `.sf3` spectrum without a charge);
/// - retention time is the stored one, 0 when none is stored (and
///   `extra["opensf3.no_retention_time"]` says so);
/// - `extra` holds `opensf3.title` (Scaffold's spectrum title) and
///   `opensf3.scan` when known.
///
/// Scaffold keeps at most 300 peaks per spectrum, so these are not the raw
/// spectra. A spectrum whose stored peaks can't be read is an `Err` item.
pub struct RunSpectra<'a> {
    project: &'a Project,
    run: usize,
    file_name: String,
}

impl RunSpectra<'_> {
    pub fn run(&self) -> &Run {
        &self.project.runs[self.run]
    }

    fn record(&self, index: usize, spectrum: &Spectrum) -> Result<SpectrumRecord, Error> {
        let native_id = format!("index={index}");
        let peaks = self
            .project
            .peaks(spectrum)
            .map_err(|problem| Error::Peaks {
                native_id: native_id.clone(),
                problem,
            })?;
        let (mz, intensity): (Vec<f64>, Vec<f32>) = peaks.unwrap_or_default().into_iter().unzip();
        let mut precursor = PrecursorInfo::default();
        if let Some(mz) = spectrum.precursor_mz() {
            precursor = precursor
                .with_selected_mz(mz)
                .with_charge(spectrum.charge.map(i32::from));
        }
        let mut record = SpectrumRecord::new(index, spectrum.scan.unwrap_or(0), native_id)
            .with_ms_level(2)
            .with_scan_mode(ScanMode::Centroid)
            .with_retention_time_sec(spectrum.rt_seconds.unwrap_or(0.0))
            .with_precursor(precursor)
            .with_peaks(mz, intensity);
        if spectrum.rt_seconds.is_none() {
            record = record.with_extra("opensf3.no_retention_time", "true");
        }
        if let Some(title) = &spectrum.title {
            record = record.with_extra("opensf3.title", title.clone());
        }
        if let Some(scan) = spectrum.scan {
            record = record.with_extra("opensf3.scan", scan.to_string());
        }
        Ok(record)
    }
}

impl SpectrumSource for RunSpectra<'_> {
    /// The source file is the Scaffold file; the run's own name, biosample
    /// and the Scaffold version are in `extra`. Scaffold files don't say
    /// which instrument recorded a run, so the instrument is the generic
    /// "instrument model" term.
    fn run_metadata(&self) -> RunMetadata {
        let run = self.run();
        let mut metadata = RunMetadata::new(
            self.file_name.clone(),
            CvTerm::new("MS:1000560", "mass spectrometer file format"),
            CvTerm::new("MS:1000774", "multiple peak list nativeID format"),
            CvTerm::new("MS:1000031", "instrument model"),
            "opensf3",
            env!("CARGO_PKG_VERSION"),
        )
        .with_extra("opensf3.format", self.project.format.as_str())
        .with_extra("opensf3.run", run.name.clone())
        .with_extra("opensf3.biosample", run.biosample.clone());
        if let Some(version) = &self.project.scaffold_version {
            metadata = metadata.with_extra("opensf3.scaffold_version", version.clone());
        }
        metadata
    }

    fn iter_spectra<'a>(&'a mut self) -> SpectrumIter<'a> {
        let this = &*self;
        Box::new(
            this.project
                .run_spectra(this.run)
                .enumerate()
                .map(move |(index, spectrum)| {
                    this.record(index, spectrum)
                        .map_err(openmassspec_core::Error::from)
                }),
        )
    }

    fn spectrum_count_hint(&self) -> Option<usize> {
        Some(self.project.run_spectra(self.run).count())
    }
}
