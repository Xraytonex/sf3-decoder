//! Scaffold `.sfdb` reading. Tiny databases in both layouts (Scaffold 5 and
//! Scaffold DDA) are written here, with only the tables and columns the
//! decoder reads; the expected values are derived by hand below. Real files
//! run only when `SFDB_SAMPLES` names a folder of them (`.sfdb` files are
//! never part of the repository).

use std::collections::{BTreeSet, HashMap};
use std::io::Write;
use std::path::{Path, PathBuf};

use openmassspec_core::conformance::assert_source_invariants;
use openmassspec_core::{write_indexed_mzml, SpectrumSource};
use opensf3::mgf::write_mgf;
use opensf3::model::{Match, Site, Spectrum};
use opensf3::mzid::write_mzid;
use opensf3::{Error, Format, Project, Reader, SfdbError};
use rusqlite::Connection;

const SCHEMA: &str = "
CREATE TABLE Version (major INTEGER, version_name TEXT);
CREATE TABLE AnalysisSoftware (ID INTEGER, SOFTWARE_NAME TEXT, SOFTWARE_VERSION TEXT);
CREATE TABLE IdentificationProtocol (ID INTEGER, ANALYSIS_SOFTWARE_ID INTEGER);
CREATE TABLE ScoreIndex (SCORE_ACCESSION TEXT, SCORE_DECREASING INTEGER, SCORE_INDEX INTEGER,
    SCORE_NAME TEXT, SCORE_TYPE TEXT);
CREATE TABLE Modification (MODIFICATION_REF TEXT, AMINO_ACID TEXT, MASS REAL,
    MODIFICATION_ACCESSION TEXT, SHORT_NAME TEXT, FULL_NAME TEXT);
CREATE TABLE MSSample (ID INTEGER, SAMPLE_NAME TEXT);
CREATE TABLE Spectrum (MSSAMPLE_ID INTEGER, SPECTRUM_ID TEXT, SPECTRUM_TITLE TEXT,
    ACQUISITION_MASS REAL, ACQUISITION_CHARGE INTEGER, SCAN_START_TIME_IN_SEC REAL,
    MASS_ARRAY TEXT, INTENSITY_ARRAY TEXT);
CREATE TABLE SpectrumData (DATA BLOB, FORMAT TEXT, SPECTRUM_ID TEXT, MSSAMPLE_ID INTEGER);
CREATE TABLE Peptide (PEPTIDE_REF TEXT, PEPTIDE_SEQUENCE TEXT, MODIFICATION_ID_ARRAY TEXT,
    PRIMARY_LOCATION_ARRAY TEXT);
CREATE TABLE Substitution (PEPTIDE_REF TEXT);
CREATE TABLE PeptideEvidence (PEPTIDE_REF TEXT, PROTEIN_ACCESSION TEXT,
    PROTEIN_START_INDEX INTEGER, PROTEIN_STOP_INDEX INTEGER, PRE_AA TEXT, POST_AA TEXT);
CREATE TABLE Protein (ACCESSION_NUMBER TEXT, ORIGINAL_PROTEIN_NAME TEXT,
    PROTEIN_SEQUENCE TEXT, IS_DECOY INTEGER);
CREATE TABLE PeptideSpectrumMatch (ID INTEGER, MSSAMPLE_ID INTEGER, SPECTRUM_ID TEXT,
    PEPTIDE_REF TEXT, ASSIGNED_CHARGE INTEGER, PEPTIDE_SCORE_ARRAY TEXT, PROTOCOL_ID INTEGER);
";

// Masses by hand (monoisotopic residues, water 18.010565, proton 1.007276):
// PEPTIDEK = P 2x97.05276 + E 2x129.04259 + T 101.04768 + I 113.08406
//   + D 115.02694 + K 128.09496 = 909.44434; + water = 927.454905;
//   2+ m/z (927.454905 + 2.014552) / 2 = 464.734729. EDITPEPK has the
//   same residues, so the same mass.
// AM+16C+57K = A 71.03711 + M 131.04049 + C 103.00919 + K 128.09496
//   = 433.18175; + water = 451.192315; + 15.994915 + 57.021464
//   = 524.208694; 2+ m/z (524.208694 + 2.014552) / 2 = 263.111623.
const SCAFFOLD_5: &str = r#"
INSERT INTO Version VALUES (1, 'Scaffold_5.1.2');
INSERT INTO AnalysisSoftware VALUES (200, 'Mascot', '2.7.0'), (201, 'Scaffold', 'Scaffold_5.1.2');
INSERT INTO IdentificationProtocol VALUES (200, 201);
INSERT INTO ScoreIndex VALUES ('MS:1001568', 0, 0, 'Scaffold:Peptide Probability', 'Peptide Score'),
    ('MS:1001171', 0, 1, 'Mascot:score', 'Peptide Score'),
    ('MS:1001371', 0, 2, 'Mascot:identity threshold', 'Peptide Score'),
    ('MS:1001579', 0, 0, 'Scaffold:Protein Probability', 'Protein Score');
INSERT INTO Modification VALUES
    ('UNIMOD:35_M', 'M', 15.9949150085449, 'UNIMOD:35', 'Oxidation', 'Oxidation'),
    ('UNIMOD:4_C', 'C', 57.0214538574219, 'UNIMOD:4', 'Carbamidomethyl', 'Carbamidomethyl');
INSERT INTO MSSample VALUES (200, 'sampleA');
INSERT INTO Spectrum VALUES
    (200, 'index=1', 'File: "C:\raw\run1.raw"; SpectrumID: "1"; scans: "101"', 464.734729, 2,
        600.5, '300.5,200.25,', '10.0,20.0,'),
    (200, 'index=2', 'File: "C:\raw\run1.raw"; SpectrumID: "2"; scans: "102"', 263.111623, 2,
        610.0, NULL, NULL);
INSERT INTO Peptide VALUES ('PEPTIDEK', 'PEPTIDEK', '', ''), ('EDITPEPK', 'EDITPEPK', '', ''),
    ('AM+16C+57K', 'AMCK', 'UNIMOD:35_M,UNIMOD:4_C', '2,3');
INSERT INTO PeptideEvidence VALUES ('PEPTIDEK', 'P1', 5, 12, 'K', 'A'),
    ('EDITPEPK', 'REV_P1', -1, 20, 'R', 'G'), ('AM+16C+57K', 'P2', 1, 4, '-', 'L');
INSERT INTO Protein VALUES ('P1', 'Protein one', 'MKKKPEPTIDEKAAA', 0),
    ('REV_P1', 'REV_P1', NULL, 0), ('P2', 'Protein two', 'XXXX', 0);
INSERT INTO PeptideSpectrumMatch VALUES
    (1, 200, 'index=1', 'EDITPEPK', 2, 'NaN,20.0,30.0', 200),
    (2, 200, 'index=1', 'PEPTIDEK', 2, 'NaN,40.0,30.0', 200),
    (3, 200, 'index=2', 'AM+16C+57K', 2, 'NaN,30.0,25.0', 200);
"#;

// Scaffold DDA: TMT6plex 229.162932 on the N-terminus and K.
// <229PEPTIDEK+229 = 927.454905 + 2 x 229.162932 = 1385.780769;
//   2+ m/z (1385.780769 + 2.014552) / 2 = 693.897661.
// <42PEPTIDEK+229 = 927.454905 + 42.010565 + 229.162932 = 1198.628402;
//   2+ m/z (1198.628402 + 2.014552) / 2 = 600.321477.
// <229MSYAPR with the M lost (as most of the real file's PSMs of that row),
//   so the mass of SYAPR + TMT:
//   S 87.03203 + Y 163.06333 + A 71.03711 + P 97.05276 + R 156.10111
//   = 574.28634; + water = 592.296905; + TMT = 821.459837; 2+ m/z
//   (821.459837 + 2.014552) / 2 = 411.737195. Its row lists
//   Met-loss+Acetyl, which doesn't fit; Met-loss alone does.
const SCAFFOLD_DDA: &str = "
INSERT INTO Version VALUES (2, '6.5.0');
INSERT INTO AnalysisSoftware VALUES (200, 'Sequest HT', '2.0.0.24, x64'),
    (201, 'Scaffold DDA', '6.5.0');
INSERT INTO IdentificationProtocol VALUES (200, 200);
INSERT INTO ScoreIndex VALUES ('MS:1001155', 0, 2, 'Sequest:xcorr', 'Peptide Score'),
    ('MS:1001493', 1, 0, 'percolator:PEP', 'Peptide Score'),
    ('MS:1002249', 0, 1, 'SEQUEST:sprank', 'Peptide Score');
INSERT INTO Modification VALUES
    ('UNIMOD:737_nY', 'n', 229.162932, 'UNIMOD:737', 'TMT6plex', 'TMT6plex'),
    ('UNIMOD:1_nY', 'n', 42.010565, 'UNIMOD:1047', 'Acetyl', 'Acetyl'),
    ('UNIMOD:737_K', 'K', 229.162932, 'UNIMOD:737', 'TMT6plex', 'TMT6plex'),
    ('UNIMOD:765_nM', 'n', -131.040485, 'UNIMOD:765', 'Met-loss', 'Met-loss'),
    ('UNIMOD:766_nM', 'n', -89.02992, 'UNIMOD:766', 'Met-loss+Acetyl', 'Met-loss+Acetyl');
INSERT INTO MSSample VALUES (200, 'Elute-1-1');
INSERT INTO Spectrum VALUES (200, '1', '1-7_7', 693.897661, 2, 2.2, NULL, NULL),
    (200, '2', '2-9_9', 600.321477, 2, 3.7, NULL, NULL),
    (200, '3', '3-11_11', 411.737195, 2, 4.0, NULL, NULL);
INSERT INTO Peptide VALUES
    ('<229PEPTIDEK+229', 'PEPTIDEK', 'UNIMOD:737_nY,UNIMOD:1_nY,UNIMOD:737_K', '0,0,8'),
    ('<42PEPTIDEK+229', 'PEPTIDEK', 'UNIMOD:1_nY,UNIMOD:737_K', '0,8'),
    ('<229MSYAPR', 'MSYAPR', 'UNIMOD:737_nY,UNIMOD:766_nM', '0,0');
INSERT INTO PeptideEvidence VALUES ('<229PEPTIDEK+229', 'P1', 5, 12, 'K', 'A'),
    ('<42PEPTIDEK+229', 'REVERSE_P9', 3, 10, 'R', '-'), ('<229MSYAPR', 'P3', 1, 6, '-', 'L');
INSERT INTO Protein VALUES ('P1', 'Protein one', NULL, 0), ('REVERSE_P9', 'decoy', NULL, 1),
    ('P3', 'Protein three', 'MSYAPRLL', 0);
INSERT INTO PeptideSpectrumMatch VALUES
    (200, 200, '1', '<229PEPTIDEK+229', 2, '0.01,1.0,3.5', 200),
    (201, 200, '2', '<42PEPTIDEK+229', 2, '0.9,1.0,0.8', 200),
    (202, 200, '3', '<229MSYAPR', 2, '0.02,1.0,2.1', 200);
";

fn write(dir: &Path, name: &str, rows: &str) -> PathBuf {
    let path = dir.join(name);
    let conn = Connection::open(&path).unwrap();
    conn.execute_batch(SCHEMA).unwrap();
    conn.execute_batch(rows).unwrap();
    path
}

fn dda(dir: &Path) -> PathBuf {
    let path = write(dir, "dda.sfdb", SCAFFOLD_DDA);
    let conn = Connection::open(&path).unwrap();
    conn.execute(
        "INSERT INTO SpectrumData VALUES (?1, 'gzip', '1', 200)",
        [blob(
            693.897661,
            &[(127.12476, 5000.0), (126.127725, 4000.0)],
        )],
    )
    .unwrap();
    path
}

fn by_scan(project: &Project, scan: u32) -> &Spectrum {
    project
        .spectra
        .iter()
        .find(|spectrum| spectrum.scan == Some(scan))
        .unwrap()
}

fn ranked(spectrum: &Spectrum, rank: u32) -> &Match {
    spectrum
        .matches
        .iter()
        .find(|found| found.rank == rank)
        .unwrap()
}

#[test]
fn scaffold_5_files_read_with_engine_scores_and_peaks() {
    let dir = tempfile::tempdir().unwrap();
    let reader = Reader::open(write(dir.path(), "s5.sfdb", SCAFFOLD_5)).unwrap();
    assert_eq!(reader.format(), Format::Sfdb);
    let project = reader.project();
    assert_eq!(project.scaffold_version.as_deref(), Some("Scaffold_5.1.2"));
    assert_eq!(project.engines.len(), 1);
    assert_eq!(project.engines[0].name, "Mascot");
    assert_eq!(project.engines[0].version.as_deref(), Some("2.7.0"));
    assert_eq!(project.runs.len(), 1);
    assert_eq!(project.runs[0].name, "run1");
    assert_eq!(project.runs[0].biosample, "sampleA");
    let scans: BTreeSet<Option<u32>> = project.spectra.iter().map(|s| s.scan).collect();
    assert_eq!(scans, BTreeSet::from([Some(101), Some(102)]));
    // Scaffold's verdicts are not read from an .sfdb.
    assert!(project.groups.is_empty() && project.thresholds.is_none());
    assert!(project
        .spectra
        .iter()
        .flat_map(|s| &s.matches)
        .all(|found| found.probability.is_none() && !found.chosen));

    // Scan 101: PEPTIDEK (40) outranks EDITPEPK (20).
    let first = by_scan(project, 101);
    let best = ranked(first, 1);
    assert_eq!(best.peptide.sequence, "PEPTIDEK");
    assert_eq!(ranked(first, 2).peptide.sequence, "EDITPEPK");
    let names: Vec<&str> = best.scores.iter().map(|s| s.name.as_ref()).collect();
    assert_eq!(
        names,
        [
            "Mascot:score",
            "Mascot:identity threshold",
            "Mascot:expectation value"
        ]
    );
    // 0.05 x 10^((30 - 40) / 10) = 0.005.
    assert!((best.scores[2].value - 0.005).abs() < 1e-12);
    assert_eq!(best.proteins[0].accession, "P1");
    assert_eq!(
        (best.proteins[0].start, best.proteins[0].end),
        (Some(5), Some(12))
    );

    // The decoy follows the accession convention; its start of -1 is
    // rebuilt from its stop: 20 - 8 + 1 = 13.
    let decoy = ranked(first, 2);
    assert_eq!(decoy.proteins[0].start, Some(13));
    assert!(project.proteins["REV_P1"].decoy);
    assert!(!project.proteins["P1"].decoy);

    let amck = ranked(by_scan(project, 102), 1);
    let modifications: Vec<(Site, Option<&str>)> = amck
        .peptide
        .modifications
        .iter()
        .map(|m| (m.site, m.name.as_deref()))
        .collect();
    assert_eq!(
        modifications,
        [
            (Site::Residue(2), Some("Oxidation")),
            (Site::Residue(3), Some("Carbamidomethyl"))
        ]
    );
    assert_eq!(amck.proteins[0].pre, Some('-'), "the protein's terminus");
    // P2's stored sequence doesn't hold AMCK, so it is not kept.
    assert_eq!(project.proteins["P2"].sequence, None);
    assert!(project
        .notes
        .iter()
        .any(|note| note.contains("1 stored protein")));
    assert!(!project
        .notes
        .iter()
        .any(|note| note.contains("no decoy proteins")));

    // Peaks: scan 101's, sorted by m/z; scan 102 stores none.
    assert_eq!(
        project.peaks(first).unwrap(),
        Some(vec![(200.25, 20.0), (300.5, 10.0)])
    );
    assert_eq!(project.peaks(by_scan(project, 102)).unwrap(), None);
}

#[test]
fn scaffold_dda_files_keep_the_named_terminal_modification() {
    let dir = tempfile::tempdir().unwrap();
    let reader = Reader::open(dda(dir.path())).unwrap();
    let project = reader.project();
    assert_eq!(project.engines[0].name, "Sequest HT");
    assert_eq!(project.runs[0].name, "Elute-1-1");
    assert_eq!(project.runs[0].biosample, "Elute-1-1");
    assert!(
        project.notes[0].contains("percolator:PEP"),
        "{:?}",
        project.notes
    );

    let tmt = ranked(by_scan(project, 7), 1);
    // The name says <229: TMT on the N-terminus, not acetyl.
    let sites: Vec<(Site, &str)> = tmt
        .peptide
        .modifications
        .iter()
        .map(|m| (m.site, m.name.as_deref().unwrap()))
        .collect();
    assert_eq!(
        sites,
        [(Site::NTerm, "TMT6plex"), (Site::Residue(8), "TMT6plex")]
    );
    let names: Vec<&str> = tmt.scores.iter().map(|s| s.name.as_ref()).collect();
    assert_eq!(names, ["SEQUEST:sprank", "Sequest:xcorr"]);
    assert_eq!(tmt.charge, Some(2));

    // Scaffold files the acetyl under UNIMOD:1047; it is Unimod's Acetyl.
    let acetyl = ranked(by_scan(project, 9), 1);
    assert_eq!(
        acetyl.peptide.modifications[0].name.as_deref(),
        Some("Acetyl")
    );
    assert!(project.proteins["REVERSE_P9"].decoy);

    // The M was lost: kept as the engine gives it, an N-terminal Met-loss
    // delta beside the TMT, on the protein where the row places it.
    let lost = ranked(by_scan(project, 11), 1);
    assert_eq!(lost.peptide.sequence, "MSYAPR");
    let masses: Vec<f64> = lost
        .peptide
        .modifications
        .iter()
        .map(|m| m.mass_delta)
        .collect();
    assert_eq!(masses, [-131.040485, 229.162932]);
    assert_eq!(lost.proteins[0].start, Some(1));
    assert!(
        project
            .notes
            .iter()
            .any(|note| note.contains("1 of them with the initial methionine")),
        "{:?}",
        project.notes
    );

    assert_eq!(
        project.peaks(by_scan(project, 7)).unwrap(),
        Some(vec![(126.127725, 4000.0), (127.12476, 5000.0)])
    );
    assert_eq!(project.peaks(by_scan(project, 9)).unwrap(), None);
}

#[test]
fn sfdb_outputs_carry_engine_scores_and_no_scaffold_verdicts() {
    let dir = tempfile::tempdir().unwrap();
    let reader = Reader::open(write(dir.path(), "s5.sfdb", SCAFFOLD_5)).unwrap();
    let project = reader.project();

    let mut mgf = Vec::new();
    let report = write_mgf(project, 0, &mut mgf).unwrap();
    assert_eq!(report.spectra, 2);
    let mgf = String::from_utf8(mgf).unwrap();
    assert!(mgf.contains("PEPMASS=464.734729\nCHARGE=2+\nSCANS=101\nRTINSECONDS=600.500\n"));
    assert!(mgf.contains("200.25 20\n300.5 10\n"), "{mgf}");

    let mut xml = Vec::new();
    let report = write_mzid(project, &mut xml).unwrap();
    assert_eq!(report.candidates, 3);
    let xml = String::from_utf8(xml).unwrap();
    for wanted in [
        r#"name="no threshold""#,
        r#"accession="MS:1001171" name="Mascot:score" value="40""#,
        r#"name="peak list scans" value="101""#,
        r#"start="13" end="20""#,
        r#"passThreshold="true""#,
    ] {
        assert!(xml.contains(wanted), "{wanted}");
    }
    for unwanted in [
        "ProteinDetection",
        "Scaffold:Peptide Probability",
        "Scaffold:Protein Probability",
        r#"passThreshold="false""#,
    ] {
        assert!(!xml.contains(unwanted), "{unwanted}");
    }
}

#[test]
fn run_spectra_meet_the_openmassspec_invariants() {
    let dir = tempfile::tempdir().unwrap();
    let dda_path = dda(dir.path());
    for path in [write(dir.path(), "s5.sfdb", SCAFFOLD_5), dda_path.clone()] {
        let reader = Reader::open(&path).unwrap();
        let mut total = 0;
        for run in 0..reader.runs().len() {
            let mut source = reader.run_spectra(run).unwrap();
            let count = assert_source_invariants(&mut source).unwrap();
            assert_eq!(source.spectrum_count_hint(), Some(count));
            total += count;
        }
        assert_eq!(total, reader.project().spectra.len());
    }

    let reader = Reader::open(&dda_path).unwrap();
    let mut source = reader.run_spectra(0).unwrap();
    let records: Vec<_> = source.iter_spectra().map(Result::unwrap).collect();
    assert_eq!(records[0].native_id, "index=0");
    assert_eq!(records[0].scan_number, 7);
    assert_eq!(records[0].ms_level, 2);
    let precursor = records[0].precursor.as_ref().unwrap();
    assert_eq!(
        (precursor.selected_mz, precursor.charge),
        (Some(693.897661), Some(2))
    );
    assert_eq!(records[0].mz, [126.127725, 127.12476]);
    assert_eq!(records[0].extra["opensf3.title"], "1-7_7");
    assert!(records[1].mz.is_empty());

    let mut mzml = Vec::new();
    let summary = write_indexed_mzml(&mut source, &mut mzml).unwrap();
    assert_eq!(summary.spectra_written, 3);
    let mzml = String::from_utf8(mzml).unwrap();
    assert!(mzml.contains(r#"<spectrum id="index=2" index="2""#));
    assert!(mzml.contains("multiple peak list nativeID format"));
    assert!(reader.run_spectra(1).is_err());
}

#[test]
fn misread_masses_refuse_the_file() {
    let dir = tempfile::tempdir().unwrap();
    // Both precursors 6 Da heavy (past the +3 isotope peak): no PSM fits.
    let rows = SCAFFOLD_5
        .replace("464.734729", "467.734729")
        .replace("263.111623", "266.111623");
    match Reader::open(write(dir.path(), "bad.sfdb", &rows)) {
        Err(Error::Sfdb(SfdbError::Layout(message))) => {
            assert!(message.contains("3 of 3 PSMs miss"), "{message}")
        }
        other => panic!("expected a refusal, got {other:?}"),
    }
}

#[test]
fn other_files_are_not_scaffold_files() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("other.sqlite");
    Connection::open(&path)
        .unwrap()
        .execute_batch("CREATE TABLE Something (x INTEGER);")
        .unwrap();
    assert!(matches!(Reader::open(&path), Err(Error::Unrecognized(_))));
    assert!(matches!(
        opensf3::sfdb::decode(&path),
        Err(SfdbError::NotSfdb(_))
    ));
    let text = dir.path().join("notes.sfdb");
    std::fs::write(&text, "not a database at all").unwrap();
    assert!(matches!(Reader::open(&text), Err(Error::Unrecognized(_))));
    let core: openmassspec_core::Error = Reader::open(&text).unwrap_err().into();
    assert!(matches!(core, openmassspec_core::Error::Format(_)));
}

/// Checks on real `.sfdb` files, run when `SFDB_SAMPLES` names a folder.
#[test]
fn real_files_read_consistently() {
    let Some(folder) = std::env::var_os("SFDB_SAMPLES") else {
        eprintln!("SFDB_SAMPLES not set; skipping");
        return;
    };
    let mut files: Vec<PathBuf> = std::fs::read_dir(folder)
        .expect("SFDB_SAMPLES is a folder")
        .map(|entry| entry.expect("readable").path())
        .filter(|path| path.extension().is_some_and(|ext| ext == "sfdb"))
        .collect();
    files.sort();
    for path in files {
        let name = path.display();
        let reader = Reader::open(&path).unwrap_or_else(|error| panic!("{name}: {error}"));
        let project = reader.project();
        assert_eq!(project.format, Format::Sfdb);
        let mut psms = 0;
        for spectrum in &project.spectra {
            let mut last: HashMap<usize, u32> = HashMap::new();
            for found in &spectrum.matches {
                psms += 1;
                let previous = last.insert(found.engine, found.rank).unwrap_or(0);
                assert_eq!(found.rank, previous + 1, "{name}: ranks run 1, 2, ...");
            }
        }
        let mut spectra = 0;
        for run in 0..project.runs.len() {
            let mut source = reader.run_spectra(run).unwrap();
            spectra += assert_source_invariants(&mut source)
                .unwrap_or_else(|error| panic!("{name}: {error:?}"));
        }
        assert_eq!(spectra, project.spectra.len(), "{name}");
        let report = write_mzid(project, &mut std::io::sink()).expect("writes");
        assert_eq!(report.candidates, psms, "{name}");
        eprintln!(
            "{name}: {} runs, {} spectra, {psms} PSMs, {} proteins",
            project.runs.len(),
            project.spectra.len(),
            project.proteins.len()
        );
    }
}

// A `SpectrumData` blob as Scaffold DDA 6.5 writes one: {1: "gzip",
// 2: gzip({1: {2: {1: scan, 2: title}}})}, the scan {1: precursor m/z,
// 3: intensity, 4: {1: peak...}, 5: 1}, each peak {1: m/z, 2: time,
// 3: intensity}.
fn blob(precursor: f64, peaks: &[(f64, f32)]) -> Vec<u8> {
    fn varint(mut value: u64) -> Vec<u8> {
        let mut out = Vec::new();
        loop {
            let byte = (value & 0x7f) as u8;
            value >>= 7;
            if value == 0 {
                out.push(byte);
                return out;
            }
            out.push(byte | 0x80);
        }
    }
    fn message(number: u64, body: &[u8]) -> Vec<u8> {
        let mut out = varint((number << 3) | 2);
        out.extend(varint(body.len() as u64));
        out.extend_from_slice(body);
        out
    }
    let mut list = Vec::new();
    for &(mz, intensity) in peaks {
        let mut peak = varint((1 << 3) | 1);
        peak.extend(mz.to_le_bytes());
        peak.extend(varint((2 << 3) | 5));
        peak.extend(2.2f32.to_le_bytes());
        peak.extend(varint((3 << 3) | 5));
        peak.extend(intensity.to_le_bytes());
        list.extend(message(1, &peak));
    }
    let mut scan = varint((1 << 3) | 1);
    scan.extend(precursor.to_le_bytes());
    scan.extend(varint((3 << 3) | 5));
    scan.extend(f32::NAN.to_le_bytes());
    scan.extend(message(4, &list));
    scan.extend(varint(5 << 3));
    scan.push(1);
    let mut group = message(1, &scan);
    group.extend(message(2, b"-7_7"));
    let top = message(1, &message(2, &group));
    let mut gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    gz.write_all(&top).unwrap();
    let mut outer = message(1, b"gzip");
    outer.extend(message(2, &gz.finish().unwrap()));
    outer
}
