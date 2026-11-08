# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.2.0] - 2026-11-08

### Added

- Scaffold `.sfdb` files (Scaffold 5.1.2 and Scaffold DDA 6.5): runs,
  spectra, stored peaks (text arrays or gzipped protobuf blobs), every
  engine's PSMs with the engine's own scores, N-terminal readings with
  Met-loss, protein evidence with start and end, decoys. Scaffold's
  probabilities, protein groups and thresholds are not read from `.sfdb`.
- `Reader::open`, which tells `.sf3` from `.sfdb` by content, and
  `Reader::run_spectra`, each run's stored peak lists as an
  `openmassspec_core::SpectrumSource`.
- `sf3decode ... mzml`: one indexed mzML 1.1 per run.
- Scan numbers from spectrum titles: `SCANS=` in MGF, `peak list scans`
  in mzIdentML, `scan_number` and `opensf3.scan` in mzML records.
- mzIdentML `PeptideEvidence` start and end when the file stores them.
- `opensf3::Error`, convertible into `openmassspec_core::Error`.
- Format notes in `docs/format/`, `CORPUS.md`, `ATTRIBUTION.md`, CI.

### Changed

- Cargo workspace: library crate `opensf3`, CLI crate `opensf3-cli`
  (binary still `sf3decode`). Edition 2021, MSRV 1.88.
- The mzIdentML writer is behind the default `mzid` feature.
- Public model types are `#[non_exhaustive]`; `Score` names and
  accessions are `Cow<'static, str>`; `Match` has the engine's `charge`;
  `ProteinRef` has `start` and `end`; `Project` has its `format`.
- MGF writes `PEPMASS` without `CHARGE` when a file stores an m/z but no
  charge (`.sfdb`); `.sf3` output is unchanged apart from `SCANS=`.

## [0.1.0] - 2026-10-12

### Added

- `sf3decode`: summary, MGF and mzIdentML 1.2 from Scaffold `.sf3`
  projects (Scaffold 4.0 to 5.3.4), with Scaffold's probabilities,
  protein groups and saved thresholds.
