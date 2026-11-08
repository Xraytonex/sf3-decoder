# opensf3 (sf3-decoder)

[![License: Apache-2.0](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](LICENSE)
[![Rust MSRV](https://img.shields.io/badge/rust-1.88%2B-orange.svg)](https://www.rust-lang.org)

> Built to fit the [OpenMassSpec](https://github.com/Sigilweaver/OpenMassSpec)
> stack: its spectra are an `openmassspec-core` `SpectrumSource`.

Rust reader for Proteome Software's Scaffold files, with no Scaffold
software or SDK: `.sf3` projects (Scaffold 4.0 to 5.3.4) and `.sfdb`
databases (Scaffold 5.1, Scaffold DDA 6.5). Only Scaffold opens these
files. opensf3 reads one and gives back what's in it: the MS runs, every
stored MS/MS peak list, every search engine's candidate peptides with
their scores, the proteins, and for `.sf3` Scaffold's own peptide and
protein probabilities, protein groups and saved thresholds.

The workspace has two crates:

- `opensf3` (library): `Reader::open` tells `.sf3` from `.sfdb` by content
  and decodes either into a `Project`; each run's peak lists are an
  `openmassspec_core::SpectrumSource`, so the stack's mzML writer can
  write them; MGF and mzIdentML 1.2 writers are included.
- `opensf3-cli`: the `sf3decode` command (summary, `mgf`, `mzid`, `mzml`).

Scaffold Viewer can export mzIdentML and MGF one file at a time through
its interface; `sf3decode` does it from the command line, for scripts and
pipelines, and includes every engine's candidates and scores.

## Known issues and gaps

Checked on the files in [CORPUS.md](CORPUS.md).

- **Scaffold versions.** `.sf3` from Scaffold 4.0 through 5.3.4 and
  `.sfdb` from Scaffold 5.1.2 and Scaffold DDA 6.5 have been checked.
  Anything else (Scaffold DDA 7.x `.sfdb`, Scaffold DDA `.sf3` if it
  exists) is untested and may be refused. A file whose structure doesn't
  hold up is refused, never half-read.
- **`.sfdb` carries no Scaffold verdicts here.** Scaffold's peptide and
  protein probabilities, protein groups and thresholds are not read from
  `.sfdb` files, because what those columns mean in this layout isn't
  settled. An `.sfdb` gives engine candidates and scores only: its
  mzIdentML has every candidate under a "no threshold" protocol (so
  `passThreshold` is true for all, as the specification asks when no
  threshold was applied) and no protein detection. Post-processor scores
  stored in the file (Percolator's PEP, PeptideProphet) are left out too.
- **Engines.** In `.sf3`, only Mascot, Sequest, X! Tandem, MaxQuant,
  Spectrum Mill and MSFragger are converted; other engines' candidates are
  left out, since where their scores sit isn't known yet. X! Tandem's
  scores are named by offset, because Scaffold's own names for them don't
  fit the values. In `.sfdb`, scores are named as the file names them;
  in mzIdentML a score is a PSI-MS cvParam only when the writer knows that
  term's exact name, else a userParam.
- **Ranks.** Candidates are ranked by each engine's own main score. In
  `.sf3` equal scores share a rank; in `.sfdb` equal scores are broken by
  the file's PSM id (a target and a decoy with the same score then rank
  in file order).
- **Peaks.** Scaffold keeps at most 300 peaks per spectrum, so MGF and
  mzML hold Scaffold's peak lists, not the raw spectra. Every spectrum is
  written as a centroided MS2 spectrum.
- **mzML.** Scaffold files don't say which instrument recorded a run, so
  the instrument is the generic "instrument model" term and the source
  file is the Scaffold file. A spectrum without a stored retention time
  gets 0 (the writer requires one) with `opensf3.no_retention_time` set.
  An `.sf3` spectrum without a charge has no m/z, so its precursor is
  empty; such spectra fail `openmassspec-core`'s precursor invariant.
- **Spectra without a charge** (in `.sf3`) go into the MGF without
  `PEPMASS`, and their candidates are left out of the mzIdentML.
- **Scan numbers** come only from titles that name exactly one scan; two
  spectra of one run naming the same scan both lose it.
- **Protein pass flags** in `.sf3` mzIdentML are rebuilt from the saved
  rule and may differ slightly from Scaffold's own list. Protein groups
  with no peptide Scaffold picked, and protein sequences with a `*` inside,
  can't be written in valid mzIdentML and are left out. The searched
  database itself isn't in the file, only the proteins it matched.
- Some Scaffold 4.11 demo files store protein sequences that don't contain
  their own peptides; those sequences are dropped.
- Quantitation (TMT/iTRAQ reporter intensities, SILAC channels) and GO
  annotations are not converted.
- Large projects need memory: 2.3 GB for a 500 MB `.sf3`, 1.1 GB for the
  1 GB PXD071927 `.sfdb`. An mzIdentML with every candidate of every rank
  can run to hundreds of MB (675 MB for PXD071927).

## Install

With Rust 1.88 or newer, the command-line tool:

```
cargo install --git https://github.com/Xraytonex/sf3-decoder opensf3-cli
```

The library, in `Cargo.toml`:

```toml
[dependencies]
opensf3 = { git = "https://github.com/Xraytonex/sf3-decoder" }
openmassspec-core = "3.0.0"
```

The mzIdentML writer is behind the default `mzid` feature.

## Quickstart

### Command line

```
sf3decode project.sf3
```

prints what the file holds: runs, spectra, each search engine's
candidates, proteins, and for `.sf3` the protein groups and saved
thresholds, plus anything the decoder had to leave out. `.sfdb` files work
the same way.

```
sf3decode project.sf3 mgf
sf3decode project.sf3 mzml
sf3decode project.sf3 mzid
```

`mgf` writes one MGF per MS run, named after the run's raw file, with
every spectrum's title, precursor m/z and charge, scan (`SCANS=`, when
the title names one), retention time and peaks. `mzml` writes one indexed
mzML 1.1 per run with the same spectra, named the same way. `mzid` writes
`project.mzid` (mzIdentML 1.2): every candidate each search engine kept
for each spectrum, ranked by that engine's score, with its scores; for
`.sf3` also Scaffold's peptide probability on each, `passThreshold` as
Scaffold's call (the candidate it picked, at or above the saved peptide
probability), and Scaffold's protein groups with their probability,
overall and per run. Spectrum `N` of a run (from 0) is `index=N` in its
MGF, its mzML, and the mzIdentML, so write them into the same folder.
Files go next to the Scaffold file unless `-o folder` says otherwise, and
existing files are only replaced with `--force`.

### Library

```rust
use openmassspec_core::write_indexed_mzml;
use opensf3::Reader;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let reader = Reader::open("project.sfdb")?;
    let project = reader.project();
    println!("{:?}: {} runs, {} spectra", reader.format(), project.runs.len(), project.spectra.len());
    for spectrum in project.spectra.iter().take(3) {
        for found in &spectrum.matches {
            println!("{} rank {}: {}", spectrum.title.as_deref().unwrap_or("?"), found.rank, found.peptide.key());
        }
    }
    for run in 0..reader.runs().len() {
        let mut source = reader.run_spectra(run)?;
        let mut out = std::fs::File::create(format!("run{run}.mzML"))?;
        write_indexed_mzml(&mut source, &mut out)?;
    }
    Ok(())
}
```

`Reader::open` detects the format; `opensf3::sf3::decode` and
`opensf3::sfdb::decode` read one kind directly. `opensf3::mgf::write_mgf`
and `opensf3::mzid::write_mzid` write MGF and mzIdentML from a `Project`.
Errors are `opensf3::Error`, which converts into
`openmassspec_core::Error`.

## Tests

```
cargo test --all
```

runs everything that needs no data files: unit tests and hand-built
`.sfdb` databases in both layouts, including an `openmassspec-core`
conformance check. The checks on real files run when `SF3_SAMPLES` or
`SFDB_SAMPLES` names a folder of them (Scaffold files are never part of
the repository):

```
SF3_SAMPLES=~/sf3-files SFDB_SAMPLES=~/sfdb-files cargo test --release
```

## How this was made

The `.sf3` format was worked out from about twenty `.sf3` files that are
publicly available, mostly from PRIDE submissions, written by Scaffold 4.0
through 5.3.4. Scaffold itself was never run, decompiled or disassembled
for this, and none of its code is in this repository. It was pretty much
brute force: unpack the file, guess what a stretch of bytes means, and
check the guess against something known independently (peptide masses
against precursor masses, scores against the search engines' ranges,
titles against retention times) across every sample until it held
everywhere. That's how open readers for most closed data formats get
written. Nothing here required accepting Scaffold's license, getting
around any protection, or looking at its code; it's reading data files for
interoperability, which is lawful. The decoder only reads; it never writes
`.sf3`.

The `.sfdb` reader was worked out the same way, from two public `.sfdb`
files (PRIDE PXD038643 and PXD071927); an `.sfdb` is a plain SQLite
database, read read-only, and never written. Its code, and the scan
numbers taken from spectrum titles, were ported from
[Tuffold](https://github.com/TUFF-Industries/tuffold) (same author,
Apache-2.0). See [ATTRIBUTION.md](ATTRIBUTION.md) and the format notes in
[docs/format/](docs/format/): [`.sf3`](docs/format/01-sf3.md),
[`.sfdb`](docs/format/02-sfdb.md).

This project is not affiliated with or endorsed by Proteome Software.
Scaffold is their trademark.

## License

Apache 2.0.
