# Credits

## Provenance

Both formats were worked out from public data files only (see
[CORPUS.md](CORPUS.md)). Scaffold was never run, decompiled or
disassembled, none of its code was read, and none is in this repository.
Nothing required accepting Scaffold's license.

## Code

The `.sfdb` reader (`crates/opensf3/src/sfdb/`: tables, engines and
scores, ranking, N-terminus readings, peak decoders), the Unimod-by-mass
naming, and the scan numbers taken from spectrum titles were ported from
[Tuffold](https://github.com/TUFF-Industries/tuffold) (`tuffold-io`), by
the same author, under the same Apache-2.0 license. The `.sf3` decoder is
this repository's own, and Tuffold's `.sf3` reader derives from the same
work.

## Prior art

`keesh0/proteo` (https://github.com/keesh0/proteo), a public extractor for
Scaffold DDA `.sfdb` files, has no license. It was read only for facts
about Scaffold DDA 7.1 files (peptide start indices written as -1,
`ProteinAssembly.PASSES_THRESHOLD` always 0), recorded in Tuffold's notes;
no code was taken from it, and those facts are not verified here (no 7.1
file has been seen).

## Standards

- mzIdentML 1.2 (HUPO-PSI): the mzIdentML writer.
- mzML 1.1 (HUPO-PSI): written by
  [openmassspec-core](https://crates.io/crates/openmassspec-core)'s
  indexed mzML writer.
- PSI-MS controlled vocabulary (psi-ms.obo): score, software and format
  terms.
- Unimod (unimod.org): modification names and accessions.
- Java Object Serialization Specification (protocol 2): `.sf3` inner
  streams.
- Protocol Buffers wire format: Scaffold DDA peak blobs.
- SQLite file format: both containers.
- Mascot Generic Format (MGF): the MGF writer.

## Validation corpus

Files were downloaded from the PRIDE Archive and other public sources
listed in [CORPUS.md](CORPUS.md), and are used read-only for checks.

## Dependencies

openmassspec-core (Apache-2.0), rusqlite with bundled SQLite (MIT; SQLite
is public domain), flate2 (MIT/Apache-2.0), thiserror (MIT/Apache-2.0),
clap (MIT/Apache-2.0).
