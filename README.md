# sf3-decoder

Binary to convert Proteome Software's .sf3 files into more commonly used formats.

Scaffold saves its projects as `.sf3` files, and only Scaffold can open them.
`sf3decode` reads one and writes out what's in it, including Scaffold's own
probabilities, so the results can be used elsewhere. Scaffold Viewer can
export mzIdentML and MGF one file at a time through its interface;
`sf3decode` does it from the command line, for scripts and pipelines, and
includes every engine's candidates and scores.

## How this was made

The format was worked out from about twenty `.sf3` files that are publicly
available, mostly from PRIDE submissions, written by Scaffold 4.0 through
5.3.4. Scaffold itself was never run, decompiled or disassembled for this,
and none of its code is in this repository. It was pretty much brute force:
unpack the file, guess what a stretch of bytes means, and check the guess
against something known independently (peptide masses against precursor
masses, scores against the search engines' ranges, titles against retention
times) across every sample until it held everywhere. That's how open readers
for most closed data formats get written. Nothing here required accepting
Scaffold's license, getting around any protection, or looking at its code;
it's reading data files for interoperability, which is lawful. The decoder
only reads; it never writes `.sf3`.

This project is not affiliated with or endorsed by Proteome Software.
Scaffold is their trademark.

## How to use

Install with Rust (1.88 or newer):

```
cargo install --git https://github.com/Xraytonex/sf3-decoder
```

Then:

```
sf3decode project.sf3
```

prints what the file holds: runs, spectra, each search engine's candidates,
proteins, protein groups, the thresholds it was saved with, and anything the
decoder had to leave out.

```
sf3decode project.sf3 mgf
```

writes one MGF per MS run, named after the run's raw file, with every
spectrum Scaffold kept: title, precursor m/z and charge, retention time and
peaks. Scaffold keeps at most 300 peaks per spectrum, so these are not the
full raw spectra. Files are written next to the `.sf3` unless `-o folder`
says otherwise, and existing files are only replaced with `--force`.

```
sf3decode project.sf3 mzid
```

writes `project.mzid` (mzIdentML 1.2, valid against the official schema).
It has every candidate each search engine kept for each spectrum, ranked by
that engine's score, with its scores, and Scaffold's peptide probability on
each one. `passThreshold` is Scaffold's call: the candidate it picked for the
spectrum, at or above the peptide probability the project was saved with.
Scaffold's protein groups come with their protein probability, overall and
per run. Their pass flag is rebuilt from the saved rule (protein probability
in a run, plus the minimum number of peptides Scaffold picked there), so it
may differ slightly from Scaffold's own list. Spectra are referenced by their
position in the MGF files above, so write both into the same folder.

The checks in `tests/samples.rs` run on real files when `SF3_SAMPLES` names
a folder of them:

```
SF3_SAMPLES=~/sf3-files cargo test --release
```

## The .sf3 format

An `.sf3` file is a SQLite database. Its `Files` table holds every part of
the project as a gzip-compressed blob, and its `ProteinGroup` table holds
Scaffold's protein probability for each protein group in each MS run. The
inner files come in three kinds.

**Java serialization streams** (starting `ac ed 00 05`). `Experiment` is the
project itself: biosamples, MS runs, protein groups (`ProteinExperimentGroup`,
with their accessions and decoy flag), the thresholds the project was saved
with (`DisplayThresholds`), and for every search engine Scaffold knows the
names of its scores. `DatabaseSearch` lists the search programs,
`ModificationTypes` the modifications, and `IdentificationCacheLookup` maps
every protein to its peptides. These are plain Java Object Serialization, so
they can be read from the public specification.

**Fixed-size records.** Most of the data is arrays of fixed-size blocks
split over numbered parts (`PeptideIdentification0`, `1`, …). Each type has
a descriptor stream giving `FILE_SHIFT`: record id `n` is block
`n & (2^FILE_SHIFT − 1)` of part `n >> FILE_SHIFT`, and every block starts
with its own id. Numbers are little-endian.

| record | size | what's in it |
|---|---|---|
| `Spectrum` | 48 | peak list id at 4, title id at 8, precursor [M+H]+ (f64) at 12, charge at 20, retention time in seconds (f32) at 32 |
| `SpectrumAnalysis` | 180 | one per spectrum and engine: spectrum id at 4, search program at 8 |
| `PeptideIdentification` | 160 | one candidate peptide: residues at 4–34 (the rest in `PeptideSequence`), flanking residues at 39 and 40, Scaffold's peptide probability (f32) at 49, its spectrum analysis at 57, engine scores as f32 from 66, three modifications at 90–114 (more in two int lists), the reported protein at 122, and at 147 a 0 if this is the candidate Scaffold picked |
| `ProteinIdentification` | 103 | the engine's protein accession at 8, cut at 22 characters |
| `PeakList` | 180 | 20 (intensity, m/z) pairs and the head of a chain holding the rest |
| `PeakListMassIntensities` | 232 | 28 more pairs per block, chained |

Strings (titles, long sequences) and peak lists are chains of blocks: each
block carries its chain's id and, in its last four bytes, the index of the
next block in the same part, −1 at the end.

**Block-based object files** (`BBProteinSequence`, `SpectrumCollation`) use
Scaffold's own variant of Java serialization: fields big-endian, lengths and
arrays little-endian, a class table in a descriptor stream. Protein sequences
and descriptions live here, and each run's `SpectrumCollation` lists its
spectra.

A few things that aren't obvious:

- Scores sit at `66 + 4k`, where `k` indexes Scaffold's own list of score
  names for that engine. The first entry is often Scaffold's discriminant
  score rather than the engine's. X! Tandem's list doesn't match its values,
  so its scores are named by offset.
- A spectrum's peaks are the 20 pairs in its `PeakList` record plus its whole
  chain. The 20 are Scaffold's picks and the chain holds the rest, isotope
  peaks included. Scaffold keeps at most 300 peaks per spectrum.
- Accessions stored with a candidate are cut at 22 characters, which turns a
  decoy like `tr|A0QSM5|A0QSM5_MYCS2-DECOY` into its target's name. They are
  completed from Scaffold's peptide → protein map, or by checking whether the
  peptide reads forwards or backwards in the stored target sequence.
- Older files keep a second, stale copy of many candidates; they are counted
  once.
- Spectrum Mill writes modified residues in lower case without saying which
  modification. The one modification the file uses on that residue elsewhere
  that fits the precursor mass is used.

## How the decoder works

`sf3decode` opens the SQLite file read-only, inflates the inner files it
needs, and parses the Java streams into plain values (nothing in the file is
ever executed). It reads the runs and their spectra from `Experiment` and
the collations, then every `PeptideIdentification`: sequence,
modifications, the engine's scores, Scaffold's probability and pick. Each
engine's candidates for a spectrum are ranked by that engine's own main
score. Every candidate's peptide mass is checked against its spectrum's
precursor mass (isotope errors allowed). One that misses by more than 1 Da
is left out, and if more than 5% miss, the file is refused as misread. The
same goes for anything structural that doesn't hold, like a block whose id
doesn't match its position or a chain that never ends: the file is refused
rather than half-read. Everything left out is counted and reported.

## Limitations

- Only files from Scaffold 4.0 through 5.3.4 have been checked. Scaffold DDA
  6.x files are untested and may be refused.
- Only Mascot, Sequest, X! Tandem, MaxQuant, Spectrum Mill and MSFragger are
  converted. Other engines' candidates are left out, since where their scores
  sit isn't known yet.
- X! Tandem's scores are named by offset, because Scaffold's own names for
  them don't fit the values. MaxQuant's delta score and MSFragger's scores
  have no PSI-MS term, so they are written as userParams.
- Scaffold keeps at most 300 peaks per spectrum, so the MGFs are not the raw
  spectra.
- Spectra Scaffold stored without a charge have no m/z: they go into the MGF
  without `PEPMASS`, and their candidates are left out of the mzIdentML.
- Protein pass flags are rebuilt from the saved rule and may differ slightly
  from Scaffold's own list.
- Protein groups with no peptide Scaffold picked, and protein sequences with
  a `*` inside, can't be written in valid mzIdentML and are left out. The
  searched database itself isn't in the file, only the proteins it matched.
- Some Scaffold 4.11 demo files store protein sequences that don't contain
  their own peptides; those sequences are dropped.
- Quantitation (TMT/iTRAQ reporter intensities, SILAC channels) and GO
  annotations are not converted.
- Large projects need a few GB of memory (2.3 GB for a 500 MB file), and an
  mzIdentML with every candidate of every rank can run to hundreds of MB.

## License

Apache 2.0.
