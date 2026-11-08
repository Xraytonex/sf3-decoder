# Scaffold `.sfdb` files

Scaffold 5 can save, and Scaffold DDA 6.x saves, a project as an `.sfdb`:
a plain SQLite database whose tables follow mzIdentML's vocabulary
(`PeptideSpectrumMatch`, `PeptideEvidence`, `ScoreIndex`, ...). opensf3
opens it read-only, never writes one, and uses no Scaffold code; the
layout below was worked out from two public files (`CORPUS.md`):

- **PXD038643** `P706_090622.sfdb` (`Version` 1, `Scaffold_5.1.2`):
  Mascot 2.7.0 results converted from an `.sf3`, two MS samples of three
  Q Exactive runs each, 28,303 PSMs on 27,627 spectra, no decoys.
- **PXD071927** `TMTsixplex_TcNMT_06-17-2025.sfdb` (`Version` 2, `6.5.0`,
  Scaffold DDA): Sequest HT 2.0 via Proteome Discoverer with Percolator,
  TMT6, 24 runs (one MS sample each), reversed decoys (`REVERSE_`).

The tables carry `DATE_CREATED`/`DATE_MODIFIED` triggers, so opening one
writably and touching a row changes it; opensf3 opens it read-only.

## Tables

| table | what opensf3 takes |
|---|---|
| `Version` | `version_name` -> the Scaffold version (`Scaffold_5.1.2`, `6.5.0`) |
| `AnalysisSoftware` | `ID`, `SOFTWARE_NAME`, `SOFTWARE_VERSION`: the engine (and Scaffold itself) |
| `IdentificationProtocol` | `ID` -> `ANALYSIS_SOFTWARE_ID`; 5.1.2 credits the protocol to Scaffold, DDA to the engine |
| `ScoreIndex` | `SCORE_TYPE` `Peptide Score` rows: `SCORE_INDEX` (position in `PEPTIDE_SCORE_ARRAY`), `SCORE_NAME`, `SCORE_ACCESSION`, `SCORE_DECREASING` |
| `PeptideSpectrumMatch` | `ID`, `MSSAMPLE_ID` + `SPECTRUM_ID` (the spectrum), `PEPTIDE_REF`, `ASSIGNED_CHARGE`, `PEPTIDE_SCORE_ARRAY` (comma-separated, `NaN` for none), `PROTOCOL_ID` |
| `Peptide` | `PEPTIDE_REF` (Scaffold's name, e.g. `VALEAC+57VK`, `<229SETPFHSR`), `PEPTIDE_SEQUENCE`, `MODIFICATION_ID_ARRAY`, `PRIMARY_LOCATION_ARRAY` |
| `Modification` | `MODIFICATION_REF` (`UNIMOD:35_M`, `Phospho_STY_S`, `UNIMOD:737_nY`), `AMINO_ACID` (a residue, or `n` for the N-terminus), `MASS` (float32 precision), `MODIFICATION_ACCESSION`, `SHORT_NAME` |
| `Substitution` | sequence variants by `PEPTIDE_REF` (none in either file; such PSMs are left out, counted) |
| `PeptideEvidence` | `PEPTIDE_REF` -> `PROTEIN_ACCESSION`, `PROTEIN_START_INDEX`/`STOP_INDEX` (1-based, inclusive), `PRE_AA`/`POST_AA` (`-` at a protein terminus) |
| `Protein` | `ACCESSION_NUMBER`, `ORIGINAL_PROTEIN_NAME`, `PROTEIN_SEQUENCE`, `IS_DECOY` |
| `MSSample` | `ID` -> `SAMPLE_NAME` (Scaffold's biosample) |
| `Spectrum` | `MSSAMPLE_ID` + `SPECTRUM_ID`, `SPECTRUM_TITLE`, `ACQUISITION_MASS` (precursor **m/z**), `ACQUISITION_CHARGE`, `SCAN_START_TIME_IN_SEC`, and in 5.1.2 the peaks (`MASS_ARRAY`, `INTENSITY_ARRAY`) |
| `SpectrumData` | DDA only: the peaks, one blob per spectrum (`FORMAT` `gzip`) |

Not read: Scaffold's verdicts (`IS_VALID`, `PASSES_THRESHOLD`,
`PRIMARY_SCORE`, `PERCOLATOR_PEP`, `IdentifiedPeptide`, `ProteinHypothesis`,
`ProteinAssembly`, `MatchToHypothesis`, `HypothesisBuilder`, thresholds),
GO tables, display settings (`Parameter`), and DDA's quant tables
(`QuantType`, `QuantChannel`, `QuantSample`, `QuantValue`, `QuantSpectra`,
`QuantEvidence`: Scaffold's reporter intensities, not read yet). Scaffold's
probabilities, protein groups and thresholds are left out on purpose: in
this layout their meaning isn't settled (`ProteinAssembly.PASSES_THRESHOLD`,
for one, is reported to be 0 everywhere in Scaffold DDA 7.1 files), and
opensf3 reports only what it can stand behind.

## Conventions (checked)

- **Scores.** 5.1.2: `Scaffold:Peptide Probability` (always `NaN`),
  `Mascot:score`, `Mascot:identity threshold`. DDA: `percolator:PEP`,
  `SEQUEST:sprank`, `Sequest:xcorr`. Scores named `Scaffold...` or
  `percolator...` (or PeptideProphet) are verdicts, not engine output, and
are left out. Mascot's expectation
  value is added from its ion score and identity threshold, as for `.sf3`.
  Rank 1 per spectrum is by the first kept score that isn't a threshold or
  a rank (`Mascot:score`, `Sequest:xcorr`), higher better unless
  `SCORE_DECREASING` is 1. `SCORE_DECREASING` is 0 for `SEQUEST:sprank`,
  whose lower values are better, so it isn't trusted beyond that. Both files hold several PSMs for some spectra (605 + 35 of 27,627 in
5.1.2), so they are ranked; equal main scores keep the file's order (PSM
id), so ranks are 1, 2, 3, ... without ties.
- **Modifications.** `PRIMARY_LOCATION_ARRAY` is 1-based on residues and 0
  for the N-terminus; the listed residue always equals `AMINO_ACID`.
  Fixed modifications are listed on every peptide (every `C` in 5.1.2
  carries `UNIMOD:4_C`; in DDA every `K` carries `UNIMOD:737_K`, checked
  over all 182,743 peptides). **A `Peptide` row doesn't pin down a PSM's
  N-terminus in Scaffold DDA.** One row stands for every PSM of Scaffold's
  peptide name, and those PSMs differ: over all 430,409 DDA PSMs,
  `<229MDYWTR` lists Met-loss+Acetyl (`UNIMOD:766_nM`) but its PSMs fit
  Met-loss (`765`) alone; 765-listed rows hold PSMs that fit 766; 8,411
  M-first rows list neither though 862 of their PSMs fit Met-loss; and
  `<42...` rows list acetyl and TMT, which their PSMs all carry, although the
  name shows only 42. So each row is read several ways: the reading whose
  N-terminal masses round to the name's `<N` (none, without one) first,
  then the other combinations of the row's listed terminal modifications,
  then each with Met-loss or Met-loss+Acetyl for an M-first peptide in a
  file that defines them. Each PSM takes the first reading its precursor
  fits (m/z x charge, within 0.05 Da, 30 ppm, or 1.5x the search's own
  precursor tolerance, on an isotope peak -1..+3). The rest are left out,
  counted; more than 5% failing refuses the file. Met-loss is kept as the engine reports it, an N-terminal delta on the M
(mzIdentML location 0), as Proteome Discoverer's tables give it. Every PSM of both files fits: 5.1.2 on its listed
  modifications, DDA with 110,574 PSMs on another reading than the name's
  (104,708 of them Met-loss).
- **Names.** `Phospho_STY` is Unimod's Phospho (mass and name prefix
  agree); DDA files terminal acetyls under the wrong Unimod numbers
  (`UNIMOD:1047`, `575`, `1163`, `1305`), so a modification whose name and mass agree with a Unimod entry gets
Unimod's own name (and so its accession in mzIdentML); otherwise the
row's Unimod accession names it only where that entry's mass is the
row's, and else the row's own name is kept.
- **Runs and scans.** 5.1.2 titles name the raw file and scan
  (`File: "Z:\...\P706_LILI_07310_WT_7T.raw"; SpectrumID: "1445"; scans:
  "5791"`), so runs come from the titles and an MS sample (biosample) holds
  several; `SPECTRUM_ID` is Scaffold's (`index=1983`). DDA titles are
  PD-style `ID-FIRST_LAST` (`1-7_7`), one scan when both match, and each MS
  sample is one raw file (`SAMPLE_NAME` `Elute-1-1`), so it is the run and
  the biosample. Two spectra of a run naming one scan both lose the scan number, counted.
- **Decoys.** 5.1.2 has none (`IS_DECOY` 0 everywhere, no decoy-looking
  accession); the decoder notes it. DDA flags reversed proteins `IS_DECOY` 1. A 1 is an explicit
  decoy; a 0 falls back to the accession convention, because the flag is
  only Scaffold's own decoy pattern (`decoyProteinRegEx` in `Parameter`),
  not the engine's statement.
- **Ranks.** Every PSM with the main score is ranked before anything is
  left out, and a spectrum whose rank-1 match is left out is left out
  whole, so a lower match never stands in as rank 1; a spectrum with a
  match lacking the main score can't be ranked and is left out too.
- **Evidence.** `PROTEIN_START_INDEX` below 1 (Scaffold DDA 7.1 is
  reported to write -1 everywhere) is rebuilt from the stop and the
  peptide's length; a peptide at two places on one protein keeps the
  first. A stored protein sequence is
  kept only if it holds every peptide placed on it (I = L): 18 of 3,505 in
  5.1.2 don't.
- **Peaks.** 5.1.2: comma-separated text on the `Spectrum` row, at most 299
  (Scaffold's 300 cap; mean 110). DDA: `SpectrumData.DATA` is a protobuf
  message `{1: "gzip", 2: payload}`, the payload gzip'd
  `{1: {2: {1: scan, 2: title fragment}}}`, the scan `{1: precursor m/z
  (double), 2: time (float), 3: precursor intensity (float, NaN),
  4: {1: peak, ...}, 5: 1}`, each peak `{1: m/z (double), 2: time (float),
  3: intensity (float)}`. The decoder follows only that path, bounds-checks
  every length, and refuses impossible peaks, messages over a million
  fields, and spectra inflating past 16 MB;
  a spectrum whose peaks don't decode is an error when its peaks are
read (MGF writes it without peaks and counts it; the mzML source yields
an `Err` for it).

## Checks on the real files

Measured with this repository's reader (`SFDB_SAMPLES=... cargo test
--release -p opensf3 --test sfdb`, and `sf3decode` on each file):

- 5.1.2 (PXD038643): 6 runs, 27,627 spectra (every one with stored
  peaks), 28,303 PSMs, 3,505 proteins, 18 stored sequences dropped, in
  0.4 s and 80 MB. Every run's mzML validates against the indexed mzML
  1.1 schema. The raw files on PRIDE are named `P706_FIRE_*` while the
  titles say `P706_LILI_07310_*`, so outputs are named after the titles.
- DDA 6.5 (PXD071927): 24 runs, 106,726 spectra (all with peaks), 430,409
  PSMs, 215,114 proteins, 110,574 PSMs on another N-terminal reading than
  their row's name (104,708 of them Met-loss), in 3.2 s and 1.1 GB. The
  PSM and spectrum counts equal those of the reader this code was ported
  from (Tuffold, same author), whose own check found the same 106,726
  spectra and scan numbers in the search's Proteome Discoverer `.msf` from
  the same dataset. The file keeps only xcorr and sprank from Sequest (no
  deltaCn).

Real-file tests run when `SFDB_SAMPLES` names a folder of `.sfdb` files
(`crates/opensf3/tests/sfdb.rs`); no `.sfdb` file is part of the
repository.
