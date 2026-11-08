# Validation corpus

Public Scaffold files the decoder was worked out from and checked on. None
of them is part of the repository; the real-file tests read them from a
folder named by `SF3_SAMPLES` or `SFDB_SAMPLES`. Only files whose source
is recorded are named by accession; for the rest the source column says
what is known.

| Accession | Source | Files | Format | Scaffold version | Notes |
|---|---|---|---|---|---|
| PXD038643 | PRIDE | `P706_090622.sfdb` | `.sfdb` | Scaffold_5.1.2 | Mascot 2.7.0; 6 runs, 27,627 spectra, 28,303 PSMs, no decoys. Checked on this release (`docs/format/02-sfdb.md`). |
| PXD071927 | PRIDE | `TMTsixplex_TcNMT_06-17-2025.sfdb` | `.sfdb` | 6.5.0 (Scaffold DDA) | Sequest HT 2.0 via Proteome Discoverer, Percolator, TMT6; 24 runs, 106,726 spectra, 430,409 PSMs. Checked on this release. |
| PXD068832 | PRIDE | an `.sf3` | `.sf3` | Scaffold 5 | Mascot. Named in the format notes (accessions cut at 22 characters); not re-run for this release. |
| PXD076594 | PRIDE | an `.sf3` | `.sf3` | Scaffold 5.2.2 | Sequest via Proteome Discoverer. Showed the 22-character accession cut; not re-run for this release. |
| PXD068782 | PRIDE | an `.sf3` | `.sf3` | Scaffold Q+S 5 | MSFragger; score layout in `docs/format/01-sf3.md`. Not re-run for this release. |
| PXD074543 | PRIDE | an `.sf3` | `.sf3` | Scaffold 5.3.4 | "Mascot (Ion Score Only)" engine name. Not re-run for this release. |
| - | eLife (file name: article 31807, figure 1, data 1) | `elife-31807-fig1-data1-v2.sf3` | `.sf3` | Scaffold_4.0.5 | Mascot 2.2.04; 4 runs, 514 spectra. Checked on this release. |
| - | not recorded | `Ctx_Cereb_Hip.sf3` | `.sf3` | Scaffold_4.4.0 | Spectrum Mill and X! Tandem; 20,303 spectra. Checked on this release. |
| - | not recorded (Scaffold Q+S SILAC demo projects) | `scaffold_q+s_silac_demo_discoverer.sf3`, `..._distiller.sf3`, `..._maxquant.sf3` | `.sf3` | Scaffold_4.11.0 | Sequest, Mascot Distiller, MaxQuant; 29,792 / 19,292 / 18,311 spectra. Checked on this release. |

The original `.sf3` work used about twenty public files from Scaffold 4.0
through 5.3.4; the ones above are those whose source is on record.
