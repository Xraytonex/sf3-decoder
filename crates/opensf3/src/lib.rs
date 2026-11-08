//! Reader for Proteome Software Scaffold files: `.sf3` projects (Scaffold
//! 4.0 to 5.3.4) and `.sfdb` databases (Scaffold 5.1, Scaffold DDA 6.5).
//!
//! [`Reader::open`] tells the two apart by content and decodes either into
//! a [`Project`]: runs, spectra, every search engine's candidates with
//! their scores, proteins, and for `.sf3` Scaffold's own probabilities,
//! protein groups and saved thresholds. Each run's stored peak lists are an
//! [`openmassspec_core::SpectrumSource`] ([`Reader::run_spectra`]), so the
//! OpenMassSpec mzML writer can write them. [`mgf`] and (with the default
//! `mzid` feature) [`mzid`] write MGF and mzIdentML 1.2.
//!
//! Nothing here runs or links Scaffold; the formats were worked out from
//! public data files (see `docs/format/`). Files are only read.

#![cfg_attr(not(test), warn(clippy::unwrap_used, clippy::expect_used))]

pub mod chemistry;
mod decoy;
mod error;
pub mod mgf;
pub mod model;
#[cfg(feature = "mzid")]
pub mod mzid;
pub mod output;
mod reader;
pub mod sf3;
pub mod sfdb;
mod titles;
mod unimod;
mod uniprot;

pub use error::Error;
pub use model::{Format, Project};
pub use reader::{detect, Reader, RunSpectra};
pub use sf3::{decode, Sf3Error};
pub use sfdb::SfdbError;
