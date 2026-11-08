//! The crate's error type.

use crate::sf3::Sf3Error;
use crate::sfdb::SfdbError;

/// Why a Scaffold file couldn't be opened or read.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum Error {
    /// The file couldn't be read at all.
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    /// Neither an `.sf3` nor an `.sfdb`.
    #[error("not a Scaffold .sf3 or .sfdb file: {0}")]
    Unrecognized(String),
    /// An `.sf3` that couldn't be decoded.
    #[error(transparent)]
    Sf3(#[from] Sf3Error),
    /// An `.sfdb` that couldn't be decoded.
    #[error(transparent)]
    Sfdb(#[from] SfdbError),
    /// A run index past the project's runs.
    #[error("the file has {runs} runs; there is no run {run}")]
    NoSuchRun { run: usize, runs: usize },
    /// A spectrum's stored peaks couldn't be read.
    #[error("spectrum {native_id}: its stored peaks can't be read: {problem}")]
    Peaks {
        native_id: String,
        problem: &'static str,
    },
}

impl From<Error> for openmassspec_core::Error {
    fn from(error: Error) -> Self {
        match error {
            Error::Io(error) => Self::Io(error),
            Error::Unrecognized(problem) => Self::Format(problem),
            other => Self::vendor(other),
        }
    }
}
