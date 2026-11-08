//! Decodes Scaffold `.sf3` projects.

#![cfg_attr(not(test), warn(clippy::unwrap_used, clippy::expect_used))]

pub mod chemistry;
mod decoy;
pub mod mgf;
pub mod model;
pub mod mzid;
pub mod output;
pub mod sf3;
mod titles;
mod unimod;
mod uniprot;

pub use model::Project;
pub use sf3::{decode, Sf3Error};
