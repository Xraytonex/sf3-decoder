//! Decodes Scaffold `.sf3` projects.

pub mod chemistry;
mod decoy;
pub mod mgf;
pub mod model;
pub mod mzid;
pub mod output;
pub mod sf3;
mod unimod;
mod uniprot;

pub use model::Project;
pub use sf3::{Sf3Error, decode};
