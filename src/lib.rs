//! Decodes Scaffold `.sf3` projects.

pub mod chemistry;
mod decoy;
pub mod model;
pub mod sf3;
mod uniprot;

pub use model::Project;
pub use sf3::{Sf3Error, decode};
