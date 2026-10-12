//! Which proteins a candidate maps to.
//!
//! Each record names the protein the engine reported (a
//! `ProteinIdentification`), whose accession is stored inline cut at 22
//! characters, so a decoy `tr|A0QSM5|A0QSM5_MYCS2-DECOY` reads as its
//! target. A name that long is completed from Scaffold's peptide → protein
//! map, as a decoy when it ends in part of `-DECOY`, else by its target's
//! stored sequence (the peptide read forwards or backwards in it).

use std::collections::{BTreeSet, HashMap};
use std::ops::Bound;

use super::project::ProteinEntry;
use super::records::Records;
use crate::decoy::is_decoy_accession;
use crate::model::ProteinRef;
use crate::uniprot::uniprot_accession;

pub(crate) const PROTEIN_BLOCK: usize = 103;
const ACCESSION: std::ops::Range<usize> = 8..31;
pub(crate) const INLINE_ACCESSION_MAX: usize = 22;

/// What the engine's reported protein can be taken as.
#[derive(Debug, PartialEq, Eq)]
enum Reported {
    Missing,
    /// Shorter than the inline maximum, so whole.
    Whole(String),
    /// Possibly cut, completed to one name.
    Completed(String),
    /// Possibly cut, and the map lists the peptide but no name starting with
    /// it: the map's proteins stand alone.
    MapOnly,
    /// Possibly cut, and the map doesn't list the peptide.
    Unknown(String),
}

fn reported(inline: &[u8], listed: Option<&[String]>) -> Reported {
    let end = inline
        .iter()
        .position(|&byte| byte == 0x0a)
        .unwrap_or(inline.len());
    let accession = String::from_utf8_lossy(&inline[..end]).into_owned();
    if accession.is_empty() {
        return Reported::Missing;
    }
    if end < INLINE_ACCESSION_MAX {
        return Reported::Whole(accession);
    }
    let Some(listed) = listed else {
        return match decoy_completion(&accession) {
            Some(decoy) => Reported::Completed(decoy),
            None => Reported::Unknown(accession),
        };
    };
    let mut completions = listed.iter().filter(|full| full.starts_with(&accession));
    match (completions.next(), completions.next()) {
        (Some(full), None) => Reported::Completed(full.clone()),
        _ => Reported::MapOnly,
    }
}

/// Scaffold's decoys are their target's accession plus `-DECOY`, so a cut
/// name ending in the start of that suffix is a decoy (`…_MYCS2-` from a
/// 21-character target).
fn decoy_completion(cut: &str) -> Option<String> {
    const SUFFIX: &str = "-DECOY";
    (1..=SUFFIX.len())
        .find(|&kept| cut.len() > kept && cut.ends_with(&SUFFIX[..kept]))
        .map(|kept| format!("{}{SUFFIX}", &cut[..cut.len() - kept]))
}

/// Names a cut accession can complete to.
pub(crate) struct Completions<'a> {
    stored: &'a HashMap<String, ProteinEntry>,
    /// Stored names with a sequence, sorted for prefix lookups.
    sequenced: BTreeSet<String>,
    /// Decoy names in Scaffold's protein map.
    decoys: BTreeSet<String>,
}

impl<'a> Completions<'a> {
    pub fn new(
        stored: &'a HashMap<String, ProteinEntry>,
        protein_map: &HashMap<String, Vec<String>>,
    ) -> Self {
        Self {
            stored,
            sequenced: stored
                .iter()
                .filter(|(_, entry)| entry.sequence.is_some())
                .map(|(name, _)| name.clone())
                .collect(),
            decoys: protein_map
                .values()
                .flatten()
                .filter(|name| is_decoy_accession(name))
                .cloned()
                .collect(),
        }
    }

    fn sequence(&self, name: &str) -> Option<&'a str> {
        self.stored
            .get(name)
            .and_then(|entry| entry.sequence.as_deref())
    }

    /// The one name in `names` starting with `cut`.
    fn only<'s>(names: &'s BTreeSet<String>, cut: &str) -> Option<&'s str> {
        let mut found = names
            .range::<str, _>((Bound::Included(cut), Bound::Unbounded))
            .map(String::as_str)
            .take_while(|name| name.starts_with(cut));
        match (found.next(), found.next()) {
            (Some(name), None) => Some(name),
            _ => None,
        }
    }
}

/// Settle a cut name by its target's stored sequence: Scaffold's decoys are
/// reversed proteins, so a peptide read forwards (I and L alike) is the
/// target's, read backwards (whole, or but its last residue) its decoy's.
/// Both, neither, or no sequence: `None`.
fn by_sequence(cut: &str, peptide: &str, names: &Completions<'_>) -> Option<String> {
    let target = if names.sequence(cut).is_some() {
        cut
    } else if let Some(uniprot) = uniprot_accession(cut).filter(|id| names.sequence(id).is_some()) {
        uniprot
    } else {
        Completions::only(&names.sequenced, cut)?
    };
    let leucine = |text: &str| text.replace('I', "L");
    let sequence = leucine(names.sequence(target)?);
    let peptide = leucine(peptide);
    let reversed: String = peptide.chars().rev().collect();
    let reversed_but_last: String = peptide.chars().rev().skip(1).collect();
    let forwards = sequence.contains(&peptide);
    let backwards = sequence.contains(&reversed)
        || (reversed_but_last.len() >= 4 && sequence.contains(&reversed_but_last));
    match (forwards, backwards) {
        (true, false) => Some(target.to_string()),
        (false, true) => Some(
            Completions::only(&names.decoys, cut)
                .map_or_else(|| format!("{cut}-DECOY"), str::to_string),
        ),
        _ => None,
    }
}

/// The candidate's proteins: the one the engine reported, then every
/// protein Scaffold maps its peptide to. `None` when the reported name may
/// be cut and nothing settles which protein it is.
pub(crate) fn proteins_of(
    peptide: &str,
    protein_record: i32,
    flanks: [u8; 2],
    protein_map: &HashMap<String, Vec<String>>,
    protein_records: &Records,
    names: &Completions<'_>,
) -> Option<Vec<ProteinRef>> {
    let flank = |byte: u8| (byte.is_ascii_alphabetic() || byte == b'-').then_some(byte as char);
    let listed = protein_map.get(peptide);
    let mut proteins: Vec<ProteinRef> = Vec::new();
    if let Some(block) = protein_records.get(protein_record) {
        let settled = match reported(&block[ACCESSION], listed.map(Vec::as_slice)) {
            Reported::Unknown(cut) => Reported::Completed(by_sequence(&cut, peptide, names)?),
            other => other,
        };
        if let Reported::Whole(accession) | Reported::Completed(accession) = settled {
            proteins.push(ProteinRef {
                accession,
                pre: flank(flanks[0]),
                post: flank(flanks[1]),
            });
        }
    }
    for accession in listed.into_iter().flatten() {
        if !proteins.iter().any(|known| &known.accession == accession) {
            proteins.push(ProteinRef {
                accession: accession.clone(),
                pre: None,
                post: None,
            });
        }
    }
    Some(proteins)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn inline(accession: &str) -> Vec<u8> {
        let mut bytes = accession.as_bytes().to_vec();
        bytes.push(0x0a);
        bytes.resize(23, 0);
        bytes
    }

    #[test]
    fn cut_accessions_complete_from_the_map_or_the_decoy_suffix() {
        assert_eq!(
            reported(&inline("sp|P1|SHORT"), None),
            Reported::Whole("sp|P1|SHORT".into())
        );
        let cut = "tr|A0QSM5|A0QSM5_MYCS2";
        let listed = vec!["tr|A0QSM5|A0QSM5_MYCS2-DECOY".to_string()];
        assert_eq!(
            reported(&inline(cut), Some(&listed)),
            Reported::Completed(listed[0].clone())
        );
        assert_eq!(
            reported(&inline("sp|A0QQ48|ECCE3_MYCS2-"), None),
            Reported::Completed("sp|A0QQ48|ECCE3_MYCS2-DECOY".into())
        );
        assert_eq!(reported(&inline(cut), None), Reported::Unknown(cut.into()));
    }

    #[test]
    fn a_cut_name_is_settled_by_reading_the_peptide_in_its_target() {
        let stored: HashMap<String, ProteinEntry> = [(
            "tr|A0QSM5|A0QSM5_MYCS2".to_string(),
            ProteinEntry {
                description: None,
                sequence: Some("MKPEPTIDERAAA".into()),
            },
        )]
        .into();
        let names = Completions::new(&stored, &HashMap::new());
        let cut = "tr|A0QSM5|A0QSM5_MYCS2";
        assert_eq!(by_sequence(cut, "PEPTLDER", &names).as_deref(), Some(cut));
        assert_eq!(
            by_sequence(cut, "REDITPEP", &names),
            Some(format!("{cut}-DECOY"))
        );
        assert_eq!(by_sequence(cut, "WWWWW", &names), None);
    }
}
