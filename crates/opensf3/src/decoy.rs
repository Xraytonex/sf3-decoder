//! Decoy recognition by accession convention.
//!
//! Target-decoy databases mark decoys with a label at the start of the
//! accession (`rev_sp|P12345|...`, `DECOY_P12345`, `XXX_...`). A label only
//! counts when a non-alphanumeric separator follows it, so accessions that
//! merely begin with those letters (`DECA1_HUMAN`, `REV3L`) stay targets.
//! Scaffold's databases put the label at the end instead (`P15975-DECOY`):
//! a trailing `decoy` after a separator counts too. Some sites label the
//! accession inside a UniProt-style name instead (`rr|REV_P12345|REV_NAME`):
//! a label at the start of the second `|` field counts as well.
//! Scaffold's own decoy flag on a protein group takes precedence over this
//! rule.

const DECOY_LABELS: &[&str] = &[
    "rev", "decoy", "xxx", "dec", "reverse", "reversed", "random", "shuffled",
];
/// Labels recognized at the end of an accession (Scaffold's convention).
const DECOY_SUFFIXES: &[&str] = &["decoy"];

pub fn is_decoy_accession(accession: &str) -> bool {
    decoy_label_rest(accession).is_some()
        || suffixed_target(accession).is_some()
        || labelled_id_field(accession)
}

/// `db|LABEL_ID|NAME`: the identifier field carries the label.
fn labelled_id_field(accession: &str) -> bool {
    let first_token = accession.split_whitespace().next().unwrap_or_default();
    let mut fields = first_token.split('|');
    fields.next();
    fields
        .next()
        .is_some_and(|id| decoy_label_rest(id).is_some_and(|rest| rest.len() > 1))
}

/// The accession before a trailing decoy label and its separator
/// (`P1-DECOY` -> `P1`); `None` when there is no such label or nothing
/// precedes it.
fn suffixed_target(accession: &str) -> Option<&str> {
    let first_token = accession.split_whitespace().next().unwrap_or_default();
    DECOY_SUFFIXES.iter().find_map(|label| {
        let cut = first_token.len().checked_sub(label.len())?;
        let (head, tail) = (first_token.get(..cut)?, first_token.get(cut..)?);
        let mut before = head.chars();
        let separator = before.next_back()?;
        let target = before.as_str();
        (tail.eq_ignore_ascii_case(label)
            && !separator.is_alphanumeric()
            && target.chars().any(char::is_alphanumeric))
        .then_some(target)
    })
}

/// What follows the decoy label in the accession's first token, separator
/// included; `None` when there is no label.
fn decoy_label_rest(accession: &str) -> Option<&str> {
    let first_token = accession.split_whitespace().next().unwrap_or_default();
    DECOY_LABELS.iter().find_map(|label| {
        let head = first_token.get(..label.len())?;
        let rest = &first_token[label.len()..];
        let separated = rest.chars().next().is_none_or(|ch| !ch.is_alphanumeric());
        (head.eq_ignore_ascii_case(label) && separated).then_some(rest)
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn recognizes_common_decoy_labels() {
        for accession in [
            "rev_sp|P1|N1",
            "DECOY_P1",
            "XXX_P1",
            "reversed_P1",
            "reverse-258", // TPP's numbered reverse decoys
            "rev",
            "##decoy",
            "Q5XJE5-DECOY",
            "rr|REV_Q3Z8Z2|REV_ATPB_DEHE1", // FGCZ's databases
        ] {
            let expected = accession != "##decoy";
            assert_eq!(is_decoy_accession(accession), expected, "{accession}");
        }
    }

    #[test]
    fn plain_accessions_that_share_letters_stay_targets() {
        for accession in [
            "revolver_P1",
            "reverser_P1",
            "DECA1_HUMAN",
            "REV3L",
            "P00001",
            "sp|P01008|ANT3_HUMAN",
            "PREDECOY1",
            "P1DECOY",
            "sp|Q9UBT6|REV3L_HUMAN",
            "zz|REV3L|contaminant",
            "sp|REV_|empty",
        ] {
            assert!(!is_decoy_accession(accession), "{accession}");
        }
    }
}
