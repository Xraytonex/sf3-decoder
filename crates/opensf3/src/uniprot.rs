//! UniProt accessions inside the accessions search engines report.

/// The UniProt accession in a protein accession: `P02768` from
/// `sp|P02768|ALBU_HUMAN`, `P02768-2`, or `P02768`; isoform suffixes
/// dropped. `None` when nothing looks like one.
pub fn uniprot_accession(accession: &str) -> Option<&str> {
    let parts: Vec<&str> = accession.split('|').collect();
    let candidate = match parts.as_slice() {
        [_, id, ..] => *id,
        [only] => *only,
        _ => return None,
    };
    let base = candidate.split('-').next().unwrap_or(candidate);
    let looks_right = (6..=10).contains(&base.len())
        && base
            .chars()
            .all(|c| c.is_ascii_uppercase() || c.is_ascii_digit())
        && base.chars().next().is_some_and(|c| c.is_ascii_uppercase())
        && base.chars().nth(1).is_some_and(|c| c.is_ascii_digit());
    looks_right.then_some(base)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn uniprot_accessions_come_out_of_fasta_style_ones() {
        assert_eq!(uniprot_accession("sp|P02768|ALBU_HUMAN"), Some("P02768"));
        assert_eq!(
            uniprot_accession("tr|A0A024R161|A0A024R161_HUMAN"),
            Some("A0A024R161")
        );
        assert_eq!(uniprot_accession("P02768-2"), Some("P02768"));
        assert_eq!(uniprot_accession("P02768"), Some("P02768"));
        assert_eq!(uniprot_accession("ENSP00000295897"), None);
        assert_eq!(uniprot_accession("DEMO0011"), None);
    }
}
