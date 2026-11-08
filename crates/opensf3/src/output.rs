//! Output file names.

use std::collections::HashSet;

use crate::model::Project;

/// Extensions of spectra files a run may be named after.
const SPECTRA_EXTENSIONS: &[&str] = &[
    "raw", "mzml", "mzxml", "mgf", "d", "wiff", "wiff2", "baf", "pkl", "dta", "mzdata", "ms2",
];

/// One MGF name per run: the run's file name without its path or spectra
/// extension, made safe for any file system, and unique.
pub fn mgf_names(project: &Project) -> Vec<String> {
    let mut taken = HashSet::new();
    project
        .runs
        .iter()
        .map(|run| {
            let base = file_stem(&run.name);
            let mut name = format!("{base}.mgf");
            let mut copy = 1;
            // Case-insensitive file systems would merge names differing in case.
            while !taken.insert(name.to_lowercase()) {
                copy += 1;
                name = format!("{base}_{copy}.mgf");
            }
            name
        })
        .collect()
}

fn file_stem(run: &str) -> String {
    let last = run.rsplit(['/', '\\']).next().unwrap_or(run);
    let mut stem = String::with_capacity(last.len());
    let mut rest = last;
    // Drop spectra extensions wherever they end a word: `x.raw (F01)`.
    while let Some(dot) = rest.find('.') {
        stem.push_str(&rest[..dot]);
        let after = &rest[dot + 1..];
        let word = after
            .find(|c: char| !c.is_ascii_alphanumeric())
            .map_or(after, |end| &after[..end]);
        if SPECTRA_EXTENSIONS.contains(&word.to_ascii_lowercase().as_str()) {
            rest = &after[word.len()..];
        } else {
            stem.push('.');
            rest = after;
        }
    }
    stem.push_str(rest);
    let mut safe = String::with_capacity(stem.len());
    for c in stem.chars() {
        let c = if c.is_ascii_alphanumeric() || "-.+".contains(c) {
            c
        } else {
            '_'
        };
        if !(c == '_' && safe.ends_with('_')) {
            safe.push(c);
        }
    }
    let safe = safe.trim_matches(['_', '.']);
    if safe.is_empty() {
        "run".to_string()
    } else {
        safe.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn run_names_become_safe_file_names() {
        assert_eq!(file_stem(r"C:\data\QE_0412.raw"), "QE_0412");
        assert_eq!(file_stem("/x/Sample 3 (rep).mzML"), "Sample_3_rep");
        assert_eq!(file_stem("QC_01.raw (F039491)"), "QC_01_F039491");
        assert_eq!(file_stem("Ctrl: fraction.1"), "Ctrl_fraction.1");
        assert_eq!(file_stem(".raw"), "run");
        assert_eq!(file_stem("///"), "run");
    }
}
