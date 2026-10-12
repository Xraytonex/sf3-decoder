//! Checks on real `.sf3` files, run when `SF3_SAMPLES` names a folder of
//! them (`.sf3` files are never part of the repository).

use std::collections::HashMap;
use std::path::PathBuf;

use sf3_decoder::chemistry::mass_error;

fn samples() -> Vec<PathBuf> {
    let Some(folder) = std::env::var_os("SF3_SAMPLES") else {
        eprintln!("SF3_SAMPLES not set; skipping");
        return Vec::new();
    };
    let mut files: Vec<PathBuf> = std::fs::read_dir(folder)
        .expect("SF3_SAMPLES is a folder")
        .map(|entry| entry.expect("readable").path())
        .filter(|path| path.extension().is_some_and(|extension| extension == "sf3"))
        .collect();
    files.sort();
    files
}

#[test]
fn every_sample_decodes_consistently() {
    for path in samples() {
        let project = sf3_decoder::decode(&path)
            .unwrap_or_else(|error| panic!("{}: {error}", path.display()));
        let name = path.display();
        assert!(!project.runs.is_empty(), "{name}");
        let mut unreadable = 0;
        let mut with_peaks = 0;
        for spectrum in &project.spectra {
            // Ranks rise within each engine (a left-out candidate leaves a gap).
            let mut last: HashMap<usize, u32> = HashMap::new();
            for found in &spectrum.matches {
                let previous = last.insert(found.engine, found.rank).unwrap_or(1);
                assert!(found.rank >= previous, "{name}");
                if let (Some(mass), Some(mh)) = (found.mass, spectrum.mh) {
                    assert!(mass_error(mass, mh).abs() <= 1.0, "{name}");
                }
            }
            match project.peaks(spectrum) {
                Ok(Some(peaks)) => {
                    with_peaks += 1;
                    assert!(peaks.len() <= 300, "{name}");
                    assert!(peaks.windows(2).all(|pair| pair[0].0 <= pair[1].0));
                }
                Ok(None) => {}
                Err(_) => unreadable += 1,
            }
        }
        assert_eq!(unreadable, 0, "{name}: unreadable peak lists");
        assert!(with_peaks > 0, "{name}: no peaks");
        let picks = project
            .spectra
            .iter()
            .flat_map(|spectrum| &spectrum.matches)
            .filter(|found| found.chosen);
        assert!(picks.clone().count() > 0, "{name}");
        assert!(
            picks.clone().all(|found| found.probability.is_some()),
            "{name}"
        );
        eprintln!(
            "{name}: {} spectra, {with_peaks} with peaks, {} groups",
            project.spectra.len(),
            project.groups.len()
        );
    }
}
