//! Scan numbers from spectrum titles.
//!
//! Scaffold keeps each spectrum's title as the search engine got it, and
//! many titles name the scan the spectrum came from. A title gives a scan
//! number only when it names exactly one: summed spectra and scan ranges
//! give none. Two spectra of one run naming the same scan would make the
//! number ambiguous, so neither keeps it.

use std::collections::HashMap;

use crate::model::Spectrum;

/// The one scan a title names, if any: `controllerType=0 ... scan=55`,
/// DTA-style `name.START.END.CHARGE[.dta|.pkl|.mgf]` with START = END, or
/// Mascot Distiller's `14262: Scan 28842 (rt=...)`.
pub(crate) fn scan_from_title(title: &str) -> Option<u32> {
    let title = title.trim();
    if title.is_empty() {
        return None;
    }
    if let Some(scan) = native_id_scan(title) {
        return Some(scan);
    }
    if let Some(scan) = dta_scan(title) {
        return Some(scan);
    }
    distiller_scan(title)
}

/// Mascot Distiller's `14262: Scan 28842` and `87: Scan 1377 (rt=14.417,
/// f=3, i=57) [file.raw]`: the digits after the last `Scan `, ending the
/// title or followed by a space.
fn distiller_scan(title: &str) -> Option<u32> {
    let (_, rest) = title.rsplit_once("Scan ")?;
    let digits = rest.split(' ').next()?;
    digits.parse::<u32>().ok()
}

/// `scan=N` in a whitespace-separated native id. Waters ids
/// (`function=... process=... scan=...`) number scans per function, so they
/// give none.
fn native_id_scan(title: &str) -> Option<u32> {
    let words = || title.split_whitespace();
    if words().any(|word| word.starts_with("function=")) {
        return None;
    }
    words().find_map(|word| {
        word.strip_prefix("scan=")
            .and_then(|value| value.parse::<u32>().ok())
    })
}

fn dta_scan(title: &str) -> Option<u32> {
    let mut parts: Vec<&str> = title.split('.').collect();
    if parts
        .last()
        .is_some_and(|ext| ["dta", "pkl", "mgf"].contains(&ext.to_ascii_lowercase().as_str()))
    {
        parts.pop();
    }
    let [.., start, end, charge] = parts.as_slice() else {
        return None;
    };
    if parts.len() < 4 || charge.parse::<u8>().is_err() {
        return None;
    }
    let (start, end) = (start.parse::<u32>().ok()?, end.parse::<u32>().ok()?);
    (start == end).then_some(start)
}

/// Clear the scan number of every spectrum that shares it with another
/// spectrum of its run; returns how many were cleared.
pub(crate) fn drop_clashing_scans(spectra: &mut [Spectrum]) -> usize {
    let mut uses: HashMap<(usize, u32), usize> = HashMap::new();
    for spectrum in spectra.iter() {
        if let Some(scan) = spectrum.scan {
            *uses.entry((spectrum.run, scan)).or_default() += 1;
        }
    }
    let mut cleared = 0;
    for spectrum in spectra.iter_mut() {
        let clashes = spectrum
            .scan
            .is_some_and(|scan| uses.get(&(spectrum.run, scan)).is_some_and(|&n| n > 1));
        if clashes {
            spectrum.scan = None;
            cleared += 1;
        }
    }
    cleared
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scan_numbers_come_only_from_plain_titles() {
        assert_eq!(scan_from_title("Control1.0850.0850.3.pkl"), Some(850));
        assert_eq!(scan_from_title("run.12.12.2.dta"), Some(12));
        assert_eq!(scan_from_title("14262: Scan 28842"), Some(28842));
        assert_eq!(
            scan_from_title("87: Scan 1377 (rt=14.417, f=3, i=57) [x.raw]"),
            Some(1377)
        );
        assert_eq!(scan_from_title("controllerType=0 scan=55"), Some(55));
        // Summed spectra, ranges, and per-function Waters scans have no
        // single scan.
        for title in [
            "run.12.14.2.dta",
            "8144: Sum of 3 scans in range 17366 to 17451 (rt=3862.31)",
            "Discoverer_Tutorial_SILAC-2-47_47",
            "function=2 process=0 scan=17",
            "",
        ] {
            assert_eq!(scan_from_title(title), None, "{title}");
        }
    }

    fn spectrum(run: usize, scan: Option<u32>) -> Spectrum {
        Spectrum::for_tests(run, scan)
    }

    #[test]
    fn a_scan_named_twice_in_one_run_is_dropped() {
        let mut spectra = vec![
            spectrum(0, Some(5)),
            spectrum(0, Some(5)),
            spectrum(1, Some(5)),
            spectrum(0, Some(6)),
            spectrum(0, None),
        ];
        assert_eq!(drop_clashing_scans(&mut spectra), 2);
        let scans: Vec<Option<u32>> = spectra.iter().map(|s| s.scan).collect();
        assert_eq!(scans, [None, None, Some(5), Some(6), None]);
    }
}
