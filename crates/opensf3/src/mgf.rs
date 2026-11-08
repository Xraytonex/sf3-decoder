//! MGF: one file per MS run, every spectrum in the run's order, with the
//! peaks Scaffold kept. An entry's position in its file is its `index=` in
//! the mzIdentML.

use std::io::{self, Write};

use crate::model::Project;

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub struct MgfReport {
    pub spectra: usize,
    /// Written without peaks: their stored peak list couldn't be read.
    pub unreadable: usize,
    /// Written without CHARGE, or without PEPMASS and CHARGE: the file
    /// stores no charge (`.sf3` stores [M+H]+, so its m/z needs one) or no
    /// precursor m/z.
    pub without_charge: usize,
}

pub fn write_mgf(project: &Project, run: usize, out: &mut impl Write) -> io::Result<MgfReport> {
    let mut report = MgfReport::default();
    for spectrum in project.run_spectra(run) {
        report.spectra += 1;
        writeln!(out, "BEGIN IONS")?;
        let title = spectrum
            .title
            .clone()
            .unwrap_or_else(|| format!("{} spectrum {}", project.format.as_str(), spectrum.id));
        writeln!(out, "TITLE={}", title.replace(['\r', '\n'], " "))?;
        // `.sf3` stores [M+H]+, so its m/z needs the charge; `.sfdb` stores
        // the m/z itself.
        match (spectrum.precursor_mz(), spectrum.charge) {
            (Some(mz), charge) => {
                writeln!(out, "PEPMASS={mz:.6}")?;
                match charge {
                    Some(charge) => writeln!(out, "CHARGE={charge}+")?,
                    None => report.without_charge += 1,
                }
            }
            (None, _) => report.without_charge += 1,
        }
        if let Some(scan) = spectrum.scan {
            writeln!(out, "SCANS={scan}")?;
        }
        if let Some(rt) = spectrum.rt_seconds {
            writeln!(out, "RTINSECONDS={rt:.3}")?;
        }
        match project.peaks(spectrum) {
            Ok(peaks) => {
                // Stored as f32; printed as such, not with f64 noise.
                for (mz, intensity) in peaks.into_iter().flatten() {
                    writeln!(out, "{} {intensity}", mz as f32)?;
                }
            }
            Err(_) => report.unreadable += 1,
        }
        writeln!(out, "END IONS")?;
    }
    Ok(report)
}
