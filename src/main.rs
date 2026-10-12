use std::fs::{File, OpenOptions};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::{Parser, ValueEnum};
use sf3_decoder::Project;
use sf3_decoder::mgf::write_mgf;
use sf3_decoder::mzid::write_mzid;
use sf3_decoder::output::mgf_names;

/// Convert Proteome Software's .sf3 files into more commonly used formats.
#[derive(Parser)]
#[command(name = "sf3decode", version)]
struct Args {
    /// The Scaffold .sf3 file.
    file: PathBuf,
    /// What to write; without it, prints what the file holds.
    format: Option<Format>,
    /// Folder to write into [default: the .sf3 file's folder].
    #[arg(short, long)]
    out: Option<PathBuf>,
    /// Replace files that already exist.
    #[arg(long)]
    force: bool,
}

#[derive(Clone, Copy, ValueEnum)]
enum Format {
    /// One MGF per MS run, with every spectrum's peaks.
    Mgf,
    /// mzIdentML 1.2: every engine's candidates, Scaffold's probabilities
    /// and protein groups.
    #[value(alias = "mzidentml")]
    Mzid,
}

fn main() -> ExitCode {
    let args = Args::parse();
    let project = match sf3_decoder::decode(&args.file) {
        Ok(project) => project,
        Err(error) => {
            eprintln!("sf3decode: {}: {error}", args.file.display());
            return ExitCode::FAILURE;
        }
    };
    let Some(format) = args.format else {
        summary(&project);
        return ExitCode::SUCCESS;
    };
    let folder = args.out.clone().unwrap_or_else(|| {
        args.file
            .parent()
            .map_or_else(|| PathBuf::from("."), Path::to_path_buf)
    });
    let written = match format {
        Format::Mgf => mgf(&project, &folder, args.force),
        Format::Mzid => {
            let stem = args.file.file_stem().unwrap_or_default().to_string_lossy();
            mzid(&project, &folder.join(format!("{stem}.mzid")), args.force)
        }
    };
    match written {
        Ok(()) => {
            for note in &project.notes {
                eprintln!("note: {note}");
            }
            ExitCode::SUCCESS
        }
        Err(error) => {
            eprintln!("sf3decode: {error}");
            ExitCode::FAILURE
        }
    }
}

fn mgf(project: &Project, folder: &Path, force: bool) -> Result<(), String> {
    for (run, name) in mgf_names(project).iter().enumerate() {
        let path = folder.join(name);
        let mut out = create(&path, force)?;
        let report = write_mgf(project, run, &mut out)
            .and_then(|report| out.flush().map(|()| report))
            .map_err(|error| format!("{}: {error}", path.display()))?;
        eprintln!("wrote {} ({} spectra)", path.display(), report.spectra);
        if report.unreadable > 0 {
            eprintln!(
                "  {} spectra written without peaks: their peak lists couldn't be read",
                report.unreadable
            );
        }
        if report.without_charge > 0 {
            eprintln!(
                "  {} spectra written without PEPMASS: Scaffold stored no charge for them",
                report.without_charge
            );
        }
    }
    Ok(())
}

fn mzid(project: &Project, path: &Path, force: bool) -> Result<(), String> {
    let mut out = create(path, force)?;
    let report = write_mzid(project, &mut out)
        .and_then(|report| out.flush().map(|()| report))
        .map_err(|error| format!("{}: {error}", path.display()))?;
    eprintln!(
        "wrote {} ({} candidates, {} protein groups)",
        path.display(),
        report.candidates,
        report.groups
    );
    if report.groups_without_peptides > 0 {
        eprintln!(
            "  {} protein groups left out: none of their peptides is a candidate Scaffold picked",
            report.groups_without_peptides
        );
    }
    if report.sequences_left_out > 0 {
        eprintln!(
            "  {} protein sequences left out: they hold a * inside, which mzIdentML doesn't allow",
            report.sequences_left_out
        );
    }
    if report.unknown_modifications > 0 {
        eprintln!(
            "  {} modifications not in the Unimod table written as \"unknown modification\" \
             with their name and mass",
            report.unknown_modifications
        );
    }
    Ok(())
}

/// A new file; an existing one only with `--force`.
fn create(path: &Path, force: bool) -> Result<BufWriter<File>, String> {
    let mut options = OpenOptions::new();
    options.write(true);
    if force {
        options.create(true).truncate(true);
    } else {
        options.create_new(true);
    }
    options.open(path).map(BufWriter::new).map_err(|error| {
        if error.kind() == std::io::ErrorKind::AlreadyExists {
            format!("{} exists; use --force to replace it", path.display())
        } else {
            format!("{}: {error}", path.display())
        }
    })
}

/// What the file holds.
fn summary(project: &Project) {
    let version = project
        .scaffold_version
        .as_deref()
        .unwrap_or("unknown version");
    println!("{} ({version})", project.name);
    println!(
        "  {} runs, {} spectra",
        project.runs.len(),
        project.spectra.len()
    );
    for (index, engine) in project.engines.iter().enumerate() {
        let matches = project
            .spectra
            .iter()
            .flat_map(|spectrum| &spectrum.matches)
            .filter(|found| found.engine == index);
        let (mut total, mut chosen) = (0, 0);
        for found in matches {
            total += 1;
            chosen += usize::from(found.chosen);
        }
        println!(
            "  {} {}: {total} candidates, {chosen} picked by Scaffold",
            engine.name,
            engine.version.as_deref().unwrap_or("")
        );
    }
    println!(
        "  {} proteins, {} protein groups",
        project.proteins.len(),
        project.groups.len()
    );
    if let Some(thresholds) = project.thresholds {
        println!(
            "  saved thresholds: peptides {:.0}%, proteins {:.0}%, at least {} peptides",
            thresholds.peptide_probability * 100.0,
            thresholds.protein_probability * 100.0,
            thresholds.min_peptides
        );
    }
    for note in &project.notes {
        println!("  note: {note}");
    }
}
