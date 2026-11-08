use std::fs::{File, OpenOptions};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::{Parser, ValueEnum};
use openmassspec_core::write_indexed_mzml;
use opensf3::mgf::write_mgf;
use opensf3::mzid::write_mzid;
use opensf3::output::mgf_names;
use opensf3::{Format as Kind, Project, Reader};

/// Convert Proteome Software's Scaffold files (.sf3, .sfdb) into more
/// commonly used formats.
#[derive(Parser)]
#[command(name = "sf3decode", version)]
struct Args {
    /// The Scaffold .sf3 or .sfdb file.
    file: PathBuf,
    /// What to write; without it, prints what the file holds.
    format: Option<Format>,
    /// Folder to write into [default: the Scaffold file's folder].
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
    /// mzIdentML 1.2: every engine's candidates with their scores, and
    /// for .sf3 Scaffold's probabilities and protein groups.
    #[value(alias = "mzidentml")]
    Mzid,
    /// One indexed mzML 1.1 per MS run, with every spectrum's peaks.
    Mzml,
}

fn main() -> ExitCode {
    let args = Args::parse();
    let reader = match Reader::open(&args.file) {
        Ok(reader) => reader,
        Err(error) => {
            eprintln!("sf3decode: {}: {error}", args.file.display());
            return ExitCode::FAILURE;
        }
    };
    let project = reader.project();
    let Some(format) = args.format else {
        summary(project);
        return ExitCode::SUCCESS;
    };
    let folder = args.out.clone().unwrap_or_else(|| {
        args.file
            .parent()
            .map_or_else(|| PathBuf::from("."), Path::to_path_buf)
    });
    let written = match format {
        Format::Mgf => mgf(project, &folder, args.force),
        Format::Mzid => {
            let stem = args.file.file_stem().unwrap_or_default().to_string_lossy();
            mzid(project, &folder.join(format!("{stem}.mzid")), args.force)
        }
        Format::Mzml => mzml(&reader, &folder, args.force),
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

/// One indexed mzML per run, named like its MGF.
fn mzml(reader: &Reader, folder: &Path, force: bool) -> Result<(), String> {
    let names = mgf_names(reader.project());
    for (run, name) in names.iter().enumerate() {
        let stem = name.strip_suffix(".mgf").unwrap_or(name);
        let path = folder.join(format!("{stem}.mzML"));
        let mut source = reader.run_spectra(run).map_err(|error| error.to_string())?;
        let mut out = create(&path, force)?;
        let written = write_indexed_mzml(&mut source, &mut out)
            .map_err(|error| error.to_string())
            .and_then(|summary| {
                out.flush()
                    .map(|()| summary.spectra_written)
                    .map_err(|error| error.to_string())
            })
            .map_err(|error| format!("{}: {error}", path.display()))?;
        eprintln!("wrote {} ({written} spectra)", path.display());
    }
    Ok(())
}

fn mzid(project: &Project, path: &Path, force: bool) -> Result<(), String> {
    let mut out = create(path, force)?;
    let report = write_mzid(project, &mut out)
        .and_then(|report| out.flush().map(|()| report))
        .map_err(|error| format!("{}: {error}", path.display()))?;
    match project.format {
        Kind::Sf3 => eprintln!(
            "wrote {} ({} candidates, {} protein groups)",
            path.display(),
            report.candidates,
            report.groups
        ),
        _ => eprintln!(
            "wrote {} ({} candidates; no protein groups: Scaffold's aren't read from an .sfdb)",
            path.display(),
            report.candidates
        ),
    }
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
    println!("{} (.{}, {version})", project.name, project.format.as_str());
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
        let version = engine.version.as_deref().unwrap_or("");
        match project.format {
            Kind::Sf3 => println!(
                "  {} {version}: {total} candidates, {chosen} picked by Scaffold",
                engine.name
            ),
            _ => println!("  {} {version}: {total} candidates", engine.name),
        }
    }
    match project.format {
        Kind::Sf3 => println!(
            "  {} proteins, {} protein groups",
            project.proteins.len(),
            project.groups.len()
        ),
        _ => println!("  {} proteins", project.proteins.len()),
    }
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
