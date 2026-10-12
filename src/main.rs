use std::path::PathBuf;
use std::process::ExitCode;

use clap::Parser;
use sf3_decoder::Project;

/// Convert Proteome Software's .sf3 files into more commonly used formats.
#[derive(Parser)]
#[command(name = "sf3decode", version)]
struct Args {
    /// The Scaffold .sf3 file.
    file: PathBuf,
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
    summary(&project);
    ExitCode::SUCCESS
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
