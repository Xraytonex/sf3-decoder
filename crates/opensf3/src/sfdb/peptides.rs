//! Peptides and their modifications in an `.sfdb` (`docs/format/02-sfdb.md`).
//!
//! A `Peptide` row lists its modifications by reference and location, but
//! one row stands for every PSM of Scaffold's peptide name, and those PSMs
//! don't all carry the same N-terminus: in Scaffold DDA a `<229MDYWTR` row
//! lists Met-loss+Acetyl while most of its PSMs fit Met-loss alone, and
//! `<42...` rows list acetyl and TMT and their PSMs carry both. So each
//! peptide has *readings*: the one its name gives first, then the other
//! combinations of its listed terminal modifications, then (for a peptide
//! starting with M, when the file defines them) each with Met-loss or
//! Met-loss+Acetyl. Each PSM takes the first reading its precursor fits.

use std::collections::{HashMap, HashSet};

use rusqlite::Connection;

use super::SfdbError;
use crate::chemistry::peptide_mass;
use crate::model::{Modification, Peptide, Site};
use crate::unimod;

/// Isotope peaks a precursor may have been picked at (-1..+3), as engines
/// search them.
const ISOTOPE_SPACING: f64 = 1.003_354_835;
/// A PSM's mass fits its precursor within this, or [`MASS_PPM`], or 1.5x
/// the search's own precursor tolerance, whichever is largest. A misread
/// modification is off by at least ~1 Da (deamidation, 0.984 Da, can hide
/// behind an isotope peak; nothing else does).
const MASS_TOLERANCE: f64 = 0.05;
const MASS_PPM: f64 = 30.0;
/// Readings per peptide at most (terminal combinations multiply).
const MAX_READINGS: usize = 64;
/// Terminal modifications one row may list at most (2^6 combinations).
const MAX_TERMINAL: usize = 6;
/// Met-loss (Unimod 765) and Met-loss+Acetyl (766): they remove the
/// peptide's first residue rather than modify it.
const MET_LOSSES: [f64; 2] = [-131.040_485, -89.029_920];

/// One `Modification` row.
#[derive(Debug, Clone)]
pub(super) struct ModDef {
    /// The residue, or `n`/`c` for a terminus.
    residue: char,
    mass: f64,
    name: Option<String>,
}

impl ModDef {
    fn met_loss(&self) -> bool {
        MET_LOSSES
            .iter()
            .any(|mass| (self.mass - mass).abs() <= unimod::SNAP_TOLERANCE)
    }

    fn at(&self, site: Site) -> Modification {
        Modification {
            site,
            mass_delta: self.mass,
            name: self.name.clone(),
        }
    }
}

pub(super) fn modifications(conn: &Connection) -> Result<HashMap<String, ModDef>, SfdbError> {
    let mut statement = conn.prepare(
        "SELECT MODIFICATION_REF, AMINO_ACID, MASS, MODIFICATION_ACCESSION, SHORT_NAME, \
         FULL_NAME FROM Modification",
    )?;
    let mut rows = statement.query([])?;
    let mut out = HashMap::new();
    while let Some(row) = rows.next()? {
        let reference: String = row.get(0)?;
        let residue: String = row.get(1)?;
        let mass: f64 = row.get(2)?;
        let accession: Option<String> = row.get(3)?;
        let short: Option<String> = row.get(4)?;
        let full: Option<String> = row.get(5)?;
        let mut letters = residue.chars();
        let (Some(letter), None) = (letters.next(), letters.next()) else {
            return Err(SfdbError::Layout(format!(
                "modification {reference} is on residue {residue:?}"
            )));
        };
        if !mass.is_finite() {
            return Err(SfdbError::Layout(format!(
                "modification {reference} has no mass"
            )));
        }
        let name = short.or(full).filter(|name| !name.is_empty());
        out.insert(
            reference,
            ModDef {
                residue: letter,
                mass,
                name: naming(name, accession.as_deref(), mass),
            },
        );
    }
    Ok(out)
}

/// Unimod's own name where the file's name and mass agree with a known
/// modification (Scaffold writes `Phospho_STY`, and files terminal acetyls
/// under other Unimod numbers), or where the file's Unimod accession has
/// this mass; else the file's name.
pub(super) fn naming(name: Option<String>, accession: Option<&str>, mass: f64) -> Option<String> {
    if let Some(known) = name
        .as_deref()
        .and_then(|name| unimod::named_by_mass(name, mass))
    {
        return Some(known.to_string());
    }
    if let Some(known) = accession.and_then(|accession| unimod::named_by_accession(accession, mass))
    {
        return Some(known.to_string());
    }
    name
}

/// One way to read a peptide row.
#[derive(Debug, Clone, PartialEq)]
pub(super) struct Reading {
    pub peptide: Peptide,
    /// The initial methionine was lost (Met-loss, with or without acetyl).
    pub met_lost: bool,
    /// The reading Scaffold's peptide name gives (counted when another is
    /// taken).
    pub named: bool,
}

pub(super) type Readings = Result<Vec<Reading>, String>;

/// Peptide ref -> its readings, or why it can't be read.
pub(super) fn peptides(
    conn: &Connection,
    modifications: &HashMap<String, ModDef>,
    tables: &HashSet<String>,
) -> Result<HashMap<String, Readings>, SfdbError> {
    let mut variants = HashSet::new();
    if tables.contains("Substitution") {
        let mut statement = conn.prepare("SELECT DISTINCT PEPTIDE_REF FROM Substitution")?;
        variants = statement
            .query_map([], |row| row.get::<_, String>(0))?
            .collect::<rusqlite::Result<_>>()?;
    }
    let mut statement = conn.prepare(
        "SELECT PEPTIDE_REF, PEPTIDE_SEQUENCE, MODIFICATION_ID_ARRAY, PRIMARY_LOCATION_ARRAY \
         FROM Peptide",
    )?;
    let mut rows = statement.query([])?;
    let mut out: HashMap<String, Readings> = HashMap::new();
    let mut seen: HashMap<String, (String, String, String)> = HashMap::new();
    while let Some(row) = rows.next()? {
        let reference: String = row.get(0)?;
        let sequence: String = row.get(1)?;
        let ids: Option<String> = row.get(2)?;
        let locations: Option<String> = row.get(3)?;
        let row = (
            sequence,
            ids.unwrap_or_default(),
            locations.unwrap_or_default(),
        );
        if let Some(other) = seen.get(&reference) {
            if *other != row {
                out.insert(
                    reference,
                    Err("its peptide is stored twice, differently".into()),
                );
            }
            continue;
        }
        let readings = if variants.contains(&reference) {
            Err("its peptide carries a sequence variant".into())
        } else {
            readings(&reference, &row.0, &row.1, &row.2, modifications)
        };
        seen.insert(reference.clone(), row);
        out.insert(reference, readings);
    }
    Ok(out)
}

pub(super) fn readings(
    reference: &str,
    sequence: &str,
    ids: &str,
    locations: &str,
    modifications: &HashMap<String, ModDef>,
) -> Readings {
    let sequence = sequence.trim().to_ascii_uppercase();
    if sequence.is_empty() || !sequence.chars().all(|residue| residue.is_ascii_uppercase()) {
        return Err("its peptide sequence isn't one".into());
    }
    let residues: Vec<char> = sequence.chars().collect();
    let length = residues.len();
    let split = |text: &str| -> Vec<String> {
        text.split(',')
            .map(str::trim)
            .filter(|part| !part.is_empty())
            .map(String::from)
            .collect()
    };
    let (ids, locations) = (split(ids), split(locations));
    if ids.len() != locations.len() {
        return Err("its modification lists differ in length".into());
    }
    let mut n_term = Vec::new();
    let mut c_term = Vec::new();
    let mut on_residues = Vec::new();
    for (id, location) in ids.iter().zip(&locations) {
        let definition = modifications
            .get(id)
            .ok_or("it names an unknown modification")?;
        let location: usize = location
            .parse()
            .map_err(|_| "a modification location isn't a number")?;
        if definition.met_loss() && location <= 1 {
            continue; // offered to every M-first peptide below
        }
        match (definition.residue, location) {
            ('n', _) | (_, 0) => n_term.push(definition.at(Site::NTerm)),
            ('c', _) => c_term.push(definition.at(Site::CTerm)),
            (residue, location) if (1..=length).contains(&location) => {
                if residues[location - 1] != residue {
                    return Err("a modification sits on a residue it isn't defined for".into());
                }
                on_residues.push(definition.at(Site::Residue(location)));
            }
            _ => return Err("a modification lies outside its peptide".into()),
        }
    }
    let (n_named, c_named) = terminal_masses(reference);
    let n_options = terminal_options(n_term, n_named)?;
    let c_options = terminal_options(c_term, c_named)?;
    // Met-loss is offered for any M-first peptide in a file that defines it,
    // listed on this row or not (rows are shared by PSMs that differ). It is
    // kept as the engine reports it, an N-terminal delta on the M.
    let mut losses: Vec<&ModDef> = modifications
        .values()
        .filter(|definition| definition.met_loss())
        .collect();
    losses.sort_by(|a, b| a.mass.total_cmp(&b.mass));
    losses.dedup_by(|a, b| a.mass == b.mass);
    let losable = residues[0] == 'M' && length > 1;
    let loss_choices: Vec<Option<&ModDef>> = std::iter::once(None)
        .chain(losses.into_iter().filter(|_| losable).map(Some))
        .collect();
    let mut out: Vec<Reading> = Vec::new();
    for loss in &loss_choices {
        for (n, n_matches) in &n_options {
            for (c, c_matches) in &c_options {
                let mut found: Vec<Modification> =
                    n.iter().chain(c).chain(&on_residues).cloned().collect();
                found.extend(loss.map(|definition| definition.at(Site::NTerm)));
                // One order for one set: by site, then mass.
                found.sort_by(|a, b| {
                    a.site
                        .cmp(&b.site)
                        .then(a.mass_delta.total_cmp(&b.mass_delta))
                });
                let peptide = Peptide::new(sequence.clone(), found);
                if out.iter().all(|reading| reading.peptide != peptide) {
                    out.push(Reading {
                        peptide,
                        met_lost: loss.is_some(),
                        named: loss.is_none() && *n_matches && *c_matches,
                    });
                }
                if out.len() > MAX_READINGS {
                    return Err("its terminal modifications can be read too many ways".into());
                }
            }
        }
    }
    Ok(out)
}

/// The nominal terminal masses in Scaffold's peptide name: `<42MELSSM...`
/// and `...K>16` (never seen). `None` when the name shows none.
pub(super) fn terminal_masses(reference: &str) -> (Option<i64>, Option<i64>) {
    let leading = reference.strip_prefix('<').and_then(|rest| {
        let digits: String = rest
            .chars()
            .take_while(|c| c.is_ascii_digit() || *c == '-' || *c == '+')
            .collect();
        digits.parse::<i64>().ok()
    });
    let trailing = reference
        .rfind('>')
        .and_then(|at| reference[at + 1..].parse::<i64>().ok());
    (leading, trailing)
}

/// Every combination of a terminus's listed modifications, those whose
/// masses round to the name's (none, when the name shows no mass) first,
/// each with whether it matches the name.
pub(super) fn terminal_options(
    candidates: Vec<Modification>,
    named: Option<i64>,
) -> Result<Vec<(Vec<Modification>, bool)>, String> {
    if candidates.len() > MAX_TERMINAL {
        return Err("it lists too many terminal modifications".into());
    }
    let target = named.unwrap_or(0);
    let mut options: Vec<(Vec<Modification>, bool)> = (0u32..1 << candidates.len())
        .map(|mask| {
            let chosen: Vec<Modification> = candidates
                .iter()
                .enumerate()
                .filter(|(i, _)| mask & (1 << i) != 0)
                .map(|(_, modification)| modification.clone())
                .collect();
            let sum: f64 = chosen.iter().map(|m| m.mass_delta).sum();
            let matches = sum.round() as i64 == target;
            (chosen, matches)
        })
        .collect();
    // Stable: the name's reading first, then by mask.
    options.sort_by_key(|(_, matches)| !matches);
    Ok(options)
}

/// Neutral monoisotopic mass of a peptide with its modifications.
pub(super) fn peptide_neutral(peptide: &Peptide) -> Option<f64> {
    let deltas: f64 = peptide.modifications.iter().map(|m| m.mass_delta).sum();
    peptide_mass(&peptide.sequence).map(|mass| mass + deltas)
}

/// The search's precursor tolerance, as the file states it.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub(super) struct Tolerance {
    pub da: f64,
    pub ppm: f64,
}

impl Tolerance {
    /// The widest precursor tolerance any identification protocol states.
    pub fn read(conn: &Connection) -> Self {
        let Ok(mut statement) =
            conn.prepare("SELECT PARENT_TOLERANCE, PARENT_TOL_PPM FROM IdentificationProtocol")
        else {
            return Self::default();
        };
        let rows = statement.query_map([], |row| {
            Ok((row.get::<_, Option<f64>>(0)?, row.get::<_, Option<i64>>(1)?))
        });
        let mut tolerance = Self::default();
        for (value, ppm) in rows.into_iter().flatten().flatten() {
            let Some(value) = value.filter(|value| value.is_finite() && *value > 0.0) else {
                continue;
            };
            if ppm == Some(1) {
                tolerance.ppm = tolerance.ppm.max(value);
            } else {
                tolerance.da = tolerance.da.max(value);
            }
        }
        tolerance
    }

    /// Whether a peptide's neutral mass lands on an isotope peak (-1..+3)
    /// of the precursor's.
    pub fn fits(self, calculated: f64, precursor: f64) -> bool {
        let tolerance = MASS_TOLERANCE
            .max(precursor.abs() * MASS_PPM * 1e-6)
            .max(1.5 * self.da)
            .max(1.5 * precursor.abs() * self.ppm * 1e-6);
        (-1..=3).any(|isotope| {
            (precursor - calculated - f64::from(isotope) * ISOTOPE_SPACING).abs() <= tolerance
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::PROTON;

    fn definition(residue: char, mass: f64, name: &str) -> ModDef {
        ModDef {
            residue,
            mass,
            name: Some(name.into()),
        }
    }

    fn table() -> HashMap<String, ModDef> {
        HashMap::from([
            ("tmt_n".into(), definition('n', 229.162932, "TMT6plex")),
            ("ac_n".into(), definition('n', 42.010565, "Acetyl")),
            ("tmt_k".into(), definition('K', 229.162932, "TMT6plex")),
            ("loss".into(), definition('n', -131.040485, "Met-loss")),
            ("ox".into(), definition('M', 15.994915, "Oxidation")),
        ])
    }

    fn masses(reading: &Reading) -> Vec<f64> {
        reading
            .peptide
            .modifications
            .iter()
            .map(|m| m.mass_delta)
            .collect()
    }

    #[test]
    fn terminal_masses_come_from_the_name() {
        assert_eq!(terminal_masses("<42MELSSM+16K+229"), (Some(42), None));
        assert_eq!(terminal_masses("<229SETPFHSR"), (Some(229), None));
        assert_eq!(terminal_masses("VALEAC+57VK"), (None, None));
    }

    #[test]
    fn the_name_gives_the_first_reading() {
        let readings = readings(
            "<229PEPTIDEK+229",
            "PEPTIDEK",
            "tmt_n,ac_n,tmt_k",
            "0,0,8",
            &table(),
        )
        .unwrap();
        let first = &readings[0];
        assert!(first.named && !first.met_lost);
        assert_eq!(masses(first), [229.162932, 229.162932]);
        // The other combinations follow: none, acetyl, both (on top of K's).
        assert_eq!(readings.len(), 4);
        assert!(readings[1..].iter().all(|reading| !reading.named));
    }

    #[test]
    fn methionine_loss_is_an_n_terminal_delta() {
        let all = readings("<229MSYAPR", "MSYAPR", "tmt_n,loss", "0,0", &table()).unwrap();
        let lost: Vec<&Reading> = all.iter().filter(|r| r.met_lost).collect();
        // Two terminal choices (TMT or not), each with Met-loss.
        assert_eq!(lost.len(), 2);
        assert!(lost.iter().all(|r| r.peptide.sequence == "MSYAPR"));
        assert_eq!(masses(lost[0]), [-131.040485, 229.162932]);
        // A peptide not starting with M has no such reading.
        let none = readings("PEPK", "PEPK", "", "", &table()).unwrap();
        assert!(none.iter().all(|r| !r.met_lost));
    }

    #[test]
    fn bad_rows_say_why() {
        assert!(readings("X", "PEPK", "nope", "1", &table())
            .unwrap_err()
            .contains("unknown modification"));
        assert!(readings("X", "PEPK", "ox", "1", &table())
            .unwrap_err()
            .contains("residue it isn't defined for"));
        assert!(readings("X", "PEPK", "ox", "9", &table())
            .unwrap_err()
            .contains("outside"));
    }

    #[test]
    fn masses_fit_on_isotope_peaks() {
        // TTEQVNAAM+16AAC+57K 2+ at m/z 705.82019 (a PSM of PXD038643). By
        // hand: T 2x101.04768, E 129.04259, Q 128.05858, V 99.06841,
        // N 114.04293, A 4x71.03711, M 131.04049, C 103.00919, K 128.09496
        // = 1318.60095; + water 18.010565 + 15.994915 + 57.021464
        // = 1409.627894; the precursor (705.82019 - 1.007276) x 2
        // = 1409.625828, 0.002 Da away.
        let peptide = Peptide::new(
            "TTEQVNAAMAACK".into(),
            vec![
                Modification {
                    site: Site::Residue(9),
                    mass_delta: 15.994915,
                    name: None,
                },
                Modification {
                    site: Site::Residue(12),
                    mass_delta: 57.021464,
                    name: None,
                },
            ],
        );
        let mass = peptide_neutral(&peptide).unwrap();
        let precursor = (705.82019 - PROTON) * 2.0;
        let tolerance = Tolerance::default();
        assert!(tolerance.fits(mass, precursor), "{mass} vs {precursor}");
        assert!(!tolerance.fits(mass + 15.994915, precursor));
        assert!(tolerance.fits(mass, precursor + ISOTOPE_SPACING));
        // An ion-trap search's 0.5 Da window widens it to 0.75 Da.
        let wide = Tolerance { da: 0.5, ppm: 0.0 };
        assert!(!tolerance.fits(mass + 0.4, precursor));
        assert!(wide.fits(mass + 0.4, precursor));
    }

    #[test]
    fn scaffold_names_map_to_unimod_where_mass_and_name_agree() {
        assert_eq!(
            naming(
                Some("Phospho_STY".into()),
                Some("Phospho_STY"),
                79.9663314819336
            ),
            Some("Phospho".into())
        );
        // Scaffold DDA files terminal acetyls under other Unimod numbers.
        assert_eq!(
            naming(Some("Acetyl".into()), Some("UNIMOD:1047"), 42.010565),
            Some("Acetyl".into())
        );
        assert_eq!(
            naming(
                Some("Carbamidomethylation".into()),
                Some("UNIMOD:4"),
                57.021464
            ),
            Some("Carbamidomethyl".into())
        );
        assert_eq!(
            naming(Some("Mystery".into()), Some("UNIMOD:35"), 12.5),
            Some("Mystery".into())
        );
    }
}
