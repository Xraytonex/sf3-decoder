//! `PeptideIdentification` records: one engine's candidate peptide for one
//! spectrum, with its modifications, scores and Scaffold's probability.

use std::collections::HashMap;

use super::project::ModificationType;
use super::records::{Records, read_f32, read_i32};
use crate::chemistry::{mass_error, peptide_mass};
use crate::model::{Modification, Peptide, Score, Site};

pub(crate) const PEPTIDE_BLOCK: usize = 160;
const SEQUENCE: std::ops::Range<usize> = 4..34;
const SEQUENCE_END: u8 = 0x0a;
const TAIL: usize = 35;
const PROBABILITY: usize = 49;
const ANALYSIS: usize = 57;
const MOD_TYPES: [usize; 3] = [90, 94, 98];
const MOD_POSITIONS: [usize; 3] = [106, 110, 114];
const MORE_MOD_TYPES: usize = 102;
const MORE_MOD_POSITIONS: usize = 118;
const PROTEIN: usize = 122;
/// 0 marks the candidate Scaffold picked for its spectrum.
const CHOICE: usize = 147;
/// Unused integer and float slots.
const NONE: i32 = 0x7fff_ffff;
/// Longer "peptides" are a misread, not a peptide.
const MAX_PEPTIDE_LENGTH: usize = 1_000;

/// (offset, name, PSI-MS accession)
type Field = (usize, &'static str, Option<&'static str>);

/// Where an engine's scores sit in a record. Scaffold names the floats at
/// `66 + 4k` itself (`AnalysisProgram.scoreTypes`); its "Discriminant
/// Score" is Scaffold's own and left out.
pub(crate) struct EngineScores {
    /// Ranks the candidates of a spectrum.
    pub primary: usize,
    pub fields: &'static [Field],
    /// Names rest on value ranges only.
    pub inferred: bool,
    pub expectation: Option<Expectation>,
}

#[derive(Clone, Copy)]
pub(crate) enum Expectation {
    /// Mascot's follows from its ion score and identity threshold.
    Mascot { score: usize, identity: usize },
    /// Stored as −log10 of the value (MSFragger).
    NegLog10 { offset: usize, name: &'static str },
}

impl Expectation {
    fn score(self, block: &[u8]) -> Option<Score> {
        let (name, accession, value) = match self {
            Self::Mascot { score, identity } => (
                "Mascot:expectation value",
                Some("MS:1001172"),
                expectation(
                    read_f32(block, score)?.into(),
                    read_f32(block, identity)?.into(),
                ),
            ),
            Self::NegLog10 { offset, name } => {
                (name, None, 10f64.powf(-f64::from(read_f32(block, offset)?)))
            }
        };
        value.is_finite().then_some(Score {
            name,
            accession,
            value,
        })
    }
}

/// Mascot's expectation value from its ion score and identity threshold:
/// 0.05 at the threshold, tenfold per 10 points.
pub(crate) fn expectation(score: f64, identity: f64) -> f64 {
    0.05 * 10f64.powf((identity - score) / 10.0)
}

pub(crate) fn engine_scores(program: &str) -> Option<EngineScores> {
    Some(match program.to_ascii_lowercase().as_str() {
        // "Mascot (Ion Score Only)" differs only at 66, which isn't read.
        "mascot" | "mascot (ion score only)" => EngineScores {
            primary: 70,
            fields: &[
                (70, "Mascot:score", Some("MS:1001171")),
                (74, "Mascot:identity threshold", Some("MS:1001371")),
                (78, "Mascot:homology threshold", Some("MS:1001370")),
            ],
            inferred: false,
            expectation: Some(Expectation::Mascot {
                score: 70,
                identity: 74,
            }),
        },
        "maxquant" => EngineScores {
            primary: 66,
            fields: &[
                (66, "Andromeda:score", Some("MS:1002338")),
                (70, "MaxQuant:delta score", None),
            ],
            inferred: false,
            expectation: None,
        },
        "spectrum mill" => EngineScores {
            primary: 70,
            fields: &[
                (70, "SpectrumMill:Score", Some("MS:1001572")),
                (74, "SpectrumMill:SPI", Some("MS:1001573")),
            ],
            inferred: false,
            expectation: None,
        },
        "sequest" => EngineScores {
            primary: 74,
            fields: &[
                (74, "SEQUEST:xcorr", Some("MS:1001155")),
                (70, "SEQUEST:deltacn", Some("MS:1001156")),
                (78, "SEQUEST:sp", Some("MS:1001157")),
            ],
            inferred: false,
            expectation: None,
        },
        "msfragger" => EngineScores {
            primary: 70,
            fields: &[
                (70, "MSFragger:hyperscore", None),
                (78, "MSFragger:nextscore", None),
            ],
            inferred: false,
            expectation: Some(Expectation::NegLog10 {
                offset: 74,
                name: "MSFragger:expectation value",
            }),
        },
        // Scaffold's list says 66 and 74 are −log expectation values and 70
        // the hyperscore, but 70 ranks targets better when lower, which no
        // hyperscore does; 66 separates targets from decoys best.
        "x! tandem" | "x!tandem" => EngineScores {
            primary: 66,
            fields: &[
                (66, "X!Tandem:sf3 offset 66", None),
                (70, "X!Tandem:sf3 offset 70", None),
                (74, "X!Tandem:sf3 offset 74", None),
            ],
            inferred: true,
            expectation: None,
        },
        _ => return None,
    })
}

/// A decoded record, not yet checked against its spectrum.
pub(crate) struct Candidate {
    pub decoded: Decoded,
    pub protein_record: i32,
    pub primary: f32,
    pub scores: Vec<Score>,
    pub probability: Option<f32>,
    pub chosen: bool,
    /// Residues before and after the peptide (ASCII).
    pub flanks: [u8; 2],
}

/// A sequence and its stated modifications. Spectrum Mill writes a
/// modified residue in lower case without naming the modification; those
/// positions are `implied`.
pub(crate) struct Decoded {
    /// Upper case.
    pub sequence: String,
    pub modifications: Vec<Modification>,
    pub implied: Vec<usize>,
}

/// Modifications a lower-case residue may stand for: the types the file
/// states on that residue elsewhere.
pub(crate) type ResidueOptions = HashMap<char, Vec<(String, f64)>>;

impl Decoded {
    fn stated_mass(&self) -> Option<f64> {
        peptide_mass(&self.sequence).map(|mass| {
            mass + self
                .modifications
                .iter()
                .map(|modification| modification.mass_delta)
                .sum::<f64>()
        })
    }

    /// Settle implied modifications: every lower-case occurrence of a
    /// residue takes the same option, and exactly one combination may fit
    /// the spectrum's [M+H]+ within `tolerance`. Returns the peptide and its
    /// neutral mass.
    pub fn resolve(
        self,
        mh: f64,
        options: &ResidueOptions,
        tolerance: f64,
    ) -> Option<(Peptide, Option<f64>)> {
        let mass = self.stated_mass();
        if self.implied.is_empty() {
            return Some((Peptide::new(self.sequence, self.modifications), mass));
        }
        let mass = mass?;
        let residues: Vec<char> = self.sequence.chars().collect();
        let mut letters: Vec<char> = self.implied.iter().map(|&index| residues[index]).collect();
        letters.sort_unstable();
        letters.dedup();
        let letter_of = |index: usize| {
            letters
                .iter()
                .position(|&letter| letter == residues[index])
                .expect("listed")
        };
        let choices: Vec<&[(String, f64)]> = letters
            .iter()
            .map(|letter| options.get(letter).map_or(&[][..], Vec::as_slice))
            .collect();
        let combinations = choices
            .iter()
            .try_fold(1_usize, |product, choice| product.checked_mul(choice.len()))
            .unwrap_or(usize::MAX);
        if combinations == 0 || combinations > 64 {
            return None;
        }
        let mut fitting = Vec::new();
        for mut combination in 0..combinations {
            let picked: Vec<&(String, f64)> = choices
                .iter()
                .map(|choice| {
                    let pick = &choice[combination % choice.len()];
                    combination /= choice.len();
                    pick
                })
                .collect();
            let extra: f64 = self
                .implied
                .iter()
                .map(|&index| picked[letter_of(index)].1)
                .sum();
            if mass_error(mass + extra, mh).abs() <= tolerance {
                fitting.push((picked, mass + extra));
            }
        }
        let [(picked, total)] = fitting.as_slice() else {
            return None;
        };
        let mut modifications = self.modifications;
        for &index in &self.implied {
            let (name, delta) = picked[letter_of(index)];
            modifications.push(Modification {
                site: Site::Residue(index + 1),
                mass_delta: *delta,
                name: Some(name.clone()).filter(|name| !name.is_empty()),
            });
        }
        Some((Peptide::new(self.sequence, modifications), Some(*total)))
    }
}

/// A record whose layout isn't the one known; the file is refused.
#[derive(Debug)]
pub(crate) struct Undecodable(pub String);

pub(crate) struct Peptides<'a> {
    pub tails: &'a Records,
    pub type_lists: &'a Records,
    pub position_lists: &'a Records,
    pub types: &'a HashMap<i32, ModificationType>,
}

impl Peptides<'_> {
    /// The sequence as stored (lower case where Spectrum Mill marks a
    /// modified residue).
    pub fn sequence(&self, block: &[u8]) -> Result<String, Undecodable> {
        let inline = &block[SEQUENCE];
        let end = inline
            .iter()
            .position(|&byte| byte == SEQUENCE_END)
            .unwrap_or(inline.len());
        let mut sequence = String::from_utf8_lossy(&inline[..end]).into_owned();
        let tail = read_i32(block, TAIL).unwrap_or(0);
        if tail != 0 {
            let rest = self
                .tails
                .string(tail)
                .ok_or_else(|| Undecodable(format!("sequence tail {tail} is missing")))?;
            sequence.push_str(&rest);
        }
        if sequence.is_empty()
            || sequence.len() > MAX_PEPTIDE_LENGTH
            || !sequence.bytes().all(|byte| byte.is_ascii_alphabetic())
        {
            return Err(Undecodable(format!(
                "sequence {sequence:?} is not residues"
            )));
        }
        Ok(sequence)
    }

    /// (0-based position, type id) pairs; a later one at the same position
    /// replaces an earlier one (e.g. carbamidomethyl, then the NEM that
    /// replaced it).
    fn modification_ids(&self, block: &[u8]) -> Result<Vec<(i32, i32)>, Undecodable> {
        let mut pairs: Vec<(i32, i32)> = MOD_TYPES
            .iter()
            .zip(MOD_POSITIONS)
            .filter_map(|(&kind, position)| {
                let kind = read_i32(block, kind)?;
                (kind != NONE).then(|| (read_i32(block, position).unwrap_or(NONE), kind))
            })
            .collect();
        let (more_types, more_positions) = (
            read_i32(block, MORE_MOD_TYPES).unwrap_or(0),
            read_i32(block, MORE_MOD_POSITIONS).unwrap_or(0),
        );
        if more_types != 0 || more_positions != 0 {
            match (
                self.type_lists.int_list(more_types),
                self.position_lists.int_list(more_positions),
            ) {
                (Some(types), Some(positions)) if types.len() == positions.len() => {
                    pairs.extend(positions.into_iter().zip(types));
                }
                _ => {
                    return Err(Undecodable(format!(
                        "modification lists {more_types}/{more_positions} are missing or unequal"
                    )));
                }
            }
        }
        let mut by_position: Vec<(i32, i32)> = Vec::with_capacity(pairs.len());
        for (position, kind) in pairs {
            match by_position.iter_mut().find(|(known, _)| *known == position) {
                Some(slot) => slot.1 = kind,
                None => by_position.push((position, kind)),
            }
        }
        Ok(by_position)
    }

    pub fn decode(&self, block: &[u8]) -> Result<Decoded, Undecodable> {
        let stored = self.sequence(block)?;
        let sequence = stored.to_ascii_uppercase();
        let residues: Vec<char> = sequence.chars().collect();
        let mut modifications = Vec::new();
        let mut stated = Vec::new();
        for (position, kind) in self.modification_ids(block)? {
            let modification = self
                .types
                .get(&kind)
                .ok_or_else(|| Undecodable(format!("unknown modification type {kind}")))?;
            let index = usize::try_from(position)
                .ok()
                .filter(|&index| index < residues.len())
                .ok_or_else(|| {
                    Undecodable(format!(
                        "{} at {position} is outside {sequence}",
                        modification.name
                    ))
                })?;
            let site = if modification.residues.contains(residues[index]) {
                stated.push(index);
                Site::Residue(index + 1)
            } else if index == 0 && modification.residues.contains('n') {
                Site::NTerm
            } else if index + 1 == residues.len() && modification.residues.contains('c') {
                Site::CTerm
            } else {
                return Err(Undecodable(format!(
                    "{} cannot sit on {} of {sequence}",
                    modification.name, residues[index]
                )));
            };
            modifications.push(Modification {
                site,
                mass_delta: modification.mono_mass,
                name: Some(modification.name.clone()).filter(|name| !name.is_empty()),
            });
        }
        let implied = stored
            .char_indices()
            .filter(|(index, residue)| residue.is_ascii_lowercase() && !stated.contains(index))
            .map(|(index, _)| index)
            .collect();
        Ok(Decoded {
            sequence,
            modifications,
            implied,
        })
    }

    /// The modification types the file states on each residue.
    pub fn residue_options(&self, records: &Records) -> ResidueOptions {
        let mut options = ResidueOptions::new();
        for block in records.ids().filter_map(|id| records.get(id)) {
            let (Ok(sequence), Ok(ids)) = (self.sequence(block), self.modification_ids(block))
            else {
                continue;
            };
            let residues: Vec<char> = sequence.to_ascii_uppercase().chars().collect();
            for (position, kind) in ids {
                let (Some(&residue), Some(modification)) = (
                    usize::try_from(position)
                        .ok()
                        .and_then(|index| residues.get(index)),
                    self.types.get(&kind),
                ) else {
                    continue;
                };
                if modification.residues.contains(residue) {
                    let known = options.entry(residue).or_default();
                    if !known.iter().any(|(name, _)| *name == modification.name) {
                        known.push((modification.name.clone(), modification.mono_mass));
                    }
                }
            }
        }
        options
    }

    pub fn candidate(&self, block: &[u8], engine: &EngineScores) -> Result<Candidate, Undecodable> {
        let scores = engine
            .fields
            .iter()
            .filter_map(|&(offset, name, accession)| {
                let value = f64::from(read_f32(block, offset)?);
                value.is_finite().then_some(Score {
                    name,
                    accession,
                    value,
                })
            })
            .chain(engine.expectation.and_then(|kind| kind.score(block)))
            .collect();
        Ok(Candidate {
            decoded: self.decode(block)?,
            protein_record: read_i32(block, PROTEIN).unwrap_or(0),
            primary: read_f32(block, engine.primary).unwrap_or(f32::NAN),
            scores,
            probability: read_f32(block, PROBABILITY)
                .filter(|probability| (0.0..=1.0).contains(probability)),
            chosen: block[CHOICE] == 0,
            flanks: [block[39], block[40]],
        })
    }
}

/// The spectrum analysis a record belongs to.
pub(crate) fn analysis_of(block: &[u8]) -> i32 {
    read_i32(block, ANALYSIS).unwrap_or(0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::PROTON;

    #[test]
    fn mascot_expectation_is_005_at_the_identity_threshold() {
        assert!((expectation(35.0, 35.0) - 0.05).abs() < 1e-12);
        assert!((expectation(45.0, 35.0) - 0.005).abs() < 1e-12);
        assert!((expectation(15.0, 35.0) - 5.0).abs() < 1e-9);
    }

    fn decoded(sequence: &str, implied: &[usize]) -> Decoded {
        Decoded {
            sequence: sequence.into(),
            modifications: Vec::new(),
            implied: implied.to_vec(),
        }
    }

    #[test]
    fn lower_case_residues_take_the_one_modification_that_fits() {
        let options: ResidueOptions = [(
            'C',
            vec![
                ("Carbamidomethyl".to_string(), 57.021_464),
                ("Nethylmaleimide".to_string(), 125.047_679),
            ],
        )]
        .into();
        let mass = peptide_mass("PEPCK").unwrap();
        let mh = mass + 125.047_679 + PROTON;
        let (peptide, total) = decoded("PEPCK", &[3]).resolve(mh, &options, 1.0).unwrap();
        assert_eq!(peptide.key(), "PEPC[Nethylmaleimide]K");
        assert!((total.unwrap() - mass - 125.047_679).abs() < 1e-9);
        assert!(
            decoded("PEPCK", &[3])
                .resolve(mh + 30.0, &options, 1.0)
                .is_none()
        );
        assert!(decoded("PEPMK", &[3]).resolve(mh, &options, 1.0).is_none());
        // Two options within the tolerance of each other are ambiguous.
        let close: ResidueOptions = [('C', vec![("A".into(), 57.0), ("B".into(), 57.5)])].into();
        let mh = mass + 57.2 + PROTON;
        assert!(decoded("PEPCK", &[3]).resolve(mh, &close, 1.0).is_none());
    }

    #[test]
    fn known_engines_have_scores_and_others_do_not() {
        assert_eq!(engine_scores("Mascot").unwrap().primary, 70);
        assert_eq!(
            engine_scores("Mascot (Ion Score Only)").unwrap().primary,
            70
        );
        assert!(engine_scores("X! Tandem").unwrap().inferred);
        // Scaffold's list for Sequest: Discriminant Score, deltaCn, XCorr,
        // Sp at 66, 70, 74, 78; the first is Scaffold's.
        let sequest = engine_scores("Sequest").unwrap();
        assert_eq!(sequest.primary, 74);
        assert_eq!(
            sequest
                .fields
                .iter()
                .map(|field| field.0)
                .collect::<Vec<_>>(),
            [74, 70, 78]
        );
        assert!(engine_scores("Byonic").is_none());
    }
}
