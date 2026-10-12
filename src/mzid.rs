//! mzIdentML 1.2.
//!
//! - One `SpectrumIdentificationList` per search engine, holding every
//!   candidate it kept for each spectrum, ranked by the engine's own score,
//!   with those scores and Scaffold's peptide probability (MS:1001568).
//! - `passThreshold` is Scaffold's verdict: the candidate it picked for the
//!   spectrum, at or above the project's saved peptide probability.
//! - Spectra are `index=N` in their run's MGF, as `sf3decode … mgf` writes it.
//! - Scaffold's protein groups make the `ProteinDetectionList`, with its
//!   protein probability (MS:1001579; per run as userParams). The schema
//!   wants peptides behind every protein, so each is backed by the
//!   candidates Scaffold picked; proteins without one are left out.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::fmt::Write as _;
use std::io::{self, Write};

use crate::chemistry::PROTON;
use crate::model::{Match, Modification, Project, Site, Spectrum, Thresholds};
use crate::output::mgf_names;
use crate::unimod;

/// Probabilities are f32 in records and f64 in the saved thresholds.
const PROBABILITY_SLACK: f64 = 1e-6;

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct MzidReport {
    pub candidates: usize,
    pub groups: usize,
    /// Groups without a candidate Scaffold picked, left out.
    pub groups_without_peptides: usize,
    /// Modifications written as "unknown modification" with their name.
    pub unknown_modifications: usize,
    /// Protein sequences with characters the schema doesn't allow (a `*`
    /// inside them).
    pub sequences_left_out: usize,
}

pub fn write_mzid(project: &Project, out: &mut impl Write) -> io::Result<MzidReport> {
    let index = Index::new(project);
    let mut writer = Writer {
        out,
        project,
        index,
        report: MzidReport::default(),
    };
    writer.document()?;
    Ok(writer.report)
}

type EvidenceKey<'a> = (usize, &'a str, Option<char>, Option<char>);

/// Ids for everything referenced across the file.
struct Index<'a> {
    peptides: Vec<&'a Match>,
    peptide_ids: HashMap<String, usize>,
    proteins: HashMap<&'a str, usize>,
    evidence: Vec<EvidenceKey<'a>>,
    evidence_ids: HashMap<EvidenceKey<'a>, usize>,
    /// Each spectrum's position in its run's MGF.
    positions: Vec<usize>,
    /// Engines with at least one candidate.
    engines: Vec<usize>,
}

impl<'a> Index<'a> {
    fn new(project: &'a Project) -> Self {
        let proteins = project
            .proteins
            .keys()
            .enumerate()
            .map(|(id, accession)| (accession.as_str(), id))
            .collect();
        let mut index = Self {
            peptides: Vec::new(),
            peptide_ids: HashMap::new(),
            proteins,
            evidence: Vec::new(),
            evidence_ids: HashMap::new(),
            positions: Vec::with_capacity(project.spectra.len()),
            engines: Vec::new(),
        };
        let mut next_position = vec![0; project.runs.len()];
        let mut engines = HashSet::new();
        for spectrum in &project.spectra {
            index.positions.push(next_position[spectrum.run]);
            next_position[spectrum.run] += 1;
            for found in &spectrum.matches {
                engines.insert(found.engine);
                let key = found.peptide.key();
                let peptide = match index.peptide_ids.get(&key) {
                    Some(&id) => id,
                    None => {
                        index.peptides.push(found);
                        index.peptide_ids.insert(key, index.peptides.len() - 1);
                        index.peptides.len() - 1
                    }
                };
                for protein in &found.proteins {
                    let key = (
                        peptide,
                        protein.accession.as_str(),
                        protein.pre,
                        protein.post,
                    );
                    if !index.evidence_ids.contains_key(&key) {
                        index.evidence.push(key);
                        index.evidence_ids.insert(key, index.evidence.len() - 1);
                    }
                }
            }
        }
        index.engines = engines.into_iter().collect();
        index.engines.sort_unstable();
        index
    }

    fn evidence_of(&self, found: &'a Match) -> Vec<usize> {
        let peptide = self.peptide_ids[&found.peptide.key()];
        found
            .proteins
            .iter()
            .map(|protein| {
                self.evidence_ids[&(
                    peptide,
                    protein.accession.as_str(),
                    protein.pre,
                    protein.post,
                )]
            })
            .collect()
    }
}

struct Writer<'a, W: Write> {
    out: &'a mut W,
    project: &'a Project,
    index: Index<'a>,
    report: MzidReport,
}

impl<'a, W: Write> Writer<'a, W> {
    fn put(&mut self, line: &str) -> io::Result<()> {
        self.out.write_all(line.as_bytes())?;
        self.out.write_all(b"\n")
    }

    fn document(&mut self) -> io::Result<()> {
        self.put(r#"<?xml version="1.0" encoding="UTF-8"?>"#)?;
        self.put(
            r#"<MzIdentML id="sf3decode" version="1.2.0" xmlns="http://psidev.info/psi/pi/mzIdentML/1.2" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://psidev.info/psi/pi/mzIdentML/1.2 https://www.psidev.info/files/mzIdentML1.2.0.xsd">"#,
        )?;
        self.put("  <cvList>")?;
        self.put(r#"    <cv id="PSI-MS" fullName="PSI-MS" uri="https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo"/>"#)?;
        self.put(
            r#"    <cv id="UNIMOD" fullName="UNIMOD" uri="http://www.unimod.org/obo/unimod.obo"/>"#,
        )?;
        self.put(r#"    <cv id="UO" fullName="UNIT-ONTOLOGY" uri="https://raw.githubusercontent.com/bio-ontology-research-group/unit-ontology/master/unit.obo"/>"#)?;
        self.put("  </cvList>")?;
        self.software()?;
        self.sequences()?;
        self.analyses()?;
        self.protocols()?;
        self.put("  <DataCollection>")?;
        self.inputs()?;
        self.put("    <AnalysisData>")?;
        for engine in self.index.engines.clone() {
            self.spectra(engine)?;
        }
        self.groups()?;
        self.put("    </AnalysisData>")?;
        self.put("  </DataCollection>")?;
        self.put("</MzIdentML>")
    }

    fn software(&mut self) -> io::Result<()> {
        self.put("  <AnalysisSoftwareList>")?;
        for engine in self.index.engines.clone() {
            let engine_info = &self.project.engines[engine];
            let version = engine_info
                .version
                .as_deref()
                .map_or_else(String::new, |v| format!(r#" version="{}""#, escape(v)));
            let name = match software_term(&engine_info.name) {
                Some((accession, name)) => {
                    format!(r#"<cvParam cvRef="PSI-MS" accession="{accession}" name="{name}"/>"#)
                }
                None => format!(r#"<userParam name="{}"/>"#, escape(&engine_info.name)),
            };
            self.put(&format!(
                r#"    <AnalysisSoftware id="AS_{engine}" name="{}"{version}>"#,
                escape(&engine_info.name)
            ))?;
            self.put(&format!("      <SoftwareName>{name}</SoftwareName>"))?;
            self.put("    </AnalysisSoftware>")?;
        }
        let version = self
            .project
            .scaffold_version
            .as_deref()
            .map_or_else(String::new, |v| format!(r#" version="{}""#, escape(v)));
        self.put(&format!(
            r#"    <AnalysisSoftware id="AS_scaffold" name="Scaffold"{version}>"#
        ))?;
        self.put(r#"      <SoftwareName><cvParam cvRef="PSI-MS" accession="MS:1001561" name="Scaffold"/></SoftwareName>"#)?;
        self.put("    </AnalysisSoftware>")?;
        self.put(&format!(
            r#"    <AnalysisSoftware id="AS_sf3decode" name="sf3decode" version="{}" uri="https://github.com/Xraytonex/sf3-decoder">"#,
            env!("CARGO_PKG_VERSION")
        ))?;
        self.put(r#"      <SoftwareName><userParam name="sf3decode"/></SoftwareName>"#)?;
        self.put("    </AnalysisSoftware>")?;
        self.put("  </AnalysisSoftwareList>")
    }

    fn sequences(&mut self) -> io::Result<()> {
        self.put("  <SequenceCollection>")?;
        for (id, (accession, protein)) in self.project.proteins.iter().enumerate() {
            // The schema takes residue letters only. A trailing `*` marks the
            // protein's end in some FASTA files; one inside can't be written.
            let sequence = protein
                .sequence
                .as_deref()
                .map(|sequence| sequence.trim_end_matches('*'))
                .filter(|sequence| {
                    let valid = sequence.bytes().all(|byte| byte.is_ascii_uppercase());
                    self.report.sequences_left_out += usize::from(!valid);
                    valid
                });
            let length = sequence.map_or_else(String::new, |s| format!(r#" length="{}""#, s.len()));
            self.put(&format!(
                r#"    <DBSequence id="DBSeq_{id}" accession="{}" searchDatabase_ref="SDB"{length}>"#,
                escape(accession)
            ))?;
            if let Some(sequence) = sequence {
                self.put(&format!("      <Seq>{sequence}</Seq>"))?;
            }
            if let Some(description) = &protein.description {
                self.put(&format!(
                    r#"      <cvParam cvRef="PSI-MS" accession="MS:1001088" name="protein description" value="{}"/>"#,
                    escape(description)
                ))?;
            }
            self.put("    </DBSequence>")?;
        }
        for id in 0..self.index.peptides.len() {
            let peptide = &self.index.peptides[id].peptide;
            self.put(&format!(r#"    <Peptide id="Pep_{id}">"#))?;
            self.put(&format!(
                "      <PeptideSequence>{}</PeptideSequence>",
                escape(&peptide.sequence)
            ))?;
            let length = peptide.sequence.len();
            for modification in &peptide.modifications {
                let (xml, known) = modification_xml(modification, length);
                self.report.unknown_modifications += usize::from(!known);
                self.put(&xml)?;
            }
            self.put("    </Peptide>")?;
        }
        for id in 0..self.index.evidence.len() {
            let (peptide, accession, pre, post) = self.index.evidence[id];
            let protein = self.index.proteins[accession];
            let mut attributes = format!(
                r#"id="PE_{id}" peptide_ref="Pep_{peptide}" dBSequence_ref="DBSeq_{protein}""#
            );
            for (name, residue) in [("pre", pre), ("post", post)] {
                if let Some(residue) = residue.filter(|r| r.is_ascii_uppercase() || *r == '-') {
                    let _ = write!(attributes, r#" {name}="{residue}""#);
                }
            }
            let decoy = self.project.proteins[accession].decoy;
            self.put(&format!(
                r#"    <PeptideEvidence {attributes} isDecoy="{decoy}"/>"#
            ))?;
        }
        self.put("  </SequenceCollection>")
    }

    fn analyses(&mut self) -> io::Result<()> {
        self.put("  <AnalysisCollection>")?;
        for engine in self.index.engines.clone() {
            self.put(&format!(
                r#"    <SpectrumIdentification id="SI_{engine}" spectrumIdentificationProtocol_ref="SIP_{engine}" spectrumIdentificationList_ref="SIL_{engine}">"#
            ))?;
            for run in 0..self.project.runs.len() {
                self.put(&format!(
                    r#"      <InputSpectra spectraData_ref="SD_{run}"/>"#
                ))?;
            }
            self.put(r#"      <SearchDatabaseRef searchDatabase_ref="SDB"/>"#)?;
            self.put("    </SpectrumIdentification>")?;
        }
        self.put(r#"    <ProteinDetection id="PD" proteinDetectionProtocol_ref="PDP" proteinDetectionList_ref="PDL">"#)?;
        for engine in self.index.engines.clone() {
            self.put(&format!(
                r#"      <InputSpectrumIdentifications spectrumIdentificationList_ref="SIL_{engine}"/>"#
            ))?;
        }
        self.put("    </ProteinDetection>")?;
        self.put("  </AnalysisCollection>")
    }

    fn protocols(&mut self) -> io::Result<()> {
        const NO_THRESHOLD: &str =
            r#"<cvParam cvRef="PSI-MS" accession="MS:1001494" name="no threshold"/>"#;
        let thresholds = self.project.thresholds;
        let peptide = thresholds.map_or_else(
            || NO_THRESHOLD.to_string(),
            |t| {
                format!(
                    r#"<cvParam cvRef="PSI-MS" accession="MS:1001577" name="Scaffold: Minimum Peptide Probability" value="{}"/>"#,
                    t.peptide_probability
                )
            },
        );
        self.put("  <AnalysisProtocolCollection>")?;
        for engine in self.index.engines.clone() {
            self.put(&format!(
                r#"    <SpectrumIdentificationProtocol id="SIP_{engine}" analysisSoftware_ref="AS_{engine}">"#
            ))?;
            self.put(r#"      <SearchType><cvParam cvRef="PSI-MS" accession="MS:1001083" name="ms-ms search"/></SearchType>"#)?;
            self.put(&format!("      <Threshold>{peptide}</Threshold>"))?;
            self.put("    </SpectrumIdentificationProtocol>")?;
        }
        self.put(r#"    <ProteinDetectionProtocol id="PDP" analysisSoftware_ref="AS_scaffold">"#)?;
        match thresholds {
            Some(t) => {
                self.put("      <Threshold>")?;
                self.put(&format!(
                    r#"        <cvParam cvRef="PSI-MS" accession="MS:1001576" name="Scaffold: Minimum Protein Probability" value="{}"/>"#,
                    t.protein_probability
                ))?;
                self.put(&format!(
                    r#"        <cvParam cvRef="PSI-MS" accession="MS:1001575" name="Scaffold: Minimum Peptide Count" value="{}"/>"#,
                    t.min_peptides
                ))?;
                self.put(&format!("        {peptide}"))?;
                self.put("      </Threshold>")?;
            }
            None => self.put(&format!("      <Threshold>{NO_THRESHOLD}</Threshold>"))?,
        }
        self.put("    </ProteinDetectionProtocol>")?;
        self.put("  </AnalysisProtocolCollection>")
    }

    fn inputs(&mut self) -> io::Result<()> {
        self.put("    <Inputs>")?;
        // The schema wants a location; the .sf3 holds the proteins, not the FASTA.
        self.put(r#"      <SearchDatabase id="SDB" location="not-recorded">"#)?;
        self.put(r#"        <FileFormat><cvParam cvRef="PSI-MS" accession="MS:1001348" name="FASTA format"/></FileFormat>"#)?;
        self.put(r#"        <DatabaseName><userParam name="the searched database (the matched proteins are in SequenceCollection)"/></DatabaseName>"#)?;
        self.put("      </SearchDatabase>")?;
        for (run, file) in mgf_names(self.project).iter().enumerate() {
            self.put(&format!(
                r#"      <SpectraData id="SD_{run}" location="{}" name="{}">"#,
                escape(&uri(file)),
                escape(&self.project.runs[run].name)
            ))?;
            self.put(r#"        <FileFormat><cvParam cvRef="PSI-MS" accession="MS:1001062" name="Mascot MGF format"/></FileFormat>"#)?;
            self.put(r#"        <SpectrumIDFormat><cvParam cvRef="PSI-MS" accession="MS:1000774" name="multiple peak list nativeID format"/></SpectrumIDFormat>"#)?;
            self.put("      </SpectraData>")?;
        }
        self.put("    </Inputs>")
    }

    fn spectra(&mut self, engine: usize) -> io::Result<()> {
        let project = self.project;
        self.put(&format!(
            r#"      <SpectrumIdentificationList id="SIL_{engine}">"#
        ))?;
        for (position, spectrum) in project.spectra.iter().enumerate() {
            if !spectrum.matches.iter().any(|found| found.engine == engine) {
                continue;
            }
            self.put(&format!(
                r#"        <SpectrumIdentificationResult id="SIR_{engine}_{}" spectrumID="index={}" spectraData_ref="SD_{}">"#,
                spectrum.id, self.index.positions[position], spectrum.run
            ))?;
            for (number, found) in spectrum.matches.iter().enumerate() {
                if found.engine == engine {
                    self.item(spectrum, number, found)?;
                }
            }
            if let Some(title) = &spectrum.title {
                self.put(&format!(
                    r#"          <cvParam cvRef="PSI-MS" accession="MS:1000796" name="spectrum title" value="{}"/>"#,
                    escape(title)
                ))?;
            }
            if let Some(rt) = spectrum.rt_seconds {
                self.put(&format!(
                    r#"          <cvParam cvRef="PSI-MS" accession="MS:1000016" name="scan start time" value="{rt}" unitCvRef="UO" unitAccession="UO:0000010" unitName="second"/>"#
                ))?;
            }
            self.put("        </SpectrumIdentificationResult>")?;
        }
        self.put("      </SpectrumIdentificationList>")
    }

    fn item(&mut self, spectrum: &Spectrum, number: usize, found: &'a Match) -> io::Result<()> {
        self.report.candidates += 1;
        let charge = spectrum.charge.unwrap_or(0);
        let mut attributes = format!(
            r#"id="{}" chargeState="{charge}" experimentalMassToCharge="{}""#,
            item_id(found.engine, spectrum, number),
            spectrum.precursor_mz().unwrap_or(0.0)
        );
        if let Some(mass) = found.mass.filter(|_| charge > 0) {
            let z = f64::from(charge);
            let _ = write!(
                attributes,
                r#" calculatedMassToCharge="{}""#,
                (mass + z * PROTON) / z
            );
        }
        let peptide = self.index.peptide_ids[&found.peptide.key()];
        let passes = passes_peptide(found, self.project.thresholds);
        self.put(&format!(
            r#"          <SpectrumIdentificationItem {attributes} peptide_ref="Pep_{peptide}" rank="{}" passThreshold="{passes}">"#,
            found.rank
        ))?;
        for evidence in self.index.evidence_of(found) {
            self.put(&format!(
                r#"            <PeptideEvidenceRef peptideEvidence_ref="PE_{evidence}"/>"#
            ))?;
        }
        for score in &found.scores {
            self.put(&format!(
                "            {}",
                param(score.name, score.accession, score.value)
            ))?;
        }
        if let Some(probability) = found.probability {
            self.put(&format!(
                "            {}",
                param(
                    "Scaffold:Peptide Probability",
                    Some("MS:1001568"),
                    f64::from(probability)
                )
            ))?;
        }
        if found.chosen {
            self.put(r#"            <userParam name="Scaffold: candidate chosen for this spectrum" value="true"/>"#)?;
        }
        self.put("          </SpectrumIdentificationItem>")
    }

    fn groups(&mut self) -> io::Result<()> {
        let project = self.project;
        // accession → evidence → items Scaffold picked
        let mut picked: HashMap<&str, BTreeMap<usize, Vec<String>>> = HashMap::new();
        // (run, accession) → sequences Scaffold picked at its threshold
        let mut confident: HashMap<(usize, &str), HashSet<&str>> = HashMap::new();
        for spectrum in &project.spectra {
            for (number, found) in spectrum.matches.iter().enumerate() {
                if !found.chosen {
                    continue;
                }
                let evidence = self.index.evidence_of(found);
                for (protein, id) in found.proteins.iter().zip(evidence) {
                    picked
                        .entry(&protein.accession)
                        .or_default()
                        .entry(id)
                        .or_default()
                        .push(item_id(found.engine, spectrum, number));
                    if passes_peptide(found, project.thresholds) {
                        confident
                            .entry((spectrum.run, &protein.accession))
                            .or_default()
                            .insert(&found.peptide.sequence);
                    }
                }
            }
        }
        self.put(r#"      <ProteinDetectionList id="PDL">"#)?;
        for (number, group) in project.groups.iter().enumerate() {
            let members: Vec<&String> = group
                .accessions
                .iter()
                .filter(|accession| picked.contains_key(accession.as_str()))
                .collect();
            if members.is_empty() {
                self.report.groups_without_peptides += 1;
                continue;
            }
            self.report.groups += 1;
            let passes = match project.thresholds {
                None => true,
                Some(t) => group.probabilities.iter().any(|&(run, probability)| {
                    let peptides: HashSet<&str> = group
                        .accessions
                        .iter()
                        .filter_map(|accession| confident.get(&(run, accession.as_str())))
                        .flatten()
                        .copied()
                        .collect();
                    probability + PROBABILITY_SLACK >= t.protein_probability
                        && peptides.len() >= t.min_peptides as usize
                }),
            };
            self.put(&format!(
                r#"        <ProteinAmbiguityGroup id="PAG_{number}">"#
            ))?;
            for (position, accession) in members.into_iter().enumerate() {
                let protein = self.index.proteins[accession.as_str()];
                self.put(&format!(
                    r#"          <ProteinDetectionHypothesis id="PDH_{number}_{protein}" dBSequence_ref="DBSeq_{protein}" passThreshold="{}">"#,
                    passes && position == 0
                ))?;
                for (evidence, items) in &picked[accession.as_str()] {
                    self.put(&format!(
                        r#"            <PeptideHypothesis peptideEvidence_ref="PE_{evidence}">"#
                    ))?;
                    for item in items {
                        self.put(&format!(
                            r#"              <SpectrumIdentificationItemRef spectrumIdentificationItem_ref="{item}"/>"#
                        ))?;
                    }
                    self.put("            </PeptideHypothesis>")?;
                }
                if position == 0 {
                    self.put(r#"            <cvParam cvRef="PSI-MS" accession="MS:1002403" name="group representative"/>"#)?;
                    self.put(r#"            <cvParam cvRef="PSI-MS" accession="MS:1002401" name="leading protein"/>"#)?;
                }
                if let Some(best) = group.best_probability() {
                    self.put(&format!(
                        "            {}",
                        param("Scaffold:Protein Probability", Some("MS:1001579"), best)
                    ))?;
                }
                for &(run, probability) in &group.probabilities {
                    self.put(&format!(
                        r#"            <userParam name="Scaffold:Protein Probability in {}" value="{probability}"/>"#,
                        escape(&project.runs[run].name)
                    ))?;
                }
                self.put("          </ProteinDetectionHypothesis>")?;
            }
            self.put(&format!(
                r#"          <cvParam cvRef="PSI-MS" accession="MS:1002415" name="protein group passes threshold" value="{passes}"/>"#
            ))?;
            self.put("        </ProteinAmbiguityGroup>")?;
        }
        self.put("      </ProteinDetectionList>")
    }
}

fn item_id(engine: usize, spectrum: &Spectrum, number: usize) -> String {
    format!("SII_{engine}_{}_{number}", spectrum.id)
}

/// Scaffold's verdict on a candidate.
fn passes_peptide(found: &Match, thresholds: Option<Thresholds>) -> bool {
    let Some(t) = thresholds else {
        return true;
    };
    found.chosen
        && found.probability.is_some_and(|probability| {
            f64::from(probability) + PROBABILITY_SLACK >= t.peptide_probability
        })
}

fn param(name: &str, accession: Option<&str>, value: f64) -> String {
    match accession {
        Some(accession) => format!(
            r#"<cvParam cvRef="PSI-MS" accession="{accession}" name="{}" value="{value}"/>"#,
            escape(name)
        ),
        None => format!(r#"<userParam name="{}" value="{value}"/>"#, escape(name)),
    }
}

/// The PSI-MS term for a search engine.
fn software_term(engine: &str) -> Option<(&'static str, &'static str)> {
    let name = engine.to_ascii_lowercase();
    Some(if name.starts_with("mascot") {
        ("MS:1001207", "Mascot")
    } else if name.starts_with("sequest") {
        ("MS:1001208", "SEQUEST")
    } else if name.starts_with("x! tandem") || name.starts_with("x!tandem") {
        ("MS:1001476", "X!Tandem")
    } else if name.starts_with("maxquant") {
        ("MS:1001583", "MaxQuant")
    } else if name.starts_with("spectrum mill") {
        ("MS:1000687", "Spectrum Mill for MassHunter Workstation")
    } else if name.starts_with("msfragger") {
        ("MS:1003014", "MSFragger")
    } else {
        return None;
    })
}

/// A `Modification` element (location 0 is the N-terminus, length + 1 the
/// C-terminus), and whether Unimod knows it.
fn modification_xml(modification: &Modification, length: usize) -> (String, bool) {
    let location = match modification.site {
        Site::NTerm => 0,
        Site::Residue(position) => position,
        Site::CTerm => length + 1,
    };
    let known = modification
        .name
        .as_deref()
        .and_then(|name| unimod::lookup(name, modification.mass_delta));
    let term = match (&known, &modification.name) {
        (Some((accession, name)), _) => {
            format!(r#"<cvParam cvRef="UNIMOD" accession="{accession}" name="{name}"/>"#)
        }
        (None, name) => format!(
            r#"<cvParam cvRef="PSI-MS" accession="MS:1001460" name="unknown modification"{}/>"#,
            name.as_ref()
                .map_or_else(String::new, |n| format!(r#" value="{}""#, escape(n)))
        ),
    };
    (
        format!(
            r#"      <Modification location="{location}" monoisotopicMassDelta="{}">{term}</Modification>"#,
            modification.mass_delta
        ),
        known.is_some(),
    )
}

/// A file name as an `xs:anyURI`.
fn uri(path: &str) -> String {
    let mut encoded = String::with_capacity(path.len());
    for byte in path.bytes() {
        if byte.is_ascii_alphanumeric() || b"-._~/:@!$&'()*+,;=".contains(&byte) {
            encoded.push(byte as char);
        } else {
            let _ = write!(encoded, "%{byte:02X}");
        }
    }
    encoded
}

fn escape(text: &str) -> String {
    let mut escaped = String::with_capacity(text.len());
    for c in text.chars() {
        match c {
            '&' => escaped.push_str("&amp;"),
            '<' => escaped.push_str("&lt;"),
            '>' => escaped.push_str("&gt;"),
            '"' => escaped.push_str("&quot;"),
            '\'' => escaped.push_str("&apos;"),
            c if (c as u32) < 0x20 && !matches!(c, '\t' | '\n' | '\r') => {}
            c => escaped.push(c),
        }
    }
    escaped
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn modifications_use_mzidentml_locations_and_unimod() {
        let oxidation = Modification {
            site: Site::Residue(3),
            mass_delta: 15.994_915,
            name: Some("Oxidation".into()),
        };
        assert_eq!(
            modification_xml(&oxidation, 8),
            (
                r#"      <Modification location="3" monoisotopicMassDelta="15.994915"><cvParam cvRef="UNIMOD" accession="UNIMOD:35" name="Oxidation"/></Modification>"#.to_string(),
                true
            )
        );
        let odd = Modification {
            site: Site::CTerm,
            mass_delta: 12.5,
            name: Some("Mystery".into()),
        };
        let (xml, known) = modification_xml(&odd, 8);
        assert!(!known);
        assert!(xml.contains(r#"location="9""#) && xml.contains(r#"value="Mystery""#));
    }

    #[test]
    fn text_is_escaped_and_locations_are_uris() {
        assert_eq!(escape(r#"a<b & "c""#), "a&lt;b &amp; &quot;c&quot;");
        assert_eq!(uri("run 1.mgf"), "run%201.mgf");
        assert_eq!(
            software_term("Mascot (Ion Score Only)").unwrap().0,
            "MS:1001207"
        );
        assert_eq!(software_term("Byonic"), None);
    }
}
