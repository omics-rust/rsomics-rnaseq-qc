use std::collections::{HashMap, HashSet};

use rsomics_bamio::raw::RawRecord;
use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

use crate::annotation::Annotation;
use crate::cli::Profile;

#[derive(Clone, Copy, Serialize)]
pub(crate) enum JunctionClass {
    #[serde(rename = "annotated")]
    Annotated,
    #[serde(rename = "partial_novel")]
    PartialNovel,
    #[serde(rename = "complete_novel")]
    CompleteNovel,
}

impl JunctionClass {
    pub fn name(self) -> &'static str {
        match self {
            Self::Annotated => "annotated",
            Self::PartialNovel => "partial_novel",
            Self::CompleteNovel => "complete_novel",
        }
    }
}

#[derive(Clone, Copy, Serialize)]
pub(crate) struct JunctionPolicy {
    pub excluded_flags: u16,
    pub minimum_mapq: u8,
    pub minimum_intron: u32,
    projection: &'static str,
    chromosome_identity: &'static str,
    chromosome_display: &'static str,
    annotation: &'static str,
    support: &'static str,
    order: &'static str,
    empty_result: &'static str,
    omitted_outputs: [&'static str; 5],
    #[serde(skip)]
    profile: Profile,
}

impl JunctionPolicy {
    pub fn new(profile: Profile, minimum_mapq: u8, minimum_intron: u32) -> Self {
        Self {
            excluded_flags: 0x704,
            minimum_mapq,
            minimum_intron,
            projection: match profile {
                Profile::Canonical => "M/D/N/=/X advance reference; I/S/H/P do not",
                Profile::Rseqc505 => {
                    "RSeQC 5.0.5 junctions: M/D/N advance reference; I/S/H/P/=/X do not"
                }
            },
            chromosome_identity: match profile {
                Profile::Canonical => "exact header identities",
                Profile::Rseqc505 => {
                    "uppercase pooled identities after strict exact annotation-header validation"
                }
            },
            chromosome_display: match profile {
                Profile::Canonical => "exact header names",
                Profile::Rseqc505 => "uppercase pooled name with every CHR replaced by chr",
            },
            annotation: "independent donor and acceptor site membership; no intron-pair or strand requirement",
            support: "one event per eligible original-record N occurrence, including mates and supplementary records; total events precede minimum-intron filtering",
            order: "first eligible accepted N encounter in BAM order",
            empty_result: "complete zero-result tables and exit 0 for valid no-N/all-short inputs; intentional RSeQC partial-output incompatibility",
            omitted_outputs: ["R", "PDF", "BED12", "Interact", "junction-saturation"],
            profile,
        }
    }
}

#[derive(Default)]
struct Endpoints {
    donors: HashSet<u64>,
    acceptors: HashSet<u64>,
}

pub(crate) struct SpliceSites {
    reference_pools: Vec<usize>,
    names: Vec<String>,
    endpoints: Vec<Endpoints>,
}

impl SpliceSites {
    pub fn new(annotation: &Annotation, references: &[String], profile: Profile) -> Result<Self> {
        let reference_ids: HashMap<_, _> = references
            .iter()
            .enumerate()
            .map(|(id, name)| (name.as_str(), id))
            .collect();
        let mut pools = HashMap::new();
        let mut names = Vec::new();
        let mut endpoints = Vec::new();
        let reference_pools = references
            .iter()
            .map(|name| {
                let key = match profile {
                    Profile::Canonical => name.clone(),
                    Profile::Rseqc505 => name.to_uppercase(),
                };
                let next = names.len();
                *pools.entry(key.clone()).or_insert_with(|| {
                    names.push(match profile {
                        Profile::Canonical => key,
                        Profile::Rseqc505 => key.replace("CHR", "chr"),
                    });
                    endpoints.push(Endpoints::default());
                    next
                })
            })
            .collect::<Vec<_>>();
        for transcript in &annotation.transcripts {
            let id = reference_ids
                .get(transcript.chrom.as_str())
                .ok_or_else(|| {
                    RsomicsError::InvalidInput(format!(
                        "annotation contig {} is absent from the BAM header",
                        transcript.chrom
                    ))
                })?;
            let sites = &mut endpoints[reference_pools[*id]];
            for pair in transcript.exons.windows(2) {
                sites.donors.insert(u64::from(pair[0].end));
                sites.acceptors.insert(u64::from(pair[1].start));
            }
        }
        Ok(Self {
            reference_pools,
            names,
            endpoints,
        })
    }

    fn classify(&self, key: JunctionKey) -> JunctionClass {
        let sites = &self.endpoints[key.pool];
        match (
            sites.donors.contains(&key.start),
            sites.acceptors.contains(&key.end),
        ) {
            (true, true) => JunctionClass::Annotated,
            (false, false) => JunctionClass::CompleteNovel,
            _ => JunctionClass::PartialNovel,
        }
    }
}

#[derive(Serialize)]
pub(crate) struct JunctionRow {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub read_count: u64,
    pub annotation: JunctionClass,
}

#[derive(Default, Serialize)]
pub(crate) struct JunctionCounts {
    total_events: u64,
    known_events: u64,
    partial_novel_events: u64,
    novel_events: u64,
    filtered_events: u64,
    total_junctions: u64,
    known_junctions: u64,
    partial_novel_junctions: u64,
    novel_junctions: u64,
}

#[derive(Serialize)]
pub(crate) struct JunctionResult {
    pub counts: JunctionCounts,
    pub rows: Vec<JunctionRow>,
}

#[derive(Clone, Copy, Eq, Hash, PartialEq)]
struct JunctionKey {
    pool: usize,
    start: u64,
    end: u64,
}

struct SupportedJunction {
    key: JunctionKey,
    read_count: u64,
    class: JunctionClass,
}

pub(crate) struct JunctionAccumulator {
    sites: SpliceSites,
    policy: JunctionPolicy,
    counts: JunctionCounts,
    row_indices: HashMap<JunctionKey, usize>,
    rows: Vec<SupportedJunction>,
}

impl JunctionAccumulator {
    pub fn new(sites: SpliceSites, policy: JunctionPolicy) -> Self {
        Self {
            sites,
            policy,
            counts: JunctionCounts::default(),
            row_indices: HashMap::new(),
            rows: Vec::new(),
        }
    }

    pub fn observe(&mut self, record: &RawRecord, cigar: &[(u8, u32)]) {
        if record.flags() & self.policy.excluded_flags != 0
            || record.mapping_quality() < self.policy.minimum_mapq
        {
            return;
        }
        let pool = self.sites.reference_pools[record.reference_sequence_id() as usize];
        let mut position = record.alignment_start() as u64;
        for &(kind, length) in cigar {
            let end = position + u64::from(length);
            match kind {
                3 => {
                    self.counts.total_events += 1;
                    if length < self.policy.minimum_intron {
                        self.counts.filtered_events += 1;
                    } else {
                        let key = JunctionKey {
                            pool,
                            start: position,
                            end,
                        };
                        let next = self.rows.len();
                        let index = *self.row_indices.entry(key).or_insert_with(|| {
                            let class = self.sites.classify(key);
                            self.counts.total_junctions += 1;
                            match class {
                                JunctionClass::Annotated => self.counts.known_junctions += 1,
                                JunctionClass::PartialNovel => {
                                    self.counts.partial_novel_junctions += 1
                                }
                                JunctionClass::CompleteNovel => self.counts.novel_junctions += 1,
                            }
                            self.rows.push(SupportedJunction {
                                key,
                                read_count: 0,
                                class,
                            });
                            next
                        });
                        let row = &mut self.rows[index];
                        row.read_count += 1;
                        match row.class {
                            JunctionClass::Annotated => self.counts.known_events += 1,
                            JunctionClass::PartialNovel => self.counts.partial_novel_events += 1,
                            JunctionClass::CompleteNovel => self.counts.novel_events += 1,
                        }
                    }
                    position = end;
                }
                0 | 2 => position = end,
                7 | 8 if self.policy.profile == Profile::Canonical => position = end,
                _ => {}
            }
        }
    }

    pub fn finish(self) -> JunctionResult {
        JunctionResult {
            counts: self.counts,
            rows: self
                .rows
                .into_iter()
                .map(|row| JunctionRow {
                    chrom: self.sites.names[row.key.pool].clone(),
                    start: row.key.start,
                    end: row.key.end,
                    read_count: row.read_count,
                    annotation: row.class,
                })
                .collect(),
        }
    }
}
