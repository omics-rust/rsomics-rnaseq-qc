use std::collections::{BTreeMap, HashSet};
use std::ops::Bound::{Excluded, Unbounded};

use rsomics_bamio::raw::RawRecord;
use rsomics_pileup::Column;
use serde::Serialize;

use crate::annotation::{Annotation, Transcript, TranscriptIndex};

struct TranscriptSamples {
    coverage: Vec<u32>,
    integrity: Vec<u32>,
    tin_denominator: usize,
    starts: HashSet<u32>,
}

pub(crate) struct SampledMetrics {
    transcripts: Vec<TranscriptSamples>,
    depths: Vec<BTreeMap<u32, u64>>,
    next_positions: Vec<Option<u32>>,
    references: Vec<usize>,
    minimum_starts: usize,
    coverage_selected: bool,
    integrity_selected: bool,
}

#[derive(Debug, Serialize)]
pub(crate) struct CoverageResult {
    pub sample: String,
    pub coverage: Vec<u64>,
    pub eligible_transcripts: usize,
}

#[derive(Debug, Serialize)]
pub(crate) struct IntegrityResult {
    pub transcripts: Vec<TranscriptIntegrity>,
    pub summary: IntegritySummary,
    pub eligible_transcripts: Vec<String>,
}

#[derive(Debug, Serialize)]
pub(crate) struct TranscriptIntegrity {
    #[serde(rename = "geneID")]
    pub name: String,
    pub chrom: String,
    pub tx_start: u32,
    pub tx_end: u32,
    #[serde(rename = "TIN")]
    pub tin: f64,
    pub eligible: bool,
    pub sampled_positions: usize,
}

#[derive(Debug, Default, Serialize)]
pub(crate) struct IntegritySummary {
    pub mean: f64,
    pub median: f64,
    pub stdev: f64,
}

impl SampledMetrics {
    pub fn new(
        annotation: &Annotation,
        reference_names: &[String],
        coverage: bool,
        integrity: bool,
        minimum_length: usize,
        tin_samples: usize,
        minimum_starts: usize,
    ) -> Self {
        let reference_ids: BTreeMap<&str, usize> = reference_names
            .iter()
            .enumerate()
            .map(|(i, name)| (name.as_str(), i))
            .collect();
        let mut depths = vec![BTreeMap::new(); reference_names.len()];
        let mut references = Vec::with_capacity(annotation.transcripts.len());
        let coverage_identities = annotation
            .transcripts
            .iter()
            .enumerate()
            .filter(|(_, tx)| tx.mrna_len() >= minimum_length)
            .map(|(index, tx)| {
                (
                    (
                        tx.chrom.as_str(),
                        tx.start,
                        tx.end,
                        tx.name.as_str(),
                        tx.reverse,
                    ),
                    index,
                )
            })
            .collect::<BTreeMap<_, _>>();
        let transcripts = annotation
            .transcripts
            .iter()
            .enumerate()
            .map(|(index, tx)| {
                let reference = reference_ids[tx.chrom.as_str()];
                references.push(reference);
                let identity = (
                    tx.chrom.as_str(),
                    tx.start,
                    tx.end,
                    tx.name.as_str(),
                    tx.reverse,
                );
                let coverage_points = if coverage
                    && tx.mrna_len() >= minimum_length
                    && coverage_identities.get(&identity) == Some(&index)
                {
                    coverage_positions(tx)
                } else {
                    Vec::new()
                };
                let mut integrity_points = if integrity {
                    integrity_positions(tx, tin_samples)
                } else {
                    Vec::new()
                };
                let tin_denominator = integrity_points.len();
                integrity_points.sort_unstable();
                integrity_points.dedup();
                for &point in coverage_points.iter().chain(&integrity_points) {
                    depths[reference].insert(point, 0);
                }
                TranscriptSamples {
                    coverage: coverage_points,
                    integrity: integrity_points,
                    tin_denominator,
                    starts: HashSet::new(),
                }
            })
            .collect();
        let next_positions = depths
            .iter()
            .map(|points| points.first_key_value().map(|(&point, _)| point))
            .collect();
        Self {
            transcripts,
            depths,
            next_positions,
            references,
            minimum_starts,
            coverage_selected: coverage,
            integrity_selected: integrity,
        }
    }

    pub fn observe_start(&mut self, record: &RawRecord, index: &TranscriptIndex) {
        if !self.integrity_selected || record.flags() & 0x304 != 0 {
            return;
        }
        let reference = record.reference_sequence_id() as usize;
        let start = record.alignment_start() as u32;
        for tx in index.overlapping(reference, start, start + 1) {
            let starts = &mut self.transcripts[tx].starts;
            if starts.len() <= self.minimum_starts {
                starts.insert(start);
            }
        }
    }

    pub fn observe_column(&mut self, column: &Column<'_>) {
        let reference = column.reference_id() as usize;
        let position = column.position() as u32;
        if self.next_positions[reference].is_none_or(|next| position < next) {
            return;
        }
        let next = self.depths[reference]
            .range(position..)
            .next()
            .map(|(&point, _)| point);
        self.next_positions[reference] = next;
        if next != Some(position) {
            return;
        }
        let Some(depth) = self.depths[reference].get_mut(&position) else {
            return;
        };
        *depth = column
            .entries()
            .filter(|entry| {
                let projection = entry.projection();
                !projection.is_deletion
                    && !projection.is_reference_skip
                    && entry
                        .record()
                        .quality_scores()
                        .get(projection.qpos)
                        .copied()
                        .unwrap_or(u8::MAX)
                        >= 13
            })
            .count() as u64;
        self.next_positions[reference] = self.depths[reference]
            .range((Excluded(position), Unbounded))
            .next()
            .map(|(&point, _)| point);
    }

    pub fn finish(
        self,
        annotation: &Annotation,
        sample: &str,
    ) -> (Option<CoverageResult>, Option<IntegrityResult>) {
        let mut coverage = Vec::<u64>::new();
        let mut coverage_eligible = 0;
        let mut rows = Vec::new();
        let mut eligible_names = Vec::new();
        let mut scores = Vec::new();
        for (i, points) in self.transcripts.into_iter().enumerate() {
            let tx = &annotation.transcripts[i];
            let depths = &self.depths[self.references[i]];
            if !points.coverage.is_empty() {
                coverage_eligible += 1;
                coverage.resize(coverage.len().max(points.coverage.len()), 0);
                let oriented = if tx.reverse {
                    points.coverage.iter().rev().copied().collect::<Vec<_>>()
                } else {
                    points.coverage
                };
                for (bin, point) in oriented.into_iter().enumerate() {
                    coverage[bin] += depths[&point];
                }
            }
            if self.integrity_selected {
                let eligible = points.starts.len() > self.minimum_starts;
                let tin = if eligible {
                    tin_score(
                        points.integrity.iter().map(|point| depths[point]),
                        points.tin_denominator,
                    )
                } else {
                    0.0
                };
                if eligible {
                    scores.push(tin);
                    eligible_names.push(tx.name.clone());
                }
                rows.push(TranscriptIntegrity {
                    name: tx.name.clone(),
                    chrom: tx.chrom.clone(),
                    tx_start: tx.start,
                    tx_end: tx.end,
                    tin,
                    eligible,
                    sampled_positions: points.tin_denominator,
                });
            }
        }
        let coverage = self.coverage_selected.then(|| CoverageResult {
            sample: sample.to_owned(),
            coverage,
            eligible_transcripts: coverage_eligible,
        });
        let integrity = self.integrity_selected.then(|| IntegrityResult {
            transcripts: rows,
            summary: summarize(&mut scores),
            eligible_transcripts: eligible_names,
        });
        (coverage, integrity)
    }
}

fn genomic_base(tx: &Transcript, mut index: usize) -> u32 {
    for exon in &tx.exons {
        let length = (exon.end - exon.start) as usize;
        if index < length {
            return exon.start + index as u32 + 1;
        }
        index -= length;
    }
    unreachable!("validated transcript index is inside its exons")
}

fn coverage_positions(tx: &Transcript) -> Vec<u32> {
    let length = tx.mrna_len();
    let mut points = (1..=100)
        .map(|percentile| {
            let k = (length - 1) as f64 * f64::from(percentile) / 100.0;
            let floor = k.floor() as usize;
            let ceil = k.ceil() as usize;
            let position = if floor == ceil {
                f64::from(genomic_base(tx, floor))
            } else {
                f64::from(genomic_base(tx, floor)) * (ceil as f64 - k)
                    + f64::from(genomic_base(tx, ceil)) * (k - floor as f64)
            };
            position.round_ties_even() as u32 - 1
        })
        .collect::<Vec<_>>();
    points.dedup();
    points
}

fn integrity_positions(tx: &Transcript, samples: usize) -> Vec<u32> {
    let length = tx.mrna_len();
    if length <= samples {
        let mut positions = vec![tx.start, tx.end - 1];
        for exon in &tx.exons {
            positions.extend(exon.clone());
        }
        return positions;
    }
    let step = length / samples;
    let mut positions = Vec::with_capacity(length / step + 2 * tx.exons.len());
    let mut offset = 0;
    let mut next = 0;
    for exon in &tx.exons {
        let exon_length = (exon.end - exon.start) as usize;
        positions.extend([exon.start, exon.end - 1]);
        while next < offset + exon_length {
            positions.push(exon.start + (next - offset) as u32);
            next += step;
        }
        offset += exon_length;
    }
    positions.sort_unstable();
    positions.dedup();
    positions
}

fn tin_score(depths: impl Iterator<Item = u64>, denominator: usize) -> f64 {
    let positive = depths.filter(|&depth| depth > 0).collect::<Vec<_>>();
    if positive.is_empty() {
        return 0.0;
    }
    let total = positive.iter().sum::<u64>() as f64;
    let entropy = positive
        .iter()
        .map(|&depth| {
            let p = depth as f64 / total;
            -p * p.ln()
        })
        .sum::<f64>();
    100.0 * entropy.exp() / denominator as f64
}

fn summarize(scores: &mut [f64]) -> IntegritySummary {
    if scores.is_empty() {
        return IntegritySummary::default();
    }
    let n = scores.len();
    let mean = scores.iter().sum::<f64>() / n as f64;
    let stdev = (scores
        .iter()
        .map(|score| (score - mean).powi(2))
        .sum::<f64>()
        / n as f64)
        .sqrt();
    scores.sort_unstable_by(f64::total_cmp);
    let median = if n.is_multiple_of(2) {
        (scores[n / 2 - 1] + scores[n / 2]) / 2.0
    } else {
        scores[n / 2]
    };
    IntegritySummary {
        mean,
        median,
        stdev,
    }
}
