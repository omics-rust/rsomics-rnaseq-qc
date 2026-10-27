use std::collections::{BTreeMap, HashMap};
use std::ops::Range;

use rsomics_bamio::raw::RawRecord;
use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

use crate::annotation::Annotation;

const CDS: usize = 0;
const UTR5: usize = 1;
const UTR3: usize = 2;
const INTRON: usize = 3;
const UP1: usize = 4;
const UP5: usize = 5;
const UP10: usize = 6;
const DOWN1: usize = 7;
const DOWN5: usize = 8;
const DOWN10: usize = 9;
const LABELS: [&str; 10] = [
    "CDS_Exons",
    "5'UTR_Exons",
    "3'UTR_Exons",
    "Introns",
    "TSS_up_1kb",
    "TSS_up_5kb",
    "TSS_up_10kb",
    "TES_down_1kb",
    "TES_down_5kb",
    "TES_down_10kb",
];

type Regions = [Vec<Range<u64>>; 10];

#[derive(Clone, Debug, Serialize)]
pub(crate) struct DistributionTotals {
    #[serde(rename = "Total Reads")]
    pub total_reads: u64,
    #[serde(rename = "Total Tags")]
    pub total_tags: u64,
    #[serde(rename = "Total Assigned Tags")]
    pub total_assigned_tags: u64,
}

#[derive(Clone, Debug, Serialize)]
pub(crate) struct DistributionGroup {
    pub bases: u64,
    pub tags: u64,
    pub tags_per_kb: f64,
}

#[derive(Clone, Debug, Serialize)]
pub(crate) struct DistributionResult {
    pub totals: DistributionTotals,
    pub groups: BTreeMap<String, DistributionGroup>,
}

pub(crate) struct Distribution {
    regions: Vec<Regions>,
    bases: [u64; 10],
    tags: [u64; 10],
    total_reads: u64,
    total_tags: u64,
    unassigned_tags: u64,
}

impl Distribution {
    pub fn new(annotation: &Annotation, reference_names: &[String]) -> Result<Self> {
        let mut reference_ids = HashMap::with_capacity(reference_names.len());
        for (reference_id, name) in reference_names.iter().enumerate() {
            if reference_ids
                .insert(name.to_uppercase(), reference_id)
                .is_some()
            {
                return Err(RsomicsError::InvalidInput(format!(
                    "distribution requires unambiguous case-normalized BAM reference names: {name}"
                )));
            }
        }
        let mut references: Vec<Regions> = reference_names
            .iter()
            .map(|_| std::array::from_fn(|_| Vec::new()))
            .collect();
        for transcript in &annotation.transcripts {
            let reference_id = reference_ids
                .get(&transcript.chrom.to_uppercase())
                .ok_or_else(|| {
                    RsomicsError::InvalidInput(format!(
                        "annotation transcript {} uses contig {} absent from the BAM header",
                        transcript.name, transcript.chrom
                    ))
                })?;
            let regions = &mut references[*reference_id];
            let cds_start = u64::from(transcript.cds_start);
            let cds_end = u64::from(transcript.cds_end);
            let (left_utr, right_utr) = if transcript.reverse {
                (UTR3, UTR5)
            } else {
                (UTR5, UTR3)
            };
            for exon in &transcript.exons {
                let start = u64::from(exon.start);
                let end = u64::from(exon.end);
                regions[CDS].push(start.max(cds_start)..end.min(cds_end));
                regions[left_utr].push(start..end.min(cds_start));
                regions[right_utr].push(start.max(cds_end)..end);
            }
            for exons in transcript.exons.windows(2) {
                regions[INTRON].push(u64::from(exons[0].end)..u64::from(exons[1].start));
            }
            let start = u64::from(transcript.start);
            let end = u64::from(transcript.end);
            for (size, upstream, downstream) in [
                (1000, UP1, DOWN1),
                (5000, UP5, DOWN5),
                (10000, UP10, DOWN10),
            ] {
                let left = start.saturating_sub(size)..start;
                let right = end..end + size;
                if transcript.reverse {
                    regions[upstream].push(right);
                    regions[downstream].push(left);
                } else {
                    regions[upstream].push(left);
                    regions[downstream].push(right);
                }
            }
        }
        let mut bases = [0; 10];
        for regions in &mut references {
            for category in regions.iter_mut() {
                *category = union(std::mem::take(category));
            }
            regions[UTR5] = subtract(&regions[UTR5], &regions[CDS]);
            regions[UTR3] = subtract(&regions[UTR3], &regions[CDS]);
            let introns = subtract(&regions[INTRON], &regions[CDS]);
            let introns = subtract(&introns, &regions[UTR5]);
            regions[INTRON] = subtract(&introns, &regions[UTR3]);
            let genic = union(
                regions[..=INTRON]
                    .iter()
                    .flat_map(|category| category.iter().cloned())
                    .collect(),
            );
            for category in &mut regions[UP1..] {
                *category = subtract(category, &genic);
            }
            for (index, category) in regions.iter().enumerate() {
                bases[index] += category
                    .iter()
                    .map(|range| range.end - range.start)
                    .sum::<u64>();
            }
        }
        Ok(Self {
            regions: references,
            bases,
            tags: [0; 10],
            total_reads: 0,
            total_tags: 0,
            unassigned_tags: 0,
        })
    }

    pub fn observe(&mut self, record: &RawRecord, cigar: &[(u8, u32)], legacy: bool) -> Result<()> {
        if record.flags() & 0x704 != 0 {
            return Ok(());
        }
        let invalid_record = |reason: &str| {
            RsomicsError::InvalidInput(format!(
                "read {}: {reason}",
                String::from_utf8_lossy(record.name())
            ))
        };
        let reference_id = usize::try_from(record.reference_sequence_id())
            .map_err(|_| invalid_record("mapped record has no reference ID"))?;
        if reference_id >= self.regions.len() {
            return Err(invalid_record("reference ID is absent from the BAM header"));
        }
        let mut position = u64::try_from(record.alignment_start())
            .map_err(|_| invalid_record("mapped record has no alignment start"))?;
        self.total_reads += 1;
        if legacy {
            for &(operation, length) in cigar {
                match operation {
                    0 => {
                        let end = position + u64::from(length);
                        self.assign_block(reference_id, position, end, true);
                        position = end;
                    }
                    2..=4 => position += u64::from(length),
                    1 | 5..=8 => {}
                    _ => return Err(invalid_record("unsupported CIGAR operation")),
                }
            }
            return Ok(());
        }
        let mut block_start = None;
        for &(operation, length) in cigar {
            match operation {
                0 | 7..=8 => {
                    block_start.get_or_insert(position);
                    position += u64::from(length);
                }
                2 | 3 => {
                    if let Some(start) = block_start.take() {
                        self.assign_block(reference_id, start, position, false);
                    }
                    position += u64::from(length);
                }
                1 | 6 => {}
                4 | 5 => {
                    if let Some(start) = block_start.take() {
                        self.assign_block(reference_id, start, position, false);
                    }
                }
                _ => return Err(invalid_record("unsupported CIGAR operation")),
            }
        }
        if let Some(start) = block_start {
            self.assign_block(reference_id, start, position, false);
        }
        Ok(())
    }

    fn assign_block(&mut self, reference_id: usize, start: u64, end: u64, legacy: bool) {
        self.total_tags += 1;
        let midpoint = start + (end - start) / 2;
        let regions = &self.regions[reference_id];
        let overlaps = |category: usize| contains(&regions[category], midpoint, legacy);
        if overlaps(CDS) {
            self.tags[CDS] += 1;
            return;
        }
        match (overlaps(UTR5), overlaps(UTR3)) {
            (true, false) => {
                self.tags[UTR5] += 1;
                return;
            }
            (false, true) => {
                self.tags[UTR3] += 1;
                return;
            }
            (true, true) => {
                self.unassigned_tags += 1;
                return;
            }
            (false, false) => {}
        }
        if overlaps(INTRON) {
            self.tags[INTRON] += 1;
            return;
        }
        if overlaps(UP10) && overlaps(DOWN10) {
            self.unassigned_tags += 1;
            return;
        }
        for category in UP1..=DOWN10 {
            if overlaps(category) {
                let last = if category <= UP10 { UP10 } else { DOWN10 };
                for tags in &mut self.tags[category..=last] {
                    *tags += 1;
                }
                return;
            }
        }
        self.unassigned_tags += 1;
    }

    pub fn finish(self) -> DistributionResult {
        DistributionResult {
            totals: DistributionTotals {
                total_reads: self.total_reads,
                total_tags: self.total_tags,
                total_assigned_tags: self.total_tags - self.unassigned_tags,
            },
            groups: LABELS
                .iter()
                .enumerate()
                .map(|(index, &label)| {
                    let bases = self.bases[index];
                    let tags = self.tags[index];
                    (
                        label.into(),
                        DistributionGroup {
                            bases,
                            tags,
                            tags_per_kb: tags as f64 * 1000.0 / (bases as f64 + 1.0),
                        },
                    )
                })
                .collect(),
        }
    }
}

fn union(mut intervals: Vec<Range<u64>>) -> Vec<Range<u64>> {
    intervals.retain(|interval| interval.start < interval.end);
    intervals.sort_unstable_by_key(|interval| (interval.start, interval.end));
    let mut merged: Vec<Range<u64>> = Vec::with_capacity(intervals.len());
    for interval in intervals {
        if let Some(last) = merged.last_mut()
            && interval.start <= last.end
        {
            last.end = last.end.max(interval.end);
        } else {
            merged.push(interval);
        }
    }
    merged
}

fn subtract(intervals: &[Range<u64>], excluded: &[Range<u64>]) -> Vec<Range<u64>> {
    let mut remaining = Vec::new();
    let mut first_excluded = 0;
    for interval in intervals {
        while first_excluded < excluded.len() && excluded[first_excluded].end <= interval.start {
            first_excluded += 1;
        }
        let mut start = interval.start;
        for mask in &excluded[first_excluded..] {
            if mask.start >= interval.end {
                break;
            }
            if start < mask.start {
                remaining.push(start..mask.start);
            }
            start = start.max(mask.end);
            if start >= interval.end {
                break;
            }
        }
        if start < interval.end {
            remaining.push(start..interval.end);
        }
    }
    remaining
}

fn contains(intervals: &[Range<u64>], point: u64, legacy: bool) -> bool {
    // RSeQC passes a zero-width query to bx, excluding both interval endpoints.
    let end = intervals.partition_point(|interval| {
        if legacy {
            interval.start < point
        } else {
            interval.start <= point
        }
    });
    end > 0 && point < intervals[end - 1].end
}

#[cfg(test)]
mod tests {
    use rsomics_bamio::raw::RawRecord;

    use super::Distribution;
    use crate::annotation::{Annotation, Transcript};

    fn record(start: i32) -> RawRecord {
        let mut record = RawRecord::default();
        record.clear_flag_bits(0xffff);
        record.set_reference_sequence_id(0);
        record.set_alignment_start(start);
        record
    }

    #[test]
    fn report_distribution_distinguishes_canonical_and_legacy_cigar() {
        let annotation = Annotation {
            transcripts: vec![Transcript {
                name: "tx".into(),
                chrom: "chr1".into(),
                start: 100,
                end: 500,
                cds_start: 200,
                cds_end: 400,
                reverse: false,
                exons: vec![100..250, 350..500],
            }],
        };
        let cigar = [
            (4, 5),
            (0, 5),
            (1, 1),
            (7, 5),
            (8, 5),
            (0, 5),
            (3, 100),
            (0, 10),
        ];
        let mut canonical = Distribution::new(&annotation, &["CHR1".into()]).unwrap();
        canonical.observe(&record(190), &cigar, false).unwrap();
        let canonical = canonical.finish();
        assert_eq!(canonical.totals.total_reads, 1);
        assert_eq!(canonical.totals.total_tags, 2);
        assert_eq!(canonical.groups["CDS_Exons"].tags, 1);
        assert_eq!(canonical.groups["Introns"].tags, 1);
        assert_eq!(canonical.groups["5'UTR_Exons"].tags, 0);
        assert_eq!(canonical.groups["CDS_Exons"].bases, 100);
        assert_eq!(canonical.groups["CDS_Exons"].tags_per_kb, 1000.0 / 101.0);

        let mut legacy = Distribution::new(&annotation, &["CHR1".into()]).unwrap();
        legacy.observe(&record(190), &cigar, true).unwrap();
        let legacy = legacy.finish();
        assert_eq!(legacy.totals.total_reads, 1);
        assert_eq!(legacy.totals.total_tags, 3);
        assert_eq!(legacy.groups["CDS_Exons"].tags, 1);
        assert_eq!(legacy.groups["Introns"].tags, 1);
        assert_eq!(legacy.groups["5'UTR_Exons"].tags, 1);
    }

    #[test]
    fn distribution_preserves_ambiguity_cumulative_flanks_and_legacy_boundary() {
        let transcript = Transcript {
            name: "plus".into(),
            chrom: "chr1".into(),
            start: 10000,
            end: 10200,
            cds_start: 10050,
            cds_end: 10150,
            reverse: false,
            exons: std::iter::once(10000..10200).collect(),
        };
        let mut opposite = transcript.clone();
        opposite.name = "minus".into();
        opposite.reverse = true;
        let mut isolated = transcript.clone();
        isolated.name = "isolated".into();
        isolated.chrom = "chr2".into();
        isolated.start = 30000;
        isolated.end = 30200;
        isolated.cds_start = 30000;
        isolated.cds_end = 30200;
        isolated.exons = std::iter::once(30000..30200).collect();
        let annotation = Annotation {
            transcripts: vec![transcript, opposite, isolated],
        };
        let names = ["chr1".into(), "chr2".into()];
        let mut canonical = Distribution::new(&annotation, &names).unwrap();
        canonical.observe(&record(10019), &[(0, 2)], false).unwrap();
        canonical.observe(&record(10049), &[(0, 2)], false).unwrap();
        canonical.observe(&record(10499), &[(0, 2)], false).unwrap();
        let mut upstream = record(29499);
        upstream.set_reference_sequence_id(1);
        upstream.set_flag_bits(0x800);
        canonical.observe(&upstream, &[(0, 2)], false).unwrap();
        let mut excluded = upstream.clone();
        excluded.set_flag_bits(0x400);
        canonical.observe(&excluded, &[(0, 2)], false).unwrap();
        let canonical = canonical.finish();
        assert_eq!(canonical.totals.total_reads, 4);
        assert_eq!(canonical.totals.total_tags, 4);
        assert_eq!(canonical.totals.total_assigned_tags, 2);
        assert_eq!(canonical.groups["CDS_Exons"].tags, 1);
        assert_eq!(canonical.groups["TSS_up_1kb"].tags, 1);
        assert_eq!(canonical.groups["TSS_up_5kb"].tags, 1);
        assert_eq!(canonical.groups["TSS_up_10kb"].tags, 1);
        assert_eq!(canonical.groups["TSS_up_1kb"].bases, 3000);

        let mut legacy = Distribution::new(&annotation, &names).unwrap();
        legacy.observe(&record(10049), &[(0, 2)], true).unwrap();
        assert_eq!(legacy.finish().totals.total_assigned_tags, 0);
        assert!(Distribution::new(&annotation, &["chr1".into(), "CHR1".into()]).is_err());
    }
}
