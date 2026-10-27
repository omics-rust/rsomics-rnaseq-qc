use std::collections::BTreeMap;

use rsomics_bamio::raw::RawRecord;
use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

use crate::annotation::{Annotation, TranscriptIndex};

#[derive(Default)]
pub(crate) struct Strandedness {
    sampled: u64,
    paired: u64,
    single: u64,
    spec1: u64,
    spec2: u64,
    ambiguous: u64,
}

#[derive(Debug, Serialize)]
pub(crate) struct StrandednessResult {
    pub layout: &'static str,
    pub fractions: BTreeMap<String, f64>,
    pub sampled: u64,
    pub paired_reads: u64,
    pub single_reads: u64,
    pub ambiguous_reads: u64,
}

impl Strandedness {
    pub fn observe(
        &mut self,
        record: &RawRecord,
        annotation: &Annotation,
        index: &TranscriptIndex,
        sample_size: u64,
        mapq_cut: u8,
    ) -> Result<()> {
        if self.sampled >= sample_size
            || record.flags() & 0x704 != 0
            || record.mapping_quality() < mapq_cut
        {
            return Ok(());
        }
        let reference = usize::try_from(record.reference_sequence_id())
            .map_err(|_| RsomicsError::InvalidInput("mapped read has no reference".into()))?;
        let start = u32::try_from(record.alignment_start())
            .map_err(|_| RsomicsError::InvalidInput("mapped read has no start".into()))?;
        let end = start
            .checked_add(u32::try_from(record.sequence_len()).map_err(|_| {
                RsomicsError::InvalidInput("read query length exceeds coordinate range".into())
            })?)
            .ok_or_else(|| RsomicsError::InvalidInput("read query interval overflows".into()))?;
        let overlaps = index.overlapping(reference, start, end);
        if overlaps.is_empty() {
            return Ok(());
        }
        let plus = overlaps.iter().any(|&i| !annotation.transcripts[i].reverse);
        let minus = overlaps.iter().any(|&i| annotation.transcripts[i].reverse);
        let paired = record.flags() & 1 != 0;
        let read1 = if paired {
            match record.flags() & 0xc0 {
                0x40 => true,
                0x80 => false,
                _ => {
                    return Err(RsomicsError::InvalidInput(format!(
                        "paired read {} must identify exactly one segment",
                        String::from_utf8_lossy(record.name())
                    )));
                }
            }
        } else {
            false
        };
        self.sampled += 1;
        if paired {
            self.paired += 1;
        } else {
            self.single += 1;
        }
        if plus && minus {
            self.ambiguous += 1;
            return Ok(());
        }
        let same = (record.flags() & 0x10 != 0) == minus;
        let spec1 = if paired { same == read1 } else { same };
        if spec1 {
            self.spec1 += 1;
        } else {
            self.spec2 += 1;
        }
        Ok(())
    }

    pub fn finish(self) -> StrandednessResult {
        let layout = match (self.paired > 0, self.single > 0) {
            (true, false) => "PairEnd",
            (false, true) => "SingleEnd",
            _ => "Mixture",
        };
        let mut fractions = BTreeMap::new();
        if matches!(layout, "PairEnd" | "SingleEnd") {
            let (first, second) = if layout == "PairEnd" {
                ("1++,1--,2+-,2-+", "1+-,1-+,2++,2--")
            } else {
                ("++,--", "+-,-+")
            };
            let spec1 = self.spec1 as f64 / self.sampled as f64;
            let spec2 = self.spec2 as f64 / self.sampled as f64;
            fractions.insert(
                "Fraction of reads failed to determine".into(),
                (1.0 - spec1 - spec2).max(0.0),
            );
            fractions.insert(format!("Fraction of reads explained by \"{first}\""), spec1);
            fractions.insert(
                format!("Fraction of reads explained by \"{second}\""),
                spec2,
            );
        }
        StrandednessResult {
            layout,
            fractions,
            sampled: self.sampled,
            paired_reads: self.paired,
            single_reads: self.single,
            ambiguous_reads: self.ambiguous,
        }
    }
}
