use rsomics_bamio::raw::RawRecord;
use serde::Serialize;

#[derive(Debug, Default, Serialize)]
pub(crate) struct MappingResult {
    #[serde(rename = "Total records:")]
    total: u64,
    #[serde(rename = "QC failed:")]
    qc_fail: u64,
    #[serde(rename = "Optical/PCR duplicate:")]
    duplicate: u64,
    #[serde(rename = "Non primary hits")]
    secondary: u64,
    #[serde(rename = "Unmapped reads:")]
    unmapped: u64,
    #[serde(rename = "mapq < mapq_cut (non-unique):")]
    low_mapq: u64,
    #[serde(rename = "mapq >= mapq_cut (unique):")]
    high_mapq: u64,
    #[serde(rename = "Read-1:")]
    read1: u64,
    #[serde(rename = "Read-2:")]
    read2: u64,
    #[serde(rename = "Reads map to '+':")]
    plus: u64,
    #[serde(rename = "Reads map to '-':")]
    minus: u64,
    #[serde(rename = "Non-splice reads:")]
    nonsplice: u64,
    #[serde(rename = "Splice reads:")]
    splice: u64,
    #[serde(rename = "Reads mapped in proper pairs:")]
    proper: u64,
    #[serde(rename = "Proper-paired reads map to different chrom:")]
    different_chromosome: u64,
}

impl MappingResult {
    pub fn observe(&mut self, record: &RawRecord, cigar: &[(u8, u32)], mapq_cut: u8) {
        self.total += 1;
        let flags = record.flags();
        if flags & 0x200 != 0 {
            self.qc_fail += 1;
        } else if flags & 0x400 != 0 {
            self.duplicate += 1;
        } else if flags & 0x100 != 0 {
            self.secondary += 1;
        } else if flags & 0x4 != 0 {
            self.unmapped += 1;
        } else if record.mapping_quality() < mapq_cut {
            self.low_mapq += 1;
        } else {
            self.high_mapq += 1;
            self.read1 += u64::from(flags & 0x40 != 0);
            self.read2 += u64::from(flags & 0x80 != 0);
            if flags & 0x10 != 0 {
                self.minus += 1;
            } else {
                self.plus += 1;
            }
            if cigar.iter().any(|&(kind, _)| kind == 3) {
                self.splice += 1;
            } else {
                self.nonsplice += 1;
            }
            if flags & 0x2 != 0 {
                self.proper += 1;
                self.different_chromosome += u64::from(
                    record.reference_sequence_id() != record.mate_reference_sequence_id(),
                );
            }
        }
    }
}
