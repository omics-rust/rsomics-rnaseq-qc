use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::ops::Range;
use std::path::Path;

use rsomics_common::{Context, Result, RsomicsError};

#[derive(Clone, Debug)]
pub(crate) struct Transcript {
    pub name: String,
    pub chrom: String,
    pub start: u32,
    pub end: u32,
    pub cds_start: u32,
    pub cds_end: u32,
    pub reverse: bool,
    pub exons: Vec<Range<u32>>,
}

impl Transcript {
    pub fn mrna_len(&self) -> usize {
        self.exons
            .iter()
            .map(|exon| (exon.end - exon.start) as usize)
            .sum()
    }
}

#[derive(Clone, Debug)]
pub(crate) struct Annotation {
    pub transcripts: Vec<Transcript>,
}

impl Annotation {
    pub fn read(path: &Path) -> Result<Self> {
        let file = File::open(path)
            .rs_with_context(|| format!("opening BED12 annotation {}", path.display()))?;
        Self::from_reader(BufReader::new(file), path)
    }

    fn from_reader(reader: impl BufRead, path: &Path) -> Result<Self> {
        let mut transcripts = Vec::new();
        for (line_index, line) in reader.lines().enumerate() {
            let context = || format!("{}:{}", path.display(), line_index + 1);
            let line = line.rs_with_context(context)?;
            let line = line.trim();
            if line.is_empty()
                || line.starts_with('#')
                || matches!(
                    line.split_ascii_whitespace().next(),
                    Some("track" | "browser")
                )
            {
                continue;
            }
            transcripts.push(parse_transcript(line).rs_with_context(context)?);
        }
        if transcripts.is_empty() {
            return Err(RsomicsError::InvalidInput(format!(
                "{}: BED12 annotation contains no transcripts",
                path.display()
            )));
        }
        Ok(Self { transcripts })
    }
}

fn parse_transcript(line: &str) -> Result<Transcript> {
    let fields: Vec<_> = line.split_ascii_whitespace().collect();
    if fields.len() < 12 {
        return Err(invalid("expected at least 12 BED columns"));
    }
    let coordinate = |index: usize, name: &str| -> Result<u32> {
        fields[index]
            .parse()
            .map_err(|error| invalid(format!("invalid {name}: {error}")))
    };
    let start = coordinate(1, "chromStart")?;
    let end = coordinate(2, "chromEnd")?;
    if start >= end {
        return Err(invalid("chromStart must be less than chromEnd"));
    }
    let reverse = match fields[5] {
        "+" => false,
        "-" => true,
        _ => return Err(invalid("strand must be + or -")),
    };
    let cds_start = coordinate(6, "thickStart")?;
    let cds_end = coordinate(7, "thickEnd")?;
    if cds_start < start || cds_start > cds_end || cds_end > end {
        return Err(invalid(
            "thickStart and thickEnd must be ordered within the transcript bounds",
        ));
    }
    let block_count = fields[9]
        .parse::<usize>()
        .map_err(|error| invalid(format!("invalid blockCount: {error}")))?;
    if block_count == 0 {
        return Err(invalid("blockCount must be positive"));
    }
    let sizes = parse_list(fields[10], "blockSizes")?;
    let offsets = parse_list(fields[11], "blockStarts")?;
    if sizes.len() != block_count || offsets.len() != block_count {
        return Err(invalid(
            "blockCount must equal both blockSizes and blockStarts lengths",
        ));
    }
    let mut exons: Vec<Range<u32>> = Vec::with_capacity(block_count);
    for (size, offset) in sizes.into_iter().zip(offsets) {
        if size == 0 {
            return Err(invalid("exon blocks must be nonempty"));
        }
        let exon_start = start
            .checked_add(offset)
            .ok_or_else(|| invalid("block start overflows the BED coordinate range"))?;
        let exon_end = exon_start
            .checked_add(size)
            .ok_or_else(|| invalid("block end overflows the BED coordinate range"))?;
        if exon_end > end {
            return Err(invalid("exon block extends beyond chromEnd"));
        }
        if exons
            .last()
            .is_some_and(|previous| previous.end > exon_start)
        {
            return Err(invalid("exon blocks must be sorted and nonoverlapping"));
        }
        exons.push(exon_start..exon_end);
    }
    if exons[0].start != start || exons[block_count - 1].end != end {
        return Err(invalid(
            "first and last exon blocks must meet the transcript bounds",
        ));
    }
    Ok(Transcript {
        name: fields[3].into(),
        chrom: fields[0].into(),
        start,
        end,
        cds_start,
        cds_end,
        reverse,
        exons,
    })
}

fn parse_list(field: &str, name: &str) -> Result<Vec<u32>> {
    field
        .strip_suffix(',')
        .unwrap_or(field)
        .split(',')
        .map(|value| {
            value
                .parse()
                .map_err(|error| invalid(format!("invalid {name}: {error}")))
        })
        .collect()
}

fn invalid(message: impl Into<String>) -> RsomicsError {
    RsomicsError::InvalidInput(message.into())
}

#[derive(Debug)]
struct IndexedTranscript {
    start: u32,
    end: u32,
    transcript: usize,
}

#[derive(Debug, Default)]
struct ReferenceTranscripts {
    entries: Vec<IndexedTranscript>,
    prefix_max_end: Vec<u32>,
}

#[derive(Debug)]
pub(crate) struct TranscriptIndex {
    references: Vec<ReferenceTranscripts>,
}

impl TranscriptIndex {
    pub fn new(annotation: &Annotation, reference_names: &[String]) -> Result<Self> {
        let mut reference_ids = HashMap::with_capacity(reference_names.len());
        for (reference_id, name) in reference_names.iter().enumerate() {
            if reference_ids.insert(name.as_str(), reference_id).is_some() {
                return Err(invalid(format!("duplicate BAM reference name: {name}")));
            }
        }
        let mut references: Vec<ReferenceTranscripts> = reference_names
            .iter()
            .map(|_| ReferenceTranscripts::default())
            .collect();
        for (transcript, row) in annotation.transcripts.iter().enumerate() {
            let reference_id = reference_ids.get(row.chrom.as_str()).ok_or_else(|| {
                invalid(format!(
                    "annotation transcript {} uses contig {} absent from the BAM header",
                    row.name, row.chrom
                ))
            })?;
            references[*reference_id].entries.push(IndexedTranscript {
                start: row.start,
                end: row.end,
                transcript,
            });
        }
        for reference in &mut references {
            reference
                .entries
                .sort_unstable_by_key(|entry| (entry.start, entry.end, entry.transcript));
            let mut max_end = 0;
            reference.prefix_max_end = reference
                .entries
                .iter()
                .map(|entry| {
                    max_end = max_end.max(entry.end);
                    max_end
                })
                .collect();
        }
        Ok(Self { references })
    }

    pub fn overlapping(&self, reference_id: usize, start: u32, end: u32) -> Vec<usize> {
        let Some(reference) = self.references.get(reference_id) else {
            return Vec::new();
        };
        if start >= end {
            return Vec::new();
        }
        let first = reference
            .prefix_max_end
            .partition_point(|&max_end| max_end <= start);
        let last = reference.entries.partition_point(|entry| entry.start < end);
        let mut transcripts: Vec<_> = reference.entries[first..last]
            .iter()
            .filter(|entry| entry.end > start)
            .map(|entry| entry.transcript)
            .collect();
        transcripts.sort_unstable();
        transcripts
    }
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;
    use std::path::Path;

    use super::{Annotation, TranscriptIndex};

    #[test]
    fn bed12_rows_retain_identity_and_index_nested_transcripts() {
        let bed = b"chr1\t100\t1000\tshared\t0\t+\t100\t1000\t0\t2\t50,50,\t0,850,\nchr1\t200\t300\tshared\t0\t-\t200\t300\t0\t1\t100,\t0,\n";
        let annotation =
            Annotation::from_reader(Cursor::new(bed), Path::new("genes.bed12")).unwrap();
        assert_eq!(annotation.transcripts.len(), 2);
        assert_eq!(annotation.transcripts[0].mrna_len(), 100);
        assert!(annotation.transcripts[1].reverse);
        let index = TranscriptIndex::new(&annotation, &["chr1".into()]).unwrap();
        assert_eq!(index.overlapping(0, 250, 251), [0, 1]);
        assert_eq!(index.overlapping(0, 999, 1000), [0]);
        assert!(index.overlapping(0, 1000, 1001).is_empty());
        assert!(TranscriptIndex::new(&annotation, &["CHR1".into()]).is_err());
    }

    #[test]
    fn malformed_bed12_fails_at_the_input_row() {
        for invalid in [
            "chr1\t100\t300\ttx\t0\t+\t100\t300\t0\t2\t100,100,\t0,50,",
            "chr1\t100\t300\ttx\t0\t+\t99\t300\t0\t1\t200,\t0,",
            "chr1\t100\t300\ttx\t0\t.\t100\t300\t0\t1\t200,\t0,",
            "chr1\t100\t300\ttx\t0\t+\t100\t300\t0\t2\t200,\t0,",
            "chr1\t100\t300\ttx\t0\t+\t100\t300\t0\t1\t0,\t0,",
        ] {
            let input = format!("# annotation\n{invalid}\n");
            let error = Annotation::from_reader(Cursor::new(input), Path::new("genes.bed12"))
                .unwrap_err()
                .to_string();
            assert!(error.contains("genes.bed12:2"), "{error}");
        }
    }
}
