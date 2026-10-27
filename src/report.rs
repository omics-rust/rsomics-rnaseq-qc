use std::fs;
use std::io::{Read, Write};
use std::num::NonZero;
use std::path::{Path, PathBuf};
use std::time::{Instant, UNIX_EPOCH};

use rsomics_bamio::raw::{self, RawRecord};
use rsomics_common::{Context, Result, RsomicsError};
use rsomics_pileup::{FlagFilter, PileupEngine, PileupOptions, RecordFilter};
use serde::Serialize;
use sha2::{Digest, Sha256};

use crate::annotation::{Annotation, TranscriptIndex};
use crate::cli::{Inputs, Profile};
use crate::metrics::distribution::{Distribution, DistributionResult};
use crate::metrics::mapping::MappingResult;
use crate::metrics::sampling::{CoverageResult, IntegrityResult, SampledMetrics};
use crate::metrics::strandedness::{Strandedness, StrandednessResult};

const METRICS: [&str; 5] = [
    "mapping",
    "strandedness",
    "distribution",
    "coverage",
    "integrity",
];

#[derive(Clone, Copy)]
pub(crate) struct Selection {
    mapping: bool,
    strandedness: bool,
    distribution: bool,
    coverage: bool,
    integrity: bool,
}

impl Selection {
    pub fn all() -> Self {
        Self {
            mapping: true,
            strandedness: true,
            distribution: true,
            coverage: true,
            integrity: true,
        }
    }
    pub fn only(metric: &str) -> Self {
        Self {
            mapping: metric == "mapping",
            strandedness: metric == "strandedness",
            distribution: metric == "distribution",
            coverage: metric == "coverage",
            integrity: metric == "integrity",
        }
    }
    fn names(self) -> Vec<&'static str> {
        [
            self.mapping,
            self.strandedness,
            self.distribution,
            self.coverage,
            self.integrity,
        ]
        .into_iter()
        .zip(METRICS)
        .filter_map(|(selected, name)| selected.then_some(name))
        .collect()
    }
    fn needs_annotation(self) -> bool {
        self.strandedness || self.distribution || self.coverage || self.integrity
    }
}

#[derive(Serialize)]
pub(crate) struct Report {
    schema_version: &'static str,
    profile: Profile,
    sample: String,
    metrics: Vec<&'static str>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub mapping: Option<MappingResult>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub strandedness: Option<StrandednessResult>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub distribution: Option<DistributionResult>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub coverage: Option<CoverageResult>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub integrity: Option<IntegrityResult>,
}

#[derive(Serialize)]
pub(crate) struct Manifest {
    schema_version: &'static str,
    product: &'static str,
    version: &'static str,
    status: &'static str,
    profile: Profile,
    input: InputIdentity,
    annotation: Option<InputIdentity>,
    references: Vec<Reference>,
    selected_metrics: Vec<&'static str>,
    omitted_metrics: Vec<&'static str>,
    bam_traversals: u8,
    original_records: u64,
    inflate_workers: usize,
    elapsed_seconds: f64,
    policy: MetricPolicy,
    known_limitations: Vec<&'static str>,
    pub outputs: Vec<OutputIdentity>,
}

#[derive(Serialize)]
struct Reference {
    name: String,
    length: u64,
}

#[derive(Serialize)]
struct InputIdentity {
    path: PathBuf,
    bytes: u64,
    modified_unix_seconds: Option<u64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    sha256: Option<String>,
}

#[derive(Serialize)]
struct MetricPolicy {
    mapping_mapq: u8,
    mapping_exclusive_precedence: [&'static str; 5],
    strandedness_sampling: &'static str,
    strandedness_sample_limit: u64,
    strandedness_overlap: &'static str,
    strandedness_excluded_flags: u16,
    distribution_excluded_flags: u16,
    distribution_mapq: u8,
    distribution_projection: &'static str,
    distribution_midpoint_query: &'static str,
    coverage_minimum_mrna_length: u32,
    coverage_percentiles: &'static str,
    tin_sample_size: usize,
    tin_eligibility: &'static str,
    tin_eligibility_excluded_flags: u16,
    tin_minimum_distinct_starts_exclusive: usize,
    tin_summary_population: &'static str,
    tin_background_subtraction: bool,
    pileup_excluded_flags: u16,
    pileup_minimum_mapq: u8,
    pileup_minimum_adjusted_base_quality: u8,
    pileup_include_anomalous_pairs: bool,
    pileup_adjust_overlaps: bool,
    pileup_maximum_depth: usize,
    pileup_baq: bool,
}

#[derive(Serialize)]
pub(crate) struct OutputIdentity {
    pub path: String,
    pub bytes: usize,
    pub sha256: String,
}

#[derive(Serialize)]
pub(crate) struct RunSummary {
    report: PathBuf,
    manifest: PathBuf,
    metrics: Vec<&'static str>,
    original_records: u64,
}

pub(crate) fn run(inputs: &Inputs, selection: Selection, emit_plain: bool) -> Result<RunSummary> {
    let started = Instant::now();
    if selection.needs_annotation() && inputs.genes.is_none() {
        return Err(RsomicsError::ConfigError(
            "--genes BED12 is required for the selected metric plan".into(),
        ));
    }
    let annotation = match &inputs.genes {
        Some(path) => Annotation::read(path)?,
        None => Annotation {
            transcripts: Vec::new(),
        },
    };
    let annotation_identity = inputs
        .genes
        .as_ref()
        .map(|path| identify(path, true))
        .transpose()?;
    let input_identity = identify(&inputs.bam, false)?;
    let workers = inputs
        .threads
        .unwrap_or_else(|| std::thread::available_parallelism().unwrap_or(NonZero::<usize>::MIN));
    let mut reader = rsomics_bamio::open_with_workers(&inputs.bam, workers)?;
    let header = reader.read_header().map_err(RsomicsError::Io)?;
    let references = header
        .reference_sequences()
        .iter()
        .map(|(name, reference)| {
            let name = std::str::from_utf8(name.as_ref()).map_err(|error| {
                RsomicsError::InvalidInput(format!("reference name is not UTF-8: {error}"))
            })?;
            Ok(Reference {
                name: name.to_owned(),
                length: usize::from(reference.length()) as u64,
            })
        })
        .collect::<Result<Vec<_>>>()?;
    let reference_names = references
        .iter()
        .map(|reference| reference.name.clone())
        .collect::<Vec<_>>();
    let index = TranscriptIndex::new(&annotation, &reference_names)?;
    {
        for tx in &annotation.transcripts {
            let reference = references
                .iter()
                .find(|reference| reference.name == tx.chrom)
                .ok_or_else(|| {
                    RsomicsError::InvalidInput(format!(
                        "annotation contig {} is absent from BAM header",
                        tx.chrom
                    ))
                })?;
            if u64::from(tx.end) > reference.length {
                return Err(RsomicsError::InvalidInput(format!(
                    "transcript {} ends beyond contig {} length {}",
                    tx.name, tx.chrom, reference.length
                )));
            }
        }
    }
    let mut mapping = selection.mapping.then(MappingResult::default);
    let mut strandedness = selection.strandedness.then(Strandedness::default);
    let mut distribution = if selection.distribution {
        Some(Distribution::new(&annotation, &reference_names)?)
    } else {
        None
    };
    let options = PileupOptions {
        filter: RecordFilter {
            flags: FlagFilter {
                skip_any_set: 0x704,
                ..FlagFilter::default()
            },
            minimum_mapping_quality: 0,
            include_anomalous_pairs: false,
        },
        adjust_overlaps: true,
        maximum_depth_per_source: Some(8000),
    };
    let mut sampled_pileup = if selection.coverage || selection.integrity {
        Some((
            PileupEngine::new(references.iter().map(|reference| reference.length), options),
            SampledMetrics::new(
                &annotation,
                &reference_names,
                selection.coverage,
                selection.integrity,
                inputs.min_mrna_length as usize,
                inputs.tin_sample_size.get(),
                inputs.tin_min_starts,
            ),
        ))
    } else {
        None
    };
    let mut record = RawRecord::default();
    let mut cigar = Vec::with_capacity(32);
    let mut original_records = 0;
    let mut last_mapped = None;
    loop {
        if raw::read_record(reader.get_mut(), &mut record)? == 0 {
            break;
        }
        original_records += 1;
        record.decode_cigar_into(&mut cigar)?;
        validate_record(&record, &cigar, &references, &mut last_mapped)?;
        if let Some(mapping) = &mut mapping {
            mapping.observe(&record, &cigar, inputs.mapq);
        }
        if let Some(strandedness) = &mut strandedness {
            strandedness.observe(
                &record,
                &annotation,
                &index,
                inputs.strand_sample_size.get(),
                inputs.mapq,
            )?;
        }
        if let Some(distribution) = &mut distribution {
            distribution.observe(&record, &cigar, inputs.profile == Profile::Rseqc505)?;
        }
        if let Some((pileup, sampled)) = &mut sampled_pileup {
            sampled.observe_start(&record, &index);
            pileup
                .push(std::mem::take(&mut record))
                .map_err(pileup_error)?;
            pileup.drain(|column| {
                sampled.observe_column(column);
                Ok::<_, RsomicsError>(())
            })?;
        }
    }
    if let Some((pileup, sampled)) = &mut sampled_pileup {
        pileup.finish().map_err(pileup_error)?;
        pileup.drain(|column| {
            sampled.observe_column(column);
            Ok::<_, RsomicsError>(())
        })?;
    }
    let sample = inputs
        .bam
        .file_stem()
        .map(|name| name.to_string_lossy().into_owned())
        .unwrap_or_else(|| "sample".into());
    let (coverage, integrity) = sampled_pileup
        .map(|(_, sampled)| sampled.finish(&annotation, &sample))
        .unwrap_or((None, None));
    let report = Report {
        schema_version: "1.0",
        profile: inputs.profile,
        sample,
        metrics: selection.names(),
        mapping,
        strandedness: strandedness.map(Strandedness::finish),
        distribution: distribution.map(Distribution::finish),
        coverage,
        integrity,
    };
    let omitted_metrics = METRICS
        .into_iter()
        .filter(|metric| !report.metrics.contains(metric))
        .chain([
            "junctions",
            "fragments",
            "bias",
            "duplication",
            "saturation",
            "bigwig-coverage",
            "picard-refflat-rrna",
        ])
        .collect();
    let manifest = Manifest {
        schema_version: "1.0",
        product: "rsomics-rnaseq-qc",
        version: env!("CARGO_PKG_VERSION"),
        status: "complete",
        profile: inputs.profile,
        input: input_identity,
        annotation: annotation_identity,
        references,
        selected_metrics: selection.names(),
        omitted_metrics,
        bam_traversals: 1,
        original_records,
        inflate_workers: workers.get(),
        elapsed_seconds: started.elapsed().as_secs_f64(),
        policy: MetricPolicy {
            mapping_mapq: inputs.mapq,
            mapping_exclusive_precedence: [
                "QC-fail",
                "duplicate",
                "secondary",
                "unmapped",
                "low-MAPQ",
            ],
            strandedness_sampling: "first-N filtered reads with transcript-span overlap; ambiguous overlaps consume the population",
            strandedness_sample_limit: inputs.strand_sample_size.get(),
            strandedness_overlap: "case-sensitive transcript spans queried by read-start plus query length",
            strandedness_excluded_flags: 0x704,
            distribution_excluded_flags: 0x704,
            distribution_mapq: 0,
            distribution_projection: if inputs.profile == Profile::Rseqc505 {
                "legacy M-only blocks; S/D/N advance reference; I/=/X ignored"
            } else {
                "full-CIGAR aligned M/=/X blocks; I/P do not split contiguous reference coverage; D/N split"
            },
            distribution_midpoint_query: if inputs.profile == Profile::Rseqc505 {
                "legacy zero-width query; region.start < midpoint < region.end"
            } else {
                "half-open point; region.start <= midpoint < region.end"
            },
            coverage_minimum_mrna_length: inputs.min_mrna_length,
            coverage_percentiles: "1..100 genomic-list interpolation with ties-to-even rounding; duplicate sampled positions collapse; negative transcripts reverse before aggregation",
            tin_sample_size: inputs.tin_sample_size.get(),
            tin_eligibility: "distinct starts inside transcript span; exclude unmapped/secondary/QC-fail, include duplicates and all MAPQ",
            tin_eligibility_excluded_flags: 0x304,
            tin_minimum_distinct_starts_exclusive: inputs.tin_min_starts,
            tin_summary_population: "all eligible transcripts including zero TIN; population standard deviation",
            tin_background_subtraction: false,
            pileup_excluded_flags: 0x704,
            pileup_minimum_mapq: 0,
            pileup_minimum_adjusted_base_quality: 13,
            pileup_include_anomalous_pairs: false,
            pileup_adjust_overlaps: true,
            pileup_maximum_depth: 8000,
            pileup_baq: false,
        },
        known_limitations: vec![
            "Global depth-8000 admission is not yet proven equivalent to RSeQC's transcript-restarted pileup at high-depth query boundaries",
            "The retained source witness does not exercise depth-8000 truncation, low base qualities, or anomalous pairs",
            "Canonical distribution intentionally differs from pinned legacy CIGAR and zero-width midpoint quirks",
            "Malformed annotation, absent annotation contigs, and malformed paired-segment identities fail instead of being silently skipped",
            "BAM identity records path, size, and mtime; BAM content is not reread solely to hash it",
            "This unpublished candidate has not passed the four-native-platform or representative-performance release gates",
        ],
        outputs: Vec::new(),
    };
    crate::output::write(inputs, &report, manifest)?;
    if emit_plain {
        writeln!(
            std::io::stdout().lock(),
            "{}: {} original records; {}",
            report.metrics.join(", "),
            original_records,
            inputs.out.join("report.json").display()
        )
        .map_err(RsomicsError::Io)?;
    }
    Ok(RunSummary {
        report: inputs.out.join("report.json"),
        manifest: inputs.out.join("manifest.json"),
        metrics: report.metrics,
        original_records,
    })
}

fn validate_record(
    record: &RawRecord,
    cigar: &[(u8, u32)],
    references: &[Reference],
    last: &mut Option<(i32, i32)>,
) -> Result<()> {
    if record.flags() & 4 != 0 {
        return Ok(());
    }
    let invalid = |reason: &str| {
        RsomicsError::InvalidInput(format!(
            "read {}: {reason}",
            String::from_utf8_lossy(record.name())
        ))
    };
    let reference = usize::try_from(record.reference_sequence_id())
        .ok()
        .and_then(|id| references.get(id))
        .ok_or_else(|| invalid("mapped reference ID is absent from header"))?;
    let start = u64::try_from(record.alignment_start())
        .map_err(|_| invalid("mapped alignment start is absent"))?;
    if start >= reference.length {
        return Err(invalid("alignment start is outside its reference"));
    }
    if !cigar.is_empty() {
        let query_span = cigar
            .iter()
            .filter(|(kind, _)| matches!(kind, 0 | 1 | 4 | 7 | 8))
            .map(|(_, len)| u64::from(*len))
            .sum::<u64>();
        if query_span != record.sequence_len() as u64 {
            return Err(invalid("CIGAR query span differs from sequence length"));
        }
        let reference_span = cigar
            .iter()
            .filter(|(kind, _)| matches!(kind, 0 | 2 | 3 | 7 | 8))
            .map(|(_, len)| u64::from(*len))
            .sum::<u64>();
        if start + reference_span > reference.length {
            return Err(invalid("CIGAR exceeds reference length"));
        }
    }
    let coordinate = (record.reference_sequence_id(), record.alignment_start());
    if last.is_some_and(|previous| coordinate < previous) {
        return Err(invalid("BAM is not coordinate sorted"));
    }
    *last = Some(coordinate);
    Ok(())
}

fn identify(path: &Path, hash: bool) -> Result<InputIdentity> {
    let metadata = fs::metadata(path)
        .rs_with_context(|| format!("reading input metadata {}", path.display()))?;
    if !metadata.is_file() {
        return Err(RsomicsError::InvalidInput(format!(
            "{} is not a regular file",
            path.display()
        )));
    }
    let sha256 = if hash {
        let mut file =
            fs::File::open(path).rs_with_context(|| format!("opening {}", path.display()))?;
        let mut digest = Sha256::new();
        let mut buffer = [0u8; 65536];
        loop {
            let n = file.read(&mut buffer).map_err(RsomicsError::Io)?;
            if n == 0 {
                break;
            }
            digest.update(&buffer[..n]);
        }
        Some(format!("{:x}", digest.finalize()))
    } else {
        None
    };
    Ok(InputIdentity {
        path: path.to_owned(),
        bytes: metadata.len(),
        modified_unix_seconds: metadata
            .modified()
            .ok()
            .and_then(|time| time.duration_since(UNIX_EPOCH).ok())
            .map(|duration| duration.as_secs()),
        sha256,
    })
}

fn pileup_error(error: rsomics_pileup::PileupError) -> RsomicsError {
    RsomicsError::InvalidInput(error.to_string())
}
