use std::num::NonZero;
use std::path::PathBuf;

use clap::{Args, Parser, Subcommand, ValueEnum};
use rsomics_common::{OutputArgs, ThreadArgs, ToolMeta};
use serde::Serialize;

#[derive(Debug, Parser)]
#[command(
    name = "rsomics-rnaseq-qc",
    version,
    about = "Coordinated RNA-seq alignment and transcript quality control"
)]
pub struct Cli {
    #[command(subcommand)]
    command: Command,
    #[command(flatten)]
    threads: ThreadArgs,
    #[command(flatten)]
    output: OutputArgs,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Run selected shared-scan metrics; defaults to the established five.
    Report(ReportInputs),
    /// Summarize mapping, flags, pairing, strands, and splice counts.
    Mapping(Inputs),
    /// Infer protocol evidence from the first eligible annotated reads.
    Strandedness(Inputs),
    /// Assign alignment-block midpoints to annotation categories.
    Distribution(Inputs),
    /// Aggregate gene-body coverage from 5-prime to 3-prime.
    Coverage(Inputs),
    /// Calculate transcript integrity and the eligible-transcript summary.
    Integrity(Inputs),
    /// Annotate supported splice junctions and summarize event and distinct counts.
    Junctions(Inputs),
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, ValueEnum)]
pub(crate) enum Metric {
    Mapping,
    Strandedness,
    Distribution,
    Coverage,
    Integrity,
    Junctions,
}

#[derive(Debug, Args)]
struct ReportInputs {
    #[command(flatten)]
    inputs: Inputs,
    /// Comma-separated metrics; defaults to mapping,strandedness,distribution,coverage,integrity.
    #[arg(
        long,
        value_enum,
        value_delimiter = ',',
        help_heading = "Metric policy"
    )]
    metrics: Option<Vec<Metric>>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, ValueEnum)]
pub(crate) enum Profile {
    #[value(name = "canonical")]
    #[serde(rename = "canonical")]
    Canonical,
    #[value(name = "rseqc-5.0.5")]
    #[serde(rename = "rseqc-5.0.5")]
    Rseqc505,
}

#[derive(Debug, Args)]
pub(crate) struct Inputs {
    /// Coordinate-sorted BAM input. No BAM index is required.
    #[arg(long, value_name = "BAM", help_heading = "Inputs")]
    pub bam: PathBuf,
    /// Validated BED12 annotation; required except for mapping.
    #[arg(long, value_name = "BED12", help_heading = "Inputs")]
    pub genes: Option<PathBuf>,
    /// Directory for the transactional JSON, TSV, and manifest output set.
    #[arg(long, value_name = "DIRECTORY", help_heading = "Outputs")]
    pub out: PathBuf,
    /// Canonical geometry or pinned RSeQC metric-specific CIGAR and pileup policies.
    #[arg(
        long,
        value_enum,
        default_value = "canonical",
        help_heading = "Compatibility"
    )]
    pub profile: Profile,
    #[arg(skip)]
    pub threads: Option<NonZero<usize>>,
    /// Mapping, strandedness and junctions MAPQ cutoff; distribution and pileup retain their own policy.
    #[arg(long, default_value_t = 30, help_heading = "Metric policy")]
    pub mapq: u8,
    /// Minimum retained intron length for junctions; shorter N events remain in totals.
    #[arg(long, default_value = "50", help_heading = "Metric policy")]
    pub min_intron: NonZero<u32>,
    /// First N eligible annotation-overlapping reads for strandedness; no random sampling.
    #[arg(long, default_value = "200000", help_heading = "Metric policy")]
    pub strand_sample_size: NonZero<u64>,
    /// TIN requires strictly more distinct starts than this cutoff.
    #[arg(long, default_value_t = 10, help_heading = "Metric policy")]
    pub tin_min_starts: usize,
    /// Target TIN sampling density; exon boundaries also contribute.
    #[arg(long, default_value = "100", help_heading = "Metric policy")]
    pub tin_sample_size: NonZero<usize>,
    /// Minimum spliced transcript length for gene-body coverage.
    #[arg(long, default_value = "100", value_parser = clap::value_parser!(u32).range(100..), help_heading = "Metric policy")]
    pub min_mrna_length: u32,
}

impl Cli {
    pub fn run(self) -> std::process::ExitCode {
        let (mut inputs, metrics) = match self.command {
            Command::Report(report) => {
                let metrics = report.metrics.as_deref().map_or_else(
                    crate::report::Selection::standard,
                    crate::report::Selection::from_metrics,
                );
                (report.inputs, metrics)
            }
            Command::Mapping(inputs) => (inputs, crate::report::Selection::only(Metric::Mapping)),
            Command::Strandedness(inputs) => {
                (inputs, crate::report::Selection::only(Metric::Strandedness))
            }
            Command::Distribution(inputs) => {
                (inputs, crate::report::Selection::only(Metric::Distribution))
            }
            Command::Coverage(inputs) => (inputs, crate::report::Selection::only(Metric::Coverage)),
            Command::Integrity(inputs) => {
                (inputs, crate::report::Selection::only(Metric::Integrity))
            }
            Command::Junctions(inputs) => {
                (inputs, crate::report::Selection::only(Metric::Junctions))
            }
        };
        inputs.threads = self.threads.requested();
        rsomics_common::run(
            &self.output,
            ToolMeta {
                name: "rsomics-rnaseq-qc",
                version: env!("CARGO_PKG_VERSION"),
            },
            || crate::report::run(&inputs, metrics, !self.output.json),
        )
    }
}
