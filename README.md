# rsomics-rnaseq-qc

Private, unpublished candidate for one coordinated RNA-seq QC workflow. The
declared slice is `report`, `mapping`, `strandedness`, `distribution`, `coverage`,
and `integrity`; later metrics are not exposed.

```sh
rsomics-rnaseq-qc report --bam sample.bam --genes genes.bed12 --out sample-qc
rsomics-rnaseq-qc report --bam sample.bam --genes genes.bed12 --out legacy-qc --profile rseqc-5.0.5
rsomics-rnaseq-qc integrity --bam sample.bam --genes genes.bed12 --out sample-tin
```

Input is one coordinate-sorted BAM and validated BED12 annotation. No BAM index
is needed. The original records are traversed once using `rsomics-bamio`. Global
`-t N` / `--threads N` controls BGZF inflation and RSeQC query-replay workers and
may appear before or after the subcommand. Omitting it retains the available-core
default. Record-level metrics observe the original records. The existing
`rsomics-pileup` engine supplies separate pileup filtering and paired-overlap
quality adjustment. Focused commands select
the same accumulators; `mapping` may omit `--genes`.

`canonical` is the default distribution policy: correct full-CIGAR geometry and
half-open point membership. Explicit `rseqc-5.0.5` retains the pinned source's
M-only blocks, reference-advancing soft clips, ignored `=`/`X`, and zero-width
midpoint query. It also restarts pileup state per source transcript query,
sharing buffered original records and replaying independent contexts in
parallel. Canonical coverage and integrity use a global pileup instead. Both
profiles read BAM bytes once; query buffering depends on live alignment and
transcript spans rather than a whole-BAM copy. Other metric definitions use
the pinned source contracts in both profiles. Depth-8000 saturation remains an
uncertified compatibility boundary, not a claim of general RSeQC equivalence.

Each invocation transactionally publishes `report.json`, selected metric TSVs,
and `manifest.json` using `rsomics-common`. The manifest records filters,
annotation SHA-256, deterministic first-N population, actual metric selection,
omitted metrics, input size/mtime, workers, elapsed time, known limitations, and
output checksums. It does not reread BAM bytes only to hash the input. Existing
unselected tables are not deleted; the manifest is authoritative for the output
set. `--json` uses the common structured command envelope on stdout; errors are
non-zero. CLI presentation is supplied by `rsomics-help`.

TIN eligibility uses strictly more than `--tin-min-starts` distinct starts,
including duplicates and all MAPQ after excluding unmapped/secondary/QC-fail.
The summary includes every eligible transcript, including zero scores, and uses
population standard deviation. All-zero eligible data and an empty eligible
population are distinguished by `eligible_transcripts`. Coverage percentiles
interpolate the genomic exon-base list, collapse repeated positions as the
source does, and orient each transcript before aggregation. For short eligible
transcripts this can produce fewer than 100 unique output bins.

Malformed annotation, contigs absent from the BAM header, invalid BAM/CIGAR,
unsorted mapped records, and malformed paired-segment identities fail loudly.
BED rows retain their order and repeated gene names; TIN retains repeated rows,
while coverage deduplicates the source's full transcript identity.

The source-query implementation passes hosted native Linux/macOS x86_64/aarch64
formatting, strict Clippy, and 11 tests in debug and release. Runs
[40590326515](https://github.com/omics-rust/rsomics-world/actions/runs/40590326515)
and [40597129172](https://github.com/omics-rust/rsomics-world/actions/runs/40597129172)
compare all five metrics with actual RSeQC single/paired semantic fixtures at
one and four workers. The latter changes only shared thread-argument plumbing.
Those synthetic report bytes agree across platforms and workers.

Whole-input performance run
[40592936910](https://github.com/omics-rust/rsomics-world/actions/runs/40592936910)
measures the earlier source-query binary on a 9,810,963-record paired BAM and
58,384 transcript rows. All five source metric gates pass. On one Linux host,
warm-cache median wall time is 1557.27 seconds for the five serial RSeQC scripts,
54.94 seconds for one native worker, and 28.95 seconds for four workers. The
serial speedup is 28.34x; observed worker gain is 1.90x. Median process peak RSS
is 21.11% lower at one worker. Native uses more output blocks; this is not a
cold-I/O, universal-dataset or linear-scaling claim. Real-input Mac/Linux TIN
rounding differs within the fixed numeric bound, unlike the synthetic fixtures'
byte equality. Original measurements and their exact binary remain separately
identified; later code does not inherit a benchmark by name alone.

Registry adoption of the tested pileup foundation and final publication review
remain required before release. The candidate stays unpublished. See
[MERGE_RECORD.md](MERGE_RECORD.md) for behavior and implementation provenance.
