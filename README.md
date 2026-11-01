# rsomics-rnaseq-qc

Private, unpublished candidate for one coordinated RNA-seq QC workflow. The
declared slice is `report`, `mapping`, `strandedness`, `distribution`, `coverage`,
`integrity`, and `junctions`; later metrics are not exposed. The default report
keeps the established five metrics. Junctions requires explicit selection.

```sh
rsomics-rnaseq-qc report --bam sample.bam --genes genes.bed12 --out sample-qc
rsomics-rnaseq-qc report --bam sample.bam --genes genes.bed12 --out legacy-qc --profile rseqc-5.0.5
rsomics-rnaseq-qc integrity --bam sample.bam --genes genes.bed12 --out sample-tin
rsomics-rnaseq-qc junctions --bam sample.bam --genes genes.bed12 --out sample-junctions --profile rseqc-5.0.5
rsomics-rnaseq-qc report --bam sample.bam --genes genes.bed12 --out sample-qc --metrics mapping,strandedness,distribution,coverage,integrity,junctions
```

Input is one coordinate-sorted BAM and validated BED12 annotation. No BAM index
is needed. The original records are traversed once using `rsomics-bamio`. Global
`-t N` / `--threads N` controls BGZF inflation and RSeQC query-replay workers and
may appear before or after the subcommand. Omitting it retains the available-core
default. Record-level metrics observe the original records. The existing
`rsomics-pileup` engine supplies separate pileup filtering and paired-overlap
quality adjustment. Focused commands select
the same accumulators; `mapping` may omit `--genes`.

Junctions observes each original record before pileup filtering or ownership
transfer, excludes flags `0x704` and MAPQ below `--mapq` (default 30), and retains
supplementary alignments. Each N occurrence is one event, including each mate;
total events include N operations shorter than `--min-intron` (positive,
default 50). Accepted junctions retain first BAM encounter order and per-event
support. Donor and acceptor knownness are independent annotation-site queries,
not intron-pair or strand membership. `junctions.tsv` contains every distinct
accepted junction; `junctions-summary.tsv` contains all nine event/distinct
counters. Both coordinates are zero-based interbase boundaries, with the
intron represented as `[intron_start,intron_end)`.

Canonical junctions uses exact contig identities and advances reference for
M/D/N/=/X, never soft clips. Explicit `rseqc-5.0.5` junctions advances only
M/D/N, pools uppercase contig identities, and displays every `CHR` substring
as `chr`. This differs from legacy distribution, which advances soft clips.
Both profiles still require strict valid BED12 and exact annotation contigs in
the BAM header. Distribution retains its separate rejection of ambiguous
case-normalized header names; select `--metrics mapping,junctions` for a
junction report on an otherwise admitted `chr1`/`CHR1` header. Junction
accumulation is serial on the shared record stream; workers accelerate BGZF
inflation, not a claimed parallel junction scan. Focused junctions initializes
no pileup or transcript-query replay context.

Valid no-N or all-short inputs produce complete atomic zero-result tables,
accurate total/filtered counters and exit 0. This intentionally differs from
RSeQC 5.0.5, whose captured no-N call leaves empty files and whose all-short
call exits 1 with partial files. Native R/PDF, BED/Interact exports, SAM/CRAM
ingress and junction saturation are not implemented or exposed.

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

The earlier five-metric source-query candidate passed hosted native Linux/macOS x86_64/aarch64
formatting, strict Clippy, and 11 tests in debug and release. Runs
[40590326515](https://github.com/omics-rust/rsomics-world/actions/runs/40590326515)
and [40597129172](https://github.com/omics-rust/rsomics-world/actions/runs/40597129172)
compare all five metrics with actual RSeQC single/paired semantic fixtures at
one and four workers. The latter changes only shared thread-argument plumbing.
Those synthetic report bytes agree across platforms and workers.
These earlier runs do not qualify the added junction workflow. Its source
fixture authority is the separately adopted six-call RSeQC 5.0.5 run
`41277337904/1` at `b4cde25cde7a7d52d4b3947004929991b897a727`;
captured data and per-file provenance are in `tests/fixtures/junctions`.
New four-platform qualification and full-input junction performance remain
separate pending gates; no earlier timings are inherited by this addition.

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
