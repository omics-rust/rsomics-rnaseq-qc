# rsomics-rnaseq-qc

One coordinated RNA-seq QC workflow. The first-release slice is
`report`, `mapping`, `strandedness`, `distribution`, `coverage`,
`integrity`, and `junctions`; later metrics are not exposed. The default report
keeps the established five metrics. Junctions requires explicit selection.

See [CI](https://github.com/omics-rust/rsomics-rnaseq-qc/actions/workflows/ci.yml)
and [GitHub releases](https://github.com/omics-rust/rsomics-rnaseq-qc/releases)
for source/native delivery. crates.io publication remains disabled.

## Install and first workflow

Native archives are target-specific for Linux and macOS on x86_64 and aarch64.
Verify the release SHA256SUMS before extracting an archive and adding its
single executable to PATH. The installed binary needs no Python or R runtime.
Project licenses and complete third-party notices accompany each archive.

The accompanying `rsomics-rnaseq-qc-VERSION-source.tar` is an uncompressed,
deterministic Git source archive, not a Cargo registry `.crate` or an offline
vendored bundle. It preserves the exact original Cargo.toml/Cargo.lock and
public pileup Git dependency. For version0.1.0:

```sh
tar -xf rsomics-rnaseq-qc-0.1.0-source.tar
cd rsomics-rnaseq-qc-0.1.0-source
cargo install --path . --locked --root /path/to/install-prefix
```

Source installation requires Rust1.91 or later and network/cache access to
the locked public dependencies, including `rsomics-pileup0.9.1` at
`93986e17f4c3a62716b0b6434737c3f46368478a`. Choose your own install prefix and
add its `bin` directory to PATH. Do not replace this Git pin with a registry
version as installation housekeeping.

From the source tree, this executable example decodes an already-qualified,
team-generated175-record paired fixture and runs the default five metrics:

```sh
python3 scripts/delivery.py fixture --out demo-input
rsomics-rnaseq-qc report --bam demo-input/paired.bam --genes demo-input/genes.bed12 \
  --out demo-qc --threads 4
rsomics-rnaseq-qc junctions --bam tests/fixtures/junctions/informative-paired/input.bam \
  --genes tests/fixtures/junctions/informative-paired/genes.bed12 \
  --out demo-junctions --profile rseqc-5.0.5 --threads 4
```

Python is used only by that optional fixture-decoding helper, not by the
installed product. The fixture JSON retains the exact original source/member/
size/hash bindings; no BAM rewriting or reindexing is performed. The existing
informative junction fixture has intentionally case-ambiguous header names,
so use its focused junction command rather than the default-five report.

For real inputs:

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
Current-owner sourcea618e6ab is separately qualified on all four native classes
by [41282068281](https://github.com/omics-rust/rsomics-world/actions/runs/41282068281):
12 tests per debug/release profile, captured ordered junction outputs and
default-five integration. This is a source/native correctness gate, not a
claim that Python ran on every target or that an earlier timing measures the
new addition.

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
identified; the measured old executable is
`48091abdb22d5d91e18a733848798b8a9ae3dd83d3d8117dd77feca08fd13707`.
Later code does not inherit a benchmark by name alone.

The separate complete-junction run
[41291319388](https://github.com/omics-rust/rsomics-world/actions/runs/41291319388)
is root-adopted for current sourcea618e6ab and measured Linux executable
`eb472af80dfa26a00819448f895ccbcceb570de3ad629ac81f9507c7534282ae`.
All9,810,963 original records and58,384 annotation rows are retained; all
66,283 ordered source junction rows/counters and native repeat/worker output
bytes match. Median source/native1/native4 wall is17.338/2.973/1.451s on that
same warm-input Linux workload. Native1/source wall ratio is.1715 and median
process peak RSS is46.3MB versus111.2MB. Four workers improve wall2.05x but
use more CPU than one worker (4.362 versus2.963 CPU-seconds); native scientific
outputs are12.73MB versus2.52MB source output. This is not lower output-volume,
cold-cache, aggregate/heap-RSS or universal/linear-scaling evidence. Do not sum
the historical five-metric timing and this separate junction timing into an
unmeasured six-metric pipeline. Fresh delivery binaries are identified in
BUILD.json and are not relabelled as either historical measured executable.

Fragments, bias, duplication, saturation/junction saturation, bigWig ingress,
Picard/refFlat/rRNA profiles, SAM/CRAM ingress, R/PDF/HTML rendering and junction
BED/Interact exports remain outside this release. GitHub source/native delivery
does not remove the separate tested-foundation/registry/publication gates.
See [MERGE_RECORD.md](MERGE_RECORD.md) for historical behavior and merge
provenance and [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) for the actual
dependency notice inventory and explicitly labelled team-owned supplements.
