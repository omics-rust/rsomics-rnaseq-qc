# Implementation and behavior provenance

This directory is the new coherent product, not a revival of a retired
operation-sized crate. It has no product or path dependencies and remains
`publish = false`.

## Team-owned implementation assets

| Asset and exact local revision | Disposition in this slice |
| --- | --- |
| `rsomics-bam-mapstat` `36a21482b2e76b03147d829469c4fd63b1e87aab` | Refactored count core into `metrics/mapping.rs`; removed reader, independent CLI, and fixed-width-only output |
| `rsomics-bam-strandedness` `6e4e20b08fcd649d6578eb9c906433476d219c2f` | Refactored first-N transcript-overlap classification into `metrics/strandedness.rs`; shared annotation/index and original-record observer |
| `rsomics-genebody-coverage` `1813f1f6e7504be12c4dcd8e8d6f6e9d180ff482` | Refactored genomic-list percentile interpolation into `metrics/sampling.rs`; no indexed per-transcript reader or generated R presentation retained |
| `rsomics-tin` `5e9fa05942c1ee999d87d475d3226638fa402c1a` | Reused exon-walk sampling idea; replaced source-disagreeing all-base dedup, entropy-zero handling, eligibility, summary population, and custom coverage scan |
| `rsomics-bam-read-dist` `2f39541572ebe2b795cba10275498f17333080ca` | Reviewed as asset only; no majority-overlap implementation retained |
| `rsomics-read-distribution` `e169e6c9b5ab8f1b6af64ea633134931102d1e58` | Reviewed as asset only; new shared annotation/midpoint logic follows the actual pinned source |

The user confirmed ownership and permission to reuse these implementations.
Their original repositories remain available on external storage. No upstream
Python source text is embedded in the Rust implementation.

## Pinned behavior source

RSeQC 5.0.5 wheel SHA-256:
`a1fca6d455886b19bc2d58aa1e13b233977a2c835f676ef93bd0c6f4099b8452`.
Git source revision: `0d2e017723c5c4b41eb3b938e1631692d187de1f`.
Five script bodies were matched after the installed shebang; seven qcmodule
files, including SAM, BED, bam_cigar, and mystat, matched byte-for-byte between
the wheel and pinned source. The retained actual source witness is run
`40520155303`, under
`/Volumes/Zane's HDD/rsomics-fixtures/evidence/rnaseq-qc-2026-10-27/run-40520155303/`.
Its single/paired cases contain 88/175 original records. This source witness is
not a performance benchmark or evidence of every default's boundary behavior.

The implemented contracts were checked against the retained script/qcmodule
logic, not inferred from retired crate names or old success claims. Source
defaults use pysam 0.24.0/samtools 1.23.1: excluded flags `0x704`, MAPQ 0,
proper-pair admission, adjusted overlap qualities, depth 8000, BQ 13, no BAQ
without a supplied FASTA. Those mechanisms are delegated to `rsomics-pileup`
rather than copied into this product. The current private candidate pins its
0.9.1 source revision while registry publication is blocked by authentication.

Canonical distribution corrects CIGAR and point semantics; explicitly selected
`rseqc-5.0.5` preserves them. The zero-width midpoint query was verified in
bx-python 0.15.1 `lib/bx/intervals/intersection.pyx` (`_intersect` requires
strict overlap). TIN's `100 * exp(entropy) / sampled_length` is positive at zero
entropy when coverage is positive; the source entropy case is `100/102`.
Eligible-zero transcripts remain in the sample summary. The RSeQC profile
replays independent pileup contexts using the actual sorted source query
bounds and original unadjusted qualities. Identical query windows share replay;
record buffering is bounded by live spans and retired batches, not a fixed
absolute memory ceiling. Canonical coverage/integrity retain one global engine.
The depth-8000 saturation boundary remains uncertified.

Upstream attribution: Wang et al., RSeQC, Bioinformatics (2012),
<https://doi.org/10.1093/bioinformatics/bts356>; Wang et al., TIN, BMC
Bioinformatics (2016), <https://doi.org/10.1186/s12864-016-2539-1>; RSeQC
<https://rseqc.sourceforge.net/> (GPL-licensed upstream; behavioral attribution,
not a relicensing claim); SAM/BAM and BED12 specifications. The product's
team-owned Rust code uses MIT OR Apache-2.0; no clean-room or certification
claim is made.

## Candidate verification ledger

Current source-query work supersedes the older global-pileup candidate below.
The complete original input (9,810,963 records, 58,384 transcript rows) matches
all five actual source metrics at one/four workers: 45,052 eligible TIN rows,
100 coverage bins, total coverage 13,624,401. Reports are byte-identical across
workers. These local boundary checks are not formal performance distributions
or fresh four-native acceptance. A whole-feature review and scoped follow-up
closed the unmapped-record query-watermark finding. The new regression keeps
placed-unmapped coordinates from advancing source query state.

Shared producer fixes retain filtered-record overlap state, reproduce HTSlib's
QName-keyed expiry, and advance emission only for admitted records. Actual
samtools comparisons and consumer regressions cover those contracts. No new
public foundation item, product dependency or copied Python implementation is
introduced. Final registry adoption, four-native acceptance and complete
CPU/RSS/I/O measurements remain open.

Historical green4 checkpoint:

- Read actual pinned behavior and successful single/paired source observations.
- Whole-workflow tests authored first; native execution intentionally deferred
  because the Mac boot APFS exceeds the 80% gate. No RED/GREEN claim is made.
- No local Cargo build/test/metadata/dependency download or installation.
- No new public foundation API or Layer B dependency.
- Standalone `rustfmt --edition 2024` on all Rust sources and workflow tests
  completed successfully on the external source directory. This is not a
  compilation or execution result.
- One independent whole-feature read-only source review found no actionable
  Critical/Important or Minor issue. It explicitly declined to certify native
  compilation/tests, complete helper equivalence, performance/depth-boundary
  equivalence, or publication readiness. No repeat review campaign is planned.
- Green4 run40546148508 at control-plane head
  1224fbf18d7af1b4b76869fc6e41b1611d7bb43b passes all four native platforms:
  format, strict Clippy, eight debug and eight release tests per platform,
  and sixteen complete single/paired one/four-worker source comparisons.
  All numeric outputs match the actual RSeQC witness; report bytes match
  across workers and platforms. Original artifact API/ZIP/log/source/binary
  evidence is retained under external fixtures run-40546148508.
- Executed archive SHA-256 is
  a451612a20af21f1d73e15fadb60f039f18ef9c613a16a1652ba3f7c9177ddb0;
  source-file manifest SHA-256 is
  a7449cb79fc889a3d2fcc4fd6fb8a56813d2c21a70d91481d2332af55f50a562;
  dependency lock SHA-256 is
  c09dea6b5397ca673d7052042d403fa55d2663e45c99672dc6f5ea6f51f308e2.
- The accepted run establishes the bounded synthetic semantic contract only.
  Representative real-input compatibility, CPU/RSS/I/O measurement, final
  metadata and release review remain open. Frozen manifest maturity wording
  predates the successful four-native gate; it is not fresh gate evidence.
