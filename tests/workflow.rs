use std::fs;
use std::io::Cursor;
use std::path::Path;
use std::process::Command;

use noodles::sam::alignment::io::Write as _;
use noodles::{bam, bgzf, sam};
use serde_json::Value;

const BED: &str = "chr1\t100\t1300\tpositive\t0\t+\t100\t1300\t0\t2\t300,300\t0,900\nchr1\t2000\t3200\teligible_zero\t0\t-\t2000\t3200\t0\t2\t300,300\t0,900\n";

fn write_bam(path: &Path, records: &str) {
    let text = format!("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:10000\n{records}");
    let mut reader = sam::io::Reader::new(Cursor::new(text.as_bytes()));
    let header = reader.read_header().unwrap();
    let file = fs::File::create(path).unwrap();
    let mut writer = bam::io::Writer::from(bgzf::io::Writer::new(file));
    writer.write_header(&header).unwrap();
    for record in reader.records() {
        writer
            .write_alignment_record(&header, &record.unwrap())
            .unwrap();
    }
    writer.try_finish().unwrap();
}

fn invoke(bam: &Path, genes: &Path, output: &Path, operation: &str) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_rsomics-rnaseq-qc"))
        .arg(operation)
        .args([
            "--bam",
            bam.to_str().unwrap(),
            "--genes",
            genes.to_str().unwrap(),
        ])
        .args([
            "--out",
            output.to_str().unwrap(),
            "--threads",
            "1",
            "--tin-min-starts",
            "0",
        ])
        .output()
        .unwrap()
}

#[test]
fn report_keeps_positive_zero_entropy_and_eligible_zero_and_focused_projection() {
    let dir = tempfile::tempdir().unwrap();
    let bam = dir.path().join("sample.bam");
    let genes = dir.path().join("genes.bed12");
    fs::write(&genes, BED).unwrap();
    write_bam(
        &bam,
        "a\t0\tchr1\t101\t60\t1M\t*\t0\t0\tA\tI\nb\t16\tchr1\t2401\t60\t1M\t*\t0\t0\tA\tI\n",
    );
    let output = dir.path().join("report");
    let result = invoke(&bam, &genes, &output, "report");
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let report: Value =
        serde_json::from_slice(&fs::read(output.join("report.json")).unwrap()).unwrap();
    assert_eq!(report["mapping"]["Total records:"], 2);
    assert_eq!(report["mapping"]["mapq >= mapq_cut (unique):"], 2);
    assert_eq!(report["strandedness"]["sampled"], 2);
    assert_eq!(
        report["integrity"]["eligible_transcripts"],
        serde_json::json!(["positive", "eligible_zero"])
    );
    let rows = report["integrity"]["transcripts"].as_array().unwrap();
    assert!((rows[0]["TIN"].as_f64().unwrap() - 100.0 / 102.0).abs() < 1e-12);
    assert_eq!(rows[1]["TIN"], 0.0);
    assert!(output.join("manifest.json").is_file());
    assert!(output.join("integrity.tsv").is_file());

    let focused = dir.path().join("focused");
    let result = invoke(&bam, &genes, &focused, "integrity");
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let focused: Value =
        serde_json::from_slice(&fs::read(focused.join("report.json")).unwrap()).unwrap();
    assert_eq!(report["integrity"], focused["integrity"]);
    assert!(focused.get("mapping").is_none());
}

#[test]
fn annotation_and_bam_errors_do_not_publish_a_partial_report() {
    let dir = tempfile::tempdir().unwrap();
    let bam = dir.path().join("bad.bam");
    let genes = dir.path().join("genes.bed12");
    fs::write(&bam, "not a BAM").unwrap();
    fs::write(&genes, BED).unwrap();
    let out = dir.path().join("bad-bam");
    assert!(!invoke(&bam, &genes, &out, "report").status.success());
    assert!(!out.join("manifest.json").exists());

    write_bam(&bam, "a\t0\tchr1\t101\t60\t1M\t*\t0\t0\tA\tI\n");
    fs::write(&genes, "chr1\t100\t50\tbroken\n").unwrap();
    let out = dir.path().join("bad-annotation");
    assert!(!invoke(&bam, &genes, &out, "report").status.success());
    assert!(!out.join("report.json").exists());
}

#[test]
fn output_set_failure_preserves_existing_report() {
    let dir = tempfile::tempdir().unwrap();
    let bam = dir.path().join("sample.bam");
    let genes = dir.path().join("genes.bed12");
    fs::write(&genes, BED).unwrap();
    write_bam(&bam, "a\t0\tchr1\t101\t60\t1M\t*\t0\t0\tA\tI\n");
    let output = dir.path().join("report");
    fs::create_dir(&output).unwrap();
    fs::write(output.join("report.json"), "original").unwrap();
    fs::create_dir(output.join("integrity.tsv")).unwrap();
    assert!(!invoke(&bam, &genes, &output, "report").status.success());
    assert_eq!(
        fs::read_to_string(output.join("report.json")).unwrap(),
        "original"
    );
    assert!(!output.join("manifest.json").exists());
}

#[test]
fn canonical_distribution_uses_equals_and_x_while_legacy_retains_pinned_projection() {
    let dir = tempfile::tempdir().unwrap();
    let bam = dir.path().join("sample.bam");
    let genes = dir.path().join("genes.bed12");
    fs::write(&genes, BED).unwrap();
    write_bam(&bam, "a\t0\tchr1\t101\t60\t1=1X\t*\t0\t0\tAA\tII\n");
    let canonical = dir.path().join("canonical");
    let result = invoke(&bam, &genes, &canonical, "distribution");
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let canonical: Value =
        serde_json::from_slice(&fs::read(canonical.join("report.json")).unwrap()).unwrap();
    assert_eq!(canonical["distribution"]["totals"]["Total Tags"], 1);
    assert_eq!(
        canonical["distribution"]["totals"]["Total Assigned Tags"],
        1
    );

    let legacy = dir.path().join("legacy");
    let result = Command::new(env!("CARGO_BIN_EXE_rsomics-rnaseq-qc"))
        .args([
            "distribution",
            "--bam",
            bam.to_str().unwrap(),
            "--genes",
            genes.to_str().unwrap(),
            "--out",
            legacy.to_str().unwrap(),
            "--threads",
            "1",
            "--profile",
            "rseqc-5.0.5",
        ])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let legacy: Value =
        serde_json::from_slice(&fs::read(legacy.join("report.json")).unwrap()).unwrap();
    assert_eq!(legacy["distribution"]["totals"]["Total Reads"], 1);
    assert_eq!(legacy["distribution"]["totals"]["Total Tags"], 0);
}

#[test]
fn filtered_alignment_does_not_cancel_valid_mate_quality_adjustment() {
    let dir = tempfile::tempdir().unwrap();
    let bam = dir.path().join("sample.bam");
    let genes = dir.path().join("genes.bed12");
    let sequence = "A".repeat(36);
    let qualities = ")".repeat(36);
    let record = |flags, pos, mapq, mate, length, hit| {
        format!(
            "pair\t{flags}\tchr1\t{pos}\t{mapq}\t36M\t=\t{mate}\t{length}\t{sequence}\t{qualities}\tNH:i:2\tHI:i:{hit}\n"
        )
    };
    write_bam(
        &bam,
        &[
            record(99, 101, 255, 111, 46, 1),
            record(97, 106, 0, 111, 41, 2),
            record(147, 111, 255, 101, -46, 1),
        ]
        .concat(),
    );
    fs::write(
        &genes,
        "chr1\t100\t300\ttx\t0\t+\t100\t300\t0\t1\t200,\t0,\n",
    )
    .unwrap();
    let output = dir.path().join("coverage");
    let result = invoke(&bam, &genes, &output, "coverage");
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let report: Value =
        serde_json::from_slice(&fs::read(output.join("report.json")).unwrap()).unwrap();
    let mut expected = vec![0; 100];
    expected[4..17].fill(1);
    assert_eq!(report["coverage"]["coverage"], serde_json::json!(expected));
}

#[test]
fn source_queries_restart_overlap_history_and_replay_original_qualities() {
    let dir = tempfile::tempdir().unwrap();
    let bam = dir.path().join("sample.bam");
    let genes = dir.path().join("genes.bed12");
    let record = |flags, pos, length, mate, template| {
        format!(
            "pair\t{flags}\tchr1\t{pos}\t60\t{length}M\t=\t{mate}\t{template}\t{}\t{}\n",
            "A".repeat(length),
            "I".repeat(length)
        )
    };
    write_bam(
        &bam,
        &[
            record(99, 51, 50, 111, 70),
            record(99, 101, 36, 111, 46),
            record(147, 111, 36, 101, -46),
        ]
        .concat(),
    );
    fs::write(
        &genes,
        "chr1\t100\t300\tpositive\t0\t+\t100\t300\t0\t1\t200,\t0,\nchr1\t100\t300\treverse\t0\t-\t100\t300\t0\t1\t200,\t0,\nchr1\t100\t102\tshort\t0\t+\t100\t102\t0\t1\t2,\t0,\nchr1\t140\t340\ttail\t0\t+\t140\t340\t0\t1\t200,\t0,\n",
    )
    .unwrap();
    let canonical = dir.path().join("canonical");
    assert!(invoke(&bam, &genes, &canonical, "report").status.success());
    let canonical: Value =
        serde_json::from_slice(&fs::read(canonical.join("report.json")).unwrap()).unwrap();
    assert_eq!(canonical["coverage"]["coverage"][4], 2);
    let mut baseline = None;
    for (workers, option, global_before) in [
        ("1", "-t", true),
        ("4", "--threads", true),
        ("1", "--threads", false),
        ("4", "-t", false),
    ] {
        let output = dir.path().join(format!("source-{workers}-{global_before}"));
        let mut command = Command::new(env!("CARGO_BIN_EXE_rsomics-rnaseq-qc"));
        if global_before {
            command.args([option, workers]);
        }
        command.args([
            "report",
            "--bam",
            bam.to_str().unwrap(),
            "--genes",
            genes.to_str().unwrap(),
            "--out",
            output.to_str().unwrap(),
            "--profile",
            "rseqc-5.0.5",
            "--tin-min-starts",
            "0",
        ]);
        if !global_before {
            command.args([option, workers]);
        }
        let result = command.output().unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        let bytes = fs::read(output.join("report.json")).unwrap();
        let report: Value = serde_json::from_slice(&bytes).unwrap();
        let mut expected = vec![0; 100];
        expected[..22].fill(1);
        expected[..2].fill(2);
        expected[78..].fill(1);
        assert_eq!(report["coverage"]["coverage"], serde_json::json!(expected));
        for row in &report["integrity"]["transcripts"].as_array().unwrap()[..2] {
            assert!((row["TIN"].as_f64().unwrap() - 2300.0 / 101.0).abs() < 1e-10);
        }
        assert_eq!(report["integrity"]["transcripts"][2]["TIN"], 50.0);
        assert_eq!(report["integrity"]["transcripts"][3]["TIN"], 0.0);
        let manifest: Value =
            serde_json::from_slice(&fs::read(output.join("manifest.json")).unwrap()).unwrap();
        assert_eq!(manifest["bam_traversals"], 1);
        assert_eq!(manifest["original_records"], 3);
        assert_eq!(
            manifest["inflate_workers"],
            workers.parse::<usize>().unwrap()
        );
        if let Some(previous) = &baseline {
            assert_eq!(&bytes, previous);
        }
        baseline = Some(bytes);
    }
}

#[test]
fn placed_unmapped_coordinates_do_not_advance_source_query_windows() {
    let dir = tempfile::tempdir().unwrap();
    let bam = dir.path().join("sample.bam");
    let genes = dir.path().join("genes.bed12");
    write_bam(
        &bam,
        "a\t0\tchr1\t101\t60\t1M\t*\t0\t0\tA\tI\nunmapped\t4\tchr1\t1001\t0\t1M\t*\t0\t0\tA\tI\nb\t0\tchr1\t111\t60\t1M\t*\t0\t0\tA\tI\n",
    );
    fs::write(
        &genes,
        "chr1\t100\t300\ttx\t0\t+\t100\t300\t0\t1\t200,\t0,\n",
    )
    .unwrap();
    let mut baseline = None;
    for workers in ["1", "4"] {
        let output = dir.path().join(format!("source-{workers}"));
        let result = Command::new(env!("CARGO_BIN_EXE_rsomics-rnaseq-qc"))
            .args([
                "report",
                "--bam",
                bam.to_str().unwrap(),
                "--genes",
                genes.to_str().unwrap(),
                "--out",
                output.to_str().unwrap(),
                "--profile",
                "rseqc-5.0.5",
                "--threads",
                workers,
                "--tin-min-starts",
                "0",
            ])
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        let bytes = fs::read(output.join("report.json")).unwrap();
        let report: Value = serde_json::from_slice(&bytes).unwrap();
        assert_eq!(report["mapping"]["Total records:"], 3);
        assert_eq!(report["mapping"]["Unmapped reads:"], 1);
        let mut expected = vec![0; 100];
        expected[4] = 1;
        assert_eq!(report["coverage"]["coverage"], serde_json::json!(expected));
        let row = &report["integrity"]["transcripts"][0];
        assert_eq!(row["eligible"], true);
        assert!((row["TIN"].as_f64().unwrap() - 200.0 / 101.0).abs() < 1e-12);
        let manifest: Value =
            serde_json::from_slice(&fs::read(output.join("manifest.json")).unwrap()).unwrap();
        assert_eq!(manifest["original_records"], 3);
        assert_eq!(manifest["bam_traversals"], 1);
        if let Some(previous) = &baseline {
            assert_eq!(&bytes, previous);
        }
        baseline = Some(bytes);
    }
}

#[test]
fn junctions_match_captured_source_and_share_the_original_record_transaction() {
    let dir = tempfile::tempdir().unwrap();
    let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/junctions");
    let run = |case: &str, name: &str, operation: &str, extra: &[&str]| {
        let out = dir.path().join(name);
        let mut command = Command::new(env!("CARGO_BIN_EXE_rsomics-rnaseq-qc"));
        command
            .arg(operation)
            .arg("--bam")
            .arg(fixtures.join(case).join("input.bam"))
            .arg("--genes")
            .arg(fixtures.join(case).join("genes.bed12"))
            .arg("--out")
            .arg(&out);
        if !extra.contains(&"--profile") {
            command.args(["--profile", "rseqc-5.0.5"]);
        }
        if !extra.contains(&"--threads") {
            command.args(["--threads", "1"]);
        }
        let result = command.args(extra).output().unwrap();
        (out, result)
    };
    let load = |out: &Path, file| -> Value {
        serde_json::from_slice(&fs::read(out.join(file)).unwrap()).unwrap()
    };
    let assert_ok = |result: &std::process::Output| {
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
    };
    let counter_names = [
        ("total_events", "Total splicing  Events"),
        ("known_events", "Known Splicing Events"),
        ("partial_novel_events", "Partial Novel Splicing Events"),
        ("novel_events", "Novel Splicing Events"),
        ("filtered_events", "Filtered Splicing Events"),
        ("total_junctions", "Total splicing  Junctions"),
        ("known_junctions", "Known Splicing Junctions"),
        (
            "partial_novel_junctions",
            "Partial Novel Splicing Junctions",
        ),
        ("novel_junctions", "Novel Splicing Junctions"),
    ];
    for (case, records) in [("informative-single", 16), ("informative-paired", 30)] {
        let expected: Value = serde_json::from_slice(
            &fs::read(fixtures.join(case).join("observation.json")).unwrap(),
        )
        .unwrap();
        assert_eq!(expected["raw_returncode"], 0);
        let stderr = fs::read_to_string(fixtures.join(case).join("process.stderr")).unwrap();
        for (_, label) in counter_names {
            let observed = stderr
                .lines()
                .find_map(|line| {
                    let (name, value) = line.split_once(':')?;
                    (name == label).then(|| value.trim().parse::<u64>().unwrap())
                })
                .unwrap();
            assert_eq!(expected["counters"][label][0], observed);
        }
        let source_rows: Vec<Vec<String>> =
            fs::read_to_string(fixtures.join(case).join("observed.junction.xls"))
                .unwrap()
                .lines()
                .skip(1)
                .map(|line| line.split_ascii_whitespace().map(str::to_owned).collect())
                .collect();
        let mut baseline = None;
        for workers in ["1", "4"] {
            let name = format!("{case}-{workers}");
            let (out, result) = run(case, &name, "junctions", &["--threads", workers]);
            assert_ok(&result);
            let report = load(&out, "report.json");
            assert_eq!(report["metrics"], serde_json::json!(["junctions"]));
            assert!(report.get("mapping").is_none());
            assert_eq!(report["junctions"]["rows"], expected["rows"]);
            for (native, source) in counter_names {
                assert_eq!(
                    report["junctions"]["counts"][native],
                    expected["counters"][source][0]
                );
            }
            let detail = fs::read_to_string(out.join("junctions.tsv")).unwrap();
            let rows: Vec<Vec<String>> = detail
                .lines()
                .skip(1)
                .map(|line| line.split_ascii_whitespace().map(str::to_owned).collect())
                .collect();
            assert_eq!(rows, source_rows);
            let summary = fs::read_to_string(out.join("junctions-summary.tsv")).unwrap();
            assert_eq!(summary.lines().count(), 10);
            for line in summary.lines().skip(1) {
                let (metric, value) = line.split_once('\t').unwrap();
                assert_eq!(
                    report["junctions"]["counts"][metric],
                    value.parse::<u64>().unwrap()
                );
            }
            let manifest = load(&out, "manifest.json");
            assert_eq!(manifest["original_records"], records);
            assert_eq!(manifest["bam_traversals"], 1);
            assert_eq!(
                manifest["inflate_workers"],
                workers.parse::<usize>().unwrap()
            );
            assert_eq!(manifest["policy"]["junctions"]["excluded_flags"], 0x704);
            assert_eq!(manifest["policy"]["junctions"]["minimum_mapq"], 30);
            assert_eq!(manifest["policy"]["junctions"]["minimum_intron"], 50);
            assert!(
                !manifest["omitted_metrics"]
                    .as_array()
                    .unwrap()
                    .contains(&serde_json::json!("junctions"))
            );
            for name in ["junctions.tsv", "junctions-summary.tsv"] {
                let identity = manifest["outputs"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .find(|row| row["path"] == name)
                    .unwrap();
                let body = fs::read(out.join(name)).unwrap();
                use sha2::{Digest, Sha256};
                assert_eq!(identity["sha256"], format!("{:x}", Sha256::digest(&body)));
                assert_eq!(identity["bytes"], body.len());
            }
            let bytes = fs::read(out.join("report.json")).unwrap();
            if let Some(previous) = &baseline {
                assert_eq!(&bytes, previous);
            }
            baseline = Some(bytes);
            let (full, result) = run(
                case,
                &format!("{case}-report-{workers}"),
                "report",
                &["--threads", workers, "--metrics", "mapping,junctions"],
            );
            assert_ok(&result);
            let full = load(&full, "report.json");
            assert_eq!(full["junctions"], report["junctions"]);
            let (default, result) = run(
                case,
                &format!("{case}-standard-{workers}"),
                "report",
                &["--threads", workers],
            );
            assert!(!result.status.success());
            assert!(
                String::from_utf8_lossy(&result.stderr).contains("unambiguous case-normalized")
            );
            assert!(!default.join("manifest.json").exists());
        }
    }
    let bam = dir.path().join("shared-stream.bam");
    let genes = dir.path().join("shared-stream.bed12");
    fs::write(&genes, BED).unwrap();
    write_bam(
        &bam,
        &format!(
            "a\t0\tchr1\t381\t60\t20M600N20M\t*\t0\t0\t{}\t{}\nb\t2048\tchr1\t381\t60\t20M600N20M\t*\t0\t0\t{}\t{}\nc\t0\tchr1\t401\t60\t10M590N20M\t*\t0\t0\t{}\t{}\n",
            "A".repeat(40),
            "I".repeat(40),
            "A".repeat(40),
            "I".repeat(40),
            "A".repeat(30),
            "I".repeat(30),
        ),
    );
    for profile in ["canonical", "rseqc-5.0.5"] {
        let mut baseline = None;
        for workers in ["1", "4"] {
            let execute = |name: &str, operation, extra: &[&str]| {
                let out = dir
                    .path()
                    .join(format!("shared-{name}-{profile}-{workers}"));
                let result = Command::new(env!("CARGO_BIN_EXE_rsomics-rnaseq-qc"))
                    .arg(operation)
                    .arg("--bam")
                    .arg(&bam)
                    .arg("--genes")
                    .arg(&genes)
                    .arg("--out")
                    .arg(&out)
                    .args([
                        "--profile",
                        profile,
                        "--threads",
                        workers,
                        "--tin-min-starts",
                        "0",
                    ])
                    .args(extra)
                    .output()
                    .unwrap();
                assert_ok(&result);
                (
                    load(&out, "report.json"),
                    fs::read(out.join("report.json")).unwrap(),
                )
            };
            let (full, bytes) = execute(
                "full",
                "report",
                &[
                    "--metrics",
                    "mapping,strandedness,distribution,coverage,integrity,junctions",
                ],
            );
            let (default, _) = execute("standard", "report", &[]);
            let (focused, _) = execute("focused", "junctions", &[]);
            assert_eq!(full["junctions"], focused["junctions"]);
            assert_eq!(full["junctions"]["counts"]["total_events"], 3);
            assert_eq!(full["junctions"]["counts"]["known_events"], 2);
            assert_eq!(full["junctions"]["counts"]["partial_novel_events"], 1);
            assert_eq!(
                full["junctions"]["rows"][1],
                serde_json::json!({"chrom":"chr1","start":410,"end":1000,"read_count":1,"annotation":"partial_novel"})
            );
            assert!(default.get("junctions").is_none());
            assert_eq!(default["metrics"].as_array().unwrap().len(), 5);
            for metric in [
                "mapping",
                "strandedness",
                "distribution",
                "coverage",
                "integrity",
            ] {
                assert_eq!(full[metric], default[metric]);
            }
            if let Some(previous) = &baseline {
                assert_eq!(&bytes, previous);
            }
            baseline = Some(bytes);
        }
    }
    let (canonical, result) = run(
        "informative-single",
        "canonical-junctions",
        "junctions",
        &["--profile", "canonical"],
    );
    assert_ok(&result);
    let canonical = load(&canonical, "report.json");
    let rows = canonical["junctions"]["rows"].as_array().unwrap();
    assert_eq!(rows.len(), 7);
    assert_eq!(
        rows[5],
        serde_json::json!({"chrom":"chr1","start":530,"end":580,"read_count":1,"annotation":"complete_novel"})
    );
    assert_eq!(
        rows[6],
        serde_json::json!({"chrom":"CHR1","start":120,"end":170,"read_count":1,"annotation":"complete_novel"})
    );
    assert_eq!(canonical["junctions"]["counts"]["known_events"], 5);
    assert_eq!(canonical["junctions"]["counts"]["novel_events"], 4);
    for profile in ["canonical", "rseqc-5.0.5"] {
        let (out, result) = run(
            "case-alias-only",
            &format!("alias-{profile}"),
            "junctions",
            &["--profile", profile],
        );
        assert!(!result.status.success());
        assert!(String::from_utf8_lossy(&result.stderr).contains("absent from the BAM header"));
        assert!(!out.join("manifest.json").exists());
        for (case, events, filtered, source_exit) in [("no-N", 0, 0, 0), ("all-short", 1, 1, 1)] {
            let expected: Value = serde_json::from_slice(
                &fs::read(fixtures.join(case).join("observation.json")).unwrap(),
            )
            .unwrap();
            assert_eq!(expected["raw_returncode"], source_exit);
            let (out, result) = run(
                case,
                &format!("{case}-{profile}"),
                "junctions",
                &["--profile", profile],
            );
            assert_ok(&result);
            let report = load(&out, "report.json");
            assert_eq!(report["junctions"]["rows"], serde_json::json!([]));
            for (name, _) in counter_names {
                let value = match name {
                    "total_events" => events,
                    "filtered_events" => filtered,
                    _ => 0,
                };
                assert_eq!(report["junctions"]["counts"][name], value);
            }
            assert_eq!(
                fs::read_to_string(out.join("junctions.tsv"))
                    .unwrap()
                    .lines()
                    .count(),
                1
            );
            assert_eq!(
                fs::read_to_string(out.join("junctions-summary.tsv"))
                    .unwrap()
                    .lines()
                    .count(),
                10
            );
            assert_eq!(load(&out, "manifest.json")["status"], "complete");
        }
    }
    let (out, result) = run(
        "informative-single",
        "selected",
        "report",
        &["--metrics", "junctions,mapping"],
    );
    assert_ok(&result);
    let report = load(&out, "report.json");
    assert_eq!(
        report["metrics"],
        serde_json::json!(["mapping", "junctions"])
    );
    assert!(report.get("coverage").is_none());
    for (name, args) in [
        ("zero-intron", vec!["--min-intron", "0"]),
        ("focused-metrics", vec!["--metrics", "junctions"]),
    ] {
        assert!(
            !run("informative-single", name, "junctions", &args)
                .1
                .status
                .success()
        );
    }
    assert!(
        !run(
            "informative-single",
            "unknown-metric",
            "report",
            &["--metrics", "imaginary"]
        )
        .1
        .status
        .success()
    );
    let (out, result) = run(
        "informative-single",
        "mapq255",
        "junctions",
        &["--mapq", "255"],
    );
    assert_ok(&result);
    assert_eq!(
        load(&out, "report.json")["junctions"]["counts"]["total_events"],
        1
    );
    let (out, result) = run(
        "informative-single",
        "min51",
        "junctions",
        &["--min-intron", "51"],
    );
    assert_ok(&result);
    assert_eq!(
        load(&out, "report.json")["junctions"]["counts"]["filtered_events"],
        9
    );
    let out = dir.path().join("atomic");
    fs::create_dir(&out).unwrap();
    fs::write(out.join("report.json"), "original-report").unwrap();
    fs::write(out.join("junctions.tsv"), "original-detail").unwrap();
    fs::write(out.join("manifest.json"), "original-manifest").unwrap();
    fs::create_dir(out.join("junctions-summary.tsv")).unwrap();
    assert!(
        !run(
            "informative-single",
            "atomic",
            "report",
            &["--metrics", "mapping,junctions"]
        )
        .1
        .status
        .success()
    );
    assert_eq!(
        fs::read_to_string(out.join("report.json")).unwrap(),
        "original-report"
    );
    assert_eq!(
        fs::read_to_string(out.join("junctions.tsv")).unwrap(),
        "original-detail"
    );
    assert_eq!(
        fs::read_to_string(out.join("manifest.json")).unwrap(),
        "original-manifest"
    );
    assert!(!out.join("mapping.tsv").exists());
}
