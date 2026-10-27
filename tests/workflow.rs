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
