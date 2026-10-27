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
