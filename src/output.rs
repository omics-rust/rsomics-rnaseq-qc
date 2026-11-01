use std::io::Write;

use rsomics_common::{AtomicFile, Context, Result, RsomicsError, reject_output_alias};
use serde::Serialize;
use serde_json::Value;
use sha2::{Digest, Sha256};

use crate::cli::Inputs;
use crate::report::{Manifest, OutputIdentity, Report};

pub(crate) fn write(inputs: &Inputs, report: &Report, mut manifest: Manifest) -> Result<()> {
    let mut outputs = vec![("report.json", json(report)?)];
    if let Some(mapping) = &report.mapping {
        outputs.push(("mapping.tsv", flat_table(mapping)?));
    }
    if let Some(strandedness) = &report.strandedness {
        let mut table = format!("metric\tvalue\nlayout\t{}\nsampled\t{}\npaired_reads\t{}\nsingle_reads\t{}\nambiguous_reads\t{}\n", strandedness.layout, strandedness.sampled, strandedness.paired_reads, strandedness.single_reads, strandedness.ambiguous_reads).into_bytes();
        for (key, value) in &strandedness.fractions {
            writeln!(table, "{key}\t{value}").map_err(RsomicsError::Io)?;
        }
        outputs.push(("strandedness.tsv", table));
    }
    if let Some(distribution) = &report.distribution {
        let mut table = b"group\tbases\ttags\ttags_per_kb\n".to_vec();
        for (group, fields) in &distribution.groups {
            writeln!(
                table,
                "{group}\t{}\t{}\t{}",
                fields.bases, fields.tags, fields.tags_per_kb
            )
            .map_err(RsomicsError::Io)?;
        }
        outputs.push(("distribution.tsv", table));
        outputs.push(("distribution-totals.tsv", flat_table(&distribution.totals)?));
    }
    if let Some(coverage) = &report.coverage {
        let mut table = b"percentile\tcoverage\n".to_vec();
        for (bin, count) in coverage.coverage.iter().enumerate() {
            writeln!(table, "{}\t{count}", bin + 1).map_err(RsomicsError::Io)?;
        }
        outputs.push(("coverage.tsv", table));
    }
    if let Some(integrity) = &report.integrity {
        let mut table =
            b"geneID\tchrom\ttx_start\ttx_end\tTIN\teligible\tsampled_positions\n".to_vec();
        for tx in &integrity.transcripts {
            writeln!(
                table,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                tx.name,
                tx.chrom,
                tx.tx_start,
                tx.tx_end,
                tx.tin,
                tx.eligible,
                tx.sampled_positions
            )
            .map_err(RsomicsError::Io)?;
        }
        outputs.push(("integrity.tsv", table));
        outputs.push(("integrity-summary.tsv", flat_table(&integrity.summary)?));
    }
    if let Some(junctions) = &report.junctions {
        let mut table = b"chrom\tintron_start\tintron_end\tread_count\tannotation\n".to_vec();
        for row in &junctions.rows {
            writeln!(
                table,
                "{}\t{}\t{}\t{}\t{}",
                row.chrom,
                row.start,
                row.end,
                row.read_count,
                row.annotation.name()
            )
            .map_err(RsomicsError::Io)?;
        }
        outputs.push(("junctions.tsv", table));
        outputs.push(("junctions-summary.tsv", flat_table(&junctions.counts)?));
    }
    manifest.outputs = outputs
        .iter()
        .map(|(name, bytes)| OutputIdentity {
            path: (*name).to_owned(),
            bytes: bytes.len(),
            sha256: format!("{:x}", Sha256::digest(bytes)),
        })
        .collect();
    outputs.push(("manifest.json", json(&manifest)?));
    std::fs::create_dir_all(&inputs.out)
        .rs_with_context(|| format!("creating report directory {}", inputs.out.display()))?;
    let mut staged = Vec::with_capacity(outputs.len());
    for (name, bytes) in outputs {
        let path = inputs.out.join(name);
        reject_output_alias(
            &path,
            std::iter::once(inputs.bam.as_path()).chain(inputs.genes.as_deref()),
        )?;
        let mut file = AtomicFile::new(&path)?;
        file.file_mut()
            .write_all(&bytes)
            .rs_with_context(|| format!("writing staged output {}", path.display()))?;
        staged.push(file);
    }
    AtomicFile::commit_all(staged)
}

fn json(value: &impl Serialize) -> Result<Vec<u8>> {
    let mut bytes = serde_json::to_vec_pretty(value).map_err(json_error)?;
    bytes.push(b'\n');
    Ok(bytes)
}

fn flat_table(value: &impl Serialize) -> Result<Vec<u8>> {
    let value = serde_json::to_value(value).map_err(json_error)?;
    let mut table = b"metric\tvalue\n".to_vec();
    if let Value::Object(fields) = value {
        for (metric, value) in fields {
            writeln!(table, "{metric}\t{value}").map_err(RsomicsError::Io)?;
        }
    }
    Ok(table)
}

fn json_error(error: serde_json::Error) -> RsomicsError {
    RsomicsError::Io(std::io::Error::other(error))
}
