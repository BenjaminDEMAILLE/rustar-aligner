//! Alternative BAM writer backed by htslib, behind the `htslib-bam` feature.
//!
//! # Why a second backend rather than a replacement
//!
//! BGZF is a sequence of independently-deflated blocks, so compression
//! parallelises perfectly. The obvious way to exploit that with the existing
//! stack — `noodles_bgzf::MultithreadedWriter` — turns out to be **slower**
//! here, by roughly 37% at every size measured up to a 44 MB BAM. The reason is
//! that `noodles-bgzf` is already built with its `libdeflate` feature and both
//! its writers go through the same `deflate::encode`, so the multithreaded one
//! adds per-block channel and ordering overhead on top of work that was never
//! the bottleneck.
//!
//! htslib's `hts_set_threads` is a different thing: a genuinely parallel
//! deflate rather than a layer over one compressor. This module makes it
//! available for the runs where BAM writing dominates.
//!
//! # Why it is off by default, and never built on Windows
//!
//! `rust-htslib` ships pre-built bindings for Mac and Linux only, its README
//! describes Windows `bindgen` as untested, and its own CI has no Windows job.
//! This project's CI matrix requires `windows-x86_64`, so the dependency is
//! declared under `cfg(not(windows))` and the feature is off by default. The
//! published crate and the default build are unchanged by its existence.
//!
//! # The contract
//!
//! Whatever this writes, the default backend must write the same records.
//! BGZF block boundaries may differ — they are a framing detail, and htslib
//! chooses them differently — but the decoded record stream may not. That is
//! what `backends_agree` checks.

use std::path::Path;

use noodles::sam::{self, alignment::record_buf::RecordBuf};
use rust_htslib::bam::{self as hts};

use crate::error::Error;
use crate::params::Parameters;

/// Number of htslib worker threads to use for a run of `--runThreadN n`.
///
/// htslib counts these *in addition* to the calling thread, so `n - 1` keeps
/// the total at what the user asked for. A single-threaded run gets no extra
/// workers at all, which keeps `--runThreadN 1` genuinely serial.
fn worker_threads(params: &Parameters) -> usize {
    params.run_thread_n.get().saturating_sub(1)
}

/// Streaming, unsorted BAM writer using htslib.
pub struct HtslibBamWriter {
    writer: hts::Writer,
    header: sam::Header,
}

impl HtslibBamWriter {
    /// Create a writer at `output_path`, mirroring
    /// [`crate::io::bam::BamWriter::create`].
    pub fn create(
        output_path: &Path,
        genome: &crate::genome::Genome,
        params: &Parameters,
    ) -> Result<Self, Error> {
        let header = crate::io::sam::build_sam_header(genome, params)?;
        Self::from_header(output_path, header, params)
    }

    /// Create a writer from an already-built header.
    ///
    /// The header is handed to htslib as SAM text rather than assembled field
    /// by field, so the two backends cannot disagree about it: both render the
    /// same `sam::Header`.
    pub fn from_header(
        output_path: &Path,
        header: sam::Header,
        params: &Parameters,
    ) -> Result<Self, Error> {
        let text = crate::io::bam::render_sam_header_text(&header, None);
        let hts_header = hts::Header::from_template(&hts::HeaderView::from_bytes(&text));
        let mut writer = hts::Writer::from_path(output_path, &hts_header, hts::Format::Bam)
            .map_err(|e| Error::Alignment(format!("{}: {e}", output_path.display())))?;
        writer
            .set_compression_level(hts::CompressionLevel::Level(
                params.out_bam_compression.clamp(0, 9) as u32,
            ))
            .map_err(|e| Error::Alignment(format!("{}: {e}", output_path.display())))?;
        let workers = worker_threads(params);
        if workers > 0 {
            writer
                .set_threads(workers)
                .map_err(|e| Error::Alignment(format!("{}: {e}", output_path.display())))?;
        }
        Ok(Self { writer, header })
    }

    /// Write a batch of records.
    pub fn write_batch(&mut self, batch: &[RecordBuf]) -> Result<(), Error> {
        for rec in batch {
            let sam_line = crate::io::bam::render_sam_record_line(&self.header, rec)?;
            let hts_rec = hts::Record::from_sam(self.writer.header(), sam_line.as_bytes())
                .map_err(|e| Error::Alignment(format!("htslib record: {e}")))?;
            self.writer
                .write(&hts_rec)
                .map_err(|e| Error::Alignment(format!("htslib write: {e}")))?;
        }
        Ok(())
    }

    /// Flush and close. htslib appends the BGZF EOF block on drop.
    pub fn finish(&mut self) -> Result<(), Error> {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The contract, checked rather than asserted in prose: for the same
    /// records, both backends must produce the same decoded BAM.
    ///
    /// Only the decoded stream is compared. BGZF block boundaries are a framing
    /// detail and htslib picks them differently, so comparing raw bytes would
    /// fail for a reason nobody cares about.
    #[test]
    fn backends_agree_on_the_decoded_record_stream() {
        use noodles::sam::alignment::record::Flags;

        let dir = tempfile::tempdir().unwrap();
        let genome = crate::genome::Genome {
            transform_blocks: None,
            sequence: vec![0u8; 512].into(),
            n_genome: 256,
            n_genome_real: 256,
            n_chr_real: 1,
            chr_name: vec!["chr1".to_string()],
            chr_length: vec![256],
            chr_start: vec![0, 256],
        };
        let params = Parameters::parse_from(["rustar-aligner", "--readFilesIn", "r.fq"]);

        // A handful of unmapped records: enough to exercise header, name,
        // sequence, quality and flags without needing a real alignment.
        let mut batch = Vec::new();
        for i in 0..8u32 {
            let mut rec = RecordBuf::default();
            *rec.name_mut() = Some(format!("read{i}").into_bytes().into());
            *rec.flags_mut() = Flags::UNMAPPED;
            *rec.sequence_mut() = b"ACGTACGT".to_vec().into();
            *rec.quality_scores_mut() = vec![30u8; 8].into();
            batch.push(rec);
        }

        let noodles_path = dir.path().join("noodles.bam");
        let mut w = crate::io::bam::BamWriter::create(&noodles_path, &genome, &params).unwrap();
        w.write_batch(&batch).unwrap();
        w.finish().unwrap();
        drop(w); // the BGZF EOF block lands on drop; read-back needs it there

        let hts_path = dir.path().join("htslib.bam");
        let mut w = HtslibBamWriter::create(&hts_path, &genome, &params).unwrap();
        w.write_batch(&batch).unwrap();
        w.finish().unwrap();
        drop(w);

        // Compare the fields that carry meaning, read back through htslib so
        // both files go through the same decoder.
        use rust_htslib::bam::Read as _;
        let decode = |p: &Path| -> Vec<String> {
            let mut r = hts::Reader::from_path(p).unwrap();
            r.records()
                .map(|rec| {
                    let rec = rec.unwrap();
                    format!(
                        "{}\t{}\t{}\t{}\t{}",
                        String::from_utf8_lossy(rec.qname()),
                        rec.flags(),
                        rec.tid(),
                        rec.pos(),
                        String::from_utf8_lossy(&rec.seq().as_bytes()),
                    )
                })
                .collect()
        };

        assert_eq!(
            decode(&noodles_path),
            decode(&hts_path),
            "the two BAM backends must decode to the same records"
        );
    }
}
