/// BAM output writer with noodles (streaming, unsorted)
use crate::error::Error;
use crate::genome::Genome;
use crate::params::Parameters;
use crate::quant::transcriptome::TranscriptomeIndex;
use byteorder::{LittleEndian, WriteBytesExt};
use noodles::sam::alignment::io::Write as SamWrite;
use noodles::sam::alignment::record_buf::RecordBuf;
use noodles::{bam, bgzf, sam};
use std::ffi::CString;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

/// Buffer for BAM records built by parallel threads
#[derive(Default)]
pub struct BufferedBamRecords {
    pub records: Vec<RecordBuf>,
}

impl BufferedBamRecords {
    /// Create new buffer with capacity
    pub fn new() -> Self {
        Self {
            records: Vec::with_capacity(10000),
        }
    }

    /// Add a record to the buffer
    pub fn push(&mut self, record: RecordBuf) {
        self.records.push(record);
    }
}

/// Convert STAR's `--outBAMcompression` integer (-1..9) to a noodles level.
///
/// STAR mapping: -1 or 0 = uncompressed, 1-9 = deflate levels, default 1.
fn bgzf_compression(level: i32) -> bgzf::io::writer::CompressionLevel {
    use bgzf::io::writer::CompressionLevel;
    match level {
        n if n <= 0 => CompressionLevel::NONE,
        n if n >= 9 => CompressionLevel::BEST,
        n => CompressionLevel::try_from(n as u8).unwrap_or_default(),
    }
}

/// Create a BGZF writer with the given STAR compression level.
///
/// Deflate runs on the global rayon pool, which `run()` sizes from
/// `--runThreadN`, while a dedicated thread writes the finished frames in
/// order. Block boundaries, frame layout and the EOF marker are the same as
/// the single-threaded `bgzf::io::Writer`: both stage into a buffer of
/// `MAX_BUF_SIZE`, both emit through `write_frame`, and libdeflate is
/// deterministic at a fixed level. The bytes are therefore identical to what
/// the serial writer produced, which
/// `multithreaded_bgzf_is_byte_identical_to_the_serial_writer` pins.
///
/// **The returned writer must not be written to from a rayon worker thread.**
/// It hands each finished block to the pool via `rayon::spawn` and blocks on a
/// channel bounded by the worker count while a dedicated thread drains it. A
/// caller that is itself occupying a worker can therefore fill the channel and
/// then wait for a compression task that has no free worker to run on: a real
/// deadlock, not a slowdown, and it is reachable at any pool size. rustar's
/// writes come from the pipeline's dispatcher thread (`run_batch_pipeline`
/// calls `consume` on its own thread, never on the pool), which is why this is
/// safe here; `write_batch` carries a `debug_assert` so a future call site that
/// breaks the rule fails loudly in tests instead of hanging in production.
fn make_bgzf_writer<W: std::io::Write + Send + 'static>(
    inner: W,
    compression: i32,
) -> bgzf::io::MultithreadedWriter<W> {
    bgzf::io::multithreaded_writer::Builder::default()
        .set_compression_level(bgzf_compression(compression))
        .build_from_writer(inner)
}

/// Below this many records in a batch, encoding on the calling thread beats
/// paying for the fan-out.
const PARALLEL_ENCODE_MIN_RECORDS: usize = 512;

/// Encode records to the raw BAM byte stream: per record, a little-endian
/// `block_size` followed by the encoded record.
///
/// This is what `bam::io::Writer::write_alignment_record` writes, and it is
/// produced here by that same call against an in-memory writer, so the bytes
/// cannot drift from the noodles encoder.
fn encode_records_into(
    header: &sam::Header,
    records: &[RecordBuf],
    out: Vec<u8>,
) -> Result<Vec<u8>, Error> {
    let mut writer = bam::io::Writer::from(out);
    for record in records {
        writer.write_alignment_record(header, record)?;
    }
    Ok(writer.into_inner())
}

/// How many records the sorted writers encode at a time on `finish()`.
///
/// The whole sorted set is already resident as `RecordBuf`s; encoding it in one
/// pass would hold the entire encoded stream alongside them. Encoding in slices
/// keeps the extra buffer bounded while still handing the pool enough records
/// per slice to be worth the fan-out.
const ENCODE_SLICE_RECORDS: usize = 65_536;

/// Guard the rule documented on `make_bgzf_writer`: writing to the
/// multithreaded BGZF writer from a rayon worker can deadlock. Cheap enough to
/// run on every batch, and only in debug builds.
fn assert_not_on_rayon_worker() {
    debug_assert!(
        rayon::current_thread_index().is_none(),
        "BGZF writes must come from outside the rayon pool: writing from a worker \
         can block the pool on its own compression tasks and deadlock"
    );
}

/// Encode `records` in bounded slices and write each to `writer`.
fn write_records_chunked<W: Write>(
    writer: &mut W,
    header: &sam::Header,
    records: &[RecordBuf],
) -> Result<(), Error> {
    for slice in records.chunks(ENCODE_SLICE_RECORDS) {
        let bytes = encode_batch(header, slice)?;
        writer.write_all(&bytes)?;
    }
    Ok(())
}

/// Encode a batch to raw BAM bytes, splitting the work across the rayon pool.
///
/// `bam::io::Writer` carries no state between records other than a scratch
/// buffer it clears at the top of every `write_alignment_record`, so encoding
/// a sub-range against a fresh in-memory writer yields exactly the bytes that
/// sub-range would have contributed to the serial stream. Concatenating the
/// sub-ranges in input order therefore reproduces the serial stream byte for
/// byte, which is why this is output-neutral rather than merely equivalent.
fn encode_batch(header: &sam::Header, batch: &[RecordBuf]) -> Result<Vec<u8>, Error> {
    let threads = rayon::current_num_threads();
    if threads < 2 || batch.len() < PARALLEL_ENCODE_MIN_RECORDS {
        return encode_records_into(header, batch, Vec::new());
    }

    use rayon::prelude::*;
    let chunk_len = batch.len().div_ceil(threads).max(64);
    // `par_chunks` is an indexed parallel iterator, so `collect` restores input
    // order regardless of which worker finished first.
    let chunks: Vec<Vec<u8>> = batch
        .par_chunks(chunk_len)
        .map(|chunk| encode_records_into(header, chunk, Vec::new()))
        .collect::<Result<_, _>>()?;

    let mut out = Vec::with_capacity(chunks.iter().map(Vec::len).sum());
    for chunk in &chunks {
        out.extend_from_slice(chunk);
    }
    Ok(out)
}

/// BAM file writer (streaming, unsorted)
///
/// This writer streams BAM records directly to disk as they're generated,
/// without buffering or sorting. The output is BGZF-compressed but unsorted.
/// Users can sort the output with `samtools sort` if needed.
pub struct BamWriter {
    writer: bgzf::io::MultithreadedWriter<BufWriter<File>>,
    header: sam::Header,
    finished: bool,
}

/// BAM file writer that collects all records in memory, sorts by coordinate,
/// then writes a single sorted BAM file on `finish()`.
///
/// The header emits `SO:coordinate`. Unmapped records sort to the end.
pub struct SortedBamWriter {
    records: Vec<RecordBuf>,
    output_path: std::path::PathBuf,
    header: sam::Header,
    compression: i32,
    limit_bam_sort_ram: u64,
}

impl BamWriter {
    /// Create a BAM writer at `output_path` with the given prepared header.
    ///
    /// Uses a lenient header writer that bypasses noodles' SAM-spec
    /// reference-name validator. STAR accepts `(` / `)` in reference names
    /// (yeast tRNA transcripts like `tP(UGG)A`), which noodles' strict
    /// validator rejects. Writing the SAM text portion of the BAM header
    /// manually sidesteps that validation while preserving every other byte.
    fn with_header(
        output_path: &Path,
        header: sam::Header,
        compression: i32,
    ) -> Result<Self, Error> {
        let buf_writer = BufWriter::new(File::create(output_path)?);
        let mut writer = make_bgzf_writer(buf_writer, compression);
        write_bam_header_lenient(&mut writer, &header, None)?;
        Ok(Self {
            writer,
            header,
            finished: false,
        })
    }

    /// Create a new BAM writer with header from genome index.
    pub fn create(output_path: &Path, genome: &Genome, params: &Parameters) -> Result<Self, Error> {
        Self::with_header(
            output_path,
            crate::io::sam::build_sam_header(genome, params)?,
            params.out_bam_compression,
        )
    }

    /// Create a BAM writer whose @SQ header lists the transcripts (not
    /// chromosomes) from `tr_idx`.  Used for
    /// `Aligned.toTranscriptome.out.bam`.
    pub fn create_transcriptome(
        output_path: &Path,
        tr_idx: &TranscriptomeIndex,
        params: &Parameters,
    ) -> Result<Self, Error> {
        let refs = tr_idx
            .tr_ids
            .iter()
            .zip(tr_idx.tr_length.iter())
            .map(|(id, len)| (id.as_str(), *len as usize));
        Self::with_header(
            output_path,
            crate::io::sam::build_sam_header_from_refs(refs, params)?,
            params.out_bam_compression,
        )
    }

    /// Write batch of buffered records (for parallel processing)
    ///
    /// # Arguments
    /// * `batch` - Slice of records to write
    pub fn write_batch(&mut self, batch: &[RecordBuf]) -> Result<(), Error> {
        assert_not_on_rayon_worker();
        let bytes = encode_batch(&self.header, batch)?;
        self.writer.write_all(&bytes)?;
        Ok(())
    }

    /// Flush and close BAM file
    pub fn finish(&mut self) -> Result<(), Error> {
        // `MultithreadedWriter::finish` shuts the workers down and panics if
        // called twice; `Drop` calls it too when it has not run yet.
        if !self.finished {
            self.writer.finish()?;
            self.finished = true;
        }
        log::info!("BAM file written successfully");
        Ok(())
    }
}

impl SortedBamWriter {
    /// Create a sorted BAM writer. Records are buffered in memory until `finish()`.
    pub fn create(
        output_path: &Path,
        genome: &crate::genome::Genome,
        params: &Parameters,
    ) -> Result<Self, Error> {
        let header = crate::io::sam::build_sam_header(genome, params)?;
        Ok(Self {
            records: Vec::new(),
            output_path: output_path.to_path_buf(),
            header,
            compression: params.out_bam_compression,
            limit_bam_sort_ram: params.limit_bam_sort_ram,
        })
    }

    /// Buffer records — no disk I/O yet.
    pub fn write_batch(&mut self, batch: &[RecordBuf]) -> Result<(), Error> {
        self.records.extend_from_slice(batch);
        Ok(())
    }

    /// Estimate memory used by buffered records (rough: 400 bytes/record for 150bp reads).
    fn estimated_ram(&self) -> u64 {
        self.records.len() as u64 * 400
    }

    fn check_ram_limit(&self) -> Result<(), Error> {
        if self.limit_bam_sort_ram > 0 {
            let est = self.estimated_ram();
            if est > self.limit_bam_sort_ram {
                return Err(Error::Alignment(format!(
                    "limitBAMsortRAM={} bytes exceeded: estimated {} bytes for {} records. \
                     Increase --limitBAMsortRAM or use --outSAMtype BAM Unsorted.",
                    self.limit_bam_sort_ram,
                    est,
                    self.records.len()
                )));
            }
        }
        Ok(())
    }

    /// Sort all buffered records by coordinate and write a single sorted BAM.
    ///
    /// Sort key: (reference_sequence_id, alignment_start).
    /// Unmapped records (no reference) sort to the end.
    pub fn finish(&mut self) -> Result<(), Error> {
        self.check_ram_limit()?;
        self.records
            .sort_by_key(|r| match (r.reference_sequence_id(), r.alignment_start()) {
                (Some(chr), Some(pos)) => (chr, pos.get()),
                _ => (usize::MAX, 0),
            });

        let buf_writer = BufWriter::new(File::create(&self.output_path)?);
        let mut writer = make_bgzf_writer(buf_writer, self.compression);
        write_bam_header_lenient(&mut writer, &self.header, Some("coordinate"))?;
        write_records_chunked(&mut writer, &self.header, &self.records)?;
        writer.finish()?;
        log::info!("Sorted BAM written ({} records)", self.records.len());
        Ok(())
    }

    /// Sort all buffered records and write to stdout (for `--outStd BAM_SortedByCoordinate`).
    pub fn finish_to_stdout(&mut self) -> Result<(), Error> {
        self.check_ram_limit()?;
        self.records
            .sort_by_key(|r| match (r.reference_sequence_id(), r.alignment_start()) {
                (Some(chr), Some(pos)) => (chr, pos.get()),
                _ => (usize::MAX, 0),
            });

        // `std::io::stdout()` rather than a lock guard: the BGZF writer owns
        // its sink on a worker thread, so the sink has to be `'static`.
        let buf_writer = BufWriter::new(std::io::stdout());
        let mut writer = make_bgzf_writer(buf_writer, self.compression);
        write_bam_header_lenient(&mut writer, &self.header, Some("coordinate"))?;
        write_records_chunked(&mut writer, &self.header, &self.records)?;
        writer.finish()?;
        log::info!(
            "Sorted BAM written to stdout ({} records)",
            self.records.len()
        );
        Ok(())
    }
}

/// Write a BAM header that tolerates reference sequence names STAR emits
/// (e.g. `tP(UGG)A`) but that noodles' SAM-spec validator rejects.
///
/// Replicates `noodles_bam::io::writer::header::write_header` byte-for-byte
/// for compliant headers; the only divergence is that the SAM text block
/// between `BAM\x01` and the binary reference list is produced via a local
/// formatter instead of `sam::io::Writer::write_header`, so forbidden-char
/// names pass through unchanged.
///
/// Binary reference list (after the text block) uses `CString::new` — the
/// only constraint there is "no interior nul", which is enforced upstream
/// via the usual UTF-8 input.
///
/// `sort_order`: if `Some("coordinate")`, injects `SO:coordinate` into the
/// @HD line. Pass `None` for unsorted output.
fn write_bam_header_lenient<W: Write>(
    writer: &mut W,
    header: &sam::Header,
    sort_order: Option<&str>,
) -> Result<(), Error> {
    const MAGIC: &[u8; 4] = b"BAM\x01";

    writer.write_all(MAGIC)?;

    // Build the SAM text block byte-for-byte identical to
    // `sam::io::Writer::write_header` minus the name validator:
    // `@HD`, `@SQ` (one per reference), `@RG` (if any), `@PG` (if any),
    // `@CO` (if any), each line terminated by `\n`.
    let text = render_sam_text_lenient(header, sort_order);
    let l_text = i32::try_from(text.len()).map_err(|_| {
        Error::Index(format!(
            "BAM SAM-text header exceeds i32::MAX bytes: {} bytes",
            text.len()
        ))
    })?;
    writer.write_i32::<LittleEndian>(l_text)?;
    writer.write_all(&text)?;

    // Binary reference list: n_ref then (l_name, name\0, l_ref) per ref.
    let refs = header.reference_sequences();
    let n_ref = i32::try_from(refs.len())
        .map_err(|_| Error::Index("BAM reference count exceeds i32::MAX".into()))?;
    writer.write_i32::<LittleEndian>(n_ref)?;
    for (name, rs) in refs {
        let c_name = CString::new(name.to_vec())
            .map_err(|e| Error::Index(format!("reference name contains interior NUL byte: {e}")))?;
        let name_bytes = c_name.as_bytes_with_nul();
        let l_name = u32::try_from(name_bytes.len()).map_err(|_| {
            Error::Index(format!(
                "reference name longer than u32::MAX: {} bytes",
                name_bytes.len()
            ))
        })?;
        writer.write_u32::<LittleEndian>(l_name)?;
        writer.write_all(name_bytes)?;
        let l_ref = i32::try_from(usize::from(rs.length()))
            .map_err(|_| Error::Index("reference length exceeds i32::MAX".into()))?;
        writer.write_i32::<LittleEndian>(l_ref)?;
    }

    Ok(())
}

fn render_sam_text_lenient(header: &sam::Header, sort_order: Option<&str>) -> Vec<u8> {
    let mut buf: Vec<u8> = Vec::new();

    // @HD line. noodles' Map<Header> in this version doesn't expose
    // sort_order/group_order via dedicated accessors; we serialize through
    // the generic `other_fields` map (which includes SO/GO when set).
    if let Some(hd) = header.header() {
        buf.extend_from_slice(b"@HD\tVN:");
        buf.extend_from_slice(hd.version().to_string().as_bytes());
        if let Some(so) = sort_order {
            buf.extend_from_slice(b"\tSO:");
            buf.extend_from_slice(so.as_bytes());
        }
        for (tag, value) in hd.other_fields() {
            buf.push(b'\t');
            buf.extend_from_slice(tag.as_ref());
            buf.push(b':');
            buf.extend_from_slice(value);
        }
        buf.push(b'\n');
    }

    // @SQ lines. Use the raw bytes so forbidden characters pass through.
    for (name, rs) in header.reference_sequences() {
        buf.extend_from_slice(b"@SQ\tSN:");
        buf.extend_from_slice(name);
        buf.extend_from_slice(b"\tLN:");
        buf.extend_from_slice(usize::from(rs.length()).to_string().as_bytes());
        // Other optional @SQ fields (AH, AN, AS, DS, M5, SP, TP, UR) —
        // rustar-aligner doesn't set any today, so skip.
        buf.push(b'\n');
    }

    // @RG lines.
    for (id, rg) in header.read_groups() {
        buf.extend_from_slice(b"@RG\tID:");
        buf.extend_from_slice(id);
        for (tag, value) in rg.other_fields() {
            buf.push(b'\t');
            buf.extend_from_slice(tag.as_ref());
            buf.push(b':');
            buf.extend_from_slice(value);
        }
        buf.push(b'\n');
    }

    // @PG lines — noodles' map doesn't guarantee insertion order; for
    // rustar-aligner we emit a single @PG with id "rustar-aligner". If more are added
    // later, pipe them in here.
    for (id, pg) in header.programs().as_ref() {
        buf.extend_from_slice(b"@PG\tID:");
        buf.extend_from_slice(id);
        for (tag, value) in pg.other_fields() {
            buf.push(b'\t');
            buf.extend_from_slice(tag.as_ref());
            buf.push(b':');
            buf.extend_from_slice(value);
        }
        buf.push(b'\n');
    }

    // @CO lines (comments).
    for comment in header.comments() {
        buf.extend_from_slice(b"@CO\t");
        buf.extend_from_slice(comment);
        buf.push(b'\n');
    }

    buf
}

/// Streaming unsorted BAM writer that writes to stdout.
pub struct BamStdoutWriter {
    writer: bgzf::io::MultithreadedWriter<BufWriter<std::io::Stdout>>,
    header: sam::Header,
    finished: bool,
}

impl BamStdoutWriter {
    pub fn create(genome: &crate::genome::Genome, params: &Parameters) -> Result<Self, Error> {
        let header = crate::io::sam::build_sam_header(genome, params)?;
        let mut writer = make_bgzf_writer(
            BufWriter::new(std::io::stdout()),
            params.out_bam_compression,
        );
        write_bam_header_lenient(&mut writer, &header, None)?;
        Ok(Self {
            writer,
            header,
            finished: false,
        })
    }

    pub fn write_batch(&mut self, batch: &[RecordBuf]) -> Result<(), Error> {
        assert_not_on_rayon_worker();
        let bytes = encode_batch(&self.header, batch)?;
        self.writer.write_all(&bytes)?;
        Ok(())
    }

    pub fn finish(&mut self) -> Result<(), Error> {
        if !self.finished {
            self.writer.finish()?;
            self.finished = true;
        }
        Ok(())
    }
}

/// Coordinate-sorted BAM writer that writes to stdout on `finish()`.
pub struct SortedBamStdoutWriter {
    records: Vec<RecordBuf>,
    header: sam::Header,
    compression: i32,
    limit_bam_sort_ram: u64,
}

impl SortedBamStdoutWriter {
    pub fn create(genome: &crate::genome::Genome, params: &Parameters) -> Result<Self, Error> {
        let header = crate::io::sam::build_sam_header(genome, params)?;
        Ok(Self {
            records: Vec::new(),
            header,
            compression: params.out_bam_compression,
            limit_bam_sort_ram: params.limit_bam_sort_ram,
        })
    }

    pub fn write_batch(&mut self, batch: &[RecordBuf]) -> Result<(), Error> {
        self.records.extend_from_slice(batch);
        Ok(())
    }

    pub fn finish(&mut self) -> Result<(), Error> {
        if self.limit_bam_sort_ram > 0 {
            let est = self.records.len() as u64 * 400;
            if est > self.limit_bam_sort_ram {
                return Err(Error::Alignment(format!(
                    "limitBAMsortRAM={} bytes exceeded: estimated {} bytes for {} records.",
                    self.limit_bam_sort_ram,
                    est,
                    self.records.len()
                )));
            }
        }
        self.records
            .sort_by_key(|r| match (r.reference_sequence_id(), r.alignment_start()) {
                (Some(chr), Some(pos)) => (chr, pos.get()),
                _ => (usize::MAX, 0),
            });
        let mut bgzf = make_bgzf_writer(BufWriter::new(std::io::stdout()), self.compression);
        write_bam_header_lenient(&mut bgzf, &self.header, Some("coordinate"))?;
        let mut bam_writer = bam::io::Writer::from(bgzf);
        for record in &self.records {
            bam_writer.write_alignment_record(&self.header, record)?;
        }
        bam_writer.finish(&self.header)?;
        log::info!(
            "Sorted BAM written to stdout ({} records)",
            self.records.len()
        );
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::transcript::{Exon, Transcript};
    use noodles::sam::alignment::record::cigar;
    use tempfile::NamedTempFile;

    fn create_test_genome() -> Genome {
        Genome {
            transform_blocks: None,
            sequence: vec![0, 1, 2, 3, 0, 1, 2, 3].into(), // ACGTACGT
            n_genome: 8,
            n_genome_real: 8,
            n_chr_real: 1,
            chr_name: vec!["chr1".to_string()],
            chr_length: vec![8],
            chr_start: vec![0, 8],
        }
    }

    fn default_params() -> Parameters {
        Parameters::parse_from(["rustar-aligner", "--readFilesIn", "test.fq"])
    }

    #[test]
    fn test_bam_writer_creation() {
        let genome = create_test_genome();
        let params = default_params();
        let temp_file = NamedTempFile::new().unwrap();

        let writer = BamWriter::create(temp_file.path(), &genome, &params);
        assert!(writer.is_ok(), "BAM writer creation should succeed");
    }

    #[test]
    fn test_bam_unmapped_write() {
        let genome = create_test_genome();
        let params = default_params();
        let temp_file = NamedTempFile::new().unwrap();

        let mut writer = BamWriter::create(temp_file.path(), &genome, &params).unwrap();

        let read_name = "read1";
        let read_seq = vec![0, 1, 2, 3]; // ACGT
        let read_qual = vec![30, 30, 30, 30];

        // Build record using SAM builder
        let record = crate::io::sam::SamWriter::build_unmapped_record(
            read_name,
            &read_seq,
            &read_qual,
            &params,
            crate::stats::UnmappedReason::Other,
        )
        .unwrap();

        let result = writer.write_batch(&[record]);
        assert!(result.is_ok(), "Writing unmapped read should succeed");

        // Finish the file
        let result = writer.finish();
        assert!(result.is_ok(), "Finishing BAM file should succeed");
    }

    #[test]
    fn test_bam_alignment_write() {
        use cigar::op::{Kind, Op};

        let genome = create_test_genome();
        let params = default_params();
        let temp_file = NamedTempFile::new().unwrap();

        let mut writer = BamWriter::create(temp_file.path(), &genome, &params).unwrap();

        // Create a simple transcript
        let transcript = Transcript {
            chr_idx: 0,
            genome_start: 100,
            genome_end: 104,
            is_reverse: false,
            exons: vec![Exon {
                genome_start: 100,
                genome_end: 104,
                read_start: 0,
                read_end: 4,
                i_frag: 0,
            }],
            cigar: vec![Op::new(Kind::Match, 4)],
            score: 0,
            n_mismatch: 0,
            n_gap: 0,
            n_junction: 0,
            junction_motifs: vec![],
            junction_annotated: vec![],
            read_seq: vec![],
        };

        let read_name = "read1";
        let read_seq = vec![0, 1, 2, 3]; // ACGT
        let read_qual = vec![30, 30, 30, 30];

        // Build records using SAM builder
        let records = crate::io::sam::SamWriter::build_alignment_records(
            read_name,
            &read_seq,
            &read_qual,
            0,
            0,
            &[transcript],
            &genome,
            &params,
            1,
        )
        .unwrap();

        let result = writer.write_batch(&records);
        assert!(result.is_ok(), "Writing alignment should succeed");

        let result = writer.finish();
        assert!(result.is_ok(), "Finishing BAM file should succeed");
    }

    #[test]
    fn test_bam_transcriptome_writer_creation() {
        use crate::junction::gtf::GtfRecord;
        use crate::quant::transcriptome::TranscriptomeIndex;
        use std::collections::HashMap;

        let genome = create_test_genome();
        // Stretch genome / exon to fit the tiny test chromosome.
        let mut attrs = HashMap::new();
        attrs.insert("gene_id".to_string(), "G1".to_string());
        attrs.insert("transcript_id".to_string(), "T1".to_string());
        let exons = vec![GtfRecord {
            seqname: "chr1".to_string(),
            feature: "exon".to_string(),
            start: 1,
            end: 8,
            strand: '+',
            attributes: attrs,
        }];
        let tr_idx = TranscriptomeIndex::from_gtf_exons(&exons, &genome).unwrap();
        assert_eq!(tr_idx.n_transcripts(), 1);

        let params = default_params();
        let temp_file = NamedTempFile::new().unwrap();
        let writer = BamWriter::create_transcriptome(temp_file.path(), &tr_idx, &params);
        assert!(
            writer.is_ok(),
            "transcriptome BAM writer creation should succeed"
        );

        // Header should contain exactly 1 @SQ entry (matching n_transcripts).
        let writer = writer.unwrap();
        assert_eq!(writer.header.reference_sequences().len(), 1);
    }

    #[test]
    fn test_bam_batch_write() {
        let genome = create_test_genome();
        let params = default_params();
        let temp_file = NamedTempFile::new().unwrap();

        let mut writer = BamWriter::create(temp_file.path(), &genome, &params).unwrap();

        // Create a batch of unmapped records
        let records = vec![
            crate::io::sam::SamWriter::build_unmapped_record(
                "read1",
                &[0, 1, 2, 3],
                &[30, 30, 30, 30],
                &params,
                crate::stats::UnmappedReason::Other,
            )
            .unwrap(),
            crate::io::sam::SamWriter::build_unmapped_record(
                "read2",
                &[0, 1, 2, 3],
                &[30, 30, 30, 30],
                &params,
                crate::stats::UnmappedReason::Other,
            )
            .unwrap(),
        ];

        let result = writer.write_batch(&records);
        assert!(result.is_ok(), "Writing batch should succeed");

        let result = writer.finish();
        assert!(result.is_ok(), "Finishing BAM file should succeed");
    }

    #[test]
    fn test_bam_compression_level_zero() {
        let genome = create_test_genome();
        let mut params = default_params();
        params.out_bam_compression = 0;
        let temp_file = NamedTempFile::new().unwrap();
        let writer = BamWriter::create(temp_file.path(), &genome, &params);
        assert!(
            writer.is_ok(),
            "BAM writer with compression=0 should succeed"
        );
    }

    /// The pre-change write path: one serial `bgzf::io::Writer`, records fed
    /// one at a time through `bam::io::Writer`. Every neutrality test below
    /// compares against these bytes, so what is being asserted is not "the two
    /// new halves agree with each other" but "the output did not move".
    fn serial_reference_bytes(header: &sam::Header, records: &[RecordBuf]) -> Vec<u8> {
        let mut bgzf = bgzf::io::writer::Builder::default()
            .set_compression_level(bgzf_compression(1))
            .build_from_writer(Vec::new());
        write_bam_header_lenient(&mut bgzf, header, None).unwrap();
        let mut writer = bam::io::Writer::from(bgzf);
        for record in records {
            writer.write_alignment_record(header, record).unwrap();
        }
        writer.try_finish().unwrap();
        writer.into_inner().into_inner()
    }

    fn unmapped_records(params: &Parameters, n: usize) -> Vec<RecordBuf> {
        (0..n)
            .map(|i| {
                crate::io::sam::SamWriter::build_unmapped_record(
                    &format!("read{i}"),
                    &[0, 1, 2, 3, 3, 2, 1, 0],
                    &[30; 8],
                    params,
                    crate::stats::UnmappedReason::Other,
                )
                .unwrap()
            })
            .collect()
    }

    /// The multithreaded BGZF writer stages into the same `MAX_BUF_SIZE`
    /// buffer, emits through the same `write_frame`, and appends the same EOF
    /// block as the serial writer, so switching to it must not move a byte.
    /// Enough records to span several BGZF blocks, so block boundaries are
    /// actually exercised rather than a single short block.
    #[test]
    fn multithreaded_bgzf_is_byte_identical_to_the_serial_writer() {
        let params = default_params();
        let header = sam::Header::default();
        let records = unmapped_records(&params, 4000);
        let expected = serial_reference_bytes(&header, &records);

        let mut writer = make_bgzf_writer(Vec::new(), 1);
        write_bam_header_lenient(&mut writer, &header, None).unwrap();
        for record in &records {
            let bytes =
                encode_records_into(&header, std::slice::from_ref(record), Vec::new()).unwrap();
            writer.write_all(&bytes).unwrap();
        }
        let got = writer.finish().unwrap();

        assert_eq!(got, expected);
    }

    /// Splitting the batch across the pool is only sound if concatenating the
    /// sub-ranges in input order reproduces the serial stream. Checked at
    /// several pool sizes, including one worker (which takes the serial branch)
    /// and a size that does not divide the batch evenly.
    #[test]
    fn parallel_encoding_is_byte_identical_at_every_worker_count() {
        let params = default_params();
        let header = sam::Header::default();
        let records = unmapped_records(&params, 3000);
        let expected = encode_records_into(&header, &records, Vec::new()).unwrap();

        for threads in [1, 3, 4, 7, 16] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap();
            let got = pool.install(|| encode_batch(&header, &records)).unwrap();
            assert_eq!(got, expected, "encoding diverged at {threads} workers");
        }
    }

    /// Batches under the fan-out threshold take the serial branch; that branch
    /// has to produce the same bytes as the parallel one, not merely a valid
    /// stream.
    #[test]
    fn a_batch_below_the_fanout_threshold_encodes_the_same_bytes() {
        let params = default_params();
        let header = sam::Header::default();
        let records = unmapped_records(&params, PARALLEL_ENCODE_MIN_RECORDS - 1);
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(8)
            .build()
            .unwrap();
        let got = pool.install(|| encode_batch(&header, &records)).unwrap();
        assert_eq!(
            got,
            encode_records_into(&header, &records, Vec::new()).unwrap()
        );
    }

    /// The sorted writers encode in bounded slices rather than in one pass, so
    /// the slice boundary must not show up in the output.
    #[test]
    fn chunked_encoding_matches_a_single_pass() {
        let params = default_params();
        let header = sam::Header::default();
        let records = unmapped_records(&params, 2500);
        let expected = encode_records_into(&header, &records, Vec::new()).unwrap();

        let mut got = Vec::new();
        for slice in records.chunks(700) {
            got.extend_from_slice(&encode_batch(&header, slice).unwrap());
        }
        assert_eq!(got, expected);
    }

    /// End to end through `BamWriter`: the file on disk is byte-identical to
    /// what the serial path wrote.
    ///
    /// Deliberately not wrapped in a `ThreadPool::install`. Doing so puts the
    /// caller on a worker and deadlocks, for the reason spelled out on
    /// `make_bgzf_writer` — this test found that the hard way. Worker-count
    /// invariance of the part that has one is covered by
    /// `parallel_encoding_is_byte_identical_at_every_worker_count`; BGZF block
    /// boundaries do not depend on the worker count at all, since the staging
    /// buffer fills to a fixed size before any block is handed to the pool.
    #[test]
    fn bam_file_is_byte_identical_to_the_serial_path() {
        let genome = create_test_genome();
        let params = default_params();
        let header = crate::io::sam::build_sam_header(&genome, &params).unwrap();
        let records = unmapped_records(&params, 2000);
        let expected = serial_reference_bytes(&header, &records);

        let temp_file = NamedTempFile::new().unwrap();
        let mut writer = BamWriter::create(temp_file.path(), &genome, &params).unwrap();
        writer.write_batch(&records).unwrap();
        writer.finish().unwrap();

        assert_eq!(std::fs::read(temp_file.path()).unwrap(), expected);
    }

    /// `finish()` shuts the BGZF workers down, and calling it twice would panic
    /// inside noodles. `Drop` also calls it. Both paths have to stay safe.
    #[test]
    fn finishing_twice_is_not_an_error() {
        let genome = create_test_genome();
        let params = default_params();
        let temp_file = NamedTempFile::new().unwrap();
        let mut writer = BamWriter::create(temp_file.path(), &genome, &params).unwrap();
        writer.finish().unwrap();
        writer.finish().unwrap();
    }

    #[test]
    fn test_sorted_bam_limit_ram_unlimited() {
        let genome = create_test_genome();
        let mut params = default_params();
        params.limit_bam_sort_ram = 0; // unlimited
        let temp_file = NamedTempFile::new().unwrap();
        let mut writer = SortedBamWriter::create(temp_file.path(), &genome, &params).unwrap();
        let result = writer.finish();
        assert!(
            result.is_ok(),
            "Sorted BAM with unlimited RAM should succeed"
        );
    }

    #[test]
    fn test_sorted_bam_limit_ram_exceeded() {
        let genome = create_test_genome();
        let mut params = default_params();
        params.limit_bam_sort_ram = 1; // 1 byte — will be exceeded by any records
        let temp_file = NamedTempFile::new().unwrap();
        let mut writer = SortedBamWriter::create(temp_file.path(), &genome, &params).unwrap();
        // Add a record to trigger the limit
        let rec = crate::io::sam::SamWriter::build_unmapped_record(
            "r1",
            &[0, 1, 2, 3],
            &[30; 4],
            &params,
            crate::stats::UnmappedReason::Other,
        )
        .unwrap();
        writer.write_batch(&[rec]).unwrap();
        let result = writer.finish();
        assert!(result.is_err(), "Should fail when RAM limit is exceeded");
    }
}
