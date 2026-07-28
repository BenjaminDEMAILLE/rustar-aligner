// Phase 6+: FASTQ reader, SAM/BAM output, SJ.out.tab

pub mod bam;
#[cfg(all(feature = "htslib-bam", not(windows)))]
pub mod bam_htslib;
pub mod fastq;
pub mod log;
pub mod sam;
