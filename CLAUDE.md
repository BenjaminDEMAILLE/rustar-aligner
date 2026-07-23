# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Important: Git Workflow

**DO NOT commit changes automatically.** The user will review, test, and commit changes themselves. Claude should:
- Make code changes as requested
- Suggest what should be committed
- Let the user handle `git add`, `git commit`, and `git push`

## Project Overview

rustar-aligner is a Rust reimplementation of [STAR](https://github.com/alexdobin/STAR) (Spliced Transcripts Alignment to a Reference), an RNA-seq aligner originally written in C++ by Alexander Dobin. Licensed under MIT to match the original STAR license.

The primary goal is a faithful port — matching the original STAR behavior as closely as possible. Extra features and divergences from the original will come in later releases/forks. When implementing, refer to the [STAR source code](https://github.com/alexdobin/STAR) to ensure correctness and behavioral parity.

## Build Commands

Rust 2024 edition. Standard Cargo commands:

```bash
cargo build                 # Debug build
cargo build --release       # Release build
cargo test                  # Run all tests
cargo test <name>           # Run a single test by name
cargo clippy --all-targets  # Lint
cargo fmt                   # Format code
```

Always run `cargo clippy --all-targets`, `cargo fmt --check`, and `cargo test` before considering a phase complete.

## Current Status

End-to-end SE + PE RNA-seq alignment with splice-junction detection, two-pass mode, 4-tier chimeric detection, GeneCounts / TranscriptomeSAM quantification, SortedByCoordinate BAM, and multi-threaded processing. 0 clippy warnings.

**Faithfulness vs STAR** (10k yeast reads, ERR12389696, 150 bp): SE **99.815% tie-adjusted** (8611/8627 non-tie reads exact; 0 STAR-only, 0 rustar-aligner-only, 1 CIGAR-only tie). PE **8390 both-mapped** (matches STAR), 0 half-mapped, **99.883% tie-adjusted** exact; 0 MAPQ inflations/deflations, 0 NH diffs, 0 proper-pair diffs. The SA is byte-for-byte identical to STAR's for the yeast genome. Remaining raw disagreements are verified genuine ties (SA-order / seeded-RNG tie-break — rustar-aligner uses `StdRng`, not STAR's mt19937).

Phase-by-phase history (what changed, when, and why) lives in **[ROADMAP.md](ROADMAP.md)** and **[CHANGELOG.md](CHANGELOG.md)**; per-phase dev notes in **[docs-old/](docs-old/)**; the published Astro Starlight docs site is in **[docs/](docs/)**. Update those, not this section, when adding phases.

## Source Layout

```
src/
  main.rs          -- Thin entry: parse CLI (clap), init logging, call lib::run()
  lib.rs           -- run() dispatches on RunMode (AlignReads | GenomeGenerate)
  cpu.rs           -- Runtime CPU/SIMD feature-detection guard; VERSION_BODY from build.rs
  params/
    mod.rs         -- STAR CLI params via clap derive, --camelCase long names, Parameters::validate()
    sam.rs         -- SAM-attribute / outSAMtype parsing helpers
  error.rs         -- Error enum with thiserror (Parameter, Io, Fasta, Index, Alignment, Gtf)
  mapq.rs          -- MAPQ calculation (STAR lookup table: n=1→255, n=2→3, n≥5→0)
  stats.rs         -- Alignment statistics, Log.final.out writer, UnmappedReason enum
  genome/
    mod.rs         -- Genome struct, padding logic, reverse complement, file writing
    fasta.rs       -- FASTA parser, base encoding (A=0,C=1,G=2,T=3,N=4)
  index/
    mod.rs         -- GenomeIndex (build + load + write)
    packed_array.rs-- Variable-width bit packing (1-64 bits per element)
    packed_stream.rs- Streaming writer for the PackedArray bit format (SA build)
    suffix_array.rs-- SA data structure + strand encoding + naive oracle for tests
    sa_build.rs    -- STAR-faithful SA construction via the caps-sa crate (sentinel transform, in-mem vs ext-mem)
    sa_index.rs    -- K-mer lookup table (35-bit entries with flags)
    io.rs          -- Load index from disk (Genome, SA, SAindex)
  align/
    mod.rs         -- Module definition
    seed.rs        -- Seed finding via hierarchical SAindex lookup + MMP search
    stitch.rs      -- Seed clustering + DP stitching + alignment extension (extendAlign)
    score.rs       -- Scoring functions (gaps, mismatches, splice junctions)
    transcript.rs  -- Transcript struct (exon coords, CIGAR, scores)
    read_align.rs  -- Per-read alignment driver
  io/
    mod.rs         -- Module exports
    fastq.rs       -- FASTQ reader (plain + gzip, noodles wrapper)
    sam.rs         -- SAM writer (header + records, noodles wrapper), BySJout disk buffering
    bam.rs         -- BAM writer (BGZF compression, streaming unsorted + in-memory sorted output)
    log.rs         -- STAR-compatible Log.out / Log.progress.out writers
  junction/
    mod.rs         -- GTF parsing, junction database, motif detection, two-pass filtering
    sj_output.rs   -- SJ.out.tab writer
    sjdb_insert.rs -- SJDB insertion into Genome+SA at genomeGenerate (ports sjdbPrepare/sjdbBuildIndex.cpp)
    gtf.rs         -- GTF parser (internal)
  quant/
    mod.rs         -- Gene-level read counting (--quantMode GeneCounts, ReadsPerGene.out.tab)
    transcriptome.rs- Transcriptome-coordinate SAM output (--quantMode TranscriptomeSAM)
  chimeric/
    mod.rs         -- Module exports
    detect.rs      -- Chimeric detection (Tier 1: transcript-pair, Tier 1b: soft-clip re-seed, Tier 2: multi-cluster, Tier 3: residual re-seed)
    segment.rs     -- ChimericSegment and ChimericAlignment data structures
    score.rs       -- Junction type classification, repeat length calculation
    output.rs      -- Chimeric.out.junction writer (14-column format)
```

## Development Philosophy — Match STAR Exactly

**Never backtrack on changes that make the method closer to STAR's logic, even if metrics temporarily regress.** A regression in alignment statistics (position agreement, splice rate, etc.) after making the code more STAR-faithful does NOT mean the change was wrong — it means the implementation is incomplete. More STAR-matching changes are needed to close the gap.

- It is better to go backwards in metrics temporarily to eventually reach 100% matching results, than to keep approximate methods that can never get there.
- Every algorithm, threshold, and code path should match STAR's C++ source exactly. Any deviation — even a small one — should be corrected.
- Do not revert STAR-faithful changes. Instead, investigate what STAR does differently in the surrounding logic and implement that too.
- When a STAR-matching change causes regression, the response is: "What else does STAR do that we don't?" — not "revert and try something different."

## Key Conventions

- **Crate/binary name is `rustar-aligner`** — lib name in Rust is `rustar_aligner`; `#![allow(non_snake_case)]` in lib.rs is retained for STAR-style camelCase identifiers
- **STAR params use `--camelCase` naming** — clap `#[arg(long = "camelCase")]` maps to snake_case Rust fields
- **Multi-value params** (genomeFastaFiles, readFilesIn, outSAMtype, outSAMattributes, chimOutType, alignSJstitchMismatchNmax, outSJfilterIntronMaxVsReadN) need explicit `num_args`
- **Negative defaults** (scoreGapNoncan=-8, readMapNumber=-1, etc.) need `allow_hyphen_values = true`
- **`outSAMtype`** is parsed as raw `Vec<String>` then structured via `Parameters::out_sam_type()` method
- **Validation** beyond clap's type checking is in `Parameters::validate()` (e.g. genomeGenerate requires FASTA files)
- **No async** — CPU-bound work; async adds complexity with zero benefit
- **Error handling** — `thiserror` for `Error` enum, `anyhow` for top-level result propagation

## Dependencies

See `Cargo.toml` for the authoritative list; the load-bearing choices:

```toml
[dependencies]
clap = { version = "4", features = ["derive"] }   # CLI parsing (derive)
anyhow = "1"                                       # top-level error propagation
thiserror = "2"                                    # Error enum
noodles = { version = "0.111", features = ["fastq", "sam", "bam", "bgzf"] }  # SAM/BAM/FASTQ IO
memmap2 = "0.9"                                    # mmap the genome/SA index
caps-sa = "0.6"                                    # suffix-array construction (see index/sa_build.rs)
rayon = "1"                                        # per-read parallel alignment
dashmap = "6"                                      # concurrent junction/gene maps
rand = "0.10"                                      # --runRNGseed tie-breaking (StdRng, NOT mt19937)
tempfile = "3"                                     # BySJout disk buffering (now a runtime dep)
bitflags = "2"                                     # SAM flag bitsets
shlex = "2"                                        # split parametersFiles / command strings
flate2 = "1"                                       # gzip FASTQ + BGZF levels
mimalloc = { version = "0.1", default-features = false }  # global allocator (bounds peak RSS, faster small allocs)
# also: log, env_logger, byteorder, bstr, chrono

[build-dependencies]                               # build.rs stamps OS/arch/SIMD/git-hash into VERSION_BODY
chrono = { version = "0.4", default-features = false, features = ["clock"] }

[dev-dependencies]
assert_cmd = "2"
predicates = "3"
```

**Release profile is deliberately slow to build**: `lto = "fat"`, `codegen-units = 1`, `strip = true` in `[profile.release]` (tuned for the seed-search / DP inner loops). For iterating, prefer `cargo build` (debug) or `cargo test`; CI overrides via `CARGO_PROFILE_RELEASE_LTO=thin` / `CODEGEN_UNITS=16`.

## Testing Pattern

- Unit tests: `#[cfg(test)]` in each module (the bulk of the suite)
- Integration tests: `tests/` directory (`alignment_features.rs`, `phase9_threading.rs`, `transcriptome_sam.rs`) — synthetic genomes, run the built binary via `assert_cmd`
- Differential testing against STAR: Python/shell harness in `test/` (NOT `tests/`) — `compare_sam.py`, `compare_pe.py`, `assess_faithfulness.py`, `compare_junctions.py`, `compare_chimeric.py`, etc. `test/ci.sh` drives them.
- Test data tiers: synthetic micro-genome → 10k yeast reads (ERR12389696) → full human genome
- **Don't hardcode a test count in prose** — it goes stale fast. Get the current number from `cargo test`. (`git grep -c '#\[test\]'` currently reports ~455 test fns.)

## Known Issues & Limitations

Full issue tracking is in **[ROADMAP.md](ROADMAP.md)**; feature status in **[docs-old/phase17_features.md](docs-old/phase17_features.md)**. The durable takeaways a future session needs:

- **Remaining SE/PE disagreements vs STAR are almost all genuine ties.** Both tools find the identical alignment set; the primary differs only because of SA-iteration order or RNG divergence (rustar-aligner tie-breaks with seeded `StdRng` via `--runRNGseed`, not STAR's mt19937). When a diff appears, **first check whether it is a tie** before treating it as a bug — the alignment sets matching is the signal.
- **One actionable SE CIGAR-only tie**: `ERR12389696.13573895` — same position (XV:218357), same AS=133, but insertion placed at read pos 100 vs STAR's 108, because the 71-base seed anchors at a different offset through a homopolymer run. Fixing it requires reproducing STAR's exact Lmapped chain.
- **4 PE AS diffs are rustar-aligner improvements, not bugs** (`.844151`, `.4972950`, …): rustar-aligner stitches a better-scoring pair than STAR's combined-window approach finds. Do not "fix" these toward STAR.
- **Not implemented**: STARsolo single-cell features (deferred).
