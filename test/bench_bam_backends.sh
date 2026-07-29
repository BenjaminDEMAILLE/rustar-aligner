#!/usr/bin/env bash
# Interleaved A/B of the two unsorted-BAM backends.
#
# The default build uses the noodles writer; --features htslib-bam swaps in the
# htslib one. Both binaries are built first, then alternated pair by pair, with
# the order flipped on even pairs so a machine that drifts (thermal, or page
# cache warming) cannot favour whichever runs first.
#
#   test/bench_bam_backends.sh <genomeDir> <reads.fastq> [threads] [pairs]
#
# Report the median of the pairs, and discard any pair whose two runs straddle
# a visible drift in the series rather than averaging over it.
set -euo pipefail

GENOME_DIR=${1:?usage: bench_bam_backends.sh <genomeDir> <reads.fastq> [threads] [pairs]}
READS=${2:?usage: bench_bam_backends.sh <genomeDir> <reads.fastq> [threads] [pairs]}
THREADS=${3:-16}
PAIRS=${4:-6}

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

# Refuse to measure on a busy machine. A single unrelated job saturating the
# cores does not just add noise, it can invert the result, and it is invisible
# in the numbers afterwards: the series simply drifts. Checking beforehand is
# the only cheap way to know the measurement means anything.
load_now() { uptime | sed -E 's/.*load averages?: *([0-9.]+).*/\1/'; }
LOAD=$(load_now)
if awk -v l="$LOAD" 'BEGIN{exit !(l > 2.0)}'; then
  echo "load average is $LOAD; other work is running and the numbers would not mean anything." >&2
  echo "wait for the machine to be idle, or set BENCH_IGNORE_LOAD=1 to override." >&2
  [ "${BENCH_IGNORE_LOAD:-0}" = "1" ] || exit 1
fi

echo "building both backends" >&2
cargo build --release --manifest-path "$ROOT/Cargo.toml" >&2
cp "$ROOT/target/release/rustar-aligner" "$WORK/noodles"
cargo build --release --features htslib-bam --manifest-path "$ROOT/Cargo.toml" >&2
cp "$ROOT/target/release/rustar-aligner" "$WORK/htslib"

run() { # $1 = binary, $2 = output dir
  rm -rf "$2"
  mkdir -p "$2"
  local t0 t1
  t0=$(python3 -c 'import time; print(time.time())')
  "$1" --genomeDir "$GENOME_DIR" --readFilesIn "$READS" --runThreadN "$THREADS" \
       --outSAMtype BAM Unsorted --outFileNamePrefix "$2/" >/dev/null 2>&1
  t1=$(python3 -c "import time; print(f'{time.time()-$t0:.2f}')")
  echo "$t1"
}

# One discarded run so the reads and the index are in page cache for pair 1.
run "$WORK/noodles" "$WORK/warm" >/dev/null

for i in $(seq 1 "$PAIRS"); do
  if (( i % 2 )); then
    a=$(run "$WORK/noodles" "$WORK/n")
    b=$(run "$WORK/htslib" "$WORK/h")
  else
    b=$(run "$WORK/htslib" "$WORK/h")
    a=$(run "$WORK/noodles" "$WORK/n")
  fi
  echo "pair$i noodles=${a}s htslib=${b}s load=$(load_now)"
done
