#!/usr/bin/env bash
# Interleaved A/B of BAM writing, before and after the parallel writer.
#
#   test/bench_bam_write.sh <genomeDir> <reads.fastq> [threads] [pairs] [baseRef]
#
# Two binaries are built: one from `baseRef` (default origin/main) and one from
# the working tree. They are then alternated pair by pair, with the order
# flipped on even pairs so a machine that drifts (thermal, page cache warming)
# cannot favour whichever runs first.
#
# Each pair times three output modes:
#
#   None                  alignment only, no BAM written
#   BAM Unsorted          streaming write
#   BAM SortedByCoordinate  sort then write
#
# `None` is not decoration. Only the difference between a BAM mode and `None`
# is the work this change touches; the total is mostly alignment and hides it.
# A previous measurement on this repo reported a BAM-writer difference that the
# workload could not have resolved, because BAM writing was 1-4% of the run.
# Report the delta, and if the delta is smaller than the run-to-run spread, say
# that the workload cannot settle the question rather than publishing a median.
set -euo pipefail

GENOME_DIR=${1:?usage: bench_bam_write.sh <genomeDir> <reads.fastq> [threads] [pairs] [baseRef]}
READS=${2:?usage: bench_bam_write.sh <genomeDir> <reads.fastq> [threads] [pairs] [baseRef]}
THREADS=${3:-16}
PAIRS=${4:-6}
BASE_REF=${5:-origin/main}

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK=$(mktemp -d)
BASE_TREE=$(mktemp -d)
trap 'rm -rf "$WORK"; git -C "$ROOT" worktree remove --force "$BASE_TREE" 2>/dev/null || rm -rf "$BASE_TREE"' EXIT

# Refuse to measure on a busy machine. A single unrelated job saturating the
# cores does not just add noise, it can invert the result, and it is invisible
# in the numbers afterwards: the series simply drifts.
#
# The check is on CPU idle, not on load average. Load average is an exponential
# average over minutes, so it stays high long after the offending job has gone
# and would refuse to measure on a machine that is now perfectly quiet.
cpu_idle() { # percent idle, sampled over one second
  top -l 2 -n 0 2>/dev/null | awk '/CPU usage/ {gsub("%","",$(NF-1)); v=$(NF-1)} END {print v+0}'
}
IDLE=$(cpu_idle)
if awk -v i="$IDLE" 'BEGIN{exit !(i < 80)}'; then
  echo "CPU is only ${IDLE}% idle; other work is running and the numbers would not mean anything." >&2
  echo "wait for the machine to be idle, or set BENCH_IGNORE_LOAD=1 to override." >&2
  [ "${BENCH_IGNORE_LOAD:-0}" = "1" ] || exit 1
fi

echo "building new (working tree)" >&2
cargo build --release --manifest-path "$ROOT/Cargo.toml" >&2
cp "$ROOT/target/release/rustar-aligner" "$WORK/new"

echo "building old ($BASE_REF)" >&2
rm -rf "$BASE_TREE"
git -C "$ROOT" worktree add --detach "$BASE_TREE" "$BASE_REF" >&2
cargo build --release --manifest-path "$BASE_TREE/Cargo.toml" >&2
cp "$BASE_TREE/target/release/rustar-aligner" "$WORK/old"

run() { # $1 = binary, $2 = outSAMtype words
  local dir="$WORK/run"
  rm -rf "$dir"
  mkdir -p "$dir"
  local t0 t1
  t0=$(python3 -c 'import time; print(time.time())')
  # shellcheck disable=SC2086 # $2 is deliberately two words for "BAM Unsorted"
  "$1" --genomeDir "$GENOME_DIR" --readFilesIn "$READS" --runThreadN "$THREADS" \
       --outSAMtype $2 --outFileNamePrefix "$dir/" >/dev/null 2>&1
  t1=$(python3 -c "import time; print(f'{time.time()-$t0:.2f}')")
  echo "$t1"
}

# One discarded run so the reads and the index are in page cache for pair 1.
run "$WORK/new" "None" >/dev/null

for mode in "None" "BAM Unsorted" "BAM SortedByCoordinate"; do
  echo "=== --outSAMtype $mode ==="
  for i in $(seq 1 "$PAIRS"); do
    if (( i % 2 )); then
      o=$(run "$WORK/old" "$mode")
      n=$(run "$WORK/new" "$mode")
    else
      n=$(run "$WORK/new" "$mode")
      o=$(run "$WORK/old" "$mode")
    fi
    echo "pair$i old=${o}s new=${n}s idle=$(cpu_idle)%"
  done
done
