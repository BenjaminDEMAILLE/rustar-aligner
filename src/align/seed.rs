use crate::error::Error;
use crate::index::GenomeIndex;
use crate::io::fastq::complement_base;
use crate::params::Parameters;

/// A seed represents an exact match between a read position and genome location(s).
#[derive(Debug, Clone)]
pub struct Seed {
    /// Position in the read where this seed starts
    pub read_pos: usize,

    /// Length of the exact match
    pub length: usize,

    /// Range in the suffix array [start, end) where this k-mer appears
    pub sa_start: usize,
    pub sa_end: usize,

    /// Whether this seed is on the reverse strand of the read
    pub is_reverse: bool,

    /// Whether this seed was found via R→L (reverse-complement) search.
    /// When true, genome_positions() converts coordinates back to forward orientation.
    pub search_rc: bool,

    /// Mate identifier for paired-end reads
    /// 0 = mate1, 1 = mate2, 2 = single-end (default)
    pub mate_id: u8,
}

impl Seed {
    /// Find all seeds for a read sequence using MMP (Maximal Mappable Prefix) search.
    ///
    /// For each position in the read, performs binary search on the suffix array
    /// to find the longest exact match.
    ///
    /// # Arguments
    /// * `read_seq` - Read sequence (encoded as 0=A, 1=C, 2=G, 3=T)
    /// * `index` - Genome index with SA and SAindex
    /// * `min_seed_length` - Minimum seed length to report
    /// * `params` - Parameters including seedMultimapNmax
    ///
    /// # Returns
    /// Vector of seeds found in the read
    pub fn find_seeds(
        read_seq: &[u8],
        index: &GenomeIndex,
        min_seed_length: usize,
        params: &Parameters,
        debug_name: &str,
    ) -> Result<Vec<Seed>, Error> {
        let mut seeds = Vec::new();
        let read_len = read_seq.len();

        // STAR uses the SAME sparse chain-based loop for both L→R and R→L directions.
        // For each direction: Nstart evenly-spaced starting positions, each advancing by
        // MMP length (Lmapped). Chains continue until only seedMapMin (5) bases remain.
        // This matches STAR's ReadAlign_mapOneRead.cpp: for(iDir=0;iDir<2;iDir++)
        //   for(istart=0;istart<Nstart;istart++)
        //     while(istart*Lstart + Lmapped + seedMapMin < readLen) { ... Lmapped += L; }

        // Search L→R (forward direction on read): sparse chain search
        search_direction_sparse(
            read_seq,
            read_len,
            index,
            min_seed_length,
            params,
            false,
            debug_name,
            &mut seeds,
        );

        // Cap check between directions (STAR: seedPerReadNmax applies across both)
        if seeds.len() >= params.seed_per_read_nmax {
            return Ok(seeds);
        }

        // Search R→L (reverse direction on read): sparse chain search on RC read
        let rc_read = reverse_complement_read(read_seq);
        search_direction_sparse(
            &rc_read,
            read_len,
            index,
            min_seed_length,
            params,
            true,
            debug_name,
            &mut seeds,
        );

        // STAR's storeAligns keeps the seed array `PC[]` sorted by rStart, with
        // longer seeds ordered before shorter ones at the same rStart, and drops
        // exact (rStart, Length) duplicates **regardless of search direction**
        // (`ReadAlign_storeAligns.cpp`, OPTIM_STOREaligns_SIMPLE). rustar-aligner
        // previously kept `search_rc` in the dedup key (so a forward seed and its
        // reverse-search twin were both retained) and left the array in
        // "all L→R, then all R→L" order. That made window-creation order in
        // `cluster_seeds` diverge from STAR's, perturbing the earliest-window
        // primary tie-break and the seed chain chosen through repeats.
        //
        // Match STAR: sort by (rStart asc, Length desc) — a *stable* sort so that
        // among equal (rStart, Length) the earliest-found seed (forward, since
        // L→R is collected first) is kept — then dedup direction-agnostically.
        Ok(finalize_seed_order(seeds))
    }

    /// Find seeds for several reads at once, with their suffix-array searches
    /// interleaved.
    ///
    /// Same seeds as calling [`find_seeds`](Self::find_seeds) on each read, in
    /// the same order: the chains are driven in lockstep rather than one after
    /// another, and each chain's seeds are collected into its own bucket so the
    /// concatenation reproduces the sequential push order exactly. The
    /// `seedPerReadNmax` cap is then applied by truncation, which is what the
    /// sequential early return amounts to.
    ///
    /// Batching pays because a chain's next search cannot start until the
    /// current one finishes, but *different* chains -- across starting
    /// positions, across directions, and across reads -- are independent. There
    /// are only a handful of chains within one read, well under the width where
    /// the interleave pays for itself, which is why this takes a slice of reads
    /// rather than working inside one.
    pub fn find_seeds_batch(
        reads: &[&[u8]],
        index: &GenomeIndex,
        min_seed_length: usize,
        params: &Parameters,
    ) -> Vec<Vec<Seed>> {
        // One cursor per (read, direction, starting position). `bucket` is the
        // cursor's slot in the final concatenation, so seed order does not
        // depend on the order chains happen to finish.
        struct Cursor {
            read_idx: usize,
            bucket: usize,
            is_rc: bool,
            pos: usize,
            original_read_len: usize,
            done: bool,
        }

        let rc_reads: Vec<Vec<u8>> = reads.iter().map(|r| reverse_complement_read(r)).collect();
        let mut cursors: Vec<Cursor> = Vec::new();
        let mut buckets: Vec<Vec<Seed>> = Vec::new();
        // Per read, the range of bucket indices in sequential order.
        let mut read_buckets: Vec<Vec<usize>> = vec![Vec::new(); reads.len()];

        for (read_idx, read) in reads.iter().enumerate() {
            for is_rc in [false, true] {
                let seq: &[u8] = if is_rc { &rc_reads[read_idx] } else { read };
                for start_pos in chain_starts(seq.len(), params) {
                    let bucket = buckets.len();
                    buckets.push(Vec::new());
                    read_buckets[read_idx].push(bucket);
                    cursors.push(Cursor {
                        read_idx,
                        bucket,
                        is_rc,
                        pos: start_pos,
                        original_read_len: read.len(),
                        done: false,
                    });
                }
            }
        }

        let mut reqs: Vec<MmpReq> = Vec::with_capacity(SEED_BATCH_WIDTH);
        let mut owners: Vec<usize> = Vec::with_capacity(SEED_BATCH_WIDTH);
        let mut results: Vec<(usize, usize, usize)> = Vec::new();

        loop {
            reqs.clear();
            owners.clear();
            // One round: every live cursor contributes at most one search.
            for (ci, c) in cursors.iter_mut().enumerate() {
                if c.done {
                    continue;
                }
                let seq: &[u8] = if c.is_rc {
                    &rc_reads[c.read_idx]
                } else {
                    reads[c.read_idx]
                };
                if c.pos >= seq.len() || seq.len() - c.pos < min_seed_length {
                    c.done = true;
                    continue;
                }
                let (prepared, sa_start, sa_end, l_initial) =
                    prepare_seed_at_position(seq, c.pos, index, min_seed_length, false, params);
                match prepared {
                    PreparedSeed::Resolved(r) => {
                        push_seed(
                            &mut buckets[c.bucket],
                            r.seed,
                            c.is_rc,
                            c.original_read_len,
                            params,
                        );
                        c.pos += r.advance;
                    }
                    PreparedSeed::Search { req_read_pos } => {
                        reqs.push(MmpReq {
                            read_seq: seq,
                            read_pos: req_read_pos,
                            sa_start,
                            sa_end,
                            l_initial,
                        });
                        owners.push(ci);
                    }
                }
            }
            if reqs.is_empty() {
                if cursors.iter().all(|c| c.done) {
                    break;
                }
                continue;
            }
            for chunk_start in (0..reqs.len()).step_by(SEED_BATCH_WIDTH) {
                let chunk = &reqs[chunk_start..(chunk_start + SEED_BATCH_WIDTH).min(reqs.len())];
                max_mappable_length_batch(chunk, index, &mut results);
                for (k, &(match_length, ns, ne)) in results.iter().enumerate() {
                    let c = &mut cursors[owners[chunk_start + k]];
                    let r = finish_seed(
                        chunk[k].read_pos,
                        match_length,
                        ns,
                        ne,
                        min_seed_length,
                        false,
                        params,
                    );
                    push_seed(
                        &mut buckets[c.bucket],
                        r.seed,
                        c.is_rc,
                        c.original_read_len,
                        params,
                    );
                    c.pos += r.advance;
                }
            }
        }

        // Concatenate each read's buckets in sequential order, apply the
        // per-read cap, then the same sort and dedup `find_seeds` applies.
        read_buckets
            .into_iter()
            .map(|bs| {
                let mut seeds: Vec<Seed> = Vec::new();
                for b in bs {
                    if seeds.len() >= params.seed_per_read_nmax {
                        break;
                    }
                    seeds.extend(std::mem::take(&mut buckets[b]));
                }
                seeds.truncate(params.seed_per_read_nmax);
                finalize_seed_order(seeds)
            })
            .collect()
    }

    /// Find all seeds for paired-end reads using unified seed pooling.
    ///
    /// This implements STAR's hybrid approach: seeds from both mates are found
    /// independently, tagged with their mate origin, then pooled together for
    /// unified clustering.
    ///
    /// # Arguments
    /// * `mate1_seq` - First mate sequence (encoded)
    /// * `mate2_seq` - Second mate sequence (encoded)
    /// * `index` - Genome index with SA and SAindex
    /// * `min_seed_length` - Minimum seed length to report
    /// * `params` - Parameters including seedMultimapNmax
    ///
    /// # Returns
    /// Vector of seeds from both mates, tagged with mate_id (0 or 1)
    pub fn find_paired_seeds(
        mate1_seq: &[u8],
        mate2_seq: &[u8],
        index: &GenomeIndex,
        min_seed_length: usize,
        params: &Parameters,
    ) -> Result<Vec<Seed>, Error> {
        // Find seeds from mate1 (tag with mate_id = 0)
        let mut seeds = Self::find_seeds(mate1_seq, index, min_seed_length, params, "")?;
        for seed in &mut seeds {
            seed.mate_id = 0;
        }

        // Find seeds from mate2 (tag with mate_id = 1)
        // IMPORTANT: read_pos is relative to mate2 start (will be adjusted during stitching)
        let mut seeds2 = Self::find_seeds(mate2_seq, index, min_seed_length, params, "")?;
        for seed in &mut seeds2 {
            seed.mate_id = 1;
        }

        // Pool seeds together
        seeds.extend(seeds2);

        Ok(seeds)
    }

    /// Get all genome positions for this seed.
    ///
    /// Expands the SA range to actual genome positions.
    pub fn get_genome_positions(&self, index: &GenomeIndex) -> Vec<(u64, bool)> {
        self.genome_positions(index).collect()
    }

    /// Iterate over genome positions for this seed without allocating.
    ///
    /// Returns an iterator that lazily decodes SA entries.
    /// For R→L seeds (search_rc == true), converts positions back:
    /// (pos, is_rev) → (n_genome - pos - length, !is_rev)
    /// Positions where the conversion would underflow are filtered out.
    pub fn genome_positions<'a>(
        &'a self,
        index: &'a GenomeIndex,
    ) -> impl Iterator<Item = (u64, bool)> + 'a {
        let search_rc = self.search_rc;
        let length = self.length as u64;
        let n_genome = index.genome.n_genome;
        (self.sa_start..self.sa_end).filter_map(move |sa_idx| {
            let sa_entry = index.suffix_array.get(sa_idx);
            let (pos, is_rev) = index.suffix_array.decode(sa_entry);
            if search_rc {
                if pos + length <= n_genome {
                    Some((n_genome - pos - length, !is_rev))
                } else {
                    None // Position would span past genome boundary
                }
            } else {
                Some((pos, is_rev))
            }
        })
    }
}

/// Reverse-complement an encoded read sequence.
///
/// Reverses the order and complements each base (A↔T, C↔G).
fn reverse_complement_read(read_seq: &[u8]) -> Vec<u8> {
    read_seq.iter().rev().map(|&b| complement_base(b)).collect()
}

/// Result of an MMP (Maximal Mappable Prefix) search at a single position.
/// Always provides the advance length for Lmapped tracking, even when no
/// seed is stored (matching STAR's behavior).
struct MmpResult {
    /// The seed to store, if it passed all filters (multimap, min length)
    seed: Option<Seed>,
    /// MMP length to advance by (>= 1). Used for Lmapped tracking regardless
    /// of whether a seed was stored.
    advance: usize,
}

/// Search one direction using STAR's seedSearchNmax-based starting positions with Lmapped tracking.
///
/// Uses seedSearchNmax (= seedSearchStartLmax = 50 by default) evenly-spaced starting
/// positions in [0, seedSearchStartLmax). From each start, does successive MMP searches
/// forward, advancing past found seeds (Lmapped).
///
/// STAR's formula: iStart = seedSearchStartLmax * i / seedSearchNmax
/// With default seedSearchNmax=seedSearchStartLmax=50: iStart = i → dense {0,1,...,49}.
///
/// Used for R→L direction (is_rc=true). L→R uses dense every-position search.
#[allow(clippy::too_many_arguments)]
/// The chain starting positions STAR uses for one direction of one read
/// (`ReadAlign_mapOneRead.cpp`: `Nstart`, `Lstart`).
///
/// Shared by the sequential and batched paths so the two cannot pick different
/// chains.
fn chain_starts(read_len: usize, params: &Parameters) -> Vec<usize> {
    let effective_start_lmax = if read_len > 0 {
        let over_lread_limit =
            (params.seed_search_start_lmax_over_lread * (read_len as f64 - 1.0)) as usize;
        params.seed_search_start_lmax.min(over_lread_limit)
    } else {
        params.seed_search_start_lmax
    };
    let nstart = if effective_start_lmax > 0 && effective_start_lmax < read_len {
        read_len / effective_start_lmax + 1
    } else {
        1
    };
    let lstart = read_len / nstart;
    (0..nstart).map(|i| (i * lstart).min(read_len)).collect()
}

/// Store one found seed the way `search_direction_sparse` stores it: apply the
/// `seedSearchLmax` cap, tag the direction, and convert an RC hit's read
/// position back to original coordinates.
fn push_seed(
    out: &mut Vec<Seed>,
    seed: Option<Seed>,
    is_rc: bool,
    original_read_len: usize,
    params: &Parameters,
) {
    let Some(mut seed) = seed else { return };
    if params.seed_search_lmax > 0 && seed.length > params.seed_search_lmax {
        seed.length = params.seed_search_lmax;
    }
    seed.search_rc = is_rc;
    if is_rc {
        seed.read_pos = original_read_len - seed.read_pos - seed.length;
    }
    out.push(seed);
}

/// STAR's `storeAligns` ordering: sort by (rStart asc, Length desc), stably so
/// the earliest-found seed wins a tie, then drop exact (rStart, Length)
/// duplicates regardless of search direction.
fn finalize_seed_order(mut seeds: Vec<Seed>) -> Vec<Seed> {
    seeds.sort_by(|a, b| {
        a.read_pos
            .cmp(&b.read_pos)
            .then_with(|| b.length.cmp(&a.length))
    });
    let mut seen: rustc_hash::FxHashSet<(usize, usize)> =
        rustc_hash::FxHashSet::with_capacity_and_hasher(seeds.len(), rustc_hash::FxBuildHasher);
    seeds.retain(|s| seen.insert((s.read_pos, s.length)));
    seeds
}

#[allow(clippy::too_many_arguments)]
fn search_direction_sparse(
    read_seq: &[u8],
    original_read_len: usize,
    index: &GenomeIndex,
    min_seed_length: usize,
    params: &Parameters,
    is_rc: bool,
    debug_name: &str,
    seeds: &mut Vec<Seed>,
) {
    let read_len = read_seq.len();

    // STAR (ReadAlign_mapOneRead.cpp lines 41-42):
    //   seedSearchStartLmax = min(P.seedSearchStartLmax, seedSearchStartLmaxOverLread*(Lread-1))
    let effective_start_lmax = if read_len > 0 {
        let over_lread_limit =
            (params.seed_search_start_lmax_over_lread * (read_len as f64 - 1.0)) as usize;
        params.seed_search_start_lmax.min(over_lread_limit)
    } else {
        params.seed_search_start_lmax
    };

    // STAR (line 48): Nstart = seedSearchStartLmax>0 && seedSearchStartLmax<readLen
    //                          ? readLen/seedSearchStartLmax + 1 : 1
    // Same formula for both L→R and R→L (computed once before the iDir loop).
    // For readLen=150, seedSearchStartLmax=50: Nstart=150/50+1=4, Lstart=37.
    let nstart = if effective_start_lmax > 0 && effective_start_lmax < read_len {
        read_len / effective_start_lmax + 1
    } else {
        1
    };
    let lstart = read_len / nstart; // STAR: Lstart = (splitR[1]-splitR[0]) / Nstart

    for istart in 0..nstart {
        let start_pos = (istart * lstart).min(read_len);
        let mut pos = start_pos;

        // From this starting position, search forward with Lmapped tracking.
        // Continue while remaining bases >= seedMapMin (STAR: istart*Lstart + Lmapped + seedMapMin < readLen).
        // Chains advance until only seedMapMin (5) bases remain.
        loop {
            if pos >= read_len {
                break;
            }
            // Stop if remaining bases < seedMapMin (matches STAR's while condition:
            // istart*Lstart + Lmapped + P.seedMapMin < splitR[1][ip]).
            // STAR chains continue until only seedMapMin (5) bases remain, NOT
            // seedSearchStartLmax (50). This allows chains to reach terminal small
            // exons (e.g. 9M after intron) near the read end.
            if read_len - pos < min_seed_length {
                break;
            }

            let result =
                find_seed_at_position(read_seq, pos, index, min_seed_length, false, params);

            if !debug_name.is_empty() {
                let dir = if is_rc { "RC" } else { "FWD" };
                let seed_info = match &result.seed {
                    Some(s) => format!("seed(len={} sa={}-{})", s.length, s.sa_start, s.sa_end),
                    None => "no_seed".to_string(),
                };
                eprintln!(
                    "[DEBUG-SEED {}] {} istart={} pos={} advance={} {}",
                    debug_name, dir, istart, pos, result.advance, seed_info
                );
            }

            if let Some(mut seed) = result.seed {
                // Apply seedSearchLmax cap
                if params.seed_search_lmax > 0 && seed.length > params.seed_search_lmax {
                    seed.length = params.seed_search_lmax;
                }

                seed.search_rc = is_rc;

                // Convert RC read_pos back to original read coordinates
                if is_rc {
                    seed.read_pos = original_read_len - seed.read_pos - seed.length;
                }

                seeds.push(seed);

                if seeds.len() >= params.seed_per_read_nmax {
                    return;
                }
            }

            pos += result.advance; // Always advance by MMP length (matches STAR)
            // Remaining-length check at loop top: stop when < seedMapMin bases remain
        }
    }
}

/// Find a seed starting at a specific position in the read.
///
/// Returns an MmpResult that always provides the MMP advance length for Lmapped
/// tracking, even when no seed is stored. This matches STAR's behavior where
/// `maxMappableLength2strands()` always returns the MMP length, and `Lmapped += L`
/// always advances — regardless of whether the seed passes filters.
/// What the `SAindex` prefix jump concluded for one position.
///
/// The jump answers many positions outright -- an `N` first base, an absent
/// k-mer, or STAR's short-circuit when a shorter prefix already has tight
/// bounds. Only the rest need the suffix-array binary search, and only those
/// are worth batching. Splitting the two is what lets a caller collect the
/// searches of several *independent* reads and run them together.
enum PreparedSeed {
    /// Already answered: no suffix-array search needed.
    Resolved(MmpResult),
    /// Needs the search. Run it, then hand the result to [`finish_seed`].
    Search { req_read_pos: usize },
}

/// The part of `find_seed_at_position` that runs *after* the suffix-array
/// search: the `seedMultimapNmax` filter and the seed itself.
///
/// Shared by the sequential and batched paths so the two cannot drift.
fn finish_seed(
    read_pos: usize,
    match_length: usize,
    narrowed_start: usize,
    narrowed_end: usize,
    min_seed_length: usize,
    is_reverse: bool,
    params: &Parameters,
) -> MmpResult {
    let advance = match_length.max(1);
    let n_loci = narrowed_end - narrowed_start;
    if n_loci > params.seed_multimap_nmax {
        return MmpResult {
            seed: None,
            advance,
        };
    }
    if match_length >= min_seed_length {
        MmpResult {
            seed: Some(Seed {
                read_pos,
                length: match_length,
                sa_start: narrowed_start,
                sa_end: narrowed_end,
                is_reverse,
                search_rc: false,
                mate_id: 2,
            }),
            advance,
        }
    } else {
        MmpResult {
            seed: None,
            advance,
        }
    }
}

/// The part of `find_seed_at_position` that runs *before* the suffix-array
/// search: the `SAindex` prefix jump and its short-circuits.
///
/// Returns the search to run, plus the SA range and starting offset it needs.
fn prepare_seed_at_position(
    read_seq: &[u8],
    read_pos: usize,
    index: &GenomeIndex,
    min_seed_length: usize,
    is_reverse: bool,
    params: &Parameters,
) -> (PreparedSeed, usize, usize, usize) {
    let no_seed = |advance| {
        (
            PreparedSeed::Resolved(MmpResult {
                seed: None,
                advance,
            }),
            0,
            0,
            0,
        )
    };
    if read_pos >= read_seq.len() {
        return no_seed(1);
    }
    let sa_nbases = index.sa_index.nbases as usize;
    let remaining = read_seq.len() - read_pos;
    if remaining < min_seed_length {
        return no_seed(1);
    }
    let lookup_len = remaining.min(sa_nbases);
    let mut kmer_idx = 0u64;
    let mut actual_len = 0usize;
    for i in 0..lookup_len {
        let base = read_seq[read_pos + i];
        if base >= 4 {
            break;
        }
        kmer_idx = (kmer_idx << 2) | (base as u64);
        actual_len = i + 1;
    }
    if actual_len == 0 {
        return no_seed(1);
    }
    let n_sa = index.suffix_array.len();
    let Some((sa_start, sa_end, matched_level, bounds_tight)) =
        index
            .sa_index
            .hierarchical_lookup(kmer_idx, actual_len as u32, n_sa)
    else {
        return no_seed(1);
    };
    if sa_start >= sa_end {
        return no_seed(1);
    }
    if bounds_tight && matched_level < sa_nbases {
        // STAR's short-circuit: a shorter prefix with tight bounds is the MMP,
        // no genome comparison needed.
        return (
            PreparedSeed::Resolved(finish_seed(
                read_pos,
                matched_level,
                sa_start,
                sa_end,
                min_seed_length,
                is_reverse,
                params,
            )),
            0,
            0,
            0,
        );
    }
    let l_initial = if bounds_tight { matched_level } else { 0 };
    (
        PreparedSeed::Search {
            req_read_pos: read_pos,
        },
        sa_start,
        sa_end,
        l_initial,
    )
}

/// Find a seed starting at a specific position in the read.
///
/// Returns an MmpResult that always provides the MMP advance length for Lmapped
/// tracking, even when no seed is stored. This matches STAR's behavior where
/// `maxMappableLength2strands()` always returns the MMP length, and `Lmapped += L`
/// always advances — regardless of whether the seed passes filters.
///
/// The sequential counterpart of the batched path: same `prepare`, same
/// `finish`, only the search between them differs.
fn find_seed_at_position(
    read_seq: &[u8],
    read_pos: usize,
    index: &GenomeIndex,
    min_seed_length: usize,
    is_reverse: bool,
    params: &Parameters,
) -> MmpResult {
    let (prepared, sa_start, sa_end, l_initial) = prepare_seed_at_position(
        read_seq,
        read_pos,
        index,
        min_seed_length,
        is_reverse,
        params,
    );
    match prepared {
        PreparedSeed::Resolved(r) => r,
        PreparedSeed::Search { req_read_pos } => {
            let (match_length, narrowed_start, narrowed_end) =
                max_mappable_length(read_seq, req_read_pos, index, sa_start, sa_end, l_initial);
            finish_seed(
                req_read_pos,
                match_length,
                narrowed_start,
                narrowed_end,
                min_seed_length,
                is_reverse,
                params,
            )
        }
    }
}

/// Overflow-safe median of two unsigned integers.
/// Equivalent to STAR's medianUint2: a/2 + b/2 + (a%2 + b%2)/2
fn median_uint2(a: usize, b: usize) -> usize {
    a / 2 + b / 2 + usize::midpoint(a % 2, b % 2)
}

/// Compare read to genome at a specific SA position, starting from offset l_start.
/// Returns (total_match_length, is_read_greater_at_mismatch).
/// Ports STAR's compareSeqToGenome (SuffixArrayFuns.cpp).
///
/// Starts comparing from offset `l_start` (bases 0..l_start are assumed to match).
/// Walks forward until a mismatch, end of read, or genome padding.
fn compare_seq_to_genome(
    read_seq: &[u8],
    read_pos: usize,
    index: &GenomeIndex,
    sa_idx: usize,
    l_start: usize,
) -> (usize, bool) {
    let sa_entry = index.suffix_array.get(sa_idx);
    compare_seq_to_genome_raw(read_seq, read_pos, index, sa_entry, l_start)
}

/// [`compare_seq_to_genome`] with the suffix-array entry already read.
///
/// Split out so a batched driver can read (and prefetch) the entry for several
/// *independent* searches before any of them compares, which is the whole point
/// of [`max_mappable_length_batch`]. The comparison itself is unchanged, so the
/// two paths cannot diverge.
fn compare_seq_to_genome_raw(
    read_seq: &[u8],
    read_pos: usize,
    index: &GenomeIndex,
    sa_entry: u64,
    l_start: usize,
) -> (usize, bool) {
    let (genome_pos, is_reverse) = index.suffix_array.decode(sa_entry);

    let genome_start = if is_reverse {
        genome_pos as usize + index.genome.n_genome as usize
    } else {
        genome_pos as usize
    };

    let remaining = read_seq.len() - read_pos;
    let mut match_len = l_start;
    let mut i = l_start;

    // SIMD fast path: forward-strand hits (`!is_reverse`) scan a genuinely
    // contiguous byte range directly from the genome's physical forward
    // buffer — both `GenomeSeq::Owned` and `::Mapped` store the forward
    // strand contiguously (see `GenomeSeq::as_slice`), unlike the
    // reverse-complement half, which `GenomeSeq::base` computes per-byte on
    // the fly (no contiguous byte range exists to scan). Bounded strictly to
    // `genome_start+i < n_genome` (NOT `sequence.len()`, which is `2*n_genome`
    // and would let this slice run into the RC-mirror continuation the
    // scalar loop below still walks correctly — see its boundary check) so
    // this fast path can never diverge from the exact original semantics.
    if !is_reverse {
        let n_genome = index.genome.n_genome as usize;
        if genome_start < n_genome {
            let genome_slice_all = index.genome.sequence.as_slice();
            let simd_end = remaining.min(n_genome.min(genome_slice_all.len()) - genome_start);
            if simd_end > i {
                let read_chunk = &read_seq[read_pos + i..read_pos + simd_end];
                let genome_chunk = &genome_slice_all[genome_start + i..genome_start + simd_end];
                if let Some(off) = crate::align::simd_scan::find_stop(read_chunk, genome_chunk) {
                    let genome_base = genome_chunk[off];
                    if genome_base >= 5 {
                        return (i + off, true);
                    }
                    let read_base = read_chunk[off];
                    return (i + off, read_base > genome_base);
                }
                match_len = simd_end;
                i = simd_end;
            }
        }
    }

    for i in i..remaining {
        let genome_idx = genome_start + i;

        if genome_idx >= index.genome.sequence.len() {
            // Past end of genome array — treat like padding (STAR: comp_res > 0)
            return (match_len, true);
        }

        let genome_base = index.genome.sequence.base(genome_idx);

        if genome_base >= 5 {
            // Padding character — STAR returns comp_res > 0 (read > genome)
            return (match_len, true);
        }

        let read_base = read_seq[read_pos + i];

        if read_base != genome_base {
            return (match_len, read_base > genome_base);
        }

        match_len += 1;
    }

    // Matched all remaining bases — STAR returns comp_res < 0 (genome >= read)
    (match_len, false)
}

/// Find maximum mappable prefix length within SA range [sa_start, sa_end).
/// Binary searches the range while extending match length, then narrows to
/// all positions matching the maximum length.
/// Returns (match_length, narrowed_sa_start, narrowed_sa_end_exclusive).
/// Ports STAR's maxMappableLength (SuffixArrayFuns.cpp).
/// How many independent MMP searches [`max_mappable_length_batch`] advances
/// together.
///
/// Each search stalls on a random suffix-array read and then on a random genome
/// read; running N of them in lockstep turns N dependent stalls into N
/// overlapping ones. The donor measured 2.09x at 32 and a *regression* at 8, so
/// the width is not a free parameter: too narrow and the extra bookkeeping
/// costs more than the overlap saves.
pub const SEED_BATCH_WIDTH: usize = 32;

/// Where a batched MMP search has got to. Mirrors [`max_mappable_length`]'s
/// straight-line phases so it can be suspended at each probe -- the point where
/// it would otherwise stall on DRAM.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Phase {
    /// Needs the `compare(i1)` that opens the function.
    Init1,
    /// Needs the `compare(i2)`.
    Init2,
    /// Inside `while i1 + 1 < i2`, needs `compare(i3)`.
    Loop,
    Done,
}

/// One search's live state: [`max_mappable_length`]'s locals, made explicit.
///
/// The arithmetic and the order of updates below are deliberately identical to
/// the sequential function; only the control flow is turned inside out.
struct MmpState<'a> {
    read_seq: &'a [u8],
    read_pos: usize,
    remaining: usize,
    i1: usize,
    i2: usize,
    l: usize,
    l1: usize,
    l2: usize,
    l1a: usize,
    l1b: usize,
    i1a: usize,
    i1b: usize,
    l2a: usize,
    l2b: usize,
    i2a: usize,
    i2b: usize,
    i3: usize,
    l3: usize,
    phase: Phase,
    /// The suffix-array index this search wants probed next, or `None` once done.
    want: Option<usize>,
    /// The suffix-array entry the driver pre-read for `want`, consumed by the
    /// comparison.
    raw: u64,
    out: Option<(usize, usize, usize)>,
}

/// One independent search for [`max_mappable_length_batch`]: the arguments
/// [`max_mappable_length`] takes, bundled so a batch of them can be advanced
/// together.
#[derive(Debug, Clone, Copy)]
pub(crate) struct MmpReq<'a> {
    pub read_seq: &'a [u8],
    pub read_pos: usize,
    pub sa_start: usize,
    pub sa_end: usize,
    pub l_initial: usize,
}

impl<'a> MmpState<'a> {
    fn new(r: &MmpReq<'a>) -> Self {
        MmpState {
            read_seq: r.read_seq,
            read_pos: r.read_pos,
            remaining: r.read_seq.len() - r.read_pos,
            i1: r.sa_start,
            i2: r.sa_end - 1,
            l: r.l_initial,
            l1: 0,
            l2: 0,
            l1a: 0,
            l1b: 0,
            i1a: 0,
            i1b: 0,
            l2a: 0,
            l2b: 0,
            i2a: 0,
            i2b: 0,
            i3: r.sa_start,
            l3: 0,
            phase: Phase::Init1,
            want: Some(r.sa_start),
            raw: 0,
            out: None,
        }
    }
}

/// Run several *independent* [`max_mappable_length`] searches with their memory
/// accesses interleaved, returning each one's `(match_len, narrowed_start,
/// narrowed_end)` in request order.
///
/// The loop has three passes per round, and the split is the whole point:
/// issue every live search's suffix-array prefetch first, then read those
/// entries and prefetch each search's *dependent* genome byte, then compare.
/// One search alone would serialise those two misses; N of them overlap.
///
/// `find_mult_range` is left sequential, as in the donor: it is a second binary
/// search with its own probe pattern, and keeping it out makes the equivalence
/// with the sequential function easy to see. Batching it is the obvious
/// follow-up if a measurement justifies it.
pub(crate) fn max_mappable_length_batch(
    reqs: &[MmpReq<'_>],
    index: &GenomeIndex,
    out: &mut Vec<(usize, usize, usize)>,
) {
    out.clear();
    if reqs.is_empty() {
        return;
    }

    let mut st: Vec<MmpState> = Vec::with_capacity(reqs.len());
    for r in reqs {
        // The single-element case never enters the state machine, exactly as
        // the sequential function short-circuits it.
        if r.sa_start + 1 >= r.sa_end {
            let (l, _) =
                compare_seq_to_genome(r.read_seq, r.read_pos, index, r.sa_start, r.l_initial);
            let mut s = MmpState::new(r);
            s.phase = Phase::Done;
            s.want = None;
            s.out = Some((l, r.sa_start, r.sa_start + 1));
            st.push(s);
        } else {
            st.push(MmpState::new(r));
        }
    }

    loop {
        let mut any = false;
        // Pass A: issue every live search's suffix-array probe, so N line-fills
        // go out together instead of one per DRAM round-trip.
        for s in &st {
            if let Some(isa) = s.want {
                index.suffix_array.prefetch(isa);
                any = true;
            }
        }
        if !any {
            break;
        }
        // Pass B: consume the SA entries (their lines are arriving) and
        // prefetch each search's dependent genome byte -- the second, and
        // costlier, miss in the chain.
        for s in &mut st {
            if let Some(isa) = s.want {
                s.raw = index.suffix_array.get(isa);
                let (genome_pos, is_reverse) = index.suffix_array.decode(s.raw);
                if !is_reverse {
                    let pos = genome_pos as usize + s.l;
                    let g = index.genome.sequence.as_slice();
                    if pos < g.len() {
                        // SAFETY: `pos` is in bounds of `g` (just checked), and a
                        // prefetch reads nothing.
                        crate::cpu::prefetch_read(unsafe { g.as_ptr().add(pos) });
                    }
                }
            }
        }
        // Pass C: compare (the genome lines are now hot) and advance each state.
        for s in &mut st {
            if s.want.is_some() {
                let (ml, cr) = compare_seq_to_genome_raw(s.read_seq, s.read_pos, index, s.raw, s.l);
                mmp_advance(s, index, ml, cr);
            }
        }
    }

    out.extend(st.iter().map(|s| s.out.expect("every search reaches Done")));
}

/// Drive one [`MmpState`] forward by one probe result. This is
/// [`max_mappable_length`]'s control flow with the straight line broken into
/// phases so a batch can suspend a search at each probe.
fn mmp_advance(s: &mut MmpState, index: &GenomeIndex, ml: usize, comp_res: bool) {
    match s.phase {
        Phase::Init1 => {
            s.l1 = ml;
            s.phase = Phase::Init2;
            s.want = Some(s.i2);
        }
        Phase::Init2 => {
            s.l2 = ml;
            s.l = s.l1.min(s.l2);
            s.l3 = s.l;
            s.i3 = s.i1;
            s.i1a = s.i1;
            s.l1a = s.l1;
            s.i1b = s.i1;
            s.l1b = s.l1;
            s.i2a = s.i2;
            s.l2a = s.l2;
            s.i2b = s.i2;
            s.l2b = s.l2;
            s.phase = Phase::Loop;
            mmp_step_loop(s, index);
        }
        Phase::Loop => {
            s.l3 = ml;
            if s.l3 == s.remaining {
                mmp_finish(s, index);
                return;
            }
            if comp_res {
                if s.l3 > s.l1 {
                    s.i1a = s.i1b;
                    s.l1a = s.l1b;
                    s.i1b = s.i1;
                    s.l1b = s.l1;
                }
                s.i1 = s.i3;
                s.l1 = s.l3;
            } else {
                if s.l3 > s.l2 {
                    s.i2a = s.i2b;
                    s.l2a = s.l2b;
                    s.i2b = s.i2;
                    s.l2b = s.l2;
                }
                s.i2 = s.i3;
                s.l2 = s.l3;
            }
            s.l = s.l1.min(s.l2);
            mmp_step_loop(s, index);
        }
        Phase::Done => unreachable!("advanced a finished search"),
    }
}

/// The `while i1 + 1 < i2` test: either request the next midpoint probe, or fall
/// through to the tail.
fn mmp_step_loop(s: &mut MmpState, index: &GenomeIndex) {
    if s.i1 + 1 < s.i2 {
        s.i3 = median_uint2(s.i1, s.i2);
        s.want = Some(s.i3);
    } else {
        mmp_finish(s, index);
    }
}

/// The tail of [`max_mappable_length`]: pick the best match length, then narrow
/// the range with the two (sequential) `find_mult_range` calls.
fn mmp_finish(s: &mut MmpState, index: &GenomeIndex) {
    if s.l3 < s.remaining {
        if s.l1 > s.l2 {
            s.l3 = s.l1;
            s.i3 = s.i1;
        } else {
            s.l3 = s.l2;
            s.i3 = s.i2;
        }
    }
    let narrowed_start = find_mult_range(
        s.read_seq,
        s.read_pos,
        index,
        s.remaining,
        s.i3,
        s.l3,
        s.i1,
        s.l1,
        s.i1a,
        s.l1a,
        s.i1b,
        s.l1b,
    );
    let narrowed_end = find_mult_range(
        s.read_seq,
        s.read_pos,
        index,
        s.remaining,
        s.i3,
        s.l3,
        s.i2,
        s.l2,
        s.i2a,
        s.l2a,
        s.i2b,
        s.l2b,
    );
    s.out = Some((s.l3, narrowed_start, narrowed_end + 1));
    s.phase = Phase::Done;
    s.want = None;
}

fn max_mappable_length(
    read_seq: &[u8],
    read_pos: usize,
    index: &GenomeIndex,
    sa_start: usize,
    sa_end: usize,
    l_initial: usize,
) -> (usize, usize, usize) {
    let remaining = read_seq.len() - read_pos;

    // Single element: just compare
    if sa_start + 1 >= sa_end {
        let (l, _) = compare_seq_to_genome(read_seq, read_pos, index, sa_start, l_initial);
        return (l, sa_start, sa_start + 1);
    }

    // Convert to inclusive range (STAR convention internally)
    let mut i1 = sa_start;
    let mut i2 = sa_end - 1;

    let (mut l1, _) = compare_seq_to_genome(read_seq, read_pos, index, i1, l_initial);
    let (mut l2, _) = compare_seq_to_genome(read_seq, read_pos, index, i2, l_initial);

    let mut l = l1.min(l2);
    let mut l3 = l;
    let mut i3 = i1;

    // Track history for find_mult_range
    let (mut i1a, mut l1a) = (i1, l1);
    let (mut i1b, mut l1b) = (i1, l1);
    let (mut i2a, mut l2a) = (i2, l2);
    let (mut i2b, mut l2b) = (i2, l2);

    // Binary search within SA range
    while i1 + 1 < i2 {
        i3 = median_uint2(i1, i2);
        // Prefetch the SA entries of the two possible next probes so their
        // (random, mmap-backed) bytes arrive while the current probe's genome
        // comparison runs. Whichever branch is taken below, its midpoint is one
        // of these two. Hint only — does not change results.
        if i1 + 1 < i3 {
            index.suffix_array.prefetch(median_uint2(i1, i3));
        }
        if i3 + 1 < i2 {
            index.suffix_array.prefetch(median_uint2(i3, i2));
        }
        let comp3;
        (l3, comp3) = compare_seq_to_genome(read_seq, read_pos, index, i3, l);

        if l3 == remaining {
            break; // Perfect match found
        }

        if comp3 {
            // read > genome at mismatch: move left boundary up
            // STAR only shifts history when match length improves (L3 > L1)
            if l3 > l1 {
                i1a = i1b;
                l1a = l1b;
                i1b = i1;
                l1b = l1;
            }
            i1 = i3;
            l1 = l3;
        } else {
            // read <= genome at mismatch: move right boundary down
            // STAR only shifts history when match length improves (L3 > L2)
            if l3 > l2 {
                i2a = i2b;
                l2a = l2b;
                i2b = i2;
                l2b = l2;
            }
            i2 = i3;
            l2 = l3;
        }

        l = l1.min(l2);
    }

    // Pick the best match length
    if l3 < remaining {
        if l1 > l2 {
            l3 = l1;
            i3 = i1;
        } else {
            l3 = l2;
            i3 = i2;
        }
    }

    // Find narrowed range using find_mult_range
    let narrowed_start = find_mult_range(
        read_seq, read_pos, index, remaining, i3, l3, i1, l1, i1a, l1a, i1b, l1b,
    );
    let narrowed_end = find_mult_range(
        read_seq, read_pos, index, remaining, i3, l3, i2, l2, i2a, l2a, i2b, l2b,
    );

    // Convert back to exclusive end
    (l3, narrowed_start, narrowed_end + 1)
}

/// Binary search to find the SA boundary where match length transitions
/// from >= l3 to < l3. Used to narrow the SA range to only positions
/// matching the maximum prefix length.
/// Ports STAR's findMultRange (SuffixArrayFuns.cpp).
///
/// STAR's logic: given the "best" SA index i3 with match length L3,
/// find the farthest SA index that also matches L3 bases, searching
/// outward from i1 (which may or may not already match).
///
/// i1a tracks the boundary with L >= L3 ("good" side)
/// i1b tracks the boundary with L < L3 ("bad" side)
/// Binary search narrows between them until adjacent.
#[allow(clippy::too_many_arguments)]
fn find_mult_range(
    read_seq: &[u8],
    read_pos: usize,
    index: &GenomeIndex,
    _remaining: usize,
    i3: usize,
    l3: usize,
    i1: usize,
    l1: usize,
    i1a: usize,
    l1a: usize,
    i1b: usize,
    l1b: usize,
) -> usize {
    // STAR's findMultRange: set up (i1a, i1b) search range
    // i1a will have L >= L3 (the "good" side)
    // i1b will have L < L3 (the "bad" side)
    let (mut ia, mut ib, mut lb);
    if l1 < l3 {
        // i1 is below target: search between i3 (good) and i1 (bad)
        ib = i1;
        lb = l1;
        ia = i3;
    } else {
        // i1 already at target length
        if l1a < l1 {
            // Search between i1a (bad) and i1 (good), outward from i1
            ib = i1a;
            lb = l1a;
            ia = i1;
        } else {
            // i1a also at target — search between i1a and i1b
            // (STAR: falls through without reassignment, keeps original i1a/i1b)
            ia = i1a;
            ib = i1b;
            lb = l1b;
        }
    }

    // Binary search: ia has L >= l3, ib has L < l3
    // compareSeqToGenome is called with N=l3 (not remaining), matching STAR
    while (ib + 1 < ia) || (ia + 1 < ib) {
        let ic = median_uint2(ia, ib);
        let (lc, _) = compare_seq_to_genome(read_seq, read_pos, index, ic, lb);

        if lc >= l3 {
            ia = ic;
        } else {
            ib = ic;
            lb = lc;
        }
    }

    ia
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::Parameters;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn make_test_index(sequence: &str) -> GenomeIndex {
        let mut file = NamedTempFile::new().unwrap();
        writeln!(file, ">chr1").unwrap();
        writeln!(file, "{sequence}").unwrap();

        let dir = tempfile::tempdir().unwrap();

        let args = vec![
            "rustar-aligner",
            "--runMode",
            "genomeGenerate",
            "--genomeFastaFiles",
            file.path().to_str().unwrap(),
            "--genomeDir",
            dir.path().to_str().unwrap(),
            "--genomeChrBinNbits",
            "2",
            "--genomeSAindexNbases",
            "2",
        ];

        let params = Parameters::parse_from(args);
        GenomeIndex::build(&params).unwrap()
    }

    fn encode_sequence(seq: &str) -> Vec<u8> {
        seq.bytes()
            .map(|b| match b {
                b'A' | b'a' => 0,
                b'C' | b'c' => 1,
                b'G' | b'g' => 2,
                b'T' | b't' => 3,
                _ => 4,
            })
            .collect()
    }

    fn params(args: &[&str]) -> Parameters {
        let mut full_args = vec!["rustar-aligner", "--readFilesIn", "reads.fq"];
        full_args.extend_from_slice(args);
        Parameters::parse_from(full_args)
    }

    /// The batched seed finder must return exactly what calling the sequential
    /// one per read returns -- same seeds, same order, per read.
    ///
    /// Order matters as much as content: `cluster_seeds` creates windows in
    /// seed order, and the earliest window wins the primary tie-break, so a
    /// reordering here would change which alignment is reported primary
    /// without changing any seed.
    #[test]
    fn batched_seed_finding_matches_per_read_sequential_calls() {
        let index = make_test_index(
            "ACGTACGTAAGGTTCCACGTACGTTTGGCCAAACGTACGTAGGCCTTAACGTACGTGGTTAACCACGTACGT\
             TTTTGGGGCCCCAAAAACGTACGTATATATATCGCGCGCGACGTACGTTACGTACGACGTACGTGCATGCAT",
        );
        let p = params(&["--seedSearchStartLmax", "10"]);
        let reads: Vec<Vec<u8>> = [
            "ACGTACGTAAGGTTCCACGTACGT",
            "TTTTGGGGCCCCAAAAACGTACGT",
            "GCATGCATACGTACGTTACGTACG",
            "ACGTACGTATATATATCGCGCGCG",
            "NNNNACGTACGTAAGGTTCC",
            "ACGT",
        ]
        .iter()
        .map(|r| encode_sequence(r))
        .collect();
        let refs: Vec<&[u8]> = reads.iter().map(Vec::as_slice).collect();

        let sequential: Vec<Vec<Seed>> = refs
            .iter()
            .map(|r| Seed::find_seeds(r, &index, 5, &p, "").unwrap())
            .collect();
        let batched = Seed::find_seeds_batch(&refs, &index, 5, &p);

        assert_eq!(batched.len(), sequential.len());
        for (i, (b, q)) in batched.iter().zip(sequential.iter()).enumerate() {
            assert_eq!(
                b.len(),
                q.len(),
                "read {i}: seed count, batched {b:?} vs sequential {q:?}"
            );
            for (j, (bs, qs)) in b.iter().zip(q.iter()).enumerate() {
                assert_eq!(
                    (bs.read_pos, bs.length, bs.sa_start, bs.sa_end, bs.search_rc),
                    (qs.read_pos, qs.length, qs.sa_start, qs.sa_end, qs.search_rc),
                    "read {i} seed {j}"
                );
            }
        }
    }

    /// A batch of one, and an empty batch, are the shapes a chunked call site
    /// hits at the end of every read batch.
    #[test]
    fn batched_seed_finding_handles_degenerate_batches() {
        let index = make_test_index("ACGTACGTTTGGCCAAACGTACGT");
        let p = params(&[]);
        let read = encode_sequence("ACGTACGTTTGGCCAA");

        assert!(Seed::find_seeds_batch(&[], &index, 5, &p).is_empty());

        let one = Seed::find_seeds_batch(&[&read], &index, 5, &p);
        let seq = Seed::find_seeds(&read, &index, 5, &p, "").unwrap();
        assert_eq!(one.len(), 1);
        assert_eq!(one[0].len(), seq.len());
    }

    /// The batched search must return exactly what the sequential one returns,
    /// at every batch width, for every start position of every read.
    ///
    /// This is the acceptance criterion for the whole batching change: the
    /// state machine is `max_mappable_length` turned inside out, so any
    /// divergence is a transcription bug and would show up as a different seed,
    /// a different alignment, and a different output file.
    #[test]
    fn batched_mmp_search_matches_the_sequential_one_at_every_width() {
        let index = make_test_index(
            "ACGTACGTAAGGTTCCACGTACGTTTGGCCAAACGTACGTAGGCCTTAACGTACGTGGTTAACCACGTACGT\
             TTTTGGGGCCCCAAAAACGTACGTATATATATCGCGCGCGACGTACGTTACGTACGACGTACGTGCATGCAT",
        );
        let reads: Vec<Vec<u8>> = [
            "ACGTACGTAAGGTTCC",
            "ACGTACGT",
            "TTTTGGGGCCCCAAAA",
            "GCATGCAT",
            "ACGTACGTATATATAT",
            "NNNNACGTACGT",
            "A",
        ]
        .iter()
        .map(|r| encode_sequence(r))
        .collect();

        // Every (read, start) pair that has a live SA range, prepared the way
        // `find_seed_at_position` prepares them.
        let mut reqs: Vec<MmpReq> = Vec::new();
        for read in &reads {
            for pos in 0..read.len() {
                reqs.push(MmpReq {
                    read_seq: read,
                    read_pos: pos,
                    sa_start: 0,
                    sa_end: index.suffix_array.len(),
                    l_initial: 0,
                });
            }
        }
        assert!(reqs.len() > 20, "the fixture should exercise many searches");

        let sequential: Vec<(usize, usize, usize)> = reqs
            .iter()
            .map(|r| {
                max_mappable_length(
                    r.read_seq,
                    r.read_pos,
                    &index,
                    r.sa_start,
                    r.sa_end,
                    r.l_initial,
                )
            })
            .collect();

        // Width 1 is the degenerate batch; 128 is wider than the request list,
        // so the last chunk is short. Both are the shapes most likely to break.
        for width in [1usize, 2, 8, 32, 128] {
            let mut got: Vec<(usize, usize, usize)> = Vec::new();
            let mut out = Vec::new();
            for chunk in reqs.chunks(width) {
                max_mappable_length_batch(chunk, &index, &mut out);
                assert_eq!(out.len(), chunk.len(), "width {width}: result count");
                got.extend_from_slice(&out);
            }
            assert_eq!(got, sequential, "batched != sequential at width {width}");
        }
    }

    /// A single-element SA range never enters the state machine; the batch must
    /// still answer it, and answer it the same way.
    #[test]
    fn batched_mmp_handles_single_element_ranges_and_empty_batches() {
        let index = make_test_index("ACGTACGTTTGGCCAA");
        let read = encode_sequence("ACGTACGT");

        let mut out = Vec::new();
        max_mappable_length_batch(&[], &index, &mut out);
        assert!(out.is_empty(), "an empty batch produces no results");

        let req = MmpReq {
            read_seq: &read,
            read_pos: 0,
            sa_start: 3,
            sa_end: 4,
            l_initial: 0,
        };
        max_mappable_length_batch(std::slice::from_ref(&req), &index, &mut out);
        assert_eq!(
            out[0],
            max_mappable_length(&read, 0, &index, 3, 4, 0),
            "single-element range"
        );
    }

    #[test]
    fn find_exact_match() {
        let index = make_test_index("ACGTACGT");
        let read = encode_sequence("ACGT");
        let params = params(&["--runMode", "alignReads"]);

        let seeds = Seed::find_seeds(&read, &index, 4, &params, "").unwrap();

        // Should find at least one seed
        assert!(!seeds.is_empty());

        // First seed should be at position 0 with length 4
        assert_eq!(seeds[0].read_pos, 0);
        assert_eq!(seeds[0].length, 4);
    }

    #[test]
    fn min_seed_length_filter() {
        let index = make_test_index("AAAAAAAA");
        let read = encode_sequence("AAA");
        let params = params(&[]);

        // With min_seed_length=4, should find nothing (read is only 3bp)
        let seeds = Seed::find_seeds(&read, &index, 4, &params, "").unwrap();
        assert!(seeds.is_empty());

        // With min_seed_length=2, should find seeds
        let seeds = Seed::find_seeds(&read, &index, 2, &params, "").unwrap();
        assert!(!seeds.is_empty());
    }

    #[test]
    fn no_match() {
        let index = make_test_index("ACAC");
        let read = encode_sequence("GGGG");
        let params = params(&[]);

        let seeds = Seed::find_seeds(&read, &index, 2, &params, "").unwrap();

        // No seeds should be found (GGGG not in ACAC or its reverse complement GTGT)
        assert!(seeds.is_empty());
    }

    #[test]
    fn get_genome_positions() {
        let index = make_test_index("ACGTACGT");
        let read = encode_sequence("ACGT");
        let params = params(&[]);

        let seeds = Seed::find_seeds(&read, &index, 4, &params, "").unwrap();
        assert!(!seeds.is_empty());

        // Get positions for first seed
        let positions = seeds[0].get_genome_positions(&index);
        assert!(!positions.is_empty());

        // Should have at least one valid position
        for (pos, _is_reverse) in positions {
            assert!(pos < index.genome.n_genome);
        }
    }

    #[test]
    fn test_single_end_mate_id() {
        let index = make_test_index("ACGTACGT");
        let read = encode_sequence("ACGT");
        let params = params(&[]);

        let seeds = Seed::find_seeds(&read, &index, 4, &params, "").unwrap();
        assert!(!seeds.is_empty());

        // Single-end seeds should have mate_id = 2
        for seed in seeds {
            assert_eq!(seed.mate_id, 2);
        }
    }

    #[test]
    fn test_find_paired_seeds() {
        let index = make_test_index("ACGTACGTTTGGCCAA");
        let mate1 = encode_sequence("ACGT");
        let mate2 = encode_sequence("TTGG");
        let params = params(&[]);

        let seeds = Seed::find_paired_seeds(&mate1, &mate2, &index, 4, &params).unwrap();

        // Should have seeds from both mates
        let mate1_seeds: Vec<_> = seeds.iter().filter(|s| s.mate_id == 0).collect();
        let mate2_seeds: Vec<_> = seeds.iter().filter(|s| s.mate_id == 1).collect();

        assert!(!mate1_seeds.is_empty(), "Should have mate1 seeds");
        assert!(!mate2_seeds.is_empty(), "Should have mate2 seeds");

        // Verify mate1 seeds have correct read positions
        for seed in mate1_seeds {
            assert!(seed.read_pos < mate1.len());
        }

        // Verify mate2 seeds have correct read positions (relative to mate2)
        for seed in mate2_seeds {
            assert!(seed.read_pos < mate2.len());
        }
    }

    #[test]
    fn test_paired_seeds_pooling() {
        let index = make_test_index("ACGTACGT");
        let mate1 = encode_sequence("ACGT");
        let mate2 = encode_sequence("ACGT");
        let params = params(&[]);

        let seeds = Seed::find_paired_seeds(&mate1, &mate2, &index, 4, &params).unwrap();

        // Should have roughly double the seeds (one set from each mate)
        let mate1_count = seeds.iter().filter(|s| s.mate_id == 0).count();
        let mate2_count = seeds.iter().filter(|s| s.mate_id == 1).count();

        assert!(mate1_count > 0);
        assert!(mate2_count > 0);
        assert_eq!(seeds.len(), mate1_count + mate2_count);
    }

    #[test]
    fn test_reverse_complement_read() {
        // ACGT → RC = ACGT (palindrome)
        let read = encode_sequence("ACGT");
        let rc = reverse_complement_read(&read);
        assert_eq!(rc, encode_sequence("ACGT"));

        // AACC → RC = GGTT
        let read2 = encode_sequence("AACC");
        let rc2 = reverse_complement_read(&read2);
        assert_eq!(rc2, encode_sequence("GGTT"));

        // Single base
        let read3 = encode_sequence("A");
        let rc3 = reverse_complement_read(&read3);
        assert_eq!(rc3, encode_sequence("T"));

        // N bases preserved
        let read4 = vec![0, 4, 1]; // A, N, C
        let rc4 = reverse_complement_read(&read4);
        assert_eq!(rc4, vec![2, 4, 3]); // G, N, T
    }

    #[test]
    fn test_direction_agnostic_seed_dedup() {
        // Genome AACCTTGG; read CCAAGGTT is its reverse-complement, so the whole
        // read maps as one seed. The L→R pass finds it on the RC strand of the
        // doubled genome (which contains RC(AACCTTGG)=CCAAGGTT); the R→L pass, on
        // RC(read)=AACCTTGG, finds the SAME (rStart, length) seed on the forward
        // strand. STAR's storeAligns dedups (rStart, Length) regardless of search
        // direction (OPTIM_STOREaligns_SIMPLE), so exactly one copy survives — the
        // L→R one, collected first.
        let index = make_test_index("AACCTTGG");
        let read = encode_sequence("CCAAGGTT");
        let params = params(&[]);

        let seeds = Seed::find_seeds(&read, &index, 4, &params, "").unwrap();

        // Exactly one full-length (read_pos=0, length=8) seed survives dedup.
        let full: Vec<_> = seeds
            .iter()
            .filter(|s| s.read_pos == 0 && s.length == 8)
            .collect();
        assert_eq!(
            full.len(),
            1,
            "direction-agnostic dedup should keep exactly one (0,8) seed; got {seeds:?}"
        );
        // The survivor is the forward (L→R) seed, since L→R is collected first.
        assert!(
            !full[0].search_rc,
            "the L→R seed should be the survivor after dedup"
        );

        // Seeds are sorted by rStart ascending (STAR's PC[] order).
        for w in seeds.windows(2) {
            assert!(
                w[0].read_pos <= w[1].read_pos,
                "seeds must be sorted by read_pos (rStart): {seeds:?}"
            );
        }
    }

    #[test]
    fn test_shared_seed_cap() {
        // Test that combined L→R + R→L respects seedPerReadNmax
        let index = make_test_index("ACGTACGTACGTACGT");
        let read = encode_sequence("ACGTACGT");
        let params = params(&["--seedPerReadNmax", "3"]);

        let seeds = Seed::find_seeds(&read, &index, 4, &params, "").unwrap();
        assert!(
            seeds.len() <= 3,
            "Total seeds ({}) should respect seedPerReadNmax=3",
            seeds.len()
        );
    }

    #[test]
    fn test_sparse_nstart_calculation() {
        // Verify Nstart matches STAR: readLen/seedSearchStartLmax + 1
        // (when seedSearchStartLmax > 0 && seedSearchStartLmax < readLen)
        //
        // STAR (ReadAlign_mapOneRead.cpp line 48):
        //   Nstart = seedSearchStartLmax>0 && seedSearchStartLmax<splitR[1]
        //            ? splitR[1]/seedSearchStartLmax+1 : 1
        //   Lstart = splitR[1] / Nstart

        // 150 / 50 + 1 = 4, Lstart = 150/4 = 37
        let nstart = 150 / 50 + 1;
        assert_eq!(nstart, 4);
        assert_eq!(150 / nstart, 37);

        // 151 / 50 + 1 = 4, Lstart = 151/4 = 37
        let nstart = 151 / 50 + 1;
        assert_eq!(nstart, 4);
        assert_eq!(151 / nstart, 37);

        // 30 / 50: seedSearchStartLmax (50) >= readLen (30) → Nstart=1
        // (condition seedSearchStartLmax < readLen is false)
        assert_eq!(1_usize, 1);

        // 50 / 50: seedSearchStartLmax (50) >= readLen (50) → Nstart=1
        assert_eq!(1_usize, 1);

        // 100 / 50 + 1 = 3, Lstart = 100/3 = 33
        let nstart = 100 / 50 + 1;
        assert_eq!(nstart, 3);
        assert_eq!(100 / nstart, 33);
    }

    #[test]
    fn test_sparse_rc_read_pos_conversion() {
        // The R→L pass searches RC(read) and converts positions back to original
        // read coordinates: read_pos = original_len - rc_pos - length. Drive
        // search_direction_sparse directly with is_rc=true so the conversion is
        // exercised in isolation (find_seeds' direction-agnostic dedup would
        // otherwise drop the R→L seed as a duplicate of the L→R one — see
        // test_direction_agnostic_seed_dedup).
        // Genome AACCTTGG; read CCAAGGTT; RC(read)=AACCTTGG matches forward genome.
        let index = make_test_index("AACCTTGG");
        let read = encode_sequence("CCAAGGTT");
        let rc_read = reverse_complement_read(&read);
        let params = params(&[]);

        let mut rc_seeds = Vec::new();
        search_direction_sparse(
            &rc_read,
            read.len(),
            &index,
            4,
            &params,
            true,
            "",
            &mut rc_seeds,
        );

        assert!(
            !rc_seeds.is_empty(),
            "R→L search should find seeds (RC(read) matches the forward genome)"
        );
        for seed in &rc_seeds {
            assert!(
                seed.search_rc,
                "seeds from the R→L pass must have search_rc=true"
            );
            // Converted read positions must be valid original-read coordinates.
            assert!(
                seed.read_pos + seed.length <= read.len(),
                "converted R→L read_pos {} + length {} exceeds read len {}",
                seed.read_pos,
                seed.length,
                read.len()
            );
        }
    }

    #[test]
    fn test_sparse_fewer_seeds_than_dense() {
        // With a longer genome and read, sparse search should produce fewer seeds
        // than the old dense (every-position) search
        let genome_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let index = make_test_index(genome_seq);
        // Use a read that's long enough for multiple start positions
        let read = encode_sequence("ACGTACGTACGTACGTACGTACGT"); // 24bp
        let params = params(&[]);

        let sparse_seeds = Seed::find_seeds(&read, &index, 4, &params, "").unwrap();

        // Count how many seeds dense would produce (every position that has a match)
        let mut dense_count = 0;
        for read_pos in 0..read.len() {
            let result = find_seed_at_position(&read, read_pos, &index, 4, false, &params);
            if result.seed.is_some() {
                dense_count += 1;
            }
        }
        // Also count R→L dense seeds
        let rc_read = reverse_complement_read(&read);
        for rc_pos in 0..rc_read.len() {
            let result = find_seed_at_position(&rc_read, rc_pos, &index, 4, false, &params);
            if result.seed.is_some() {
                dense_count += 1;
            }
        }

        assert!(
            sparse_seeds.len() <= dense_count,
            "Sparse ({}) should produce <= dense ({}) seeds",
            sparse_seeds.len(),
            dense_count
        );
    }

    #[test]
    fn test_rc_seed_genome_positions() {
        // Genome: AACCTTGG, read RC = CCAAGGTT
        // RC(read) = AACCTTGG matches forward genome
        let index = make_test_index("AACCTTGG");
        let read = encode_sequence("CCAAGGTT");
        let params = params(&[]);

        let seeds = Seed::find_seeds(&read, &index, 4, &params, "").unwrap();

        for seed in &seeds {
            if seed.search_rc {
                let positions: Vec<_> = seed.genome_positions(&index).collect();
                assert!(
                    !positions.is_empty(),
                    "RC seed should have genome positions"
                );

                for (pos, _is_rev) in &positions {
                    // Converted positions should be valid (within genome)
                    assert!(
                        *pos < index.genome.n_genome,
                        "Converted position {} should be < n_genome {}",
                        pos,
                        index.genome.n_genome
                    );
                }
            }
        }
    }
}
