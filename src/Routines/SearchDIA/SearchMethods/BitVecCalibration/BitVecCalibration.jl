#==========================================================
BitVecCalibration: Learn a 256-entry lookup table from fragment
index bitmask patterns to discriminate targets from decoys.

Runs between QuadTuning and MainSearch. Cross-file pooled:

process_file! (per file):
  1. Sample top-TIC MS2 scans
  2. Run fragment index in accumulator mode (PatternAccumulator)
     — tallies 256-bin target/decoy counts directly in the scoring
       loop with zero tuple allocation
  3. Merge per-thread counts into global accumulators on results

summarize_results! (once, after all files):
  4. Lock-then-merge: compute CI for all 256 bins, lock bins
     with clear verdicts (confident pass or fail), then greedily
     merge only undecided bins using (b1,k3,k5,k8) hierarchy.
     Recompute CI after each merge; lock when verdict is clear.
  5. High-signal sparse bins (e.g. Tc=100, Dc=2) lock as PASS
     immediately — never merged or contaminated
  6. Store single Vector{Bool}(256) filter for ALL files
==========================================================#

"""
    _ci_test(t, d, z, s) -> (excess_rate, se, lower_ci, upper_ci)

Compute library-corrected excess rate and confidence interval for a single bin.
`s` is the FDR scale factor (n_library_targets / n_library_decoys).
Pseudocount α=1 added to both Tc and Dc to avoid divide-by-zero and
shrink extreme rates from sparse bins toward zero.
"""
@inline function _ci_test(t::Int, d::Int, z::Float64, s::Float64=1.0)
    α = 1.0
    tf = Float64(t) + α
    df = (Float64(d) + α) * s            # expected decoys at target library scale
    n_var = tf + (Float64(d) + α) * s * s # Var(Tc - Dc*s) ≈ Tc + Dc*s²
    r = (tf - df) / df                    # corrected excess rate
    se = sqrt(n_var) / df
    return (r, se, r - z * se, r + z * se)
end

"""
    _adaptive_merge(tc, dc, z, min_excess_rate) -> NamedTuple

Lock-then-merge algorithm for adaptive bin resolution:

1. Compute CI test for all 256 bins
2. Lock bins with clear verdicts:
   - PASS: lower_CI ≥ min_excess_rate (confident signal, even sparse bins)
   - FAIL: upper_CI < min_excess_rate AND total ≥ N_min (confident noise)
3. Among remaining undecided bins, greedily merge the sparsest with its
   nearest undecided neighbor (hierarchy: (b1,k3,k5,k8) > (b1,k3,k5) > (b1,k3) > (b1) > any)
4. Recompute CI on merged bin; lock if verdict is now clear
5. Repeat until all bins are locked or no valid merges remain

Returns NamedTuple with fields:
- partition, tc_merged, dc_merged: bin groups and pooled counts
- passes: Bool per bin (true = pass)
- merge_levels: count of merges at each hierarchy level (1-2, within coarse group only)
- n_locked_early: bins locked before any merging
- n_min: computed N_min threshold
"""
function _adaptive_merge(tc::Vector{Int}, dc::Vector{Int}, z::Float64, min_r::Float64, fdr_scale::Float64=1.0)
    @assert length(tc) == 256 && length(dc) == 256
    n_min = ceil(Int, (2.0 * z / min_r)^2)

    # Pre-compute similarity keys for each pattern (0-indexed pattern value)
    pat_b1 = Vector{Int}(undef, 256)
    pat_k3 = Vector{Int}(undef, 256)
    pat_k5 = Vector{Int}(undef, 256)
    pat_k8 = Vector{Int}(undef, 256)
    for p in 0:255
        bits = ntuple(i -> (p >> (i-1)) & 1, 8)
        pat_b1[p+1] = bits[1]
        pat_k3[p+1] = sum(bits[i] for i in 1:3)
        pat_k5[p+1] = sum(bits[i] for i in 1:5)
        pat_k8[p+1] = count_ones(UInt8(p))
    end

    # Initialize bins
    BinType = NamedTuple{(:patterns, :tc, :dc), Tuple{Vector{Int}, Int, Int}}
    bins = Vector{Union{Nothing, BinType}}(undef, 256)
    for i in 1:256
        bins[i] = (patterns=[i], tc=tc[i], dc=dc[i])
    end

    # Lock status: :undecided, :pass, :fail
    status = fill(:undecided, 256)

    # ── Step 1: Lock bins that already have clear verdicts ──
    n_locked_early = 0
    for i in 1:256
        t, d = bins[i].tc, bins[i].dc
        # Patterns that can never be observed (k < min_count_ones) are locked as FAIL
        if t + d == 0
            status[i] = :fail
            n_locked_early += 1
            continue
        end
        _, _, lo, hi = _ci_test(t, d, z, fdr_scale)
        if lo >= min_r
            status[i] = :pass
            n_locked_early += 1
        elseif hi < min_r && (t + d) >= n_min
            status[i] = :fail
            n_locked_early += 1
        end
    end

    merge_levels = zeros(Int, 2)

    # ── Step 2: Greedy merge of undecided bins (within coarse group only) ──
    # Merges are restricted to the same (b1, k3, k5) coarse group (CG18).
    # Level 1: same (b1, k3, k5, k8) — differ only in which bits 6-8 are set
    # Level 2: same (b1, k3, k5) — differ in total popcount
    # No cross-group merges: if no undecided neighbor in the same group, lock as-is.
    while true
        # Find sparsest undecided bin
        min_total = typemax(Int)
        min_idx = 0
        for i in 1:length(bins)
            (isnothing(bins[i]) || status[i] != :undecided) && continue
            total = bins[i].tc + bins[i].dc
            if total < min_total
                min_total = total
                min_idx = i
            end
        end
        min_idx == 0 && break  # no undecided bins left

        # Find nearest UNDECIDED neighbor within same coarse group
        ref_pat = bins[min_idx].patterns[1]
        ref_b1 = pat_b1[ref_pat]
        ref_k3 = pat_k3[ref_pat]
        ref_k5 = pat_k5[ref_pat]
        ref_k8 = pat_k8[ref_pat]

        merged = false
        for level in 1:2
            best_idx = 0
            best_total = -1
            for j in 1:length(bins)
                (isnothing(bins[j]) || j == min_idx || status[j] != :undecided) && continue
                rep = bins[j].patterns[1]
                is_match = if level == 1
                    pat_b1[rep] == ref_b1 && pat_k3[rep] == ref_k3 &&
                    pat_k5[rep] == ref_k5 && pat_k8[rep] == ref_k8
                else
                    pat_b1[rep] == ref_b1 && pat_k3[rep] == ref_k3 && pat_k5[rep] == ref_k5
                end
                if is_match && (bins[j].tc + bins[j].dc) > best_total
                    best_total = bins[j].tc + bins[j].dc
                    best_idx = j
                end
            end
            if best_idx > 0
                new_pats = vcat(bins[best_idx].patterns, bins[min_idx].patterns)
                new_tc = bins[best_idx].tc + bins[min_idx].tc
                new_dc = bins[best_idx].dc + bins[min_idx].dc
                bins[best_idx] = (patterns=new_pats, tc=new_tc, dc=new_dc)
                bins[min_idx] = nothing
                status[min_idx] = :fail  # consumed
                merge_levels[level] += 1
                merged = true

                # Recompute CI on merged bin — lock if verdict is now clear
                _, _, lo, hi = _ci_test(new_tc, new_dc, z, fdr_scale)
                if lo >= min_r
                    status[best_idx] = :pass
                elseif hi < min_r && (new_tc + new_dc) >= n_min
                    status[best_idx] = :fail
                end
                break
            end
        end

        if !merged
            # No undecided neighbor — lock with current verdict
            t, d = bins[min_idx].tc, bins[min_idx].dc
            _, _, _, hi = _ci_test(t, d, z, fdr_scale)
            status[min_idx] = hi >= min_r ? :pass : :fail
        end
    end

    # ── Step 3: Lock any remaining undecided bins ──
    for i in 1:length(bins)
        (isnothing(bins[i]) || status[i] != :undecided) && continue
        t, d = bins[i].tc, bins[i].dc
        _, _, _, hi = _ci_test(t, d, z, fdr_scale)
        status[i] = hi >= min_r ? :pass : :fail
    end

    # Compact into output arrays
    partition = Vector{Vector{Int}}()
    tc_merged = Vector{Int}()
    dc_merged = Vector{Int}()
    passes = Vector{Bool}()
    for i in 1:length(bins)
        isnothing(bins[i]) && continue
        push!(partition, bins[i].patterns)
        push!(tc_merged, bins[i].tc)
        push!(dc_merged, bins[i].dc)
        push!(passes, status[i] == :pass)
    end

    return (partition=partition, tc_merged=tc_merged, dc_merged=dc_merged,
            passes=passes, merge_levels=merge_levels, n_locked_early=n_locked_early,
            n_min=n_min)
end

# Format a pattern as 8-bit string, rank 1 (top fragment) on left: "11010010"
function _pat_str(pat_idx::Int)
    p = UInt8(pat_idx - 1)
    return String([((p >> i) & 1) == 1 ? '1' : '0' for i in 0:7])
end

struct BitVecCalibrationSearch <: SearchMethod end

struct BitVecCalibrationParameters <: FragmentIndexSearchParameters
    isotope_err_bounds::Tuple{UInt8, UInt8}
    min_index_search_score::UInt8
    spec_order::Set{Int64}
    min_excess_rate::Float32
    diagnostic_stratify::Bool   # DIAGNOSTIC ONLY: stratify pattern counts by scan peak count
    n_stratify_bins::Int        # number of peak-count quantile bins for the diagnostic
    diagnostic_all_scans::Bool  # DIAGNOSTIC: stratify over ALL MS2 scans, not just the
                                # TIC-priority adaptive sample (removes the dense-scan bias)
    diagnostic_strata_target::Int  # DIAGNOSTIC: if >0, RANDOM-sample scans within each
                                   # peak-count stratum until it reaches this many counts
                                   # (instead of TIC-priority), so sparse strata are sampled

    function BitVecCalibrationParameters(params::PioneerParameters)
        # Per-bitvec-pattern minimum target-excess rate for admission. Lower =>
        # more candidate (scan, precursor) cells admitted to scoring (more IDs,
        # more runtime). Surfaced as a config knob
        # (optimization.bitvec_calibration.min_excess_rate) so different
        # thresholds can be swept in one REPL session without recompiling;
        # defaults to BITVEC_MIN_EXCESS_RATE (0.03) when absent.
        bitvec_cfg = get(params.optimization, :bitvec_calibration, (;))
        min_excess_rate = Float32(get(bitvec_cfg, :min_excess_rate, BITVEC_MIN_EXCESS_RATE))
        # DIAGNOSTIC: when on, additionally accumulate per-pattern target/decoy
        # counts stratified into peak-count quantile bins and write
        # bitvec_stratified_counts.arrow. Does NOT change the live filter.
        diagnostic_stratify = Bool(get(bitvec_cfg, :diagnostic_stratify, false))
        n_stratify_bins = Int(get(bitvec_cfg, :n_stratify_bins, 10))
        diagnostic_all_scans = Bool(get(bitvec_cfg, :diagnostic_all_scans, false))
        diagnostic_strata_target = Int(get(bitvec_cfg, :diagnostic_strata_target, 0))
        new(
            (UInt8(1), UInt8(0)),  # isotope_err_bounds
            UInt8(3),              # min_score for CountFilter fallback (normal search path)
            Set{Int64}([2]),       # MS2 only
            min_excess_rate,
            diagnostic_stratify,
            n_stratify_bins,
            diagnostic_all_scans,
            diagnostic_strata_target,
        )
    end
end

# z-score used only for the diagnostic Arrow files (lower_ci/upper_ci columns
# in bitvec_model_global.arrow). Live filter decisions don't use it.
# BITVEC_MIN_EXCESS_RATE (defined a few lines below) was formerly surfaced as
# fragment_index_search.bitvec_calibration.min_excess_rate but no shipping
# config ever overrode the constant.
const BITVEC_DIAGNOSTIC_Z = Float32(1.28)

struct BitVecCalibrationResults <: SearchResults
    global_tc::Vector{Int}  # length 256, accumulated across all files
    global_dc::Vector{Int}  # length 256, accumulated across all files
    # DIAGNOSTIC (sized 256×B, B=0 when off): per-pattern target/decoy counts
    # stratified into B peak-count quantile bins, summed across files (bin = the
    # b-th per-file crowding decile). strat_peak_sum/strat_nscans give the
    # count-weighted mean peak count per bin for the report's x-axis.
    strat_tc::Matrix{Int}
    strat_dc::Matrix{Int}
    strat_peak_sum::Vector{Float64}
    strat_nscans::Vector{Int}
end

#==========================================================
Interface
==========================================================#

get_parameters(::BitVecCalibrationSearch, params::Any) = BitVecCalibrationParameters(params)

function init_search_results(::BitVecCalibrationSearch, params::BitVecCalibrationParameters, ::SearchContext)
    B = params.diagnostic_stratify ? params.n_stratify_bins : 0
    return BitVecCalibrationResults(
        zeros(Int, 256), zeros(Int, 256),
        zeros(Int, 256, B), zeros(Int, 256, B), zeros(Float64, B), zeros(Int, B),
    )
end

reset_results!(::BitVecCalibrationResults) = nothing

const BITVEC_MIN_TOTAL_COUNTS = Int64(10_000_000)
const BITVEC_INITIAL_SCANS = Int64(2000)
const BITVEC_COARSE_GROUP_THRESHOLD = Int64(1_000_000)
const BITVEC_MIN_EXCESS_RATE = Float32(0.03)

function _bitvec_excess_rank_table(
    target_counts::AbstractVector{<:Integer},
    decoy_counts::AbstractVector{<:Integer},
    fdr_scale::Real = 1.0,
)
    length(target_counts) == 256 || throw(ArgumentError("target count vector must have length 256"))
    length(decoy_counts) == 256 || throw(ArgumentError("decoy count vector must have length 256"))

    α = 1.0
    s = Float64(fdr_scale)
    rates = Vector{Float64}(undef, 256)
    excess = Vector{Float64}(undef, 256)
    @inbounds for i in 1:256
        target = Float64(target_counts[i]) + α
        decoy = s * (Float64(decoy_counts[i]) + α)
        ex = target - decoy
        excess[i] = ex
        rates[i] = ex / decoy
    end

    order = sortperm(1:256, by = i -> (-rates[i], -excess[i], i))
    ranks = Vector{UInt16}(undef, 256)
    @inbounds for (rank, i) in pairs(order)
        ranks[i] = UInt16(rank)
    end
    return ranks
end

function process_file!(
    results::BitVecCalibrationResults,
    params::BitVecCalibrationParameters,
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::MassSpecData
)
    t_start = time()
    file_name = getParsedFileName(search_context, ms_file_idx)

    spec_lib = getSpecLib(search_context)
    precursors = getPrecursors(spec_lib)
    is_decoy = getIsDecoy(precursors)
    partitioned_index = get_fragment_index(spec_lib, params)
    qtm = getQuadTransmissionModel(search_context, ms_file_idx)
    mem = getMassErrorModel(search_context, ms_file_idx)
    rt_to_irt = getRtIrtModel(search_context, ms_file_idx)
    irt_tol = get_irt_tolerance(search_context, params, ms_file_idx)
    prec_mzs = getMz(precursors)
    prec_irts = getIrt(precursors)

    scan_priority = get_ms2_scan_priority_order(spectra)
    total_ms2 = length(scan_priority)
    scan_to_prec_idx = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra))

    # Adaptive accumulation: collect scans until total counts ≥ 1M or all scans used
    file_tc = zeros(Int, 256)
    file_dc = zeros(Int, 256)
    prev = 0
    scan_target = min(BITVEC_INITIAL_SCANS, total_ms2)
    total_counts = 0

    while prev < total_ms2
        batch_scans = Int.(scan_priority[(prev+1):scan_target])
        acc = PatternAccumulator(Threads.nthreads(), is_decoy; min_score=UInt8(2))

        searchFragmentIndexPartitionMajorHinted(
            scan_to_prec_idx, partitioned_index, spectra, batch_scans,
            Threads.nthreads(), params, qtm, mem, rt_to_irt, irt_tol, prec_mzs;
            pattern_accumulator=acc, precursor_irts=prec_irts,
            scratch=getFragIndexScratch(search_context))

        batch_tc, batch_dc = merge_accumulator(acc)
        file_tc .+= batch_tc
        file_dc .+= batch_dc
        prev = scan_target
        total_counts = sum(file_tc) + sum(file_dc)

        total_counts >= BITVEC_MIN_TOTAL_COUNTS && break

        # Estimate additional scans from count rate
        scans_remaining = total_ms2 - prev
        scans_remaining <= 0 && break
        rate = total_counts / max(prev, 1)
        if rate > 0
            additional = clamp(ceil(Int, (BITVEC_MIN_TOTAL_COUNTS - total_counts) / rate * 1.5),
                               1, scans_remaining)
        else
            additional = min(prev, scans_remaining)
        end
        scan_target = min(prev + additional, total_ms2)
        scan_target <= prev && break
    end

    # Also accumulate into global counts for diagnostics
    results.global_tc .+= file_tc
    results.global_dc .+= file_dc

    # Compute per-file filter
    fdr_scale = Float64(getLibraryFdrScaleFactor(search_context))
    min_excess_rate = Float64(params.min_excess_rate)
    α = 1.0  # pseudocount

    filter_table = Vector{Bool}(undef, 256)
    setBitVecExcessRanks!(search_context, ms_file_idx,
        _bitvec_excess_rank_table(file_tc, file_dc, fdr_scale))

    if total_counts >= BITVEC_COARSE_GROUP_THRESHOLD
        # Sufficient data: per-bin excess rate decision
        for i in 1:256
            tc_i = Float64(file_tc[i]) + α
            dc_i = fdr_scale * (Float64(file_dc[i]) + α)
            filter_table[i] = (tc_i - dc_i) / dc_i > min_excess_rate
        end
        mode = "per-bin"
    else
        # Sparse data: coarse group (b1, k3, k5) — 18 groups
        # All patterns in a group pass or fail together
        group_key = Vector{NTuple{3, Int}}(undef, 256)
        for p in 0:255
            b1 = (p >> 0) & 1
            k3 = count_ones(UInt8(p) & 0x07)
            k5 = count_ones(UInt8(p) & 0x1f)
            group_key[p + 1] = (b1, k3, k5)
        end

        group_tc = Dict{NTuple{3, Int}, Float64}()
        group_dc = Dict{NTuple{3, Int}, Float64}()
        for p in 0:255
            key = group_key[p + 1]
            group_tc[key] = get(group_tc, key, 0.0) + Float64(file_tc[p + 1])
            group_dc[key] = get(group_dc, key, 0.0) + Float64(file_dc[p + 1])
        end

        group_pass = Dict{NTuple{3, Int}, Bool}()
        for key in keys(group_tc)
            tc_g = group_tc[key] + α
            dc_g = fdr_scale * (group_dc[key] + α)
            group_pass[key] = (tc_g - dc_g) / dc_g > min_excess_rate
        end
        for p in 0:255
            filter_table[p + 1] = group_pass[group_key[p + 1]]
        end
        mode = "coarse(b1,k3,k5)"
    end

    setBitVecFilter!(search_context, ms_file_idx, filter_table)
    n_pass = count(filter_table)
    elapsed = round(time() - t_start, digits=2)
    @debug_l1 "BitVec: $file_name — $(prev) scans, $(total_counts) counts, $(n_pass)/256 pass ($mode, min_excess=$(round(min_excess_rate*100, digits=1))%, $(elapsed)s)"

    # DIAGNOSTIC: re-accumulate the same sampled scans, stratified by peak count.
    # Live filter above is untouched.
    if params.diagnostic_stratify
        if params.diagnostic_strata_target > 0
            # Random-sample within each peak-count stratum until the count target.
            accumulate_stratified_random!(results, params, search_context, ms_file_idx, spectra,
                Int.(scan_priority), scan_to_prec_idx, partitioned_index,
                qtm, mem, rt_to_irt, irt_tol, prec_mzs, prec_irts, is_decoy)
        else
            # All MS2 scans (unbiased) vs the TIC-priority sample the live filter used.
            diag_scans = params.diagnostic_all_scans ? Int.(scan_priority) : Int.(scan_priority[1:prev])
            accumulate_stratified!(results, params, search_context, ms_file_idx, spectra,
                diag_scans, scan_to_prec_idx, partitioned_index,
                qtm, mem, rt_to_irt, irt_tol, prec_mzs, prec_irts, is_decoy)
        end
    end
end

"""
    accumulate_stratified!(results, ...)

DIAGNOSTIC ONLY. Bins the file's sampled MS2 scans into `n_stratify_bins`
peak-count quantile bins, re-runs the fragment-index search once per bin (each
scan processed exactly once, in its bin), and folds the per-(pattern, bin)
target/decoy counts into the global stratified table. Reuses the existing
`PatternAccumulator` + search — no change to the live gate or the hot loop.
"""
function accumulate_stratified!(
    results::BitVecCalibrationResults,
    params::BitVecCalibrationParameters,
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::MassSpecData,
    sampled_scans::Vector{Int},
    scan_to_prec_idx::Vector{Union{Missing, UnitRange{Int64}}},
    partitioned_index,
    qtm, mem, rt_to_irt, irt_tol,
    prec_mzs::AbstractVector{Float32},
    prec_irts::AbstractVector{Float32},
    is_decoy::AbstractVector{Bool},
)
    B = params.n_stratify_bins
    n = length(sampled_scans)
    n < B && return nothing
    pc = Float64[length(getMzArray(spectra, s)) for s in sampled_scans]
    # Nearest-rank quantile edges (pure Base; avoids a Statistics dependency).
    sp = sort(pc); np = length(sp)
    edges = Float64[sp[clamp(floor(Int, q * (np - 1)) + 1, 1, np)] for q in range(0.0, 1.0; length = B + 1)]
    binof(x) = clamp(searchsortedlast(edges, x), 1, B)
    nthreads = Threads.nthreads()
    for b in 1:B
        binscans = Int[sampled_scans[i] for i in 1:n if binof(pc[i]) == b]
        isempty(binscans) && continue
        acc = PatternAccumulator(nthreads, is_decoy; min_score = UInt8(2))
        searchFragmentIndexPartitionMajorHinted(
            scan_to_prec_idx, partitioned_index, spectra, binscans, nthreads,
            params, qtm, mem, rt_to_irt, irt_tol, prec_mzs;
            pattern_accumulator = acc, precursor_irts = prec_irts,
            scratch = getFragIndexScratch(search_context))
        tc_b, dc_b = merge_accumulator(acc)
        @inbounds for p in 1:256
            results.strat_tc[p, b] += tc_b[p]
            results.strat_dc[p, b] += dc_b[p]
        end
        results.strat_nscans[b] += length(binscans)
        results.strat_peak_sum[b] += sum(pc[i] for i in 1:n if binof(pc[i]) == b)
    end
    return nothing
end

"""
    accumulate_stratified_random!(results, ...)

DIAGNOSTIC ONLY. Bins ALL MS2 scans into `n_stratify_bins` peak-count quantile
strata, then **randomly** samples scans within each stratum (not TIC-priority,
which skips sparse spectra) until that stratum reaches `diagnostic_strata_target`
counts (or its scans are exhausted). Random subsampling preserves the per-pattern
excess-rate estimate; the targets just guarantee adequate counts in every
stratum — including the sparse ones the live sampler under-samples.
"""
function accumulate_stratified_random!(
    results::BitVecCalibrationResults,
    params::BitVecCalibrationParameters,
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::MassSpecData,
    all_scans::Vector{Int},
    scan_to_prec_idx::Vector{Union{Missing, UnitRange{Int64}}},
    partitioned_index,
    qtm, mem, rt_to_irt, irt_tol,
    prec_mzs::AbstractVector{Float32},
    prec_irts::AbstractVector{Float32},
    is_decoy::AbstractVector{Bool},
)
    B = params.n_stratify_bins
    target = params.diagnostic_strata_target
    n = length(all_scans)
    n < B && return nothing
    pc = Float64[length(getMzArray(spectra, s)) for s in all_scans]
    sp = sort(pc); np = length(sp)
    edges = Float64[sp[clamp(floor(Int, q * (np - 1)) + 1, 1, np)] for q in range(0.0, 1.0; length = B + 1)]
    binof(x) = clamp(searchsortedlast(edges, x), 1, B)
    bin_scans = [Int[] for _ in 1:B]
    bin_pc    = [Float64[] for _ in 1:B]
    @inbounds for i in 1:n
        b = binof(pc[i]); push!(bin_scans[b], all_scans[i]); push!(bin_pc[b], pc[i])
    end
    rng = Random.MersenneTwister(1844 + ms_file_idx)
    nthreads = Threads.nthreads()
    for b in 1:B
        m = length(bin_scans[b]); m == 0 && continue
        order = Random.shuffle(rng, collect(1:m))
        acc = PatternAccumulator(nthreads, is_decoy; min_score = UInt8(2))
        processed = 0; pcsum = 0.0; chunk = 2000
        while processed < m
            idxs = order[(processed + 1):min(processed + chunk, m)]
            batch = Int[bin_scans[b][j] for j in idxs]
            searchFragmentIndexPartitionMajorHinted(
                scan_to_prec_idx, partitioned_index, spectra, batch, nthreads,
                params, qtm, mem, rt_to_irt, irt_tol, prec_mzs;
                pattern_accumulator = acc, precursor_irts = prec_irts,
                scratch = getFragIndexScratch(search_context))
            pcsum += sum(bin_pc[b][j] for j in idxs)
            processed += length(idxs)
            tc, dc = merge_accumulator(acc)
            (sum(tc) + sum(dc)) >= target && break
        end
        tc, dc = merge_accumulator(acc)
        @inbounds for p in 1:256
            results.strat_tc[p, b] += tc[p]; results.strat_dc[p, b] += dc[p]
        end
        results.strat_nscans[b] += processed
        results.strat_peak_sum[b] += pcsum
    end
    return nothing
end

function process_search_results!(
    ::BitVecCalibrationResults,
    ::BitVecCalibrationParameters,
    ::SearchContext,
    ::Int64,
    ::MassSpecData
)
    return nothing
end

function summarize_results!(
    results::BitVecCalibrationResults,
    params::BitVecCalibrationParameters,
    search_context::SearchContext
)
    n_files = length(getMassSpecData(search_context).file_paths)

    # Per-file filters already installed in process_file!.
    # Write global diagnostic Arrow files from accumulated counts.
    tc = results.global_tc
    dc = results.global_dc
    Nt = sum(tc); Nd = sum(dc)
    if Nd == 0 || Nt == 0
        return
    end

    z = Float64(BITVEC_DIAGNOSTIC_Z)
    min_r = Float64(params.min_excess_rate)
    fdr_scale = Float64(getLibraryFdrScaleFactor(search_context))
    result = _adaptive_merge(tc, dc, z, min_r, fdr_scale)
    partition = result.partition
    tc_m = result.tc_merged
    dc_m = result.dc_merged
    group_passes = result.passes
    n_bins = length(partition)

    excess_rate = zeros(Float64, n_bins)
    se_rate = zeros(Float64, n_bins)
    upper_ci = zeros(Float64, n_bins)
    lower_ci = zeros(Float64, n_bins)
    for b in 1:n_bins
        r, se, lo, hi = _ci_test(tc_m[b], dc_m[b], z, fdr_scale)
        excess_rate[b] = r; se_rate[b] = se; lower_ci[b] = lo; upper_ci[b] = hi
    end

    # ── 4. Write Arrow diagnostics ──
    out_dir = getDataOutDir(search_context)

    # Raw 256-bin counts
    Arrow.write(joinpath(out_dir, "bitvec_counts_256.arrow"), DataFrame(
        pattern = collect(UInt8(0):UInt8(255)),
        target_count = tc,
        decoy_count = dc,
        count_ones = UInt8[count_ones(UInt8(p)) for p in 0:255],
    ))

    # Per-merged-bin model
    bin_ids = Vector{Int}()
    bin_pat_strs = Vector{String}()
    bin_tc = Vector{Int}()
    bin_dc = Vector{Int}()
    bin_excess = Vector{Int}()
    bin_rate = Vector{Float32}()
    bin_lo = Vector{Float32}()
    bin_hi = Vector{Float32}()
    bin_pass = Vector{Bool}()
    bin_n_pats = Vector{Int}()

    for b in 1:n_bins
        pats = partition[b]
        push!(bin_ids, b)
        push!(bin_pat_strs, join([_pat_str(p) for p in pats], ","))
        push!(bin_tc, tc_m[b])
        push!(bin_dc, dc_m[b])
        push!(bin_excess, tc_m[b] - dc_m[b])
        push!(bin_rate, Float32(excess_rate[b]))
        push!(bin_lo, Float32(lower_ci[b]))
        push!(bin_hi, Float32(upper_ci[b]))
        push!(bin_pass, group_passes[b])
        push!(bin_n_pats, length(pats))
    end

    Arrow.write(joinpath(out_dir, "bitvec_model_global.arrow"), DataFrame(
        bin = bin_ids,
        n_patterns = bin_n_pats,
        patterns = bin_pat_strs,
        target_count = bin_tc,
        decoy_count = bin_dc,
        excess_count = bin_excess,
        excess_rate = bin_rate,
        lower_ci = bin_lo,
        upper_ci = bin_hi,
        passes = bin_pass,
    ))

    # ── 5. Diagnostic log (global pooled counts, for reference only) ──
    n_passing_bins = count(group_passes)
    @debug_l1 "BitVecCalibration global: $n_files files, $(Nt + Nd) total counts ($Nt target, $Nd decoy)\n"

    bin_order = sortperm(1:n_bins, by=b -> (-group_passes[b], -excess_rate[b]))
    for b in bin_order
        (tc_m[b] + dc_m[b] == 0) && continue
        pats = partition[b]
        n_p = length(pats)
        r_str = @sprintf("%.4f", excess_rate[b])
        lo_str = @sprintf("%.4f", lower_ci[b])
        hi_str = @sprintf("%.4f", upper_ci[b])
        pass_str = group_passes[b] ? "  ✓ " : "  ✗ "
        @debug_l2 "│ $(lpad(b, 3)) │ $(lpad(n_p, 4)) │ $(lpad(tc_m[b], 12)) │ $(lpad(dc_m[b], 12)) │ $(rpad(r_str, 8)) │ $(lo_str) … $(hi_str) │$(pass_str)│\n"

        # Show individual patterns in this bin, sorted by total count (most populated first)
        sorted_pats = sort(pats, by=p -> tc[p] + dc[p], rev=true)
        n_show = n_p <= 8 ? n_p : 6
        for p in sorted_pats[1:n_show]
            pat_val = UInt8(p - 1)
            @debug_l2 "│     │      │  $(_pat_str(p)) k=$(count_ones(pat_val)) │ $(lpad(tc[p], 8)) tc │ $(lpad(dc[p], 8)) dc │                      │      │\n"
        end
        if n_p > 8
            @debug_l2 "│     │      │  ... $(n_p - n_show) more patterns                                                │      │\n"
        end
    end
    @debug_l2 "└─────┴──────┴──────────────┴──────────────┴──────────┴──────────────────────┴──────┘\n"

    # Per-count_ones summary (using global counts + simple excess rate)
    α = 1.0
    @debug_l1 "\ncount_ones summary (global):\n"
    for k in 0:8
        mask_pats = [p for p in 0:255 if count_ones(UInt8(p)) == k]
        t_k = sum(tc[p+1] for p in mask_pats)
        d_k = sum(dc[p+1] for p in mask_pats)
        n_pass_k = count(((Float64(tc[p+1]) + α) - fdr_scale * (Float64(dc[p+1]) + α)) /
                         (fdr_scale * (Float64(dc[p+1]) + α)) > min_r for p in mask_pats)
        n_pats_k = length(mask_pats)
        (t_k + d_k > 0) && @debug_l1 "  k=$k: $(n_pats_k) patterns, $(t_k) target, $(d_k) decoy, $(n_pass_k)/$(n_pats_k) pass\n"
    end

    # ── 6. DIAGNOSTIC: stratified (pattern × peak-count bin) counts ──
    if params.diagnostic_stratify && size(results.strat_tc, 2) > 0
        B = size(results.strat_tc, 2)
        mean_pc = [results.strat_nscans[b] > 0 ? results.strat_peak_sum[b]/results.strat_nscans[b] : 0.0 for b in 1:B]
        αs = 1.0
        rows_pat=Int[]; rows_bin=Int[]; rows_k=Int[]; rows_mpc=Float64[]
        rows_tc=Int[]; rows_dc=Int[]; rows_excess=Float64[]; rows_rate=Float64[]; rows_pass=Bool[]; rows_ns=Int[]
        for p in 1:256, b in 1:B
            t=Float64(results.strat_tc[p,b]); d=fdr_scale*Float64(results.strat_dc[p,b])
            rate=((t+αs)-(d+αs))/(d+αs)
            push!(rows_pat,p-1); push!(rows_bin,b); push!(rows_k,count_ones(UInt8(p-1))); push!(rows_mpc,mean_pc[b])
            push!(rows_tc,results.strat_tc[p,b]); push!(rows_dc,results.strat_dc[p,b])
            push!(rows_excess, t - fdr_scale*Float64(results.strat_dc[p,b]))
            push!(rows_rate, rate); push!(rows_pass, rate > min_r); push!(rows_ns, results.strat_nscans[b])
        end
        Arrow.write(joinpath(out_dir, "bitvec_stratified_counts.arrow"), DataFrame(
            pattern=rows_pat, bin=rows_bin, count_ones=rows_k, mean_peak_count=rows_mpc,
            target_count=rows_tc, decoy_count=rows_dc, excess_count=rows_excess,
            excess_rate=rows_rate, would_pass=rows_pass, n_scans_bin=rows_ns))
        @user_info "BitVec stratified diagnostic: wrote bitvec_stratified_counts.arrow " *
                   "($B bins, mean peak counts $(round.(mean_pc, digits=0)))"
    end
end
