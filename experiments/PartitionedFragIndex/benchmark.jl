#!/usr/bin/env julia
#
# Benchmark: Partitioned Fragment Index vs. NativeFragmentIndex
#
# Usage:
#   julia --threads=auto --project=. experiments/PartitionedFragIndex/benchmark.jl /path/to/params.json
#

using Pioneer
using Random: shuffle!

# ── Load experiment modules ──────────────────────────────────────────────────
include(joinpath(@__DIR__, "partitioned_types.jl"))
include(joinpath(@__DIR__, "build_partitioned_index.jl"))
include(joinpath(@__DIR__, "search_partitioned_index.jl"))

# ── Parse CLI argument ──────────────────────────────────────────────────────
if isempty(ARGS)
    error("Usage: julia --threads=auto --project=. experiments/PartitionedFragIndex/benchmark.jl <params.json>")
end
const PARAMS_PATH = ARGS[1]
isfile(PARAMS_PATH) || error("Params file not found: $PARAMS_PATH")

# ── Load parameters ─────────────────────────────────────────────────────────
println("Loading parameters from $PARAMS_PATH ...")
Pioneer.checkParams(PARAMS_PATH)
params = Pioneer.parse_pioneer_parameters(PARAMS_PATH)
mkpath(params.paths[:results])

# ── Resolve data paths ──────────────────────────────────────────────────────
MS_DATA_DIR = params.paths[:ms_data]
if !isabspath(MS_DATA_DIR)
    MS_DATA_DIR = joinpath(@__DIR__, "../../", MS_DATA_DIR)
end
isdir(MS_DATA_DIR) || error("ms_data directory does not exist: $MS_DATA_DIR")

MS_TABLE_PATHS = [joinpath(MS_DATA_DIR, f)
                  for f in readdir(MS_DATA_DIR)
                  if isfile(joinpath(MS_DATA_DIR, f)) && endswith(f, ".arrow")]
isempty(MS_TABLE_PATHS) && error("No .arrow files found in: $MS_DATA_DIR")

SPEC_LIB_DIR = params.paths[:library]
if !isabspath(SPEC_LIB_DIR)
    SPEC_LIB_DIR = joinpath(@__DIR__, "../../", SPEC_LIB_DIR)
end

# ── Load spectral library ──────────────────────────────────────────────────
println("Loading spectral library ...")
SPEC_LIB = Pioneer.loadSpectralLibrary(SPEC_LIB_DIR, params)

# ── Initialize search context ──────────────────────────────────────────────
println("Initializing search context ...")
SEARCH_CONTEXT = Pioneer.initSearchContext(
    SPEC_LIB,
    Pioneer.parseIsoXML(Pioneer.isotope_spline_path()),
    Pioneer.ArrowTableReference(MS_TABLE_PATHS),
    Threads.nthreads(),
    250000
)
Pioneer.setDataOutDir!(SEARCH_CONTEXT, params.paths[:results])

# ── Run calibration ─────────────────────────────────────────────────────────
println("Running ParameterTuningSearch ...")
Pioneer.execute_search(Pioneer.ParameterTuningSearch(), SEARCH_CONTEXT, params)

println("Running NceTuningSearch ...")
Pioneer.execute_search(Pioneer.NceTuningSearch(), SEARCH_CONTEXT, params)

println("Running QuadTuningSearch ...")
Pioneer.execute_search(Pioneer.QuadTuningSearch(), SEARCH_CONTEXT, params)

# ── Prepare search arguments ───────────────────────────────────────────────
println("\nPreparing search arguments ...")
search_parameters = Pioneer.FirstPassSearchParameters(params)
fragment_index     = Pioneer.getFragmentIndex(Pioneer.getSpecLib(SEARCH_CONTEXT))
search_data        = Pioneer.getSearchData(SEARCH_CONTEXT)

msdr = Pioneer.getMassSpecData(SEARCH_CONTEXT)
ms_file_idx, spectra = first(enumerate(msdr))

qtm              = Pioneer.getQuadTransmissionModel(SEARCH_CONTEXT, ms_file_idx)
mem              = Pioneer.getMassErrorModel(SEARCH_CONTEXT, ms_file_idx)
rt_to_irt_spline = Pioneer.getRtIrtModel(SEARCH_CONTEXT, ms_file_idx)
irt_tol          = Pioneer.getIrtErrors(SEARCH_CONTEXT)[ms_file_idx]

thread_tasks = Pioneer.partition_scans(spectra, Threads.nthreads())

# Precursor m/z vector for post-filter
precursor_mzs = Pioneer.getMz(Pioneer.getPrecursors(Pioneer.getSpecLib(SEARCH_CONTEXT)))

# ── Materialize native index ───────────────────────────────────────────────
println("Materializing fragment index ...")
materialize_time = @elapsed begin
    native_frag_index = Pioneer.materialize(fragment_index)
end
n_frags = length(native_frag_index.fragments)
n_frag_bins = length(native_frag_index.fragment_bins)
n_rt_bins = length(native_frag_index.rt_bins)
native_mem_mb = (sizeof(native_frag_index.fragment_bins) +
                 sizeof(native_frag_index.rt_bins) +
                 sizeof(native_frag_index.fragments)) / 1024^2
println("  Materialization: $(round(materialize_time, digits=3))s")
println("  fragment_bins: $n_frag_bins  rt_bins: $n_rt_bins  fragments: $n_frags")
println("  Native memory: $(round(native_mem_mb, digits=1)) MB")

# ── Build partitioned index from library ──────────────────────────────────
println("\nBuilding partitioned index from library ...")
build_time = @elapsed begin
    partitioned_index = build_partitioned_index_from_lib(
        Pioneer.getSpecLib(SEARCH_CONTEXT);
        partition_width=5.0f0,
        frag_bin_tol_ppm=10.0f0,
        rt_bin_tol=1.0f0,
    )
end
total_part_frags = sum(length(Pioneer.getFragments(p)) for p in partitioned_index.partitions)
total_part_fbins = sum(length(Pioneer.getFragBins(p)) for p in partitioned_index.partitions)
part_mem_mb = sum(
    sizeof(Pioneer.getFragBins(p)) + sizeof(Pioneer.getRTBins(p)) + sizeof(Pioneer.getFragments(p))
    for p in partitioned_index.partitions
) / 1024^2
println("  Build time: $(round(build_time, digits=3))s")
println("  Partitions: $(partitioned_index.n_partitions)")
println("  Total fragments across partitions: $total_part_frags (original: $n_frags)")
println("  Total frag_bins across partitions: $total_part_fbins (original: $n_frag_bins)")
println("  Partitioned memory: $(round(part_mem_mb, digits=1)) MB")

# ── Correctness validation (score-level) ──────────────────────────────────
println("\n" * "="^80)
println("CORRECTNESS VALIDATION — score-level comparison")
println("="^80)

# Sample 500 MS2 scans stratified by peak count (100 per quintile)
let
    # Collect all MS2 scan indices and their peak counts
    candidates = Tuple{Int, Int}[]  # (scan_idx, n_peaks)
    for (_, task_scans) in thread_tasks
        for s in task_scans
            (s <= 0 || s > length(spectra)) && continue
            Pioneer.getMsOrder(spectra, s) ∉ Pioneer.getSpecOrder(search_parameters) && continue
            push!(candidates, (s, length(Pioneer.getMzArray(spectra, s))))
        end
    end
    sort!(candidates, by=last)

    n_quantiles = 5
    per_quantile = 100
    chunk_size = max(1, length(candidates) ÷ n_quantiles)
    sampled = Int[]
    for q in 1:n_quantiles
        lo = (q - 1) * chunk_size + 1
        hi = q == n_quantiles ? length(candidates) : q * chunk_size
        pool = [c[1] for c in candidates[lo:hi]]
        shuffle!(pool)
        append!(sampled, pool[1:min(per_quantile, length(pool))])
    end
    sort!(sampled)
    global all_ms2_scans = sampled

    peak_counts = [length(Pioneer.getMzArray(spectra, s)) for s in sampled]
    println("Validating on $(length(sampled)) MS2 scans (stratified by peak count)")
    println("  Peak count range: $(minimum(peak_counts)) – $(maximum(peak_counts)), median: $(peak_counts[length(peak_counts)÷2])")
end

# Helper: snapshot counter scores after searchScan! / searchScanPartitioned! but before
# filterPrecursorMatches! and reset!.  Returns Dict{UInt32, UInt8} of {prec_id => score}.
function snapshot_counter(counter::Pioneer.Counter{UInt32, UInt8})
    scores = Dict{UInt32, UInt8}()
    @inbounds for i in 1:(Pioneer.getSize(counter) - 1)
        pid = Pioneer.getID(counter, i)
        pid == 0 && continue
        sc = counter.counts[pid]
        sc == 0 && continue
        scores[pid] = sc
    end
    return scores
end

# Run scan-by-scan, capturing raw scores before filter/reset
validation_task_id = first(first(thread_tasks))
counter = Pioneer.getPrecursorScores(search_data[validation_task_id])

# Baseline: mirrors searchFragmentIndex scan loop
baseline_scores = Dict{Int, Dict{UInt32, UInt8}}()
let rt_bin_idx = 1
    rt_bins_ref = Pioneer.getRTBins(native_frag_index)
    frag_bins_ref = Pioneer.getFragBins(native_frag_index)
    fragments_ref = Pioneer.getFragments(native_frag_index)
    for scan_idx in all_ms2_scans
        irt_lo, irt_hi = Pioneer.getRTWindow(rt_to_irt_spline(Pioneer.getRetentionTime(spectra, scan_idx)), irt_tol)
        while rt_bin_idx < length(rt_bins_ref) && Pioneer.getHigh(rt_bins_ref[rt_bin_idx]) < irt_lo
            rt_bin_idx += 1
        end
        while rt_bin_idx > 1 && Pioneer.getLow(rt_bins_ref[rt_bin_idx]) > irt_lo
            rt_bin_idx -= 1
        end

        Pioneer.searchScan!(
            counter, rt_bins_ref, frag_bins_ref, fragments_ref,
            Pioneer.getMzArray(spectra, scan_idx),
            Pioneer.getIntensityArray(spectra, scan_idx),
            rt_bin_idx, irt_hi, mem,
            Pioneer.getQuadTransmissionFunction(qtm, Pioneer.getCenterMz(spectra, scan_idx), Pioneer.getIsolationWidthMz(spectra, scan_idx)),
            Pioneer.getIsotopeErrBounds(search_parameters)
        )

        baseline_scores[scan_idx] = snapshot_counter(counter)
        Pioneer.reset!(counter)
    end
end

# Partitioned: mirrors searchFragmentIndexPartitioned scan loop
partitioned_scores = Dict{Int, Dict{UInt32, UInt8}}()
let
    for scan_idx in all_ms2_scans
        irt_lo, irt_hi = Pioneer.getRTWindow(rt_to_irt_spline(Pioneer.getRetentionTime(spectra, scan_idx)), irt_tol)

        quad_func = Pioneer.getQuadTransmissionFunction(qtm, Pioneer.getCenterMz(spectra, scan_idx), Pioneer.getIsolationWidthMz(spectra, scan_idx))
        iso_bounds = Pioneer.getIsotopeErrBounds(search_parameters)

        searchScanPartitioned!(
            counter, partitioned_index, irt_lo, irt_hi,
            Pioneer.getMzArray(spectra, scan_idx),
            mem, quad_func, iso_bounds
        )

        # Post-filter: zero out false positives outside actual quad window
        actual_prec_min = Float32(Pioneer.getPrecMinBound(quad_func) - Pioneer.NEUTRON * first(iso_bounds) / 2)
        actual_prec_max = Float32(Pioneer.getPrecMaxBound(quad_func) + Pioneer.NEUTRON * last(iso_bounds) / 2)
        @inbounds for idx in 1:(Pioneer.getSize(counter) - 1)
            pid = Pioneer.getID(counter, idx)
            pid == 0 && continue
            pmz = precursor_mzs[pid]
            if pmz < actual_prec_min || pmz > actual_prec_max
                counter.counts[pid] = zero(UInt8)
            end
        end

        partitioned_scores[scan_idx] = snapshot_counter(counter)
        Pioneer.reset!(counter)
    end
end

# Compare scores
let n_perfect = 0, n_score_mismatch = 0, n_id_mismatch = 0, n_empty_both = 0,
    total_precs_checked = 0, max_score_diff = UInt8(0), printed = 0,
    n_part_lower = 0, n_part_higher = 0, total_score_mismatches = 0
    for scan_idx in all_ms2_scans
        base = baseline_scores[scan_idx]
        part = partitioned_scores[scan_idx]

        if isempty(base) && isempty(part)
            n_empty_both += 1
            n_perfect += 1
            continue
        end

        base_keys = Set(keys(base))
        part_keys = Set(keys(part))

        if base_keys != part_keys
            n_id_mismatch += 1
            if printed < 5
                only_b = setdiff(base_keys, part_keys)
                only_p = setdiff(part_keys, base_keys)
                println("  Scan $scan_idx ID mismatch: baseline_only=$(length(only_b)), partitioned_only=$(length(only_p))")
                if length(only_b) <= 3
                    for pid in only_b
                        println("    baseline-only: prec_id=$pid score=$(base[pid]) prec_mz=$(precursor_mzs[pid])")
                    end
                end
                if length(only_p) <= 3
                    for pid in only_p
                        println("    partitioned-only: prec_id=$pid score=$(part[pid]) prec_mz=$(precursor_mzs[pid])")
                    end
                end
                printed += 1
            end
            continue
        end

        # Same ID set — compare scores
        scores_match = true
        for pid in base_keys
            total_precs_checked += 1
            if base[pid] != part[pid]
                scores_match = false
                total_score_mismatches += 1
                diff = abs(Int(base[pid]) - Int(part[pid]))
                max_score_diff = max(max_score_diff, UInt8(min(diff, 255)))
                if part[pid] < base[pid]
                    n_part_lower += 1
                else
                    n_part_higher += 1
                end
                if printed < 5
                    println("  Scan $scan_idx prec_id=$pid SCORE MISMATCH: baseline=$(base[pid]) vs partitioned=$(part[pid]) prec_mz=$(precursor_mzs[pid])")
                    printed += 1
                end
            end
        end

        if scores_match
            n_perfect += 1
        else
            n_score_mismatch += 1
        end
    end

    n_total = length(all_ms2_scans)
    println("\nScore-level validation summary:")
    println("  Total scans:       $n_total")
    println("  Perfect match:     $n_perfect ($(round(100*n_perfect/n_total, digits=1))%)")
    println("  ID set mismatch:   $n_id_mismatch")
    println("  Score mismatch:    $n_score_mismatch")
    println("  Empty both:        $n_empty_both")
    println("  Precursors checked: $total_precs_checked")
    if n_score_mismatch > 0
        println("  Max score diff:    $max_score_diff")
        println("  Score mismatches:  $total_score_mismatches precursors total")
        println("    partitioned < baseline: $n_part_lower")
        println("    partitioned > baseline: $n_part_higher")
    end
    if n_id_mismatch == 0 && n_score_mismatch == 0
        println("  PASS — all scans produce identical precursor ID sets and scores")
    end
end

# ── Benchmark ──────────────────────────────────────────────────────────────
println("\n" * "="^80)
println("BENCHMARK")
println("="^80)

# Warmup both paths
println("Warmup ...")
let
    s2p = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra))
    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin
            tid = first(thread_task)
            Pioneer.searchFragmentIndex(
                s2p, native_frag_index, spectra, last(thread_task),
                search_data[tid], search_parameters, qtm, mem, rt_to_irt_spline, irt_tol
            )
        end
    end
    fetch.(tasks)

    s2p2 = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra))
    tasks2 = map(thread_tasks) do thread_task
        Threads.@spawn begin
            tid = first(thread_task)
            searchFragmentIndexPartitioned(
                s2p2, partitioned_index, spectra, last(thread_task),
                search_data[tid], search_parameters, qtm, mem, rt_to_irt_spline, irt_tol,
                precursor_mzs
            )
        end
    end
    fetch.(tasks2)
end

# Timed baseline
println("\nRunning baseline (NativeFragmentIndex) ...")
baseline_time = @elapsed begin
    s2p_base = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra))
    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin
            tid = first(thread_task)
            Pioneer.searchFragmentIndex(
                s2p_base, native_frag_index, spectra, last(thread_task),
                search_data[tid], search_parameters, qtm, mem, rt_to_irt_spline, irt_tol
            )
        end
    end
    fetch.(tasks)
end
println("  Baseline time: $(round(baseline_time, digits=3))s")

# Timed partitioned
println("Running partitioned (PartitionedFragmentIndex) ...")
partitioned_time = @elapsed begin
    s2p_part = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra))
    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin
            tid = first(thread_task)
            searchFragmentIndexPartitioned(
                s2p_part, partitioned_index, spectra, last(thread_task),
                search_data[tid], search_parameters, qtm, mem, rt_to_irt_spline, irt_tol,
                precursor_mzs
            )
        end
    end
    fetch.(tasks)
end
println("  Partitioned time: $(round(partitioned_time, digits=3))s")

speedup = baseline_time / partitioned_time
println("\n" * "-"^40)
println("  Baseline:    $(round(baseline_time, digits=3))s")
println("  Partitioned: $(round(partitioned_time, digits=3))s")
println("  Speedup:     $(round(speedup, digits=2))x")
println("  Build cost:  $(round(build_time, digits=3))s ($(round(build_time/baseline_time*100, digits=1))% of one baseline search)")
println("-"^40)

# ── Serialize diagnostic data ─────────────────────────────────────────────
using Serialization

diag_path = joinpath(dirname(PARAMS_PATH), "partitioned_diag.jls")
println("\nSerializing diagnostic data to $diag_path ...")

# Collect a few mismatched scans with full context
diag_scans = Dict{Int, Any}()
let
    rt_bins_ref  = Pioneer.getRTBins(native_frag_index)
    frag_bins_ref = Pioneer.getFragBins(native_frag_index)
    fragments_ref = Pioneer.getFragments(native_frag_index)

    rt_bin_idx_b = 1
    n_collected = 0

    for scan_idx in all_ms2_scans
        n_collected >= 10 && break

        irt_lo, irt_hi = Pioneer.getRTWindow(rt_to_irt_spline(Pioneer.getRetentionTime(spectra, scan_idx)), irt_tol)

        while rt_bin_idx_b < length(rt_bins_ref) && Pioneer.getHigh(rt_bins_ref[rt_bin_idx_b]) < irt_lo
            rt_bin_idx_b += 1
        end
        while rt_bin_idx_b > 1 && Pioneer.getLow(rt_bins_ref[rt_bin_idx_b]) > irt_lo
            rt_bin_idx_b -= 1
        end

        quad_func = Pioneer.getQuadTransmissionFunction(qtm, Pioneer.getCenterMz(spectra, scan_idx), Pioneer.getIsolationWidthMz(spectra, scan_idx))
        iso_bounds = Pioneer.getIsotopeErrBounds(search_parameters)

        # Baseline
        Pioneer.searchScan!(counter, rt_bins_ref, frag_bins_ref, fragments_ref,
            Pioneer.getMzArray(spectra, scan_idx), Pioneer.getIntensityArray(spectra, scan_idx),
            rt_bin_idx_b, irt_hi, mem, quad_func, iso_bounds)
        base_snap = snapshot_counter(counter)
        Pioneer.reset!(counter)

        # Partitioned
        searchScanPartitioned!(counter, partitioned_index, irt_lo, irt_hi,
            Pioneer.getMzArray(spectra, scan_idx), mem, quad_func, iso_bounds)
        actual_prec_min = Float32(Pioneer.getPrecMinBound(quad_func) - Pioneer.NEUTRON * first(iso_bounds) / 2)
        actual_prec_max = Float32(Pioneer.getPrecMaxBound(quad_func) + Pioneer.NEUTRON * last(iso_bounds) / 2)
        @inbounds for idx in 1:(Pioneer.getSize(counter) - 1)
            pid = Pioneer.getID(counter, idx)
            pid == 0 && continue
            pmz = precursor_mzs[pid]
            if pmz < actual_prec_min || pmz > actual_prec_max
                counter.counts[pid] = zero(UInt8)
            end
        end
        part_snap = snapshot_counter(counter)
        Pioneer.reset!(counter)

        # Check for score mismatch
        has_mismatch = false
        for pid in keys(base_snap)
            if haskey(part_snap, pid) && base_snap[pid] != part_snap[pid]
                has_mismatch = true
                break
            end
        end
        !has_mismatch && continue

        n_collected += 1
        masses_vec = collect(Float32, skipmissing(Pioneer.getMzArray(spectra, scan_idx)))
        diag_scans[scan_idx] = (
            scan_idx      = scan_idx,
            baseline_scores = base_snap,
            partitioned_scores = part_snap,
            rt_bin_idx    = rt_bin_idx_b,
            irt_lo        = irt_lo,
            irt_hi        = irt_hi,
            prec_min      = actual_prec_min,
            prec_max      = actual_prec_max,
            center_mz     = Pioneer.getCenterMz(spectra, scan_idx),
            iso_width     = Pioneer.getIsolationWidthMz(spectra, scan_idx),
            iso_bounds    = iso_bounds,
            masses        = masses_vec,
        )
    end
end

diag_data = (
    # Mismatched scan data
    diag_scans        = diag_scans,
    # Native index (plain vectors, serializable)
    native_frag_bins  = native_frag_index.fragment_bins,
    native_rt_bins    = native_frag_index.rt_bins,
    native_fragments  = native_frag_index.fragments,
    # Partitioned index
    partitioned_index = partitioned_index,
    # Precursor mz lookup
    precursor_mzs     = collect(Float32, precursor_mzs),
    # Mass error model parameters (for getCorrectedMz / getMzBoundsReverse)
    mass_err_model    = mem,
)

serialize(diag_path, diag_data)
println("  Saved $(length(diag_scans)) mismatched scans + both indices")
println("  File size: $(round(filesize(diag_path) / 1024^2, digits=1)) MB")
println("  Analyze with: julia --project=. experiments/PartitionedFragIndex/diagnose.jl $diag_path")
