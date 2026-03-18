#!/usr/bin/env julia
#
# Benchmark: SoA + hybrid binary→SIMD hint-based search
#
# Compares:
#   - Baseline (native monolithic index)
#   - Hybrid binary→SIMD with threshold sweep (binary narrows to ≤threshold, SIMD finishes)
#
# Usage:
#   julia --threads=auto --project=. experiments/PartitionedFragIndex/benchmark_step_and_hints.jl /path/to/params.json

using Pioneer

include(joinpath(@__DIR__, "partitioned_types.jl"))
include(joinpath(@__DIR__, "build_partitioned_index.jl"))
include(joinpath(@__DIR__, "search_partitioned_index.jl"))

if isempty(ARGS)
    error("Usage: julia --threads=auto --project=. benchmark_step_and_hints.jl <params.json>")
end
const PARAMS_PATH = ARGS[1]
isfile(PARAMS_PATH) || error("Params file not found: $PARAMS_PATH")

# ── Load & calibrate ──────────────────────────────────────────────────────
println("Loading parameters from $PARAMS_PATH ...")
Pioneer.checkParams(PARAMS_PATH)
params = Pioneer.parse_pioneer_parameters(PARAMS_PATH)
mkpath(params.paths[:results])

MS_DATA_DIR = params.paths[:ms_data]
MS_TABLE_PATHS = [joinpath(MS_DATA_DIR, f) for f in readdir(MS_DATA_DIR)
                  if isfile(joinpath(MS_DATA_DIR, f)) && endswith(f, ".arrow")]
SPEC_LIB_DIR = params.paths[:library]
SPEC_LIB = Pioneer.loadSpectralLibrary(SPEC_LIB_DIR, params)

SC = Pioneer.initSearchContext(SPEC_LIB, Pioneer.parseIsoXML(Pioneer.isotope_spline_path()),
    Pioneer.ArrowTableReference(MS_TABLE_PATHS), Threads.nthreads(), 250000)
Pioneer.setDataOutDir!(SC, params.paths[:results])

Pioneer.execute_search(Pioneer.ParameterTuningSearch(), SC, params)
Pioneer.execute_search(Pioneer.NceTuningSearch(), SC, params)
Pioneer.execute_search(Pioneer.QuadTuningSearch(), SC, params)

search_parameters = Pioneer.FirstPassSearchParameters(params)
search_data       = Pioneer.getSearchData(SC)
msdr = Pioneer.getMassSpecData(SC)
ms_file_idx, spectra = first(enumerate(msdr))
qtm              = Pioneer.getQuadTransmissionModel(SC, ms_file_idx)
mem              = Pioneer.getMassErrorModel(SC, ms_file_idx)
rt_to_irt_spline = Pioneer.getRtIrtModel(SC, ms_file_idx)
irt_tol          = Pioneer.getIrtErrors(SC)[ms_file_idx]
thread_tasks     = Pioneer.partition_scans(spectra, Threads.nthreads())
precursor_mzs    = Pioneer.getMz(Pioneer.getPrecursors(Pioneer.getSpecLib(SC)))

# ── Build indices ─────────────────────────────────────────────────────────
println("\nMaterializing native fragment index ...")
native_frag_index = Pioneer.materialize(Pioneer.getFragmentIndex(Pioneer.getSpecLib(SC)))

println("Building partitioned index (7 frags, weighted, SoA layout) ...")
pi = build_partitioned_index_from_lib(Pioneer.getSpecLib(SC);
    partition_width=5.0f0, frag_bin_tol_ppm=10.0f0, rt_bin_tol=1.0f0,
    rank_to_score=UInt8[8, 4, 4, 2, 2, 1, 1])

# Collect all valid MS2 scan indices
all_scans = Int[]
for (_, ts) in thread_tasks
    for s in ts
        (s <= 0 || s > length(spectra)) && continue
        Pioneer.getMsOrder(spectra, s) ∉ Pioneer.getSpecOrder(search_parameters) && continue
        push!(all_scans, s)
    end
end
sort!(all_scans)
println("$(length(all_scans)) MS2 scans, $(Threads.nthreads()) threads\n")

# ── Warmup all code paths ────────────────────────────────────────────────
println("Warmup ...")
let
    # Baseline
    s2p = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra))
    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin
            tid = first(thread_task)
            Pioneer.searchFragmentIndex(s2p, native_frag_index, spectra, last(thread_task),
                search_data[tid], search_parameters, qtm, mem, rt_to_irt_spline, irt_tol)
        end
    end
    fetch.(tasks)

    # Hybrid binary→SIMD (warmup with several threshold values)
    for lt in [1, 8, 32, 128, 1000000]
        s2p3 = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra))
        searchFragmentIndexPartitionMajorHinted(s2p3, pi, spectra, all_scans,
            Threads.nthreads(), search_parameters, qtm, mem, rt_to_irt_spline, irt_tol,
            precursor_mzs; linear_threshold=lt)
    end
end

# ── Timed runs ───────────────────────────────────────────────────────────
println("\n" * "="^70)
println("TIMED BENCHMARK — Hybrid Binary->SIMD Search")
println("="^70)

# Baseline (native monolithic)
println("\nRunning baseline (native, monolithic) ...")
t_baseline = @elapsed begin
    s2p = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra))
    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin
            tid = first(thread_task)
            Pioneer.searchFragmentIndex(s2p, native_frag_index, spectra, last(thread_task),
                search_data[tid], search_parameters, qtm, mem, rt_to_irt_spline, irt_tol)
        end
    end
    fetch.(tasks)
end

# Hinted search (default threshold=32)
println("Running hybrid search (default threshold=$HINT_LINEAR_THRESHOLD) ...")
t_hinted = @elapsed begin
    s2p = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra))
    searchFragmentIndexPartitionMajorHinted(s2p, pi, spectra, all_scans,
        Threads.nthreads(), search_parameters, qtm, mem, rt_to_irt_spline, irt_tol,
        precursor_mzs)
end

# Threshold sweep: binary narrows to ≤threshold, then SIMD finishes
# threshold=1 → almost all binary (SIMD scans 1 element)
# threshold=1000000 → almost all SIMD (binary never runs)
linear_thresholds = [1, 8, 16, 32, 64, 128, 256, 1000000]
hint_threshold_times = Dict{Int, Float64}()
for lt in linear_thresholds
    label = lt >= 1000000 ? "all-SIMD" : "$lt"
    println("Running hybrid threshold=$label ...")
    t = @elapsed begin
        local s2p = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra))
        searchFragmentIndexPartitionMajorHinted(s2p, pi, spectra, all_scans,
            Threads.nthreads(), search_parameters, qtm, mem, rt_to_irt_spline, irt_tol,
            precursor_mzs; linear_threshold=lt)
    end
    hint_threshold_times[lt] = t
end

# ── Results ──────────────────────────────────────────────────────────────
println("\n" * "="^70)
println("RESULTS")
println("="^70)
println()
println("  Baseline (native monolithic):  $(round(t_baseline, digits=3))s   1.00x")
println()
println("  Hybrid binary->SIMD (default threshold=$HINT_LINEAR_THRESHOLD): $(round(t_hinted, digits=3))s   $(round(t_baseline/t_hinted, digits=2))x")
println()
println("  Hybrid threshold sweep (binary narrows to <=N, SIMD finishes):")
println("    threshold=1 is ~pure binary, threshold=1000000 is ~pure SIMD")
for lt in linear_thresholds
    t = hint_threshold_times[lt]
    speedup = t_baseline / t
    label = lt >= 1000000 ? "     all-SIMD" : lpad(lt, 13)
    best = (t == minimum(values(hint_threshold_times))) ? " <- best" : ""
    println("    threshold=$label:  $(lpad(round(t, digits=3), 7))s   $(lpad(round(speedup, digits=2), 5))x$best")
end

best_lt = argmin(hint_threshold_times)
best_lt_t = hint_threshold_times[best_lt]
println()
println("-"^70)
println("  Best hybrid threshold: $best_lt  ($(round(best_lt_t, digits=3))s, $(round(t_baseline/best_lt_t, digits=2))x)")
println("-"^70)

# ── Correctness validation ──────────────────────────────────────────────
println("\n" * "="^70)
println("CORRECTNESS VALIDATION — threshold=1 (mostly binary) vs threshold=1000000 (all SIMD)")
println("="^70)

# threshold=1 does full binary search then SIMD scans 1 element
# threshold=1000000 does pure SIMD linear scan
# Both must produce identical results.
n_validate = min(500, length(all_scans))
validate_scans = all_scans[1:n_validate]
println("Comparing on $n_validate scans ...")

max_local = maximum(p -> Int(p.n_local_precs), getPartitions(pi); init=0)
const VALIDATION_MIN_SCORE = UInt8(1)

function score_scan_hinted(pfi, scan_idx, spectra, qtm, mem, rt_to_irt_spline, irt_tol,
                            precursor_mzs; linear_threshold::Int=HINT_LINEAR_THRESHOLD)
    irt_lo, irt_hi = Pioneer.getRTWindow(
        rt_to_irt_spline(Pioneer.getRetentionTime(spectra, scan_idx)), irt_tol)
    quad_func = Pioneer.getQuadTransmissionFunction(
        qtm, Pioneer.getCenterMz(spectra, scan_idx),
        Pioneer.getIsolationWidthMz(spectra, scan_idx))
    iso_bounds = (UInt8(1), UInt8(0))
    prec_min = Float32(Pioneer.getPrecMinBound(quad_func) - Pioneer.NEUTRON * first(iso_bounds) / 2)
    prec_max = Float32(Pioneer.getPrecMaxBound(quad_func) + Pioneer.NEUTRON * last(iso_bounds) / 2)
    masses = Pioneer.getMzArray(spectra, scan_idx)

    lc = LocalCounter(UInt16, UInt8, max_local + 1)
    result = Dict{UInt32, UInt8}()

    first_k, last_k = get_partition_range(pfi, prec_min, prec_max)
    for k in first_k:last_k
        partition = getPartition(pfi, k)
        _score_partition_hinted!(lc, partition, irt_lo, irt_hi, masses, mem;
                                  linear_threshold=linear_threshold)

        l2g = partition.local_to_global
        @inbounds for i in 1:(lc.size - 1)
            lid = lc.ids[i]
            score = lc.counts[lid]
            score < VALIDATION_MIN_SCORE && continue
            global_pid = l2g[lid]
            pmz = precursor_mzs[global_pid]
            (pmz < prec_min || pmz > prec_max) && continue
            prev = get(result, global_pid, zero(UInt8))
            if score > prev
                result[global_pid] = score
            end
        end
        reset!(lc)
    end
    return result
end

# Compare threshold=1 (mostly binary) vs threshold=1000000 (all SIMD)
n_mismatch_ids, n_mismatch_scores, n_total_ids, n_scans_with_hits = let
    _n_mismatch_ids = 0
    _n_mismatch_scores = 0
    _n_total_ids = 0
    _n_scans_with_hits = 0
    for scan_idx in validate_scans
        binary = score_scan_hinted(pi, scan_idx, spectra, qtm, mem, rt_to_irt_spline, irt_tol,
                                    precursor_mzs; linear_threshold=1)
        linear = score_scan_hinted(pi, scan_idx, spectra, qtm, mem, rt_to_irt_spline, irt_tol,
                                    precursor_mzs; linear_threshold=1000000)

        bin_ids = Set(keys(binary))
        lin_ids = Set(keys(linear))
        _n_total_ids += length(bin_ids)
        if !isempty(bin_ids) || !isempty(lin_ids)
            _n_scans_with_hits += 1
        end

        diff = symdiff(bin_ids, lin_ids)
        _n_mismatch_ids += length(diff)

        for pid in intersect(bin_ids, lin_ids)
            if binary[pid] != linear[pid]
                _n_mismatch_scores += 1
            end
        end
    end
    (_n_mismatch_ids, _n_mismatch_scores, _n_total_ids, _n_scans_with_hits)
end

println()
println("  Scans with hits: $n_scans_with_hits / $n_validate")
if n_mismatch_ids == 0 && n_mismatch_scores == 0
    println("  PASS: threshold=1 (binary) == threshold=1000000 (SIMD) on all $n_validate scans")
    println("    ($n_total_ids total precursor IDs, 0 mismatches)")
else
    println("  FAIL:")
    println("    ID set mismatches:    $n_mismatch_ids / $n_total_ids")
    println("    Score mismatches:     $n_mismatch_scores")
end
println("="^70)
