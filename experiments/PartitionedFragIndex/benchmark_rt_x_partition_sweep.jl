#!/usr/bin/env julia
#
# Benchmark: RT bin tolerance × partition width cross-sweep
#
# Tests wider RT bins (fewer per-scan RT bin iterations) with different
# partition widths (prec m/z). The idea: wider RT bins → fewer loops,
# then post-filter by exact precursor iRT on the backend.
#
# Usage:
#   julia --threads=auto --project=. experiments/PartitionedFragIndex/benchmark_rt_x_partition_sweep.jl /path/to/params.json

using Pioneer

include(joinpath(@__DIR__, "partitioned_types.jl"))
include(joinpath(@__DIR__, "build_partitioned_index.jl"))
include(joinpath(@__DIR__, "search_partitioned_index.jl"))

if isempty(ARGS)
    error("Usage: julia --threads=auto --project=. benchmark_rt_x_partition_sweep.jl <params.json>")
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
SPEC_LIB = Pioneer.loadSpectralLibrary(params.paths[:library], params)

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

println("\nCalibrated iRT tolerance: $irt_tol")
println("Search window: $(round(2*irt_tol, digits=2)) iRT units")

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

# ── Materialize baseline ─────────────────────────────────────────────────
println("\nMaterializing native fragment index ...")
native_frag_index = Pioneer.materialize(Pioneer.getFragmentIndex(Pioneer.getSpecLib(SC)))

# ── Configuration ─────────────────────────────────────────────────────────
const FRAG_BIN_PPM = 10.0f0  # Hold constant
const LINEAR_THRESHOLD = UInt32(128)
const RT_BIN_TOLS = Float32[1.0, 2.0, 3.0]
const PARTITION_WIDTHS = Float32[2.5, 5.0, 10.0]

# ── Build all index combinations ─────────────────────────────────────────
indices = Dict{Tuple{Float32, Float32}, Any}()

for rt_tol in RT_BIN_TOLS
    for pw in PARTITION_WIDTHS
        println("\nBuilding index (rt_bin_tol=$(rt_tol), partition_width=$(pw)) ...")
        t_build = @elapsed begin
            indices[(rt_tol, pw)] = build_partitioned_index_from_lib(Pioneer.getSpecLib(SC);
                partition_width=pw, frag_bin_tol_ppm=FRAG_BIN_PPM, rt_bin_tol=rt_tol,
                rank_to_score=UInt8[8, 4, 4, 2, 2, 1, 1])
        end
        pi_cur = indices[(rt_tol, pw)]
        n_parts = pi_cur.n_partitions
        n_rt_bins = sum(length(getRTBins(p)) for p in getPartitions(pi_cur))
        n_frags = sum(length(getFragments(p)) for p in getPartitions(pi_cur))
        total_mem = sum(
            sizeof(getFragments(p)) +
            sizeof(p.fragment_bins.lows) + sizeof(p.fragment_bins.highs) +
            sizeof(p.fragment_bins.first_bins) + sizeof(p.fragment_bins.last_bins) +
            sizeof(p.skip_hints) + sizeof(p.local_to_global)
            for p in getPartitions(pi_cur)
        ) / 1024^2
        avg_rt = n_rt_bins / n_parts
        println("  Build: $(round(t_build, digits=1))s  partitions=$n_parts  rt_bins=$n_rt_bins (avg $(round(avg_rt, digits=1))/part)  frags=$n_frags  mem=$(round(total_mem, digits=0))MB")
        rt_bins_per_scan = 2 * irt_tol / rt_tol
        println("  Expected RT bins/scan: ~$(round(rt_bins_per_scan, digits=1))")
    end
end

println("\n$(length(all_scans)) MS2 scans, $(Threads.nthreads()) threads")

# ── Helper ────────────────────────────────────────────────────────────────
function count_precursor_hits(scan_to_prec_idx)
    total = 0
    for entry in scan_to_prec_idx
        if !ismissing(entry)
            total += length(entry)
        end
    end
    return total
end

# ── Warmup ────────────────────────────────────────────────────────────────
println("\nWarmup ...")
let
    # Baseline
    local s2p = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra))
    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin
            tid = first(thread_task)
            Pioneer.searchFragmentIndex(s2p, native_frag_index, spectra, last(thread_task),
                search_data[tid], search_parameters, qtm, mem, rt_to_irt_spline, irt_tol)
        end
    end
    fetch.(tasks)

    # Warmup each combination
    for rt_tol in RT_BIN_TOLS
        for pw in PARTITION_WIDTHS
            local s2p2 = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra))
            searchFragmentIndexPartitionMajorHinted(s2p2, indices[(rt_tol, pw)], spectra, all_scans,
                Threads.nthreads(), search_parameters, qtm, mem, rt_to_irt_spline, irt_tol,
                precursor_mzs; linear_threshold=LINEAR_THRESHOLD)
        end
    end
end

# ── Timed runs ────────────────────────────────────────────────────────────
println("\n" * "="^80)
println("TIMED BENCHMARK — RT bin tol × Partition width (frag_ppm=$FRAG_BIN_PPM, threshold=$LINEAR_THRESHOLD)")
println("="^80)

# Baseline
println("\nRunning baseline (native, monolithic) ...")
s2p_base = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra))
t_baseline = @elapsed begin
    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin
            tid = first(thread_task)
            Pioneer.searchFragmentIndex(s2p_base, native_frag_index, spectra, last(thread_task),
                search_data[tid], search_parameters, qtm, mem, rt_to_irt_spline, irt_tol)
        end
    end
    fetch.(tasks)
end
baseline_hits = count_precursor_hits(s2p_base)
println("  $(round(t_baseline, digits=2))s, $baseline_hits hits")

# Each combination
results = Dict{Tuple{Float32, Float32}, NamedTuple{(:time, :hits), Tuple{Float64, Int}}}()
for rt_tol in RT_BIN_TOLS
    for pw in PARTITION_WIDTHS
        println("Running rt_bin_tol=$(rt_tol), partition_width=$(pw) ...")
        local s2p = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra))
        t = @elapsed begin
            searchFragmentIndexPartitionMajorHinted(s2p, indices[(rt_tol, pw)], spectra, all_scans,
                Threads.nthreads(), search_parameters, qtm, mem, rt_to_irt_spline, irt_tol,
                precursor_mzs; linear_threshold=LINEAR_THRESHOLD)
        end
        hits = count_precursor_hits(s2p)
        results[(rt_tol, pw)] = (time=t, hits=hits)
        println("  $(round(t, digits=2))s, $hits hits, $(round(t_baseline/t, digits=2))x")
    end
end

# ── Results ───────────────────────────────────────────────────────────────
println("\n" * "="^80)
println("RESULTS — RT bin tol × Partition width")
println("="^80)
println("\nBaseline: $(round(t_baseline, digits=2))s, $baseline_hits hits")
println("iRT tolerance (calibrated): $(round(irt_tol, digits=3)), search window: $(round(2*irt_tol, digits=2)) iRT")
println()

# Table: rows = rt_bin_tol, columns = partition_width
local header = "  rt_bin \\ part_w |"
for pw in PARTITION_WIDTHS
    header *= "   $(pw) Da           |"
end
println(header)
println("  " * "-"^(length(header)-2))

for rt_tol in RT_BIN_TOLS
    row = "  rt=$(lpad(rt_tol, 4))       |"
    for pw in PARTITION_WIDTHS
        r = results[(rt_tol, pw)]
        speedup = t_baseline / r.time
        row *= " $(lpad(round(r.time, digits=2), 6))s $(lpad(round(speedup, digits=2), 5))x $(lpad(r.hits, 8)) |"
    end
    println(row)
end

# Summary: best combination
println()
best_key = argmin(Dict(k => v.time for (k, v) in results))
best_r = results[best_key]
println("-"^80)
println("  Best: rt_bin_tol=$(best_key[1]), partition_width=$(best_key[2])")
println("    $(round(best_r.time, digits=2))s, $(round(t_baseline/best_r.time, digits=2))x speedup, $(best_r.hits) hits")
println("-"^80)

# Also show hit ratio relative to baseline
println("\nHit ratios (partitioned / baseline):")
row = "  rt_bin \\ part_w |"
for pw in PARTITION_WIDTHS
    row *= "   $(pw) Da  |"
end
println(row)
println("  " * "-"^(length(row)-2))
for rt_tol in RT_BIN_TOLS
    row = "  rt=$(lpad(rt_tol, 4))       |"
    for pw in PARTITION_WIDTHS
        r = results[(rt_tol, pw)]
        ratio = r.hits / baseline_hits
        row *= " $(lpad(round(ratio, digits=3), 8)) |"
    end
    println(row)
end
println("="^80)
