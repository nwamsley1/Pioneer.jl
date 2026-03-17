#!/usr/bin/env julia
#
# Compare search optimization approaches:
#   Baseline: exponential + binary search (current)
#   Opt A:    m/z lookup table (reciprocal multiply → bucket → approx index)
#   Opt B:    delta-based skip (reciprocal multiply → estimated jump)
#
# Usage:
#   julia --threads=auto --project=. experiments/PartitionedFragIndex/benchmark_search_opts.jl /path/to/params.json

using Pioneer

include(joinpath(@__DIR__, "partitioned_types.jl"))
include(joinpath(@__DIR__, "build_partitioned_index.jl"))
include(joinpath(@__DIR__, "search_partitioned_index.jl"))
include(joinpath(@__DIR__, "SearchOptA", "search_opt_a.jl"))
include(joinpath(@__DIR__, "SearchOptB", "search_opt_b.jl"))

if isempty(ARGS)
    error("Usage: julia --threads=auto --project=. benchmark_search_opts.jl <params.json>")
end
const PARAMS_PATH = ARGS[1]
isfile(PARAMS_PATH) || error("Params file not found: $PARAMS_PATH")

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
min_score        = Pioneer.getMinIndexSearchScore(search_parameters)

# Build index
println("\nBuilding partitioned index (7 frags, weighted) ...")
pi = build_partitioned_index_from_lib(Pioneer.getSpecLib(SC);
    partition_width=5.0f0, frag_bin_tol_ppm=10.0f0, rt_bin_tol=1.0f0,
    rank_to_score=UInt8[8, 4, 4, 2, 2, 1, 1])

# Also materialize native index for baseline comparison
println("Materializing native fragment index ...")
native_frag_index = Pioneer.materialize(Pioneer.getFragmentIndex(Pioneer.getSpecLib(SC)))

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

# ── Warmup ────────────────────────────────────────────────────────────────
println("Warmup ...")
let
    # Baseline (native)
    s2p = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra))
    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin
            tid = first(thread_task)
            Pioneer.searchFragmentIndex(s2p, native_frag_index, spectra, last(thread_task),
                search_data[tid], search_parameters, qtm, mem, rt_to_irt_spline, irt_tol)
        end
    end
    fetch.(tasks)

    # Current partitioned (exp+binary)
    s2p2 = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra))
    searchFragmentIndexPartitionMajor(s2p2, pi, spectra, all_scans,
        Threads.nthreads(), search_parameters, qtm, mem, rt_to_irt_spline, irt_tol, precursor_mzs)

    # Opt A (lookup)
    s2p3 = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra))
    searchPartitionMajorLookup(s2p3, pi, spectra, all_scans,
        Threads.nthreads(), min_score, qtm, mem, rt_to_irt_spline, irt_tol, precursor_mzs)

    # Opt B (delta)
    s2p4 = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra))
    searchPartitionMajorDelta(s2p4, pi, spectra, all_scans,
        Threads.nthreads(), min_score, qtm, mem, rt_to_irt_spline, irt_tol, precursor_mzs)
end

# ── Timed runs ────────────────────────────────────────────────────────────
println("\n" * "="^70)
println("TIMED BENCHMARK")
println("="^70)

# Baseline (native, monolithic)
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

# Current partitioned (exp+binary)
println("Running partitioned (exp+binary search) ...")
t_current = @elapsed begin
    s2p = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra))
    searchFragmentIndexPartitionMajor(s2p, pi, spectra, all_scans,
        Threads.nthreads(), search_parameters, qtm, mem, rt_to_irt_spline, irt_tol, precursor_mzs)
end

# Opt A (lookup table)
println("Running Opt A (lookup table) ...")
t_lookup = @elapsed begin
    s2p = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra))
    searchPartitionMajorLookup(s2p, pi, spectra, all_scans,
        Threads.nthreads(), min_score, qtm, mem, rt_to_irt_spline, irt_tol, precursor_mzs)
end

# Opt B (delta skip)
println("Running Opt B (delta skip) ...")
t_delta = @elapsed begin
    s2p = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra))
    searchPartitionMajorDelta(s2p, pi, spectra, all_scans,
        Threads.nthreads(), min_score, qtm, mem, rt_to_irt_spline, irt_tol, precursor_mzs)
end

println("\n" * "-"^70)
println("  Baseline (native):           $(round(t_baseline, digits=2))s   1.00x")
println("  Partitioned (exp+binary):    $(round(t_current, digits=2))s   $(round(t_baseline/t_current, digits=2))x")
println("  Opt A (lookup table):        $(round(t_lookup, digits=2))s   $(round(t_baseline/t_lookup, digits=2))x")
println("  Opt B (delta skip):          $(round(t_delta, digits=2))s   $(round(t_baseline/t_delta, digits=2))x")
println("-"^70)
