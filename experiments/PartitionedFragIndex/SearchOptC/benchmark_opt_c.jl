#!/usr/bin/env julia
#
# Benchmark: Fixed-width bin O(1) lookup vs baseline and current partitioned
#
# Usage:
#   julia --threads=auto --project=. experiments/PartitionedFragIndex/SearchOptC/benchmark_opt_c.jl /path/to/params.json

using Pioneer

include(joinpath(@__DIR__, "..", "partitioned_types.jl"))
include(joinpath(@__DIR__, "..", "build_partitioned_index.jl"))
include(joinpath(@__DIR__, "..", "search_partitioned_index.jl"))
include(joinpath(@__DIR__, "types_opt_c.jl"))
include(joinpath(@__DIR__, "build_opt_c.jl"))
include(joinpath(@__DIR__, "search_opt_c.jl"))

if isempty(ARGS)
    error("Usage: julia --threads=auto --project=. benchmark_opt_c.jl <params.json>")
end
const PARAMS_PATH = ARGS[1]
isfile(PARAMS_PATH) || error("Params file not found: $PARAMS_PATH")

println("Loading and calibrating ...")
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

sp = Pioneer.FirstPassSearchParameters(params)
search_data = Pioneer.getSearchData(SC)
msdr = Pioneer.getMassSpecData(SC)
ms_file_idx, spectra = first(enumerate(msdr))
qtm = Pioneer.getQuadTransmissionModel(SC, ms_file_idx)
mem = Pioneer.getMassErrorModel(SC, ms_file_idx)
rts = Pioneer.getRtIrtModel(SC, ms_file_idx)
irt_tol = Pioneer.getIrtErrors(SC)[ms_file_idx]
tt = Pioneer.partition_scans(spectra, Threads.nthreads())
pmzs = Pioneer.getMz(Pioneer.getPrecursors(Pioneer.getSpecLib(SC)))
min_score = Pioneer.getMinIndexSearchScore(sp)

# Materialize native index
println("\nMaterializing native index ...")
native_idx = Pioneer.materialize(Pioneer.getFragmentIndex(Pioneer.getSpecLib(SC)))

# Build current partitioned index
println("Building partitioned index (current, variable bins) ...")
t_build_var = @elapsed pi_var = build_partitioned_index_from_lib(Pioneer.getSpecLib(SC);
    partition_width=5.0f0, frag_bin_tol_ppm=10.0f0, rt_bin_tol=1.0f0,
    rank_to_score=UInt8[8,4,4,2,2,1,1])

# Build fixed-bin index
println("\nBuilding fixed-bin index (0.005 Da bins) ...")
t_build_fix = @elapsed pi_fix = build_fixed_bin_index(Pioneer.getSpecLib(SC);
    partition_width=5.0f0, rt_bin_tol=1.0f0,
    rank_to_score=UInt8[8,4,4,2,2,1,1])
println("  Build time: $(round(t_build_fix, digits=2))s")

# Collect MS2 scans
all_scans = Int[]
for (_, ts) in tt; for s in ts
    (s <= 0 || s > length(spectra)) && continue
    Pioneer.getMsOrder(spectra, s) ∉ Pioneer.getSpecOrder(sp) && continue
    push!(all_scans, s)
end; end
sort!(all_scans)
println("\n$(length(all_scans)) MS2 scans, $(Threads.nthreads()) threads")

# Warmup
println("\nWarmup ...")
let
    s2p = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra))
    tasks = map(tt) do thread_task
        Threads.@spawn begin
            tid = first(thread_task)
            Pioneer.searchFragmentIndex(s2p, native_idx, spectra, last(thread_task),
                search_data[tid], sp, qtm, mem, rts, irt_tol)
        end
    end
    fetch.(tasks)

    s2p2 = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra))
    searchFragmentIndexPartitionMajor(s2p2, pi_var, spectra, all_scans,
        Threads.nthreads(), sp, qtm, mem, rts, irt_tol, pmzs)

    s2p3 = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra))
    searchPartitionMajorFixed(s2p3, pi_fix, spectra, all_scans,
        Threads.nthreads(), min_score, qtm, mem, rts, irt_tol, pmzs)
end

# Timed runs
println("\n" * "="^70)
println("TIMED BENCHMARK")
println("="^70)

println("\nBaseline (native, monolithic) ...")
t_base = @elapsed begin
    s2p = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra))
    tasks = map(tt) do thread_task
        Threads.@spawn begin
            tid = first(thread_task)
            Pioneer.searchFragmentIndex(s2p, native_idx, spectra, last(thread_task),
                search_data[tid], sp, qtm, mem, rts, irt_tol)
        end
    end
    fetch.(tasks)
end

println("Partitioned (variable bins, exp+binary search) ...")
t_var = @elapsed begin
    s2p = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra))
    searchFragmentIndexPartitionMajor(s2p, pi_var, spectra, all_scans,
        Threads.nthreads(), sp, qtm, mem, rts, irt_tol, pmzs)
end

println("Fixed-bin (0.005 Da, O(1) lookup) ...")
t_fix = @elapsed begin
    s2p = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra))
    searchPartitionMajorFixed(s2p, pi_fix, spectra, all_scans,
        Threads.nthreads(), min_score, qtm, mem, rts, irt_tol, pmzs)
end

println("\n" * "-"^70)
println("  Baseline (native):              $(round(t_base, digits=2))s   1.00x")
println("  Partitioned (variable bins):    $(round(t_var, digits=2))s   $(round(t_base/t_var, digits=2))x")
println("  Fixed-bin (O(1) lookup):        $(round(t_fix, digits=2))s   $(round(t_base/t_fix, digits=2))x")
println("  Fixed-bin improvement over var: $(round(t_var/t_fix, digits=2))x faster")
println("-"^70)
