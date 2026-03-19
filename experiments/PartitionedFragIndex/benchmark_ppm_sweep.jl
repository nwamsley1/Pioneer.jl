#!/usr/bin/env julia
#
# Benchmark: Fragment bin tolerance (PPM) sweep
#
# Builds partitioned indices at different frag_bin_tol_ppm values and benchmarks
# the hinted SoA+SIMD search (threshold=128) for each.
#
# Usage:
#   julia --threads=auto --project=. experiments/PartitionedFragIndex/benchmark_ppm_sweep.jl /path/to/params.json

using Pioneer

include(joinpath(@__DIR__, "partitioned_types.jl"))
include(joinpath(@__DIR__, "build_partitioned_index.jl"))
include(joinpath(@__DIR__, "search_partitioned_index.jl"))

if isempty(ARGS)
    error("Usage: julia --threads=auto --project=. benchmark_ppm_sweep.jl <params.json>")
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

# ── Build partitioned indices at each PPM tolerance ───────────────────────
const PPM_VALUES = Float32[2.5, 5.0, 10.0, 15.0, 20.0]
const LINEAR_THRESHOLD = UInt32(128)

indices = Dict{Float32, Any}()
for ppm in PPM_VALUES
    println("\nBuilding partitioned index (frag_bin_tol_ppm=$(ppm)) ...")
    t_build = @elapsed begin
        indices[ppm] = build_partitioned_index_from_lib(Pioneer.getSpecLib(SC);
            partition_width=5.0f0, frag_bin_tol_ppm=ppm, rt_bin_tol=1.0f0,
            rank_to_score=UInt8[8, 4, 4, 2, 2, 1, 1])
    end
    pi_cur = indices[ppm]
    n_frag_bins = sum(length(p.fragment_bins.highs) for p in getPartitions(pi_cur))
    n_frags = sum(length(getFragments(p)) for p in getPartitions(pi_cur))
    mem_frags = sum(sizeof(getFragments(p)) for p in getPartitions(pi_cur)) / 1024^2
    mem_bins = sum(sizeof(p.fragment_bins.lows) + sizeof(p.fragment_bins.highs) +
                   sizeof(p.fragment_bins.first_bins) + sizeof(p.fragment_bins.last_bins)
                   for p in getPartitions(pi_cur)) / 1024^2
    mem_hints = sum(sizeof(p.skip_hints) for p in getPartitions(pi_cur)) / 1024^2
    mem_l2g = sum(sizeof(p.local_to_global) for p in getPartitions(pi_cur)) / 1024^2
    println("  Build time: $(round(t_build, digits=2))s")
    println("  frag_bins: $n_frag_bins  fragments: $n_frags")
    println("  Memory: frags=$(round(mem_frags, digits=1))MB  bins=$(round(mem_bins, digits=1))MB  hints=$(round(mem_hints, digits=1))MB  l2g=$(round(mem_l2g, digits=1))MB  total=$(round(mem_frags+mem_bins+mem_hints+mem_l2g, digits=1))MB")
end

println("\n$(length(all_scans)) MS2 scans, $(Threads.nthreads()) threads")

# ── Warmup ────────────────────────────────────────────────────────────────
println("\nWarmup ...")
let
    # Baseline warmup
    s2p = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra))
    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin
            tid = first(thread_task)
            Pioneer.searchFragmentIndex(s2p, native_frag_index, spectra, last(thread_task),
                search_data[tid], search_parameters, qtm, mem, rt_to_irt_spline, irt_tol)
        end
    end
    fetch.(tasks)

    # Partitioned warmup (each PPM)
    for ppm in PPM_VALUES
        s2p2 = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra))
        searchFragmentIndexPartitionMajorHinted(s2p2, indices[ppm], spectra, all_scans,
            Threads.nthreads(), search_parameters, qtm, mem, rt_to_irt_spline, irt_tol,
            precursor_mzs; linear_threshold=LINEAR_THRESHOLD)
    end
end

# ── Timed runs ────────────────────────────────────────────────────────────
println("\n" * "="^70)
println("TIMED BENCHMARK — PPM Tolerance Sweep (threshold=$LINEAR_THRESHOLD)")
println("="^70)

# Baseline
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

# Each PPM value
ppm_times = Dict{Float32, Float64}()
for ppm in PPM_VALUES
    println("Running partitioned (ppm=$(ppm), threshold=$LINEAR_THRESHOLD) ...")
    t = @elapsed begin
        s2p = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra))
        searchFragmentIndexPartitionMajorHinted(s2p, indices[ppm], spectra, all_scans,
            Threads.nthreads(), search_parameters, qtm, mem, rt_to_irt_spline, irt_tol,
            precursor_mzs; linear_threshold=LINEAR_THRESHOLD)
    end
    ppm_times[ppm] = t
end

# ── Results ───────────────────────────────────────────────────────────────
println("\n" * "="^70)
println("RESULTS — PPM Tolerance Sweep (threshold=$LINEAR_THRESHOLD)")
println("="^70)
println()
println("  Baseline (native monolithic):  $(round(t_baseline, digits=3))s   1.00x")
println()

# Collect stats for table
println("  PPM     | Time (s) | Speedup | Frag Bins    | Fragments    | Memory (MB)")
println("  --------|----------|---------|--------------|--------------|------------")
for ppm in PPM_VALUES
    t = ppm_times[ppm]
    speedup = t_baseline / t
    pi_cur = indices[ppm]
    n_frag_bins = sum(length(p.fragment_bins.highs) for p in getPartitions(pi_cur))
    n_frags = sum(length(getFragments(p)) for p in getPartitions(pi_cur))
    total_mem = sum(
        sizeof(getFragments(p)) +
        sizeof(p.fragment_bins.lows) + sizeof(p.fragment_bins.highs) +
        sizeof(p.fragment_bins.first_bins) + sizeof(p.fragment_bins.last_bins) +
        sizeof(p.skip_hints) + sizeof(p.local_to_global)
        for p in getPartitions(pi_cur)
    ) / 1024^2
    best = (t == minimum(values(ppm_times))) ? " <- best" : ""
    println("  $(lpad(ppm, 7)) | $(lpad(round(t, digits=3), 8)) | $(lpad(round(speedup, digits=2), 7))x | $(lpad(n_frag_bins, 12)) | $(lpad(n_frags, 12)) | $(lpad(round(total_mem, digits=1), 10))$best")
end

best_ppm = argmin(ppm_times)
best_t = ppm_times[best_ppm]
println()
println("-"^70)
println("  Best PPM: $(best_ppm)  ($(round(best_t, digits=3))s, $(round(t_baseline/best_t, digits=2))x)")
println("-"^70)
println()
