#!/usr/bin/env julia
#
# Benchmark: Fragment bin PPM × Search tolerance cross-sweep
#
# For each search tolerance (the MassErrorModel left/right ppm), builds partitioned
# indices at different frag_bin_tol_ppm and benchmarks the hinted SoA+SIMD search.
# Reports time AND precursor hit count since wider tolerances produce more matches.
#
# Usage:
#   julia --threads=auto --project=. experiments/PartitionedFragIndex/benchmark_ppm_x_tol_sweep.jl /path/to/params.json

using Pioneer

include(joinpath(@__DIR__, "partitioned_types.jl"))
include(joinpath(@__DIR__, "build_partitioned_index.jl"))
include(joinpath(@__DIR__, "search_partitioned_index.jl"))

if isempty(ARGS)
    error("Usage: julia --threads=auto --project=. benchmark_ppm_x_tol_sweep.jl <params.json>")
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
mem_calibrated   = Pioneer.getMassErrorModel(SC, ms_file_idx)
rt_to_irt_spline = Pioneer.getRtIrtModel(SC, ms_file_idx)
irt_tol          = Pioneer.getIrtErrors(SC)[ms_file_idx]
thread_tasks     = Pioneer.partition_scans(spectra, Threads.nthreads())
precursor_mzs    = Pioneer.getMz(Pioneer.getPrecursors(Pioneer.getSpecLib(SC)))

println("\nCalibrated MassErrorModel:")
println("  mass_offset: $(mem_calibrated.mass_offset)")
println("  left_tol:  $(Pioneer.getLeftTol(mem_calibrated)) ppm")
println("  right_tol: $(Pioneer.getRightTol(mem_calibrated)) ppm")
println("  total window: ~$(round(Pioneer.getLeftTol(mem_calibrated) + Pioneer.getRightTol(mem_calibrated), digits=1)) ppm")

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
const PPM_VALUES = Float32[2.5, 5.0, 10.0, 15.0, 20.0]
const LINEAR_THRESHOLD = UInt32(128)

# Search tolerances to test: (left_ppm, right_ppm)
# Keep the calibrated offset, only change the tolerance window
const SEARCH_TOLS = [
    (5.0f0,  5.0f0,  "10 ppm (±5)"),
    (7.5f0,  7.5f0,  "15 ppm (±7.5)"),
    (10.0f0, 10.0f0, "20 ppm (±10)"),
    (15.0f0, 15.0f0, "30 ppm (±15)"),
    # Calibrated (asymmetric)
    (Pioneer.getLeftTol(mem_calibrated), Pioneer.getRightTol(mem_calibrated),
     "calibrated (~$(round(Int, Pioneer.getLeftTol(mem_calibrated) + Pioneer.getRightTol(mem_calibrated))) ppm)"),
    (20.0f0, 20.0f0, "40 ppm (±20)"),
    (25.0f0, 25.0f0, "50 ppm (±25)"),
]

function make_mem(left_tol::Float32, right_tol::Float32)
    Pioneer.MassErrorModel(mem_calibrated.mass_offset, (left_tol, right_tol))
end

# ── Build partitioned indices (one per PPM — shared across search tols) ───
indices = Dict{Float32, Any}()
for ppm in PPM_VALUES
    println("\nBuilding partitioned index (frag_bin_tol_ppm=$(ppm)) ...")
    t_build = @elapsed begin
        indices[ppm] = build_partitioned_index_from_lib(Pioneer.getSpecLib(SC);
            partition_width=5.0f0, frag_bin_tol_ppm=ppm, rt_bin_tol=1.0f0,
            rank_to_score=UInt8[8, 4, 4, 2, 2, 1, 1])
    end
    println("  Build: $(round(t_build, digits=1))s")
end

println("\n$(length(all_scans)) MS2 scans, $(Threads.nthreads()) threads")

# ── Helper: count precursor hits from a search run ────────────────────────
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
    # Baseline with calibrated mem
    local s2p = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra))
    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin
            tid = first(thread_task)
            Pioneer.searchFragmentIndex(s2p, native_frag_index, spectra, last(thread_task),
                search_data[tid], search_parameters, qtm, mem_calibrated, rt_to_irt_spline, irt_tol)
        end
    end
    fetch.(tasks)

    # Partitioned warmup (just one PPM, several tols)
    for mem_test in [mem_calibrated, make_mem(10.0f0, 10.0f0), make_mem(25.0f0, 25.0f0)]
        local s2p2 = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra))
        searchFragmentIndexPartitionMajorHinted(s2p2, indices[10.0f0], spectra, all_scans,
            Threads.nthreads(), search_parameters, qtm, mem_test, rt_to_irt_spline, irt_tol,
            precursor_mzs; linear_threshold=LINEAR_THRESHOLD)
    end
end

# ── Timed runs ────────────────────────────────────────────────────────────
println("\n" * "="^80)
println("TIMED BENCHMARK — PPM × Search Tolerance Cross-Sweep (threshold=$LINEAR_THRESHOLD)")
println("="^80)

# Baseline for each search tolerance
baseline_times = Dict{String, Float64}()
baseline_hits  = Dict{String, Int}()

for (l_tol, r_tol, label) in SEARCH_TOLS
    test_mem = make_mem(l_tol, r_tol)
    println("\nRunning baseline (native) with search tol=$label ...")
    local s2p = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra))
    t = @elapsed begin
        tasks = map(thread_tasks) do thread_task
            Threads.@spawn begin
                tid = first(thread_task)
                Pioneer.searchFragmentIndex(s2p, native_frag_index, spectra, last(thread_task),
                    search_data[tid], search_parameters, qtm, test_mem, rt_to_irt_spline, irt_tol)
            end
        end
        fetch.(tasks)
    end
    baseline_times[label] = t
    baseline_hits[label] = count_precursor_hits(s2p)
    println("  $(round(t, digits=2))s, $(baseline_hits[label]) precursor hits")
end

# Partitioned for each (search_tol × bin_ppm)
results = Dict{Tuple{String, Float32}, NamedTuple{(:time, :hits), Tuple{Float64, Int}}}()

for (l_tol, r_tol, label) in SEARCH_TOLS
    test_mem = make_mem(l_tol, r_tol)
    println("\n--- Search tolerance: $label ---")
    for ppm in PPM_VALUES
        println("  Running ppm=$ppm ...")
        local s2p = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra))
        t = @elapsed begin
            searchFragmentIndexPartitionMajorHinted(s2p, indices[ppm], spectra, all_scans,
                Threads.nthreads(), search_parameters, qtm, test_mem, rt_to_irt_spline, irt_tol,
                precursor_mzs; linear_threshold=LINEAR_THRESHOLD)
        end
        hits = count_precursor_hits(s2p)
        results[(label, ppm)] = (time=t, hits=hits)
        println("    $(round(t, digits=2))s, $hits hits")
    end
end

# ── Results table ─────────────────────────────────────────────────────────
println("\n" * "="^80)
println("RESULTS — PPM × Search Tolerance Cross-Sweep")
println("="^80)

for (l_tol, r_tol, label) in SEARCH_TOLS
    bt = baseline_times[label]
    bh = baseline_hits[label]
    println("\n  Search tolerance: $label")
    println("  Baseline: $(round(bt, digits=2))s, $bh hits")
    println("  Bin PPM |  Time (s) | Speedup | Hits       | Hit ratio vs baseline")
    println("  --------|-----------|---------|------------|----------------------")
    for ppm in PPM_VALUES
        r = results[(label, ppm)]
        speedup = bt / r.time
        hit_ratio = bh > 0 ? r.hits / bh : NaN
        best = (r.time == minimum(results[(label, p)].time for p in PPM_VALUES)) ? " *" : ""
        println("  $(lpad(ppm, 7)) | $(lpad(round(r.time, digits=3), 9)) | $(lpad(round(speedup, digits=2), 7))x | $(lpad(r.hits, 10)) | $(lpad(round(hit_ratio, digits=4), 20))$best")
    end
end

# ── Summary: best PPM per tolerance ───────────────────────────────────────
println("\n" * "="^80)
println("SUMMARY — Best bin PPM for each search tolerance")
println("="^80)
println("  Search Tol              | Best PPM | Time (s) | Speedup | Baseline (s)")
println("  ------------------------|----------|----------|---------|-------------")
for (l_tol, r_tol, label) in SEARCH_TOLS
    bt = baseline_times[label]
    best_ppm = PPM_VALUES[argmin([results[(label, p)].time for p in PPM_VALUES])]
    best_r = results[(label, best_ppm)]
    speedup = bt / best_r.time
    println("  $(rpad(label, 24))| $(lpad(best_ppm, 8)) | $(lpad(round(best_r.time, digits=2), 8)) | $(lpad(round(speedup, digits=2), 7))x | $(lpad(round(bt, digits=2), 11))")
end
println("="^80)
