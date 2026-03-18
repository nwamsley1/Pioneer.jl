#!/usr/bin/env julia
#
# Profile the fragment index search (Phase 1 of LibrarySearch) during FirstPassSearch.
#
# Usage:
#   julia --threads=auto --project=. src/Routines/SearchDIA/profile_fragment_index.jl /path/to/params.json
#
# This runs the calibration pipeline (ParameterTuning, NceTuning, QuadTuning),
# then profiles ONLY the `searchFragmentIndex` calls using the full fragment index
# (not the smaller presearch index).

using Pioneer
using Profile
using PProf

# ── Parse CLI argument ──────────────────────────────────────────────────────
if isempty(ARGS)
    error("Usage: julia --threads=auto --project=. profile_fragment_index.jl <params.json>")
end
const PARAMS_PATH = ARGS[1]
isfile(PARAMS_PATH) || error("Params file not found: $PARAMS_PATH")

# ── Load parameters ─────────────────────────────────────────────────────────
println("Loading parameters from $PARAMS_PATH ...")
Pioneer.checkParams(PARAMS_PATH)
params = Pioneer.parse_pioneer_parameters(PARAMS_PATH)
mkpath(params.paths[:results])

# ── Resolve data paths (mirrors SearchDIA.jl) ───────────────────────────────
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

# ── Load spectral library ───────────────────────────────────────────────────
println("Loading spectral library ...")
SPEC_LIB = Pioneer.loadSpectralLibrary(SPEC_LIB_DIR, params)

# ── Initialize search context ───────────────────────────────────────────────
println("Initializing search context ...")
SEARCH_CONTEXT = Pioneer.initSearchContext(
    SPEC_LIB,
    Pioneer.parseIsoXML(Pioneer.isotope_spline_path()),
    Pioneer.ArrowTableReference(MS_TABLE_PATHS),
    Threads.nthreads(),
    250000
)
Pioneer.setDataOutDir!(SEARCH_CONTEXT, params.paths[:results])

# ── Run calibration steps (required before FirstPassSearch) ──────────────────
println("Running ParameterTuningSearch ...")
Pioneer.execute_search(Pioneer.ParameterTuningSearch(), SEARCH_CONTEXT, params)

println("Running NceTuningSearch ...")
Pioneer.execute_search(Pioneer.NceTuningSearch(), SEARCH_CONTEXT, params)

println("Running QuadTuningSearch ...")
Pioneer.execute_search(Pioneer.QuadTuningSearch(), SEARCH_CONTEXT, params)

# ── Prepare FirstPass search arguments ───────────────────────────────────────
println("Preparing fragment index search arguments ...")

search_parameters = Pioneer.FirstPassSearchParameters(params)
fragment_index     = Pioneer.getFragmentIndex(Pioneer.getSpecLib(SEARCH_CONTEXT))
search_data        = Pioneer.getSearchData(SEARCH_CONTEXT)

# Load the first MS file's spectra
msdr = Pioneer.getMassSpecData(SEARCH_CONTEXT)
ms_file_idx, spectra = first(enumerate(msdr))

qtm              = Pioneer.getQuadTransmissionModel(SEARCH_CONTEXT, ms_file_idx)
mem              = Pioneer.getMassErrorModel(SEARCH_CONTEXT, ms_file_idx)
rt_to_irt_spline = Pioneer.getRtIrtModel(SEARCH_CONTEXT, ms_file_idx)
irt_tol          = Pioneer.getIrtErrors(SEARCH_CONTEXT)[ms_file_idx]

thread_tasks = Pioneer.partition_scans(spectra, Threads.nthreads())

# ── Materialize Arrow arrays to native Julia vectors ─────────────────────────
println("Materializing fragment index to native vectors ...")
materialize_time = @elapsed begin
    native_frag_index = Pioneer.materialize(fragment_index)
end
n_frag_bins = length(native_frag_index.fragment_bins)
n_rt_bins   = length(native_frag_index.rt_bins)
n_frags     = length(native_frag_index.fragments)
mem_mb = (sizeof(native_frag_index.fragment_bins) +
          sizeof(native_frag_index.rt_bins) +
          sizeof(native_frag_index.fragments)) / 1024^2
println("  Materialization: $(round(materialize_time, digits=3))s")
println("  fragment_bins: $n_frag_bins  rt_bins: $n_rt_bins  fragments: $n_frags")
println("  Native memory: $(round(mem_mb, digits=1)) MB")

# ── Warmup (compile everything, exclude JIT from profile) ────────────────────
println("Warmup run (native vectors) ...")
let
    scan_to_prec_idx = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra))
    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin
            thread_id = first(thread_task)
            Pioneer.searchFragmentIndex(
                scan_to_prec_idx,
                native_frag_index,
                spectra,
                last(thread_task),
                search_data[thread_id],
                search_parameters,
                qtm,
                mem,
                rt_to_irt_spline,
                irt_tol
            )
        end
    end
    fetch.(tasks)
end

# ── Benchmark: timed runs for comparison ─────────────────────────────────────
println("\nBenchmarking searchFragmentIndex ...")

# Native vectors
native_time = @elapsed begin
    scan_to_prec_idx_native = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra))
    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin
            thread_id = first(thread_task)
            Pioneer.searchFragmentIndex(
                scan_to_prec_idx_native,
                native_frag_index,
                spectra,
                last(thread_task),
                search_data[thread_id],
                search_parameters,
                qtm,
                mem,
                rt_to_irt_spline,
                irt_tol
            )
        end
    end
    fetch.(tasks)
end

# Arrow (original)
arrow_time = @elapsed begin
    scan_to_prec_idx_arrow = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra))
    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin
            thread_id = first(thread_task)
            Pioneer.searchFragmentIndex(
                scan_to_prec_idx_arrow,
                fragment_index,
                spectra,
                last(thread_task),
                search_data[thread_id],
                search_parameters,
                qtm,
                mem,
                rt_to_irt_spline,
                irt_tol
            )
        end
    end
    fetch.(tasks)
end

println("  Arrow (original): $(round(arrow_time, digits=3))s")
println("  Native vectors:   $(round(native_time, digits=3))s")
println("  Speedup:          $(round(arrow_time / native_time, digits=2))x")
println("  Materialize cost: $(round(materialize_time, digits=3))s ($(round(materialize_time/native_time*100, digits=1))% of one search)")

# ── Profile Phase 1 (searchFragmentIndex) with native vectors ────────────────
println("\nProfiling searchFragmentIndex with native vectors ($(Threads.nthreads()) threads) ...")
Profile.clear()
scan_to_prec_idx = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra))
@profile begin
    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin
            thread_id = first(thread_task)
            Pioneer.searchFragmentIndex(
                scan_to_prec_idx,
                native_frag_index,
                spectra,
                last(thread_task),
                search_data[thread_id],
                search_parameters,
                qtm,
                mem,
                rt_to_irt_spline,
                irt_tol
            )
        end
    end
    fetch.(tasks)
end

# ── Print flat profile summary ───────────────────────────────────────────────
println("\n", "="^80)
println("FLAT PROFILE — NATIVE VECTORS (top 40 by count)")
println("="^80)
Profile.print(IOContext(stdout, :displaysize => (200, 200)), maxdepth=40, noisefloor=1.0, sortedby=:count)

# ── Save profile & generate flame graph ──────────────────────────────────────
prof_path = joinpath(dirname(PARAMS_PATH), "fragment_index_profile_native.pb.gz")
println("\nSaving profile to $prof_path ...")
pprof(out=prof_path, web=false)
println("Done. Open with:  go tool pprof -http=:8080 $prof_path")
