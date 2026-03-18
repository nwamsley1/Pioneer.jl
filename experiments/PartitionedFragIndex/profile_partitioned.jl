#!/usr/bin/env julia
#
# Profile the partitioned fragment index search (7 frags, weighted).
#
# Usage:
#   julia --threads=auto --project=. experiments/PartitionedFragIndex/profile_partitioned.jl /path/to/params.json
#

using Pioneer
using Profile
using PProf

# ── Load experiment modules ──────────────────────────────────────────────────
include(joinpath(@__DIR__, "partitioned_types.jl"))
include(joinpath(@__DIR__, "build_partitioned_index.jl"))
include(joinpath(@__DIR__, "search_partitioned_index.jl"))

# ── Parse CLI argument ──────────────────────────────────────────────────────
if isempty(ARGS)
    error("Usage: julia --threads=auto --project=. experiments/PartitionedFragIndex/profile_partitioned.jl <params.json>")
end
const PARAMS_PATH = ARGS[1]
isfile(PARAMS_PATH) || error("Params file not found: $PARAMS_PATH")

# ── Load parameters & calibrate ──────────────────────────────────────────────
println("Loading parameters from $PARAMS_PATH ...")
Pioneer.checkParams(PARAMS_PATH)
params = Pioneer.parse_pioneer_parameters(PARAMS_PATH)
mkpath(params.paths[:results])

MS_DATA_DIR = params.paths[:ms_data]
!isabspath(MS_DATA_DIR) && (MS_DATA_DIR = joinpath(@__DIR__, "../../", MS_DATA_DIR))
MS_TABLE_PATHS = [joinpath(MS_DATA_DIR, f)
                  for f in readdir(MS_DATA_DIR)
                  if isfile(joinpath(MS_DATA_DIR, f)) && endswith(f, ".arrow")]

SPEC_LIB_DIR = params.paths[:library]
!isabspath(SPEC_LIB_DIR) && (SPEC_LIB_DIR = joinpath(@__DIR__, "../../", SPEC_LIB_DIR))

println("Loading spectral library ...")
SPEC_LIB = Pioneer.loadSpectralLibrary(SPEC_LIB_DIR, params)

println("Initializing search context ...")
SEARCH_CONTEXT = Pioneer.initSearchContext(
    SPEC_LIB,
    Pioneer.parseIsoXML(Pioneer.isotope_spline_path()),
    Pioneer.ArrowTableReference(MS_TABLE_PATHS),
    Threads.nthreads(),
    250000
)
Pioneer.setDataOutDir!(SEARCH_CONTEXT, params.paths[:results])

println("Running calibration ...")
Pioneer.execute_search(Pioneer.ParameterTuningSearch(), SEARCH_CONTEXT, params)
Pioneer.execute_search(Pioneer.NceTuningSearch(), SEARCH_CONTEXT, params)
Pioneer.execute_search(Pioneer.QuadTuningSearch(), SEARCH_CONTEXT, params)

# ── Prepare search arguments ────────────────────────────────────────────────
println("\nPreparing search arguments ...")
search_parameters = Pioneer.FirstPassSearchParameters(params)
search_data       = Pioneer.getSearchData(SEARCH_CONTEXT)
msdr = Pioneer.getMassSpecData(SEARCH_CONTEXT)
ms_file_idx, spectra = first(enumerate(msdr))
qtm              = Pioneer.getQuadTransmissionModel(SEARCH_CONTEXT, ms_file_idx)
mem              = Pioneer.getMassErrorModel(SEARCH_CONTEXT, ms_file_idx)
rt_to_irt_spline = Pioneer.getRtIrtModel(SEARCH_CONTEXT, ms_file_idx)
irt_tol          = Pioneer.getIrtErrors(SEARCH_CONTEXT)[ms_file_idx]
thread_tasks     = Pioneer.partition_scans(spectra, Threads.nthreads())
precursor_mzs    = Pioneer.getMz(Pioneer.getPrecursors(Pioneer.getSpecLib(SEARCH_CONTEXT)))

# ── Build partitioned index (7 frags, weighted) ─────────────────────────────
println("\nBuilding partitioned index (7 frags, weighted) ...")
build_time = @elapsed begin
    partitioned_index = build_partitioned_index_from_lib(
        Pioneer.getSpecLib(SEARCH_CONTEXT);
        partition_width=5.0f0,
        frag_bin_tol_ppm=10.0f0,
        rt_bin_tol=1.0f0,
        rank_to_score=UInt8[8, 4, 4, 2, 2, 1, 1],
    )
end
println("  Build time: $(round(build_time, digits=3))s")

# ── Collect all valid MS2 scan indices ───────────────────────────────────────
all_ms2_scan_idxs = Int[]
for (_, task_scans) in thread_tasks
    for s in task_scans
        (s <= 0 || s > length(spectra)) && continue
        Pioneer.getMsOrder(spectra, s) ∉ Pioneer.getSpecOrder(search_parameters) && continue
        push!(all_ms2_scan_idxs, s)
    end
end
sort!(all_ms2_scan_idxs)
println("Total MS2 scans: $(length(all_ms2_scan_idxs))")

# ── Warmup ───────────────────────────────────────────────────────────────────
println("\nWarmup (hinted SoA+SIMD) ...")
let
    s2p = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra))
    searchFragmentIndexPartitionMajorHinted(
        s2p, partitioned_index, spectra, all_ms2_scan_idxs,
        Threads.nthreads(), search_parameters, qtm, mem, rt_to_irt_spline, irt_tol,
        precursor_mzs
    )
end

# ── Timed run ────────────────────────────────────────────────────────────────
println("Timed run ...")
timed = @elapsed begin
    s2p = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra))
    searchFragmentIndexPartitionMajorHinted(
        s2p, partitioned_index, spectra, all_ms2_scan_idxs,
        Threads.nthreads(), search_parameters, qtm, mem, rt_to_irt_spline, irt_tol,
        precursor_mzs
    )
end
println("  Search time: $(round(timed, digits=3))s")

# ── Profile ──────────────────────────────────────────────────────────────────
println("\nProfiling (3 iterations) ...")
Profile.clear()
@profile for _ in 1:3
    s2p = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra))
    searchFragmentIndexPartitionMajorHinted(
        s2p, partitioned_index, spectra, all_ms2_scan_idxs,
        Threads.nthreads(), search_parameters, qtm, mem, rt_to_irt_spline, irt_tol,
        precursor_mzs
    )
end

# ── Save profile ─────────────────────────────────────────────────────────────
out_path = joinpath(dirname(PARAMS_PATH), "partitioned_profile.pb.gz")
pprof(out=out_path, web=false)
println("\nProfile saved to $out_path")
println("View with: go tool pprof -http=:8080 $out_path")
