# Shared types and dispatch helpers for the per-precursor scan loop.
#
# Defines `PrecursorIndex` / `PerScanPrecursorIndex` (used by `process_scans_fused!`),
# `get_irt_tolerance`, `resize_if_needed!`, `post_design_matrix!`,
# `compute_distance_metrics!`, and `score_psms!` — each dispatches on the
# concrete `SearchParameters` subtype to vary behavior across MainSearch,
# ParameterTuning, and QuadTuning.
#
# The historical top-level `process_scans!` function (the classic multi-pass
# pipeline) is gone; all live search methods route through
# `process_scans_fused!` in `process_scans_fused.jl`.

#==========================================================
Precursor index types — control transition list reuse via dispatch
==========================================================#

"""
    PrecursorIndex

Abstract type for precursor index variants. Wraps `scan_to_prec_idx` and
`precursors_passed` with dispatch-based behavior for transition list reuse.
"""
abstract type PrecursorIndex end

"""
    PerScanPrecursorIndex

Standard per-scan precursor ranges. Used by tuning stages.
`selectTransitions!` is called for every scan.
"""
struct PerScanPrecursorIndex <: PrecursorIndex
    scan_to_prec_idx::Vector{Union{Missing, UnitRange{Int64}}}
    precursors_passed::Vector{UInt32}
end


# Accessor interface
@inline get_prec_range(pi::PrecursorIndex, scan_idx::Int) = pi.scan_to_prec_idx[scan_idx]
@inline get_precursors(pi::PrecursorIndex) = pi.precursors_passed

# Always rebuild transitions for per-scan index
@inline should_rebuild_transitions(::PerScanPrecursorIndex, ::Int, ::Float32, ::Float32) = true

#==========================================================
Dispatch helpers for library_search (LibrarySearch.jl)
==========================================================#
# These three helpers let a single `library_search` function serve every
# FragmentIndexSearchParameters subtype. Each returns a different value for
# ParameterTuning vs the default path:
#
#   get_fragment_index   — which partitioned fragment index to search against
#   get_irt_tolerance    — where to read the iRT tolerance from
#   get_nce_models       — how many NCE models to iterate (1, or N for NCE grid)

# ParameterTuning searches the smaller presearch index; everything else uses the full index.
get_fragment_index(spec_lib::SpectralLibrary, ::ParameterTuningSearchParameters) = getPresearchPartitionedIndex(spec_lib)
get_fragment_index(spec_lib::SpectralLibrary, ::FragmentIndexSearchParameters) = getPartitionedIndex(spec_lib)

# ParameterTuning carries its own irt_tol on the params struct (set before any
# per-file calibration exists). All later methods read the calibrated value from
# SearchContext, falling back to Inf if the file hasn't been calibrated yet.
get_irt_tolerance(ctx::SearchContext, params::ParameterTuningSearchParameters, ::Int64) = getIRTTol(params)
function get_irt_tolerance(ctx::SearchContext, ::FragmentIndexSearchParameters, ms_file_idx::Int64)
    haskey(getIrtErrors(ctx), ms_file_idx) ? getIrtErrors(ctx)[ms_file_idx] : Float32(Inf)
end

# Returns a vector of (nce_model, tag) pairs that library_search iterates over.
# Most methods return a single calibrated model (tag=nothing → no :nce column).
get_nce_models(ctx::SearchContext, ::FragmentIndexSearchParameters, ms_file_idx::Int64) =
    [(getNceModel(ctx, ms_file_idx), nothing)]

#==========================================================
Dispatched array accessors
==========================================================#

get_scored_psms(sd::SearchDataStructures, ::MainSearchParameters) = getMainSearchScoredPsms(sd)
get_scored_psms(sd::SearchDataStructures, ::ParameterTuningSearchParameters) = getTuningScoredPsms(sd)
get_unscored_psms(sd::SearchDataStructures, ::MainSearchParameters) = getMainUnscoredPsms(sd)
get_unscored_psms(sd::SearchDataStructures, ::ParameterTuningSearchParameters) = getTuningUnscoredPsms(sd)
# QuadTuningSearchParameters dispatches live in QuadTuningSearch/utils.jl
# (loaded after this file — QuadTuning types aren't available yet here).

#==========================================================
Dispatched helpers
==========================================================#

"""
    resize_if_needed!(search_data, params)

Grow pre-allocated arrays when the number of active precursors exceeds capacity.
Dispatches on params type to resize the correct set of arrays.
"""
function resize_if_needed!(search_data::SearchDataStructures, params::MainSearchParameters)
    weights = getTempWeights(search_data)
    if n_active(getIdToCol(search_data)) > length(weights)
        resize_arrays!(search_data, weights)
    end
end

function resize_if_needed!(search_data::SearchDataStructures, params::ParameterTuningSearchParameters)
    weights = getTempWeights(search_data)
    if n_active(getIdToCol(search_data)) > length(weights)
        resize_arrays!(search_data, weights)
    end
end



"""
    post_design_matrix!(search_data, Hs, params) -> Bool

Post-design-matrix processing. Returns `true` if scoring should proceed.
- Simple path: no-op, always returns true.
- MainSearch: initialize weights, solve deconvolution, update precursor weights.
"""
function post_design_matrix!(search_data::SearchDataStructures, Hs::AbstractSparseDesignMatrix, params::MainSearchParameters)
    weights = getTempWeights(search_data)
    initialize_weights!(getIdToCol(search_data), weights, getPrecursorWeights(search_data))
    # DIAGNOSTIC (PIONEER_PMM_DUMP_DIR, PIONEER_PMM_DUMP_EVERY): serialize every n-th solve's
    # inputs (design matrix, observed y, warm-start weights) as a solver test problem, in the
    # Dict format test/UnitTests/test_poissonMM.jl reads. Inert unless set.
    if PMM_DUMP_EVERY[] > 0
        _c = Threads.atomic_add!(PMM_DUMP_COUNTER, 1)
        _c % PMM_DUMP_EVERY[] == 0 && _pmm_dump_problem(Hs, search_data, weights, _c)
    end
    converged, n_iter = solve_deconvolution!(
        params.deconvolution_solver,
        Hs, getResiduals(search_data), weights, getColNorm2(search_data),
        getMu(search_data), getObserved(search_data),
        params.max_iter_outer, params.max_diff
    )
    # Diagnostic tally (main search): solves and outer iterations, reported per library_search.
    Threads.atomic_add!(DECONV_SOLVES, 1)
    Threads.atomic_add!(DECONV_ITERS, Int(n_iter))
    if converged
        update_precursor_weights!(getIdToCol(search_data), weights, getPrecursorWeights(search_data))
        zero_negligible_weights!(weights, Hs.n)
    end
    return converged
end

"""Main-search deconvolution tallies (solves, outer iterations); reset and logged by library_search."""
const DECONV_SOLVES = Threads.Atomic{Int}(0)
const DECONV_ITERS  = Threads.Atomic{Int}(0)

function post_design_matrix!(search_data::SearchDataStructures, Hs::AbstractSparseDesignMatrix, params::ParameterTuningSearchParameters)
    weights = getTempWeights(search_data)
    initialize_weights!(getIdToCol(search_data), weights, getPrecursorWeights(search_data))
    converged = first(solve_deconvolution!(
        OLSSolver(),
        Hs, getResiduals(search_data), weights, getColNorm2(search_data),
        getMu(search_data), getObserved(search_data),
        DECONV_MAX_ITER, DECONV_CONVERGENCE_TOL
    ))
    if converged
        update_precursor_weights!(getIdToCol(search_data), weights, getPrecursorWeights(search_data))
        zero_negligible_weights!(weights, Hs.n)
    end
    return converged
end


"""
    compute_distance_metrics!(Hs, search_data, params)

Compute spectral distance metrics between observed and predicted spectra.
- Simple path: iterative peak removal with relative improvement threshold.
- MainSearch: single-pass metrics using deconvolution weights and residuals.
"""
function compute_distance_metrics!(Hs::AbstractSparseDesignMatrix, search_data::SearchDataStructures, params::MainSearchParameters)
    getDistanceMetrics(getTempWeights(search_data), getResiduals(search_data),
        Hs, getMainSearchSpectralScores(search_data))
end

function compute_distance_metrics!(Hs::AbstractSparseDesignMatrix, search_data::SearchDataStructures, params::ParameterTuningSearchParameters)
    getDistanceMetrics(getTempWeights(search_data), getResiduals(search_data),
        Hs, getMainSearchSpectralScores(search_data))
end



"""
    score_psms!(search_data, params, Hs, scan_idx, nmatches, nmisses,
                spectra, last_val, cycle_idx) -> Int64

Score PSMs for the current scan. Returns updated `last_val`.
- Simple path: Score! with spectral contrast, matched ratio, and rank filters.
- MainSearch: Score! with deconvolution weights and cycle index.
"""
function score_psms!(
    search_data::SearchDataStructures,
    params::MainSearchParameters,
    Hs::AbstractSparseDesignMatrix,
    scan_idx::Int64,
    nmatches::Int64,
    nmisses::Int64,
    spectra::MassSpecData,
    last_val::Int64,
    ms_file_idx::Int64,
    cycle_idx::Int64;
    mem::AbstractMassErrorModel = SimpleMassErrorModel(0f0, (0f0, 0f0))
)
    score_result = Score!(
        getMainSearchScoredPsms(search_data),
        getMainUnscoredPsms(search_data),
        getMainSearchSpectralScores(search_data),
        getTempWeights(search_data),
        getIdToCol(search_data),
        ms_file_idx,
        cycle_idx,
        nmatches / (nmatches + nmisses),
        last_val,
        Hs.n,
        Float32(sum(getIntensityArray(spectra, scan_idx))),
        scan_idx;
        block_size = 500000,
        default_top3_ll = get_default_top3_ll(mem)
    )
    return score_result.last_val
end

function score_psms!(
    search_data::SearchDataStructures,
    params::ParameterTuningSearchParameters,
    Hs::AbstractSparseDesignMatrix,
    scan_idx::Int64,
    nmatches::Int64,
    nmisses::Int64,
    spectra::MassSpecData,
    last_val::Int64,
    ms_file_idx::Int64,
    cycle_idx::Int64;
    mem::AbstractMassErrorModel = SimpleMassErrorModel(0f0, (0f0, 0f0))
)
    score_result = Score!(
        getTuningScoredPsms(search_data),
        getTuningUnscoredPsms(search_data),
        getMainSearchSpectralScores(search_data),
        getTempWeights(search_data),
        getIdToCol(search_data),
        ms_file_idx,
        cycle_idx,
        nmatches / (nmatches + nmisses),
        last_val,
        Hs.n,
        Float32(sum(getIntensityArray(spectra, scan_idx))),
        scan_idx;
        block_size = 500000,
        default_top3_ll = get_default_top3_ll(mem)
    )
    return score_result.last_val
end


const PMM_DUMP_EVERY = Ref(0)
const PMM_DUMP_COUNTER = Threads.Atomic{Int}(0)
const PMM_DUMP_DIR = Ref("")
function _pmm_dump_problem(Hs, search_data, weights, counter)
    n_rows = Int(Hs.m); n_cols = Int(Hs.n); n_vals = Int(Hs.n_vals)
    n_cols == 0 && return
    # observed intensities: rows are peaks; Hs.x holds the observed value per nonzero entry
    y = zeros(Float32, n_rows)
    @inbounds for i in 1:n_vals
        r = Int(Hs.rowval[i]); r >= 1 && r <= n_rows && (y[r] = Hs.x[i])
    end
    d = Dict{Symbol,Any}(:n_rows => n_rows, :n_cols => n_cols, :n_vals => n_vals,
        :colptr => Int64.(Hs.colptr[1:n_cols + 1]), :rowval => Int64.(Hs.rowval[1:n_vals]),
        :nzval => Float32.(Hs.nzval[1:n_vals]), :y => y,
        :w_before_all => Float32.(weights[1:n_cols]), :counter => counter)
    path = joinpath(PMM_DUMP_DIR[], "pmm_problem_$(lpad(counter, 8, '0')).dat")
    open(io -> serialize(io, d), path, "w")
    return nothing
end
