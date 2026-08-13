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
    converged = first(solve_deconvolution!(
        params.deconvolution_solver,
        Hs, getResiduals(search_data), weights, getColNorm2(search_data),
        getMu(search_data), getObserved(search_data),
        params.max_iter_outer, params.max_diff
    ))
    if converged
        update_precursor_weights!(getIdToCol(search_data), weights, getPrecursorWeights(search_data))
        zero_negligible_weights!(weights, Hs.n)
    end
    return converged
end

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
