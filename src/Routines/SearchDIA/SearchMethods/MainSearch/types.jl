"""
    DEFAULT_INDEX_SEARCH_MIN_SCORE

CountFilter fallback used by `searchFragmentIndexPartitionMajorHinted` only
when BitVecCalibration hasn't written a learned per-file LUT for this MS
file. Every production run goes through BitVec which sets a per-file LUT,
so this fallback rarely (if ever) fires. Formerly user-tunable as
`fragment_index_search.min_score`.
"""
const DEFAULT_INDEX_SEARCH_MIN_SCORE = UInt8(15)

"""
    _resolve_n_isotopes(search_section) -> Int

Return `search.n_isotopes` (preferred) or fall back to the legacy nested
`search.fragment_settings.n_isotopes` location for old configs. Defaults
to 2 if neither is present (matches `defaultSearchParams.json`).
"""
function _resolve_n_isotopes(search_section)
    if hasproperty(search_section, :n_isotopes)
        return Int(search_section.n_isotopes)
    elseif hasproperty(search_section, :fragment_settings) &&
           hasproperty(search_section.fragment_settings, :n_isotopes)
        return Int(search_section.fragment_settings.n_isotopes)
    else
        return 2
    end
end

"""
    _resolve_q_value_threshold(global_section) -> Float32

Return `global.q_value_threshold` (preferred) or fall back to the legacy
nested `global.scoring.q_value_threshold` location for old configs.
Defaults to 0.01 if neither is present (matches `defaultSearchParams.json`).
"""
function _resolve_q_value_threshold(global_section)
    if hasproperty(global_section, :q_value_threshold)
        return Float32(global_section.q_value_threshold)
    elseif hasproperty(global_section, :scoring) &&
           hasproperty(global_section.scoring, :q_value_threshold)
        return Float32(global_section.scoring.q_value_threshold)
    else
        return 0.01f0
    end
end

"""
Parameters for main search (deconvolution, scoring, and prescore aggregation).
"""
struct MainSearchParameters{P<:PrecEstimation, I<:IsotopeTraceType} <: FragmentIndexSearchParameters
    # Core parameters
    isotope_err_bounds::Tuple{UInt8, UInt8}
    min_fraction_transmitted::Float32
    n_frag_isotopes::Int64
    max_frag_rank::UInt8
    spec_order::Set{Int64}

    # Deconvolution parameters (MS2)
    lambda::Float32
    reg_type::RegularizationType
    deconvolution_solver::DeconvolutionSolver
    max_iter_outer::Int64
    max_diff::Float32

    # PSM filtering
    min_y_count::Int64
    min_spectral_contrast::Float32
    min_log2_matched_ratio::Float32
    min_topn_of_m::Tuple{Int64, Int64}
    max_best_rank::Int64

    # Precursor estimation strategy
    isotope_tracetype::I
    prec_estimation::P

    # Prescore fragment settings (may differ from full deconvolution)
    prescore_n_frag_isotopes::Int64
    prescore_max_frag_rank::UInt8
    # Fragment index search parameter (used when MainSearch runs fragment index inline)
    min_index_search_score::UInt8

    # Pre-filter: require marginal candidates to appear in ≥ N scans
    prefilter_min_scan_count::Int64

    function MainSearchParameters(params::PioneerParameters)
        # Extract relevant parameter groups
        quant_params = params.search

        # Hardcoded constants — formerly user-tunable via global.isotope_settings
        # but every shipping config used the same values. The combine-traces
        # branch was already a no-op (both arms produced SeperateTraces).
        isotope_trace_type = SeperateTraces()
        isotope_bounds = (UInt8(1), UInt8(0))
        min_fraction_transmitted = 0.25f0
        prec_estimation = PartialPrecCapture()

        # max_rank = typemax(UInt8) → no runtime rank cap. The library already
        # enforces a per-precursor rank cap at build time (BuildSpecLib's
        # max_frag_rank=10), so a separate runtime knob is redundant.
        max_frag_rank = typemax(UInt8)

        # n_isotopes lives at search.n_isotopes (flattened). Fall back to
        # the legacy nested location search.fragment_settings.n_isotopes
        # so old configs keep working.
        n_isotopes_val = _resolve_n_isotopes(quant_params)

        new{typeof(prec_estimation), typeof(isotope_trace_type)}(
            isotope_bounds,
            Float32(min_fraction_transmitted),
            Int64(n_isotopes_val),
            max_frag_rank,
            Set{Int64}([2]),

            Float32(0.0),     # lambda (no regularization)
            NoNorm(),         # reg_type
            PoissonMMSolver(),  # OLS / Lasso / AdaptiveLasso paths retained in git history
            DECONV_MAX_ITER,          # max_iter_outer
            DECONV_CONVERGENCE_TOL,   # max_diff

            Int64(0),   # min_y_count hardcoded
            Float32(0), # min_spectral_contrast hardcoded
            Float32(-10),  # min_log2_matched_ratio hardcoded
            (0, 3),     # min_topn_of_m hardcoded
            Int64(max_frag_rank),

            isotope_trace_type,
            prec_estimation,

            Int64(n_isotopes_val),    # prescore_n_frag_isotopes
            max_frag_rank,             # prescore_max_frag_rank

            # min_index_search_score is the CountFilter fallback used only
            # when BitVecCalibration didn't write a per-file LUT (every
            # production run goes through BitVec, so this falls back rarely).
            DEFAULT_INDEX_SEARCH_MIN_SCORE,

            0,                  # prefilter_min_scan_count (formerly fragment_index_search.prefilter_min_scan_count; never overridden)
        )
    end
end

getIsotopeTraceType(p::MainSearchParameters) = p.isotope_tracetype
getPrefilterMinScanCount(p::MainSearchParameters) = p.prefilter_min_scan_count
