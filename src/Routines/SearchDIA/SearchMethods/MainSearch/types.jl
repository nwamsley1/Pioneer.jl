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

# Deconvolution regularization from config (optimization.deconvolution.{lambda,reg_type}).
# Defaults to lambda=0 / NoNorm (unregularized PoissonMM) so existing configs are unchanged.
# reg_type string: "none"->NoNorm, "l1"->L1Norm, "l2"->L2Norm. Enables PMM L1/L2 penalty sweeps.
function _resolve_deconv_reg(params::PioneerParameters)
    lam = 0.0f0
    reg_str = "none"
    if hasproperty(params, :optimization) && hasproperty(params.optimization, :deconvolution)
        dec = params.optimization.deconvolution
        hasproperty(dec, :lambda)   && (lam = Float32(dec.lambda))
        hasproperty(dec, :reg_type) && (reg_str = lowercase(String(dec.reg_type)))
    end
    reg = reg_str in ("l1", "l1norm") ? L1Norm() :
          reg_str in ("l2", "l2norm") ? L2Norm() : NoNorm()
    return lam, reg
end

# Phospho-localization isomer-competition feature toggle
# (optimization.localization.enabled). Default false so non-phospho runs pay
# nothing (the per-scan isomer sweep is skipped entirely). See
# add_isomer_competition_features! (Idea-1, Phase A).
function _resolve_localization_enabled(params::PioneerParameters)
    if hasproperty(params, :optimization) && hasproperty(params.optimization, :localization)
        loc = params.optimization.localization
        hasproperty(loc, :enabled) && return Bool(loc.enabled)
    end
    return false
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

    # Phospho-localization: compute per-scan isomer-competition features
    localization_enabled::Bool

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

        # Deconvolution regularization (config-driven; default lambda=0/NoNorm)
        deconv_lambda, deconv_reg = _resolve_deconv_reg(params)

        new{typeof(prec_estimation), typeof(isotope_trace_type)}(
            isotope_bounds,
            Float32(min_fraction_transmitted),
            Int64(n_isotopes_val),
            max_frag_rank,
            Set{Int64}([2]),

            deconv_lambda,    # lambda (config: optimization.deconvolution.lambda)
            deconv_reg,       # reg_type (config: optimization.deconvolution.reg_type)
            PoissonMMSolver(deconv_lambda, deconv_reg),  # PMM with optional L1/L2 penalty
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

            _resolve_localization_enabled(params),  # localization_enabled (config: optimization.localization.enabled)
        )
    end
end

getIsotopeTraceType(p::MainSearchParameters) = p.isotope_tracetype
getPrefilterMinScanCount(p::MainSearchParameters) = p.prefilter_min_scan_count
getLocalizationEnabled(p::MainSearchParameters) = p.localization_enabled
