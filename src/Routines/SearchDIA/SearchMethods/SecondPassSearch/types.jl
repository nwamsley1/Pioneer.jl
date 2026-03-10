"""
    SecondPassSearch

Second pass search: full-feature deconvolution on globally-filtered precursors.
"""
struct SecondPassSearch <: SearchMethod end

#==========================================================
Type Definitions
==========================================================#

"""
Results container for second pass search.
"""
struct SecondPassSearchResults <: SearchResults
    psms::Base.Ref{DataFrame}
    ms1_psms::Base.Ref{DataFrame}
    file_fwhms::Dict{Int, @NamedTuple{median_fwhm::Float32, mad_fwhm::Float32}}
end

"""
Parameters for second pass search.
"""
struct SecondPassSearchParameters{P<:PrecEstimation, I<:IsotopeTraceType, A<:PrescoreAggregationStrategy} <: FragmentIndexSearchParameters
    # Core parameters
    isotope_err_bounds::Tuple{UInt8, UInt8}
    min_fraction_transmitted::Float32
    n_frag_isotopes::Int64
    max_frag_rank::UInt8
    spec_order::Set{Int64}
    match_between_runs::Bool

    # Deconvolution parameters (MS2)
    lambda::Float32
    reg_type::RegularizationType
    max_iter_newton::Int64
    max_iter_bisection::Int64
    max_iter_outer::Int64
    accuracy_newton::Float32
    accuracy_bisection::Float32
    max_diff::Float32

    # MS1 deconvolution parameters
    ms1_lambda::Float32
    ms1_reg_type::RegularizationType
    ms1_huber_delta::Float32

    # PSM filtering
    min_y_count::Int64
    min_frag_count::Int64
    min_spectral_contrast::Float32
    min_log2_matched_ratio::Float32
    min_topn_of_m::Tuple{Int64, Int64}
    max_best_rank::Int64

    # Precursor estimation strategy
    isotope_tracetype::I
    prec_estimation::P

    # Collect MS1 data?
    ms1_scoring::Bool

    # Global prescore q-value threshold
    global_prescore_qvalue_threshold::Float32

    # Phase 1 prescore fragment settings (may differ from Phase 2)
    prescore_n_frag_isotopes::Int64
    prescore_max_frag_rank::UInt8

    # Prescore aggregation strategy
    prescore_aggregation::A

    function SecondPassSearchParameters(params::PioneerParameters)
        # Extract relevant parameter groups
        global_params = params.global_settings
        quant_params = params.quant_search
        frag_params = quant_params.fragment_settings
        deconv_params = params.optimization.deconvolution

        # Determine isotope trace type based on global settings
        isotope_trace_type = if haskey(global_params.isotope_settings, :combine_traces) &&
                               global_params.isotope_settings.combine_traces
            SeperateTraces() #CombineTraces(0.0f0)  # Default min_fraction_transmitted
        else
            SeperateTraces()
        end

        isotope_bounds = global_params.isotope_settings.err_bounds_quant_search
        min_fraction_transmitted = global_params.isotope_settings.min_fraction_transmitted
        prec_estimation = global_params.isotope_settings.partial_capture ? PartialPrecCapture() : FullPrecCapture()

        # Parse MS2 regularization type
        reg_type = deconv_params.ms2.reg_type
        if reg_type == "none"
            reg_type = NoNorm()
        elseif reg_type == "l1"
            reg_type = L1Norm()
        elseif reg_type == "l2"
            reg_type = L2Norm()
        else
            reg_type = NoNorm()
            @user_warn "Warning. MS2 reg type `$reg_type` not recognized. Using NoNorm. Accepted types are `none`, `l1`, `l2`"
        end

        # Parse MS1 regularization type
        ms1_reg_type = deconv_params.ms1.reg_type
        if ms1_reg_type == "none"
            ms1_reg_type = NoNorm()
        elseif ms1_reg_type == "l1"
            ms1_reg_type = L1Norm()
        elseif ms1_reg_type == "l2"
            ms1_reg_type = L2Norm()
        else
            ms1_reg_type = NoNorm()
            @user_warn "Warning. MS1 reg type `$ms1_reg_type` not recognized. Using NoNorm. Accepted types are `none`, `l1`, `l2`"
        end

        ms1_scoring = Bool(global_params.ms1_scoring)

        # Global prescore q-value threshold (default 0.05 if not in JSON)
        global_prescore_qval = if haskey(quant_params, :global_prescore_qvalue_threshold)
            Float32(quant_params.global_prescore_qvalue_threshold)
        else
            0.05f0
        end

        # Phase 1 prescore fragment settings (fallback to main fragment_settings)
        prescore_n_frag_isotopes = if haskey(quant_params, :prescore_fragment_settings) &&
                                      haskey(quant_params.prescore_fragment_settings, :n_isotopes)
            Int64(quant_params.prescore_fragment_settings.n_isotopes)
        else
            Int64(frag_params.n_isotopes)
        end

        prescore_max_frag_rank = if haskey(quant_params, :prescore_fragment_settings) &&
                                    haskey(quant_params.prescore_fragment_settings, :max_rank)
            UInt8(quant_params.prescore_fragment_settings.max_rank)
        else
            UInt8(frag_params.max_rank)
        end

        # Prescore aggregation strategy (default: PEPCalibratedAggregation)
        prescore_aggregation = if haskey(quant_params, :prescore_aggregation)
            agg_str = string(quant_params.prescore_aggregation)
            if agg_str == "raw_logodds"
                RawLogOddsAggregation()
            else
                PEPCalibratedAggregation()
            end
        else
            PEPCalibratedAggregation()
        end

        new{typeof(prec_estimation), typeof(isotope_trace_type), typeof(prescore_aggregation)}(
            (UInt8(first(isotope_bounds)), UInt8(last(isotope_bounds))),
            Float32(min_fraction_transmitted),
            Int64(frag_params.n_isotopes),
            UInt8(frag_params.max_rank),
            Set{Int64}([2]),
            Bool(global_params.match_between_runs),

            Float32(deconv_params.ms2.lambda),
            reg_type,
            Int64(deconv_params.newton_iters),
            Int64(deconv_params.bisection_iters),
            Int64(deconv_params.outer_iters),
            Float32(deconv_params.newton_accuracy),
            Float32(deconv_params.newton_accuracy),
            Float32(deconv_params.max_diff),

            Float32(deconv_params.ms1.lambda),
            ms1_reg_type,
            Float32(deconv_params.ms1.huber_delta),

            Int64(frag_params.min_y_count),
            Int64(frag_params.min_count),
            Float32(frag_params.min_spectral_contrast),
            Float32(frag_params.min_log2_ratio),
            (Int64(first(frag_params.min_top_n)), Int64(last(frag_params.min_top_n))),
            Int64(frag_params.max_rank),

            isotope_trace_type,
            prec_estimation,

            ms1_scoring,

            global_prescore_qval,

            prescore_n_frag_isotopes,
            prescore_max_frag_rank,

            prescore_aggregation
        )
    end
end

getIsotopeTraceType(p::SecondPassSearchParameters) = p.isotope_tracetype
