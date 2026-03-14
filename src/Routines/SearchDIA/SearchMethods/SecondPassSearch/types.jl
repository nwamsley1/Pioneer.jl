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
struct SecondPassSearchParameters{P<:PrecEstimation, I<:IsotopeTraceType} <: FragmentIndexSearchParameters
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
    prescore_min_frag_count::Int64
    prescore_min_spectral_contrast::Float32
    prescore_min_log2_matched_ratio::Float32
    prescore_min_topn_of_m::Tuple{Int64, Int64}

    # Dynamic range filter
    dynamic_range::Float32
    prescore_dynamic_range::Float32

    function SecondPassSearchParameters(params::PioneerParameters)
        # Extract relevant parameter groups
        global_params = params.global_settings
        quant_params = params.search
        frag_params = quant_params.fragment_settings
        deconv_params = params.optimization.deconvolution

        # Determine isotope trace type based on global settings
        isotope_trace_type = if haskey(global_params.isotope_settings, :combine_traces) &&
                               global_params.isotope_settings.combine_traces
            SeperateTraces()
        else
            SeperateTraces()
        end

        # Hardcoded isotope error bounds (always (1,0))
        isotope_bounds = (UInt8(1), UInt8(0))
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

        # Hardcoded MS1 deconvolution parameters
        ms1_reg_type = L2Norm()

        # Hardcoded: ms1_scoring always true, match_between_runs always false
        ms1_scoring = true

        # Global prescore q-value threshold (default 0.05 if not in JSON)
        global_prescore_qval = if haskey(quant_params, :global_prescore_qvalue_threshold)
            Float32(quant_params.global_prescore_qvalue_threshold)
        else
            0.05f0
        end

        # Prescore fragment settings default to same as search fragment_settings
        prescore_n_frag_isotopes = Int64(frag_params.n_isotopes)
        prescore_max_frag_rank = UInt8(frag_params.max_rank)
        prescore_min_frag_count = Int64(frag_params.min_count)
        # Hardcoded prescore filter values
        prescore_min_spectral_contrast = 0.0f0
        prescore_min_log2_matched_ratio = -10.0f0
        prescore_min_topn_of_m = (0, 3)

        # Dynamic range filter (backward compatible defaults)
        dynamic_range = if haskey(frag_params, :dynamic_range)
            Float32(frag_params.dynamic_range)
        else
            Float32(1e-3)
        end
        prescore_dynamic_range = dynamic_range

        new{typeof(prec_estimation), typeof(isotope_trace_type)}(
            isotope_bounds,
            Float32(min_fraction_transmitted),
            Int64(frag_params.n_isotopes),
            UInt8(frag_params.max_rank),
            Set{Int64}([2]),
            false,  # match_between_runs hardcoded false

            Float32(deconv_params.ms2.lambda),
            reg_type,
            Int64(50),   # max_iter_newton hardcoded
            Int64(100),  # max_iter_bisection hardcoded
            Int64(deconv_params.outer_iters),
            Float32(10),  # accuracy_newton hardcoded
            Float32(10),  # accuracy_bisection hardcoded
            Float32(deconv_params.max_diff),

            Float32(0.0001),  # ms1_lambda hardcoded
            ms1_reg_type,
            Float32(1e9),     # ms1_huber_delta hardcoded

            Int64(0),   # min_y_count hardcoded
            Int64(frag_params.min_count),
            Float32(0), # min_spectral_contrast hardcoded
            Float32(-10),  # min_log2_matched_ratio hardcoded
            (0, 3),     # min_topn_of_m hardcoded
            Int64(frag_params.max_rank),

            isotope_trace_type,
            prec_estimation,

            ms1_scoring,

            global_prescore_qval,

            prescore_n_frag_isotopes,
            prescore_max_frag_rank,
            prescore_min_frag_count,
            prescore_min_spectral_contrast,
            prescore_min_log2_matched_ratio,
            prescore_min_topn_of_m,

            dynamic_range,
            prescore_dynamic_range
        )
    end
end

getIsotopeTraceType(p::SecondPassSearchParameters) = p.isotope_tracetype
