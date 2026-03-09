# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

"""
    SecondPassSearch

Second pass search method using optimized parameters from initial searches.

This search:
1. Uses optimized parameters from first pass search
2. Performs PSM identification with previously calculated Huber delta
3. Tracks retention time windows for efficient searching
4. Records chromatogram information for later analysis

# Example Implementation
```julia
# Define search parameters
params = Dict(
    :second_pass_params => Dict(
        "min_y_count" => 1,
        "min_frag_count" => 3,
        "min_spectral_contrast" => 0.1,
        "min_log2_matched_ratio" => -3.0,
        "min_topn_of_m" => (3, 5),
        "max_best_rank" => 3,
        "n_frag_isotopes" => 2,
        "max_frag_rank" => 10
    )
)

# Execute search
results = execute_search(SecondPassSearch(), search_context, params)
```
"""
struct SecondPassSearch <: SearchMethod end

#==========================================================
Type Definitions
==========================================================#

"""
Results container for second pass search.
"""
struct SecondPassSearchResults <: SearchResults
    psms::Base.Ref{DataFrame}          # PSMs for each file
    ms1_psms::Base.Ref{DataFrame}
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

        new{typeof(prec_estimation), typeof(isotope_trace_type)}(
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

            ms1_scoring
        )
    end
end

getIsotopeTraceType(p::SecondPassSearchParameters) = p.isotope_tracetype
#==========================================================
Interface Implementation
==========================================================#

get_parameters(::SecondPassSearch, params::Any) = SecondPassSearchParameters(params)

function init_search_results(::P, search_context::SearchContext) where {P<:SecondPassSearchParameters}
    second_pass_psms = joinpath(getDataOutDir(search_context), "temp_data", "second_pass_psms")
    !isdir(second_pass_psms) && mkdir(second_pass_psms)
    return SecondPassSearchResults(
        DataFrame(),
        DataFrame()
    )
end


#==========================================================
Core Processing Methods
==========================================================#

"""
Process a single file for second pass search.
"""
function process_file!(
    results::SecondPassSearchResults,
    params::P, 
    search_context::SearchContext,    
    ms_file_idx::Int64,
    spectra::MassSpecData
) where {P<:SecondPassSearchParameters}

    # Check if file should be skipped due to previous failure
    if check_and_skip_failed_file(search_context, ms_file_idx, "SecondPassSearch")
        return results  # Return early with unchanged results
    end

    try
        # Load fragment index mapping instead of RT index
        frag_match_path = getFragmentIndexMatches(getMSData(search_context), ms_file_idx)
        scan_to_prec_idx, precursors_passed = load_fragment_index_matches(
            frag_match_path, length(spectra)
        )

        # Perform second pass search using fragment index matches
        psms = perform_second_pass_search(
            spectra,
            scan_to_prec_idx,
            precursors_passed,
            search_context,
            params,
            ms_file_idx,
            MS2CHROM()
        )

        results.psms[] = psms

    catch e
        # Handle failures gracefully using helper function (logs full stacktrace)
        handle_search_error!(search_context, ms_file_idx, "SecondPassSearch", e, createFallbackResults!, results)
    end

    return results
end

"""
Create empty PSM results for a failed file in SecondPassSearch.
"""
function createFallbackResults!(results::SecondPassSearchResults, ms_file_idx::Int64)
    # Create empty PSM DataFrame with proper schema for regular PSMs
    empty_psms = DataFrame(
        ms_file_idx = UInt32[],
        scan_idx = UInt32[], 
        precursor_idx = UInt32[],
        rt = Float32[],
        q_value = Float32[],
        score = Float32[], 
        prob = Float32[]
    )
    
    # Create empty MS1 PSM DataFrame with proper schema
    empty_ms1_psms = DataFrame(
        ms_file_idx = UInt32[],
        scan_idx = UInt32[], 
        precursor_idx = UInt32[],
        rt = Float32[],
        q_value = Float32[],
        score = Float32[], 
        prob = Float32[]
    )
    
    # Set empty results (don't append since this file failed)
    results.psms[] = empty_psms
    results.ms1_psms[] = empty_ms1_psms
end


function process_search_results!(
    results::SecondPassSearchResults,
    params::P,
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::MassSpecData
) where {P<:SecondPassSearchParameters}

    # Check if file should be skipped due to previous failure
    if check_and_skip_failed_file(search_context, ms_file_idx, "SecondPassSearch results processing")
        return nothing
    end

    try
        psms = results.psms[]
        file_name = getParsedFileName(search_context, ms_file_idx)

        # Phase 1: Compute 29 features on ALL PSMs (all scans, not best-per-precursor)
        prepare_psm_features!(psms, params, search_context, ms_file_idx, spectra)

        if nrow(psms) == 0
            @debug_l2 "No PSMs for file $ms_file_idx after feature computation"
            setSecondPassPsms!(getMSData(search_context), ms_file_idx, "")
            return nothing
        end

        # Train LightGBM on ALL PSMs, select best scan per precursor
        best_psms, scores, q_values = train_lgbm_and_select_best(psms)

        # DIA-NN diagnostics at various q-value thresholds
        diann_file, diann_global = load_diann_reference(file_name)
        if !isempty(diann_file) || !isempty(diann_global)
            best_targets = best_psms[!, :target]
            for q_thresh in [0.01, 0.05, 0.10]
                passing = Set{UInt32}(best_psms[i, :precursor_idx]
                    for i in 1:nrow(best_psms) if (q_values[i] <= q_thresh) && best_targets[i])
                n_file = isempty(diann_file) ? 0 : length(intersect(passing, diann_file))
                n_global = isempty(diann_global) ? 0 : length(intersect(passing, diann_global))
                @info "    Prescore q≤$q_thresh: $(length(passing)) targets, DIA-NN file=$n_file global=$n_global"
            end
        end

        # Write prescore table (scores only, NO fold Arrow files yet)
        prescore_dir = joinpath(getDataOutDir(search_context), "temp_data", "prescore_scores")
        mkpath(prescore_dir)
        score_df = DataFrame(
            precursor_idx = best_psms[!, :precursor_idx],
            lgbm_prob = scores,
            target = best_psms[!, :target]
        )
        writeArrow(joinpath(prescore_dir, "$(file_name).arrow"), score_df)

        # Register base path for summarize_results! (fold files written in Phase 2)
        base_dir = joinpath(getDataOutDir(search_context), "temp_data", "second_pass_psms")
        setSecondPassPsms!(getMSData(search_context), ms_file_idx,
                          joinpath(base_dir, file_name))

    catch e
        file_name = try
            getMassSpecData(search_context).file_id_to_name[ms_file_idx]
        catch
            "file_$ms_file_idx"
        end
        reason = "SecondPassSearch failed: $(typeof(e))"
        markFileFailed!(search_context, ms_file_idx, reason)
        setFailedIndicator!(getMSData(search_context), ms_file_idx, true)
        @user_warn "Second pass search failed for MS data file: $file_name. Error: $(sprint(showerror, e))\n$(sprint(Base.show_backtrace, catch_backtrace()))"
        setSecondPassPsms!(getMSData(search_context), ms_file_idx, "")
    end

    return nothing
end


"""
Reset results containers.
"""
function reset_results!(results::SecondPassSearchResults)
    results.psms[] = DataFrame()
end

function summarize_results!(
    results::SecondPassSearchResults,
    params::P,
    search_context::SearchContext
) where {P<:SecondPassSearchParameters}

    @info "=== SecondPassSearch: Phase 2 — Global aggregation + final search ==="

    # Step 1: Global prescore aggregation → passing precursor set
    passing_precs = aggregate_prescore_globally!(search_context)

    # Step 2: Re-search each file with only passing precursors
    msdr = getMassSpecData(search_context)
    ms_data = getMSData(search_context)
    n_files = length(ms_data)
    n_rerun = 0
    n_psms_total = 0

    @info "Re-deconvolving $(n_files) files with $(length(passing_precs)) globally-passing precursors"

    for ms_file_idx in 1:n_files
        getFailedIndicator(ms_data, ms_file_idx) && continue
        base_path = getSecondPassPsms(ms_data, ms_file_idx)
        isempty(base_path) && continue

        file_name = getParsedFileName(search_context, ms_file_idx)
        @info "  Phase 2 file $ms_file_idx ($file_name)..."

        try
            spectra = getMSData(msdr, ms_file_idx)

            # Re-deconvolve with only globally-passing precursors (reduced competition)
            psms = rerun_search_with_precursor_filter(
                spectra, search_context, params, ms_file_idx, passing_precs
            )

            if nrow(psms) == 0
                @info "    No PSMs after global re-search, skipping"
                continue
            end

            # Compute 29 features on ALL PSMs
            prepare_psm_features!(psms, params, search_context, ms_file_idx, spectra)

            if nrow(psms) == 0
                @info "    No PSMs after feature computation, skipping"
                continue
            end

            # Train LightGBM → best scan per precursor
            best_psms, _, _ = train_lgbm_and_select_best(psms)

            # Add ScoringSearch-required columns
            initialize_prob_group_features!(best_psms, params.match_between_runs)

            # Write fold-split Arrow files for ScoringSearch
            base_dir = joinpath(getDataOutDir(search_context), "temp_data", "second_pass_psms")
            for fold in UInt8[0, 1]
                fold_path = joinpath(base_dir, "$(file_name)_fold$(fold).arrow")
                fold_mask = best_psms.cv_fold .== fold
                if any(fold_mask)
                    writeArrow(fold_path, best_psms[fold_mask, :])
                elseif isfile(fold_path)
                    rm(fold_path)
                end
            end

            n_rerun += 1
            n_psms_total += nrow(best_psms)
            @info "    $(nrow(best_psms)) precursors written to fold Arrow files"

        catch e
            @user_warn "Phase 2 failed for file $ms_file_idx ($file_name): $(sprint(showerror, e))"
        end
    end

    @info "=== Phase 2 complete: $n_rerun files, $n_psms_total total precursors ==="

    return nothing
end
