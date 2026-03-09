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

Second pass search: full-feature deconvolution on globally-filtered precursors.

Pipeline:
1. process_file!: Load filtered_fragment_matches (from FirstPassSearch), deconvolve with
   full fragment settings, compute all 29 features, select best scan per precursor using
   Phase 1 LightGBM scan preference, add multi-scan aggregates, write fold-split Arrow files
2. process_search_results!: No-op (all work done in process_file!)
3. summarize_results!: Logging only
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

            ms1_scoring,

            global_prescore_qval,

            prescore_n_frag_isotopes,
            prescore_max_frag_rank
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
Loads filtered fragment matches, deconvolves with full settings, computes all features,
selects best scan per precursor, and writes fold-split Arrow files.
"""
function process_file!(
    results::SecondPassSearchResults,
    params::P,
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::MassSpecData
) where {P<:SecondPassSearchParameters}

    if check_and_skip_failed_file(search_context, ms_file_idx, "SecondPassSearch")
        return results
    end

    try
        t_start = time()
        file_name = getParsedFileName(search_context, ms_file_idx)

        # Load filtered fragment matches (produced by FirstPassSearch::summarize_results!)
        filtered_path = getFilteredFragmentMatches(getMSData(search_context), ms_file_idx)
        if isempty(filtered_path) || !isfile(filtered_path)
            @info "  No filtered fragment matches for file $ms_file_idx ($file_name), skipping"
            return results
        end

        scan_to_prec_idx, precursors_passed = load_fragment_index_matches(
            filtered_path, length(spectra)
        )
        t_load = time()

        # Deconvolve with full Phase 2 fragment settings
        psms = perform_second_pass_search(
            spectra,
            scan_to_prec_idx,
            precursors_passed,
            search_context,
            params,
            ms_file_idx,
            MS2CHROM()  # Uses params.n_frag_isotopes and params.max_frag_rank by default
        )
        t_deconv = time()

        if nrow(psms) == 0
            @info "    No PSMs after deconvolution, skipping"
            return results
        end

        # Compute all 29 features on ALL PSMs
        prepare_psm_features!(psms, params, search_context, ms_file_idx, spectra)
        t_features = time()

        if nrow(psms) == 0
            @info "    No PSMs after feature computation, skipping"
            return results
        end

        # Retrieve Phase 1 best scans from FirstPassSearch
        first_pass_results = get_results(search_context, FirstPassSearch)
        prec_best_scan = first_pass_results.prec_best_scan

        # Build file-specific scan lookup: precursor_idx → scan_idx from Phase 1
        file_best_scans = Dictionary{UInt32, UInt32}()
        for (pid, file_scans) in pairs(prec_best_scan)
            if haskey(file_scans, ms_file_idx)
                insert!(file_best_scans, pid, file_scans[ms_file_idx])
            end
        end

        # Select best scan per precursor using Phase 1 LightGBM scan preference
        best_psms = get_best_psm_per_precursor_by_prescore_scan(psms, file_best_scans)
        t_select = time()

        n_matched = count(row -> haskey(file_best_scans, row.precursor_idx) &&
            row.scan_idx == file_best_scans[row.precursor_idx], eachrow(best_psms))
        @info "    $(nrow(best_psms)) precursors ($(n_matched) matched Phase 1 scan, $(nrow(best_psms) - n_matched) fell back to weight)"

        # Add multi-scan aggregate features for ScoringSearch
        add_multi_scan_aggregates!(best_psms, psms)

        # Add ScoringSearch-required columns
        initialize_prob_group_features!(best_psms, params.match_between_runs)

        # Write fold-split Arrow files for ScoringSearch
        base_dir = joinpath(getDataOutDir(search_context), "temp_data", "second_pass_psms")
        setSecondPassPsms!(getMSData(search_context), ms_file_idx,
                          joinpath(base_dir, file_name))

        for fold in UInt8[0, 1]
            fold_path = joinpath(base_dir, "$(file_name)_fold$(fold).arrow")
            fold_mask = best_psms.cv_fold .== fold
            if any(fold_mask)
                writeArrow(fold_path, best_psms[fold_mask, :])
            elseif isfile(fold_path)
                rm(fold_path)
            end
        end
        t_write = time()

        r = s -> round(s, digits=2)
        @info "    SecondPassSearch: $(nrow(best_psms)) precursors written to fold Arrow files\n" *
              "      load=$(r(t_load - t_start))s, deconv=$(r(t_deconv - t_load))s, features=$(r(t_features - t_deconv))s, " *
              "select=$(r(t_select - t_features))s, write=$(r(t_write - t_select))s, total=$(r(t_write - t_start))s"

    catch e
        file_name = try
            getMassSpecData(search_context).file_id_to_name[ms_file_idx]
        catch
            "file_$ms_file_idx"
        end
        @user_warn "SecondPassSearch failed for file $ms_file_idx ($file_name): $(sprint(showerror, e))\n$(sprint(Base.show_backtrace, catch_backtrace()))"
        markFileFailed!(search_context, ms_file_idx, "SecondPassSearch failed: $(typeof(e))")
        setFailedIndicator!(getMSData(search_context), ms_file_idx, true)
    end

    return results
end

"""
Create empty PSM results for a failed file in SecondPassSearch.
"""
function createFallbackResults!(results::SecondPassSearchResults, ms_file_idx::Int64)
    empty_psms = DataFrame(
        ms_file_idx = UInt32[],
        scan_idx = UInt32[],
        precursor_idx = UInt32[],
        rt = Float32[],
        q_value = Float32[],
        score = Float32[],
        prob = Float32[]
    )
    empty_ms1_psms = DataFrame(
        ms_file_idx = UInt32[],
        scan_idx = UInt32[],
        precursor_idx = UInt32[],
        rt = Float32[],
        q_value = Float32[],
        score = Float32[],
        prob = Float32[]
    )
    results.psms[] = empty_psms
    results.ms1_psms[] = empty_ms1_psms
end

"""
No additional per-file processing needed (all work done in process_file!).
"""
function process_search_results!(
    results::SecondPassSearchResults,
    params::P,
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::MassSpecData
) where {P<:SecondPassSearchParameters}
    return nothing
end

"""
Reset results containers.
"""
function reset_results!(results::SecondPassSearchResults)
    results.psms[] = DataFrame()
end

"""
Summarize results — logging only, all per-file work done in process_file!.
"""
function summarize_results!(
    results::SecondPassSearchResults,
    params::P,
    search_context::SearchContext
) where {P<:SecondPassSearchParameters}

    @info "=== SecondPassSearch complete ==="
    @info "Phase 2 settings: n_frag_isotopes=$(params.n_frag_isotopes), max_frag_rank=$(params.max_frag_rank)"

    # Count total precursors written
    ms_data = getMSData(search_context)
    n_files = length(ms_data)
    n_processed = 0
    n_psms_total = 0

    for ms_file_idx in 1:n_files
        getFailedIndicator(ms_data, ms_file_idx) && continue
        base_path = getSecondPassPsms(ms_data, ms_file_idx)
        isempty(base_path) && continue

        for fold in UInt8[0, 1]
            fold_path = "$(base_path)_fold$(fold).arrow"
            isfile(fold_path) || continue
            tbl = Arrow.Table(fold_path)
            n_psms_total += length(tbl[:precursor_idx])
        end
        n_processed += 1
    end

    @info "  $n_processed files processed, $n_psms_total total precursors in fold files"

    return nothing
end
