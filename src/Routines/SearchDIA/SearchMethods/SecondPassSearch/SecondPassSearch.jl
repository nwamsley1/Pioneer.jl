# Types (SecondPassSearch, SecondPassSearchResults, SecondPassSearchParameters)
# are defined in types.jl, loaded before utils.jl via importScripts.jl

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
