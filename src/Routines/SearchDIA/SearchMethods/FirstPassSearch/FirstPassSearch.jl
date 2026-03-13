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
    FirstPassSearch

First pass deconvolution and LightGBM prescore search.

Pipeline:
1. process_file!: Load fragment_index_matches, deconvolve with prescore fragment settings
2. process_search_results!: Compute prescore features, train LightGBM, select best scan per precursor
3. summarize_results!: Global prescore aggregation, filter fragment_index_matches to passing precursors,
   write filtered_fragment_matches for SecondPassSearch
"""
struct FirstPassSearch <: SearchMethod end

#==========================================================
Type Definitions
==========================================================#

"""
Results container for first pass prescore search.
"""
struct FirstPassSearchResults <: SearchResults
    psms::Base.Ref{DataFrame}
    file_fwhms::Dict{Int, @NamedTuple{median_fwhm::Float32, mad_fwhm::Float32}}
end

#==========================================================
Interface Implementation
==========================================================#

get_parameters(::FirstPassSearch, params::Any) = SecondPassSearchParameters(params)

function init_search_results(::FirstPassSearch, ::P, search_context::SearchContext) where {P<:SecondPassSearchParameters}
    prescore_dir = joinpath(getDataOutDir(search_context), "temp_data", "prescore_scores")
    mkpath(prescore_dir)
    filtered_dir = joinpath(getDataOutDir(search_context), "temp_data", "filtered_fragment_matches")
    mkpath(filtered_dir)
    first_pass_psms_dir = joinpath(getDataOutDir(search_context), "temp_data", "first_pass_psms")
    mkpath(first_pass_psms_dir)
    second_pass_psms_dir = joinpath(getDataOutDir(search_context), "temp_data", "second_pass_psms")
    mkpath(second_pass_psms_dir)
    return FirstPassSearchResults(
        DataFrame(),
        Dict{Int, @NamedTuple{median_fwhm::Float32, mad_fwhm::Float32}}()
    )
end

#==========================================================
Core Processing Methods
==========================================================#

"""
Process a single file: load fragment index matches and deconvolve with prescore settings.
"""
function process_file!(
    results::FirstPassSearchResults,
    params::P,
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::MassSpecData
) where {P<:SecondPassSearchParameters}

    if check_and_skip_failed_file(search_context, ms_file_idx, "FirstPassSearch")
        return results
    end

    try
        t_start = time()

        # Load fragment index mapping from FragmentIndexSearch
        frag_match_path = getFragmentIndexMatches(getMSData(search_context), ms_file_idx)
        scan_to_prec_idx, precursors_passed = load_fragment_index_matches(
            frag_match_path, length(spectra)
        )
        n_input_pairs = sum(ismissing(r) ? 0 : length(r) for r in scan_to_prec_idx)
        file_name = getParsedFileName(search_context, ms_file_idx)
        @info "FirstPassSearch input: $file_name — $(n_input_pairs) (scan, precursor) pairs" *
              " | isotope_err_bounds=$(params.isotope_err_bounds), n_frag_isotopes=$(params.prescore_n_frag_isotopes), max_frag_rank=$(params.prescore_max_frag_rank)"
        t_load = time()

        # Deconvolve with Phase 1 prescore fragment settings
        search_result = perform_second_pass_search(
            spectra,
            scan_to_prec_idx,
            precursors_passed,
            search_context,
            params,
            ms_file_idx,
            MS2CHROM();
            n_frag_isotopes = params.prescore_n_frag_isotopes,
            max_frag_rank = params.prescore_max_frag_rank,
            min_frag_count = params.prescore_min_frag_count,
            min_spectral_contrast = params.prescore_min_spectral_contrast,
            min_log2_matched_ratio = params.prescore_min_log2_matched_ratio,
            min_topn_of_m = params.prescore_min_topn_of_m,
            dynamic_range = params.prescore_dynamic_range,
            first_pass = true
        )

        psms = search_result.psms
        t_deconv = time()

        results.psms[] = psms

        @info "FirstPassSearch deconvolution: $(nrow(psms)) PSMs, load=$(round(t_load - t_start, digits=2))s, deconv=$(round(t_deconv - t_load, digits=2))s"

    catch e
        handle_search_error!(search_context, ms_file_idx, "FirstPassSearch", e, createFirstPassFallbackResults!, results)
    end

    return results
end

"""
Create empty PSM results for a failed file.
"""
function createFirstPassFallbackResults!(results::FirstPassSearchResults, ms_file_idx::Int64)
    results.psms[] = DataFrame(
        ms_file_idx = UInt32[],
        scan_idx = UInt32[],
        precursor_idx = UInt32[],
        rt = Float32[],
        q_value = Float32[],
        score = Float32[],
        prob = Float32[]
    )
end

"""
Per-file scoring: compute prescore features, train LightGBM, select best scan per precursor.
"""
function process_search_results!(
    results::FirstPassSearchResults,
    params::P,
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::MassSpecData
) where {P<:SecondPassSearchParameters}

    if check_and_skip_failed_file(search_context, ms_file_idx, "FirstPassSearch results processing")
        return nothing
    end

    try
        t_start = time()
        psms = results.psms[]
        file_name = getParsedFileName(search_context, ms_file_idx)

        # Compute only prescore features (skip columns recomputed in Phase 2)
        prepare_psm_features!(psms, params, search_context, ms_file_idx, spectra, prescore_only=true)
        t_features = time()

        if nrow(psms) == 0
            @debug_l2 "No PSMs for file $ms_file_idx after feature computation"
            return nothing
        end

        # Train LightGBM on ALL PSMs, select best scan per precursor
        n_total_psms = nrow(psms)
        best_psms, scores, q_values, lgbm_timings = train_lgbm_and_select_best(psms)
        best_psms[!, :lgbm_prob] = scores
        best_psms[!, :prescore_q_value] = q_values
        t_lgbm = time()

        # RT recalibration: refit iRT spline from high-confidence PSMs
        recalibrate_rt!(search_context, ms_file_idx, best_psms, best_psms[!, :lgbm_prob])
        t_recal = time()

        # Write per-fold prescore tables (scores only — needed by aggregate_prescore_globally!)
        prescore_dir = joinpath(getDataOutDir(search_context), "temp_data", "prescore_scores")
        mkpath(prescore_dir)
        precursors = getPrecursors(getSpecLib(search_context))
        score_df = DataFrame(
            precursor_idx = best_psms[!, :precursor_idx],
            lgbm_prob = best_psms[!, :lgbm_prob],
            target = best_psms[!, :target],
            irt_obs = best_psms[!, :irt_obs],
            irt_pred = best_psms[!, :irt_pred],
            rt = best_psms[!, :rt],
            scan_idx = best_psms[!, :scan_idx]
        )
        score_cv_fold = UInt8[getCvFold(precursors, pid) for pid in score_df.precursor_idx]
        for fold in UInt8[0, 1]
            fold_mask = score_cv_fold .== fold
            any(fold_mask) && writeArrow(joinpath(prescore_dir, "$(file_name)_fold$(fold).arrow"), score_df[fold_mask, :])
        end
        t_prescore_write = time()

        # Phase 2: compute isotopes_captured (needs spectra data)
        get_isotopes_captured!(
            best_psms,
            getIsotopeTraceType(params),
            getQuadTransmissionModel(search_context, ms_file_idx),
            getSearchData(search_context),
            best_psms[!, :scan_idx],
            getCharge(getPrecursors(getSpecLib(search_context))),
            getMz(getPrecursors(getSpecLib(search_context))),
            getSulfurCount(getPrecursors(getSpecLib(search_context))),
            getCenterMzs(spectra),
            getIsolationWidthMzs(spectra)
        )
        # Filter by precursor_fraction_transmitted
        to_remove = findall(best_psms[!, :precursor_fraction_transmitted] .< params.min_fraction_transmitted)
        deleteat!(best_psms, to_remove)
        best_psms[!, :ms_file_idx] .= UInt32(ms_file_idx)
        t_phase2 = time()

        # Hardcoded FWHM — replace with real estimation when bypassing SecondPass
        results.file_fwhms[ms_file_idx] = (median_fwhm=0.2f0, mad_fwhm=0.2f0)

        # Drop vector columns that can't be serialized to Arrow
        dropVectorColumns!(best_psms)

        # Write per-fold first_pass_psms for summarize_results! to filter
        first_pass_psms_dir = joinpath(getDataOutDir(search_context), "temp_data", "first_pass_psms")
        best_cv_fold = UInt8[getCvFold(precursors, pid) for pid in best_psms.precursor_idx]
        for fold in UInt8[0, 1]
            fold_mask = best_cv_fold .== fold
            any(fold_mask) && writeArrow(joinpath(first_pass_psms_dir, "$(file_name)_fold$(fold).arrow"), best_psms[fold_mask, :])
        end
        t_write = time()

        # Timing summary
        n_pass_1 = count((best_psms[!, :prescore_q_value] .<= 0.01) .& best_psms[!, :target])
        n_pass_5 = count((best_psms[!, :prescore_q_value] .<= 0.05) .& best_psms[!, :target])
        t_total = t_write - t_start
        r = s -> round(s, digits=2)
        println()
        @info "FirstPassSearch scoring: $(n_total_psms) PSMs → $(nrow(best_psms)) precursors ($n_pass_1 @ 1%, $n_pass_5 @ 5% FDR)\n" *
              "  features=$(r(t_features - t_start))s, recal=$(r(t_recal - t_lgbm))s, phase2=$(r(t_phase2 - t_prescore_write))s, write=$(r(t_write - t_phase2))s\n" *
              "  lgbm: matrix=$(r(lgbm_timings.matrix))s, train_cv=$(r(lgbm_timings.train_cv))s, best=$(r(lgbm_timings.best))s, qval=$(r(lgbm_timings.qval))s\n" *
              "  total=$(r(t_total))s"
        println()

    catch e
        file_name = try
            getMassSpecData(search_context).file_id_to_name[ms_file_idx]
        catch
            "file_$ms_file_idx"
        end
        reason = "FirstPassSearch scoring failed: $(typeof(e))"
        markFileFailed!(search_context, ms_file_idx, reason)
        setFailedIndicator!(getMSData(search_context), ms_file_idx, true)
        @user_warn "First pass scoring failed for MS data file: $file_name. Error: $(sprint(showerror, e))\n$(sprint(Base.show_backtrace, catch_backtrace()))"
    end

    return nothing
end

"""
Reset results containers.
"""
function reset_results!(results::FirstPassSearchResults)
    results.psms[] = DataFrame()
end

"""
Global aggregation: aggregate prescore across files, filter fragment_index_matches
to passing precursors, write filtered_fragment_matches for SecondPassSearch.
"""

"""
Typed inner function for Phase 2 library-lookup columns in `summarize_results!`.
Accepts DataFrame columns (already `Vector{T}`) and library arrays as `AbstractVector`
so the compiler specializes on concrete types, eliminating dynamic dispatch.
"""
function _compute_phase2_columns!(
        prec_idx_col::AbstractVector{UInt32},
        irt_obs_col::AbstractVector{Float32},
        prec_irt, prec_mz_arr, prec_pair_idxs, entrap_group_ids,
        irt_diff_col::Vector{Float32},
        prec_mz_col::Vector{Float32},
        pair_id_col::Vector{UInt32},
        entrap_col::Vector{UInt8})
    @inbounds for i in eachindex(prec_idx_col)
        pid = prec_idx_col[i]
        irt_diff_col[i] = abs(irt_obs_col[i] - prec_irt[pid])
        prec_mz_col[i] = prec_mz_arr[pid]
        pair_id_col[i] = extract_pair_idx(prec_pair_idxs, pid)
        entrap_col[i] = entrap_group_ids[pid]
    end
    return nothing
end

function summarize_results!(
    results::FirstPassSearchResults,
    params::P,
    search_context::SearchContext
) where {P<:SecondPassSearchParameters}

    r(t) = round(t; digits=2)
    t_total_start = time()

    @info "=== FirstPassSearch: Global prescore aggregation + fragment index filtering ==="

    # Step 1: Per-fold global prescore aggregation → passing precursor sets + Phase 1 scan lookup
    t1_start = time()
    passing_fold0, best_scan_fold0 = aggregate_prescore_globally!(
        search_context, params.global_prescore_qvalue_threshold, params.prescore_aggregation;
        fold_suffix="_fold0")
    passing_fold1, best_scan_fold1 = aggregate_prescore_globally!(
        search_context, params.global_prescore_qvalue_threshold, params.prescore_aggregation;
        fold_suffix="_fold1")
    passing_precs = union(passing_fold0, passing_fold1)
    prec_best_scan = merge(best_scan_fold0, best_scan_fold1)
    t1 = time() - t1_start

    # Store results for SecondPassSearch to retrieve Phase 1 best scans
    store_results!(search_context, FirstPassSearch, (passing_precs=passing_precs, prec_best_scan=prec_best_scan))

    # Step 2: Load first_pass_psms, filter to passing precursors, add deferred columns, write fold-split second_pass_psms
    t2_start = time()
    ms_data = getMSData(search_context)
    n_files = length(ms_data)
    first_pass_psms_dir = joinpath(getDataOutDir(search_context), "temp_data", "first_pass_psms")
    second_pass_dir = joinpath(getDataOutDir(search_context), "temp_data", "second_pass_psms")

    # Library lookups (shared across files)
    precursors = getPrecursors(getSpecLib(search_context))
    prec_irt = getIrt(precursors)
    prec_mz_arr = getMz(precursors)
    prec_pair_idxs = getPairIdx(precursors)
    entrap_group_ids = getEntrapmentGroupId(precursors)

    n_processed_files = 0
    n_total_precs = 0
    n_kept_precs = 0

    for ms_file_idx in 1:n_files
        getFailedIndicator(ms_data, ms_file_idx) && continue

        file_name = getParsedFileName(ms_data, ms_file_idx)

        # Register base path for ScoringSearch to find fold files
        base_path = joinpath(second_pass_dir, file_name)
        setSecondPassPsms!(ms_data, ms_file_idx, base_path)

        file_has_data = false
        n_before_file = 0
        n_after_file = 0

        for fold in UInt8[0, 1]
            passing = fold == 0 ? passing_fold0 : passing_fold1
            psm_path = joinpath(first_pass_psms_dir, "$(file_name)_fold$(fold).arrow")
            fold_out_path = "$(base_path)_fold$(fold).arrow"

            if !isfile(psm_path)
                isfile(fold_out_path) && rm(fold_out_path)
                continue
            end

            # Load this fold's first pass PSMs
            tbl = DataFrame(Tables.columntable(Arrow.Table(psm_path)))
            n_before = nrow(tbl)
            n_before_file += n_before
            n_total_precs += n_before

            # Filter to this fold's passing precursors
            mask = in.(tbl[!, :precursor_idx], Ref(passing))
            tbl = tbl[mask, :]
            n_after = nrow(tbl)
            n_after_file += n_after
            n_kept_precs += n_after

            if n_after == 0
                isfile(fold_out_path) && rm(fold_out_path)
                continue
            end

            # Add Phase 2 library-lookup columns (deferred from process_search_results!)
            N = nrow(tbl)
            irt_diff_col = Vector{Float32}(undef, N)
            prec_mz_col = Vector{Float32}(undef, N)
            pair_id_col = Vector{UInt32}(undef, N)
            entrap_col = Vector{UInt8}(undef, N)
            _compute_phase2_columns!(
                tbl[!, :precursor_idx], tbl[!, :irt_obs],
                prec_irt, prec_mz_arr, prec_pair_idxs, entrap_group_ids,
                irt_diff_col, prec_mz_col, pair_id_col, entrap_col
            )
            tbl[!, :irt_diff] = irt_diff_col
            tbl[!, :prec_mz] = prec_mz_col
            tbl[!, :pair_id] = pair_id_col
            tbl[!, :entrapment_group_id] = entrap_col

            sort!(tbl, :rt)
            initialize_prob_group_features!(tbl, params.match_between_runs)
            dropVectorColumns!(tbl)

            writeArrow(fold_out_path, tbl)
            file_has_data = true
        end

        if file_has_data
            n_processed_files += 1
            pct = round(100.0 * n_after_file / max(1, n_before_file), digits=1)
            @debug "  $file_name: $n_after_file / $n_before_file precursors kept ($pct%)"
        end
    end
    t2 = time() - t2_start

    overall_pct = round(100.0 * n_kept_precs / max(1, n_total_precs), digits=1)
    @info "Second pass PSMs (from first pass): $n_kept_precs / $n_total_precs precursors ($overall_pct%) across $n_processed_files files"

    # Step 3: Compute chromatographic tolerance from fold files
    t3_start = time()
    if n_processed_files > 0
        compute_chromatographic_tolerance!(search_context, results.file_fwhms, ms_data, n_files)
    else
        @warn "No files processed in FirstPassSearch — skipping chromatographic tolerance computation"
    end
    t3 = time() - t3_start

    # Step 4: Filter fragment_index_matches to passing precursors (still needed by downstream)
    t4_start = time()
    filtered_dir = joinpath(getDataOutDir(search_context), "temp_data", "filtered_fragment_matches")
    mkpath(filtered_dir)

    n_filtered_files = 0
    n_total_entries = 0
    n_kept_entries = 0

    for ms_file_idx in 1:n_files
        getFailedIndicator(ms_data, ms_file_idx) && continue

        frag_match_path = getFragmentIndexMatches(ms_data, ms_file_idx)
        isempty(frag_match_path) && continue
        !isfile(frag_match_path) && continue

        file_name = getParsedFileName(ms_data, ms_file_idx)

        tbl = DataFrame(Tables.columntable(Arrow.Table(frag_match_path)))
        n_before = nrow(tbl)
        n_total_entries += n_before

        mask = in.(tbl[!, :precursor_idx], Ref(passing_precs))
        tbl = tbl[mask, :]
        n_after = nrow(tbl)
        n_kept_entries += n_after

        filtered_path = joinpath(filtered_dir, "$(file_name).arrow")
        writeArrow(filtered_path, tbl)
        setFilteredFragmentMatches!(ms_data, ms_file_idx, filtered_path)

        n_filtered_files += 1
        pct = round(100.0 * n_after / max(1, n_before), digits=1)
        @debug "  $file_name: $n_after / $n_before entries kept ($pct%)"
    end
    t4 = time() - t4_start

    overall_pct_frag = round(100.0 * n_kept_entries / max(1, n_total_entries), digits=1)
    @info "Filtered fragment matches: $n_kept_entries / $n_total_entries entries ($overall_pct_frag%) across $n_filtered_files files"

    t_total = time() - t_total_start
    @info "FirstPassSearch summarize: aggregation=$(r(t1))s, fold_write=$(r(t2))s, chrom_tol=$(r(t3))s, frag_filter=$(r(t4))s, total=$(r(t_total))s"

    return nothing
end
