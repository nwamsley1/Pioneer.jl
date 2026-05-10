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
    MainSearch

Fragment index search, spectral deconvolution, and LightGBM prescore search.

Pipeline:
1. process_file!: Run fragment index search, deconvolve with prescore fragment settings
2. process_search_results!: Compute prescore features, train LightGBM, select best scan per precursor
3. summarize_results!: Global prescore aggregation, write fold-split second_pass_psms for ScoringSearch
"""
struct MainSearch <: SearchMethod end

#==========================================================
Type Definitions
==========================================================#

"""
Parameters for fragment index search phase within MainSearch.
"""
struct MainSearchFragIndexParameters <: FragmentIndexSearchParameters
    isotope_err_bounds::Tuple{UInt8, UInt8}
    min_index_search_score::UInt8
    spec_order::Set{Int64}

    function MainSearchFragIndexParameters(::PioneerParameters)
        new(
            (UInt8(1), UInt8(0)),  # isotope_err_bounds hardcoded
            DEFAULT_INDEX_SEARCH_MIN_SCORE,
            Set{Int64}([2])
        )
    end
end

"""
Results container for main search.
"""
struct MainSearchResults <: SearchResults
    psms::Base.Ref{DataFrame}
end

#==========================================================
Interface Implementation
==========================================================#

get_parameters(::MainSearch, params::Any) = MainSearchParameters(params)

function init_search_results(::MainSearch, params::P, search_context::SearchContext) where {P<:MainSearchParameters}
    main_search_psms_dir = joinpath(getDataOutDir(search_context), "temp_data", "main_search_psms")
    mkpath(main_search_psms_dir)
    second_pass_psms_dir = joinpath(getDataOutDir(search_context), "temp_data", "second_pass_psms")
    mkpath(second_pass_psms_dir)

    # One-shot sortedness check (cheaper than per-thread/per-scan asserts;
    # fails fast with an actionable error if the library is rank-sorted
    # rather than m/z-sorted within each precursor).
    let t = time()
        lookup = getFragmentLookupTable(getSpecLib(search_context))
        verify_mz_sorted(getFragments(lookup), lookup.prec_frag_ranges)
        @debug_l1 "  MainSearch: verify_mz_sorted OK " *
                   "($(round(time() - t, digits=2))s)"
    end

    return MainSearchResults(
        DataFrame()
    )
end

#==========================================================
Core Processing Methods
==========================================================#

"""
Process a single file: load fragment index matches and deconvolve with prescore settings.
"""
function process_file!(
    results::MainSearchResults,
    params::P,
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::MassSpecData
) where {P<:MainSearchParameters}

    t_file_start = time()
    file_name = getParsedFileName(search_context, ms_file_idx)

    psms = library_search(spectra, search_context, params, ms_file_idx)
    results.psms[] = psms
    t_lib_search = time() - t_file_start

    @debug_l1 "MainSearch: $file_name — $(nrow(psms)) PSMs, library_search=$(round(t_lib_search, digits=2))s"

    return results
end

"""
Per-file scoring: compute prescore features, train LightGBM, select best scan per precursor.
"""
function process_search_results!(
    results::MainSearchResults,
    params::P,
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::MassSpecData
) where {P<:MainSearchParameters}

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

    # Add per-scan competition features: weight_ratio_at_scan, weight_rank_at_scan
    add_scan_competition_features!(psms)

    # Add chromatographic neighborhood features (3-scan window for each precursor)
    add_neighborhood_features!(psms)

    # Add iRT distance from current scan to this precursor's deconv-apex scan
    add_apex_distance_feature!(psms)

    # Phase 1 MS1 features: per-PSM isotope-match against nearest MS1 scan
    add_ms1_features!(psms, spectra, search_context, ms_file_idx)

    # Pre-LGBM diagnostic dump (independent of post-LGBM dump). Captures
    # PSMs *before* LGBM training so single-class runs (target-only or
    # decoy-only) still produce usable PSM tables.
    if get(ENV, "PIONEER_DUMP_PRE_LGBM", "0") == "1"
        diag_dir = joinpath(getDataOutDir(search_context), "temp_data", "pre_lgbm_psms_diag")
        isdir(diag_dir) || mkpath(diag_dir)
        cols = [:precursor_idx, :scan_idx, :target, :weight, :gof, :fitted_manhattan_distance,
                :fitted_hellinger, :max_matched_residual, :max_unmatched_residual,
                :total_ions, :total_ions_iso, :y_count, :poisson, :err_norm,
                :missed_cleavage, :Mox, :spectrum_peak_count, :sequence_length,
                :weight_ratio_at_scan, :weight_rank_at_scan, :irt_error,
                :best_gof_3scan, :best_manhattan_3scan, :best_max_residual_3scan,
                :irt_dist_best_gof_3scan, :irt_dist_best_manhattan_3scan,
                :irt_dist_best_max_residual_3scan, :irt_dist_to_weight_apex,
                :ms1_iso_count, :ms1_m0_matched, :ms1_m0_mass_err_ppm, :ms1_log_iso_obs_pred,
                :ms1_m2_matched, :ms1_m3_matched,
                :ms1_log_iso_obs_pred_m2, :ms1_log_iso_obs_pred_m3,
                :ms1_corr_weight_m0, :ms1_corr_m0_m1,
                :ms1_apex_offset_irt, :ms1_weight_apex_to_m0_apex_irt,
                :ms1_m0_intensity, :ms1_m1_intensity,
                :irt_obs]
        keep = intersect(cols, propertynames(psms))
        Arrow.write(joinpath(diag_dir, "$(ms_file_idx).arrow"), psms[!, keep])
    end

    # Train LightGBM on ALL PSMs, select best scan per precursor
    n_total_psms = nrow(psms)
    # Diagnostic: dump pre-reduction PSMs (one row per (precursor_idx, scan_idx))
    # Diagnostic: stash file_idx + dir so train_lgbm_and_select_best can dump pre-reduction
    if get(ENV, "PIONEER_DUMP_ALL_PSMS", "0") == "1"
        Pioneer.DIAG_DUMP_FILE_IDX[] = ms_file_idx
        Pioneer.DIAG_DUMP_DIR[]      = joinpath(getDataOutDir(search_context), "temp_data", "all_psms_diag")
    else
        Pioneer.DIAG_DUMP_FILE_IDX[] = 0
    end
    best_psms, scores, lgbm_timings = train_lgbm_and_select_best(psms)
    best_psms[!, :lgbm_prob] = scores
    t_lgbm = time()

    # RT recalibration: refit iRT spline from high-confidence PSMs
    recalibrate_rt!(search_context, ms_file_idx, best_psms, best_psms[!, :lgbm_prob])

    # Update irt_obs/irt_error with post-recalibration model so downstream
    # steps (ScoringSearch features, RT index, chromatogram extraction) are consistent
    new_rt_model = getRtIrtModel(search_context, ms_file_idx)
    best_psms[!, :irt_obs] .= new_rt_model.(best_psms[!, :rt])
    best_psms[!, :irt_error] .= abs.(best_psms[!, :irt_obs] .- best_psms[!, :irt_pred])
    t_recal = time()

    # Compute isotopes_captured and filter by quad transmission
    get_isotopes_captured!(
        best_psms,
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

    # Drop vector columns that can't be serialized to Arrow
    dropVectorColumns!(best_psms)

    # Write per-fold main_search_psms (used by both aggregate_prescore_globally! and summarize_results!)
    precursors = getPrecursors(getSpecLib(search_context))
    main_search_psms_dir = joinpath(getDataOutDir(search_context), "temp_data", "main_search_psms")
    best_cv_fold = UInt8[getCvFold(precursors, pid) for pid in best_psms.precursor_idx]
    for fold in UInt8[0, 1]
        fold_mask = best_cv_fold .== fold
        any(fold_mask) && writeArrow(joinpath(main_search_psms_dir, "$(file_name)_fold$(fold).arrow"), best_psms[fold_mask, :])
    end
    t_write = time()

    # Timing summary
    t_total = t_write - t_start
    r = s -> round(s, digits=2)
    @debug_l1 "MainSearch scoring: $file_name — $(n_total_psms) PSMs → $(nrow(best_psms)) precursors " *
               "(feat=$(r(t_features-t_start))s, lgbm=$(r(lgbm_timings.train_cv))s, recal=$(r(t_recal-t_lgbm))s, " *
               "phase2=$(r(t_phase2-t_recal))s, write=$(r(t_write-t_phase2))s, total=$(r(t_total))s)"

    return nothing
end

"""
Reset results containers.
"""
function reset_results!(results::MainSearchResults)
    results.psms[] = DataFrame()
end

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
    results::MainSearchResults,
    params::P,
    search_context::SearchContext
) where {P<:MainSearchParameters}

    r(t) = round(t; digits=2)
    t_total_start = time()

    ms_data = getMSData(search_context)
    n_files = length(ms_data)
    main_search_psms_dir = joinpath(getDataOutDir(search_context), "temp_data", "main_search_psms")
    precursors = getPrecursors(getSpecLib(search_context))
    lib_irt = getIrt(precursors)

    # Step 1: Per-fold global prescore aggregation → passing precursor sets + RT-binned tolerance
    t1_start = time()
    fold0_result = aggregate_prescore_globally!(search_context; fold_suffix="_fold0")
    fold1_result = aggregate_prescore_globally!(search_context; fold_suffix="_fold1")
    passing_fold0 = fold0_result.passing
    passing_fold1 = fold1_result.passing
    # Use whichever fold's RT-binned tolerance is available (they should be similar)
    rt_binned_tol = fold0_result.rt_binned_tol !== nothing ?
                    fold0_result.rt_binned_tol : fold1_result.rt_binned_tol
    passing_precs = union(passing_fold0, passing_fold1)
    # One-line summary (combined across folds; per-fold detail is at @debug_l1)
    n_targets_total = fold0_result.n_targets_pass + fold1_result.n_targets_pass
    n_decoys_total  = fold0_result.n_decoys_pass  + fold1_result.n_decoys_pass
    n_pass_total    = n_targets_total + n_decoys_total
    @user_info "  Global prescore: $n_pass_total precursors pass q≤$(fold0_result.qvalue_threshold) ($n_targets_total targets + $n_decoys_total decoys)"
    t1 = time() - t1_start

    store_results!(search_context, MainSearch, (passing_precs=passing_precs,))

    # Step 2: Load main_search_psms, filter to passing precursors, add deferred columns, write fold-split second_pass_psms
    t2_start = time()
    second_pass_dir = joinpath(getDataOutDir(search_context), "temp_data", "second_pass_psms")

    # Library lookups (shared across files)
    prec_irt = lib_irt
    prec_mz_arr = getMz(precursors)
    prec_pair_idxs = getPairIdx(precursors)
    entrap_group_ids = getEntrapmentGroupId(precursors)

    n_processed_files = 0
    n_total_precs = 0
    n_kept_precs = 0

    for ms_file_idx in 1:n_files
        file_name = getParsedFileName(ms_data, ms_file_idx)

        base_path = joinpath(second_pass_dir, file_name)

        file_has_data = false
        n_before_file = 0
        n_after_file = 0

        for fold in UInt8[0, 1]
            passing = fold == 0 ? passing_fold0 : passing_fold1
            psm_path = joinpath(main_search_psms_dir, "$(file_name)_fold$(fold).arrow")
            fold_out_path = "$(base_path)_fold$(fold).arrow"

            if !isfile(psm_path)
                isfile(fold_out_path) && rm(fold_out_path)
                continue
            end

            # Load this fold's main search PSMs
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
            initialize_prob_group_features!(tbl)
            dropVectorColumns!(tbl)

            writeArrow(fold_out_path, tbl)
            file_has_data = true
        end

        if file_has_data
            # Only register a second-pass path for files that actually
            # produced fold output. Empty files leave the registered path
            # untouched (typically empty string), so downstream code that
            # iterates `getSecondPassPsms` paths skips them naturally
            # without needing a separate failed-files set.
            setSecondPassPsms!(ms_data, ms_file_idx, base_path)
            n_processed_files += 1
            pct = round(100.0 * n_after_file / max(1, n_before_file), digits=1)
            @debug "  $file_name: $n_after_file / $n_before_file precursors kept ($pct%)"
        end
    end
    t2 = time() - t2_start

    overall_pct = round(100.0 * n_kept_precs / max(1, n_total_precs), digits=1)
    @debug_l1 "MainSearch passing PSMs: $n_kept_precs / $n_total_precs precursors ($overall_pct%) across $n_processed_files files"

    # Step 3: Compute chromatographic tolerance from fold files
    t3_start = time()
    if n_processed_files > 0
        if rt_binned_tol !== nothing
            compute_rt_binned_tolerance!(search_context, rt_binned_tol, ms_data, n_files)
        else
            @debug_l1 "No RT-binned tolerance available — chromatogram integration will use per-file iRT tolerance from recalibrate_rt!"
        end
    else
        @warn "No files processed in MainSearch — skipping chromatographic tolerance computation"
    end
    t3 = time() - t3_start

    t_total = time() - t_total_start
    @debug_l1 "MainSearch summarize: aggregation=$(r(t1))s, fold_write=$(r(t2))s, chrom_tol=$(r(t3))s, total=$(r(t_total))s"

    return nothing
end
