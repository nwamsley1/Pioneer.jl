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

"""
    _summarize_psm_counts(best_psms, stage_label, ms_file_idx, file_name)

Compute per-file q-values and PEPs from `best_psms[!, :lgbm_prob]` + `:target`
and emit a one-line @user_info summary with target counts at q≤.001, q≤.01,
PEP≤.01, PEP≤.05. Diagnostic only — purely transparent (no mutation of
best_psms columns). Promoted from @debug_l2 during the tuning phase.
"""
function _summarize_psm_counts(best_psms::DataFrame, stage_label::AbstractString,
                                ms_file_idx::Integer, file_name::AbstractString)
    nrow(best_psms) == 0 && return
    probs = Float32.(best_psms[!, :lgbm_prob])
    is_t  = Vector{Bool}(best_psms[!, :target])
    qv = zeros(Float32, length(probs))
    get_qvalues!(probs, is_t, qv; doSort=true, fdr_scale_factor=1.0f0)
    peps = Vector{Float32}(undef, length(probs))
    get_PEP!(probs, is_t, peps; doSort=true, fdr_scale_factor=1.0f0)
    n_t_q001  = count((qv .<= 0.001f0) .& is_t)
    n_t_q01   = count((qv .<= 0.01f0)  .& is_t)
    n_t_pep01 = count((peps .<= 0.01f0) .& is_t)
    n_t_pep05 = count((peps .<= 0.05f0) .& is_t)
    # Trailing "\n" pushes the next ProgressBar tick to a fresh line so this
    # message doesn't get overwritten when the per-file iteration advances.
    @user_info "  [$stage_label] (file_idx=$ms_file_idx, $file_name): " *
               "$(nrow(best_psms)) best-per-precursor PSMs; " *
               "targets q≤.001=$n_t_q001  q≤.01=$n_t_q01  PEP≤.01=$n_t_pep01  PEP≤.05=$n_t_pep05\n"
end

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

    @user_info "  MainSearch process_file! (file_idx=$ms_file_idx, $file_name): " *
               "$(nrow(psms)) PSMs from deconv; library_search elapsed: $(round(t_lib_search, digits=2))s"

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
    t_prepare = @elapsed prepare_psm_features!(psms, params, search_context, ms_file_idx, spectra, prescore_only=true)
    t_features = time()

    if nrow(psms) == 0
        @debug_l2 "No PSMs for file $ms_file_idx after feature computation"
        return nothing
    end

    t_competition = @elapsed add_scan_competition_features!(psms)
    t_neighborhood = 0.0  # neighborhood feature pass removed 2026-05-18 (all 9
                          # consumers — best_*_5scan, irt_dist_best_*_5scan,
                          # worst_*_11scan, best_max_residual_3scan — were
                          # net-harmful on Astral/Exploris/MTAC per overnight
                          # sweep; dropped from PRESCORE_FEATURES & ADVANCED).
    # add_apex_distance_feature! was deleted 2026-05-18 — all 19 columns it
    # produced (irt_dist_to_weight_apex, relative_position, dist_to_relative_center,
    # shape_m2/m1/0/p1/p2, weight_apex_relative_pos, weight_chrom_skewness/kurtosis,
    # apex_to_edge_weight_ratio, n_above_hm_{left,right}_of_apex, hm_asymmetry,
    # weight_chrom_gaussian_r2/sigma/apex_irt, gof_minus_max_gof_precursor,
    # gof_fraction_of_total_precursor) were dead by grep across src/. Pass-1 was
    # already skipping the call; Pass-2 silently dropped the columns via the
    # hasproperty filter in the LGBM trainer. See 1afebf9d for historical
    # individual-feature LGBM gains if a future campaign wants to resurrect.
    t_apex = 0.0
    t_ms1 = @elapsed add_ms1_features!(psms, spectra, search_context, ms_file_idx)

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
                :irt_dist_to_weight_apex,
                :ms1_m0_mass_err_ppm,
                :ms1_corr_weight_m0, :ms1_corr_m0_m1,
                :ms1_apex_offset_irt, :ms1_weight_apex_to_m0_apex_irt,
                :ms1_m0_intensity, :ms1_m1_intensity,
                :frag1_int, :frag2_int, :frag3_int, :frag4_int, :frag5_int, :frag6_int,
                :frag_corr_top1_top2, :frag_corr_top1_top3, :frag_corr_top1_weight,
                :frag_corr_mean_pairwise, :frag_corr_min_pairwise, :frag_corr_top3_weight,
                :frag_apex_dispersion_irt, :n_correlated_fragments,
                :irt_obs, :rt]
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
    _summarize_psm_counts(best_psms, "after best-per-precursor", ms_file_idx, file_name)
    t_lgbm = time()

    # ============================================================
    # OPTIONAL 2-PASS REFINEMENT (gated by PIONEER_PEP_FILTER_REDECONV)
    #
    # Compute PEP per precursor from pass-1 LGBM scores. Remove precursors
    # with PEP > threshold (default 0.9). Re-run library_search restricted
    # to surviving precursors → cleaner deconv (less peak competition) →
    # re-train LGBM on pass-2 PSMs. Replaces best_psms with pass-2 output.
    # ============================================================
    if get(ENV, "PIONEER_PEP_FILTER_REDECONV", "0") == "1"
        pass1_psms_n = nrow(psms)
        pass1_precs_n = nrow(best_psms)

        # Filter selection: qval if PIONEER_QVAL_FILTER_THR set, else PEP
        use_qval = haskey(ENV, "PIONEER_QVAL_FILTER_THR")
        thr = use_qval ? Float32(parse(Float64, ENV["PIONEER_QVAL_FILTER_THR"])) :
                          Float32(parse(Float64, get(ENV, "PIONEER_PEP_FILTER_THR", "0.9")))

        # Compute the discriminant
        probs = Float32.(best_psms[!, :lgbm_prob])
        targets_bool = Vector{Bool}(best_psms[!, :target])

        if use_qval
            # Compute per-file q-values on best-per-precursor pass-1 scores
            sort_idx = sortperm(probs; rev=true)
            qv = zeros(Float64, length(probs))
            let nT = 0, nD = 0
                for k in eachindex(sort_idx)
                    i = sort_idx[k]
                    if targets_bool[i]; nT += 1; else; nD += 1; end
                    qv[i] = nD / max(nT, 1)
                end
                # Monotonize from end → start (revisit in score-descending order)
                fdr = Inf
                for k in length(sort_idx):-1:1
                    i = sort_idx[k]
                    if qv[i] > fdr; qv[i] = fdr; else; fdr = qv[i]; end
                end
            end
            allowed_pids = Set{UInt32}()
            for i in eachindex(qv)
                if qv[i] <= thr
                    push!(allowed_pids, UInt32(best_psms.precursor_idx[i]))
                end
            end
            filter_desc = "qval≤$(thr)"
        else
            peps = Vector{Float32}(undef, length(probs))
            get_PEP!(probs, targets_bool, peps; doSort=true, fdr_scale_factor=1.0f0)
            allowed_pids = Set{UInt32}()
            for i in 1:length(peps)
                if peps[i] <= thr
                    push!(allowed_pids, UInt32(best_psms.precursor_idx[i]))
                end
            end
            filter_desc = "PEP≤$(thr)"
        end
        pass1_kept = length(allowed_pids)
        pass1_removed = pass1_precs_n - pass1_kept
        @user_info "  PEP-filter 2-pass (file_idx=$ms_file_idx, $file_name): " *
                   "pass1 PSMs=$pass1_psms_n → best-per-precursor=$pass1_precs_n; " *
                   "$filter_desc keeps $pass1_kept precursors, removes $pass1_removed " *
                   "($(round(100*pass1_removed/max(pass1_precs_n,1), digits=1))%)"

        if !isempty(allowed_pids)
            t_pass2_lib_start = time()
            psms2 = library_search(spectra, search_context, params, ms_file_idx;
                                    allowed_precursors=allowed_pids)
            t_pass2_lib = time() - t_pass2_lib_start
            pass2_psms_n = nrow(psms2)
            @user_info "  PEP-filter 2-pass (file_idx=$ms_file_idx): pass2 library_search returned " *
                       "$pass2_psms_n (precursor,scan) PSMs (was $pass1_psms_n in pass-1) " *
                       "in $(round(t_pass2_lib, digits=2))s"

            if nrow(psms2) > 0
                # Re-do feature prep on pass-2 PSMs (same flow as pass-1)
                prepare_psm_features!(psms2, params, search_context, ms_file_idx, spectra, prescore_only=true)
                add_scan_competition_features!(psms2)
                add_ms1_features!(psms2, spectra, search_context, ms_file_idx)

                # Re-train LGBM and re-select best per precursor on pass-2 PSMs
                best_psms2, scores2, _ = train_lgbm_and_select_best(psms2)
                best_psms2[!, :lgbm_prob] = scores2
                pass2_precs_n = nrow(best_psms2)
                @user_info "  PEP-filter 2-pass (file_idx=$ms_file_idx): pass2 LGBM → " *
                           "$pass2_precs_n best-per-precursor PSMs (was $pass1_precs_n in pass-1, " *
                           "$(pass2_precs_n - pass1_kept) added beyond surviving allowlist)"

                # Replace pass-1 best_psms with pass-2
                best_psms = best_psms2
            else
                @user_info "  PEP-filter 2-pass (file_idx=$ms_file_idx): pass2 produced 0 PSMs, keeping pass-1 results"
            end
        else
            @user_info "  PEP-filter 2-pass (file_idx=$ms_file_idx): no precursors survived PEP filter, keeping pass-1 results"
        end
    end
    # ============================================================

    # ============================================================
    # PAIR COMPETITION (default ON; disable with PIONEER_PAIR_COMPETITION=0)
    # For each (precursor, partner_decoy) pair where both rows are present
    # in best_psms, drop the lower-scoring partner. Operates on lgbm_prob
    # (post-pass-2 if 2-pass enabled, else pass-1). Equivalent to 1-vs-1
    # target-decoy competition before q-value computation.
    #
    # Empirically (Olsen Exploris one-file, non-entrap library): +2,556 q≤.001,
    # +1,131 q≤.01, +89 protein groups vs single-pass without paircomp.
    # ============================================================
    if get(ENV, "PIONEER_PAIR_COMPETITION", "1") != "0"
        precursors = getPrecursors(getSpecLib(search_context))
        partner_col = precursors.data[:partner_precursor_idx]
        n0 = nrow(best_psms)
        # Build score + row-index lookup keyed by precursor_idx (file is fixed)
        score_by_pid = Dict{UInt32, Float32}()
        row_by_pid   = Dict{UInt32, Int}()
        for i in 1:n0
            pid = UInt32(best_psms.precursor_idx[i])
            score_by_pid[pid] = Float32(best_psms.lgbm_prob[i])
            row_by_pid[pid]   = i
        end
        keep_mask = trues(n0)
        seen = Set{Tuple{UInt32, UInt32}}()
        n_dropped_t = 0; n_dropped_d = 0; n_pairs = 0
        for i in 1:n0
            pid = UInt32(best_psms.precursor_idx[i])
            ptr_raw = partner_col[Int(pid)]
            partner = ismissing(ptr_raw) ? UInt32(0) : UInt32(ptr_raw)
            (partner == 0 || partner == pid) && continue
            haskey(score_by_pid, partner) || continue
            pair_key = (min(pid, partner), max(pid, partner))
            pair_key in seen && continue
            push!(seen, pair_key)
            my_score = score_by_pid[pid]
            ptr_score = score_by_pid[partner]
            # Drop the lower one. Ties: drop the decoy (target wins).
            loser = if my_score < ptr_score
                pid
            elseif my_score > ptr_score
                partner
            else
                best_psms.target[i] ? partner : pid
            end
            j = row_by_pid[loser]
            keep_mask[j] = false
            n_pairs += 1
            if best_psms.target[j]; n_dropped_t += 1; else; n_dropped_d += 1; end
        end
        deleteat!(best_psms, .!keep_mask)
        @user_info "  Pair competition (file_idx=$ms_file_idx, $file_name): " *
                   "$n_pairs pairs found; dropped $n_dropped_t targets + $n_dropped_d decoys " *
                   "($(n0) → $(nrow(best_psms)) best-per-precursor PSMs)"
        _summarize_psm_counts(best_psms, "after paircomp", ms_file_idx, file_name)
    end
    # ============================================================

    # ============================================================
    # PER-FILE PEP FILTER (default ON at PEP ≤ 0.9)
    # PIONEER_MAIN_PEP_FILTER_THR sets the threshold; set ≥ 1.0 to disable.
    # Computed on best-per-precursor PSMs (post-paircomp) so each precursor
    # gets one PEP value. Precursors with PEP > threshold never enter the
    # second_pass arrow files and therefore can't reach ScoringSearch.
    # ============================================================
    pep_filter_thr = parse(Float32, get(ENV, "PIONEER_MAIN_PEP_FILTER_THR", "0.9"))
    if pep_filter_thr < 1.0f0 && nrow(best_psms) > 0
        probs_filt = Float32.(best_psms[!, :lgbm_prob])
        is_t_filt  = Vector{Bool}(best_psms[!, :target])
        peps_filt  = Vector{Float32}(undef, length(probs_filt))
        get_PEP!(probs_filt, is_t_filt, peps_filt; doSort=true, fdr_scale_factor=1.0f0)
        keep = peps_filt .<= pep_filter_thr
        n_before_pep = nrow(best_psms)
        n_drop_t = count(.!keep .& is_t_filt)
        n_drop_d = count(.!keep .& .!is_t_filt)
        deleteat!(best_psms, .!keep)
        @user_info "  PEP filter (file_idx=$ms_file_idx, $file_name): PEP > $pep_filter_thr drops " *
                   "$n_drop_t targets + $n_drop_d decoys " *
                   "($n_before_pep → $(nrow(best_psms)) best-per-precursor PSMs)"
        _summarize_psm_counts(best_psms, "after PEP filter", ms_file_idx, file_name)
    end
    # ============================================================

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
    @user_info "  MainSearch process_search_results! (file_idx=$ms_file_idx, $file_name): " *
               "$n_total_psms PSMs → $(nrow(best_psms)) precursors  total=$(r(t_total))s\n" *
               "    feature pass: prepare=$(r(t_prepare))s  competition=$(r(t_competition))s  " *
               "neighborhood=$(r(t_neighborhood))s  apex=$(r(t_apex))s  ms1=$(r(t_ms1))s\n" *
               "    lgbm=$(r(lgbm_timings.train_cv))s  recal=$(r(t_recal-t_lgbm))s  " *
               "phase2=$(r(t_phase2-t_recal))s  write=$(r(t_write-t_phase2))s"

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
        prec_mz_arr, entrap_group_ids,
        prec_mz_col::Vector{Float32},
        entrap_col::Vector{UInt8})
    @inbounds for i in eachindex(prec_idx_col)
        pid = prec_idx_col[i]
        prec_mz_col[i] = prec_mz_arr[pid]
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
    # PIONEER_GLOBAL_PRESCORE_FILTER (default "0" = report only, no filter).
    # When "1", the passing-set is enforced and PSMs of non-passing precursors are
    # dropped before ScoringSearch reads them.
    enforce_filter = get(ENV, "PIONEER_GLOBAL_PRESCORE_FILTER", "0") == "1"
    @user_info "  Global prescore: $n_pass_total precursors pass q≤$(fold0_result.qvalue_threshold) ($n_targets_total targets + $n_decoys_total decoys)" *
               (enforce_filter ? " [FILTER ENFORCED]" : " [report only, no filter — set PIONEER_GLOBAL_PRESCORE_FILTER=1 to enforce]")
    t1 = time() - t1_start

    store_results!(search_context, MainSearch, (passing_precs=passing_precs,))

    # Step 2: Load main_search_psms, filter to passing precursors, add deferred columns, write fold-split second_pass_psms
    t2_start = time()
    second_pass_dir = joinpath(getDataOutDir(search_context), "temp_data", "second_pass_psms")

    # Library lookups (shared across files)
    prec_irt = lib_irt
    prec_mz_arr = getMz(precursors)
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

            # Filter to this fold's passing precursors — only when enforce_filter is set.
            if enforce_filter
                mask = in.(tbl[!, :precursor_idx], Ref(passing))
                tbl = tbl[mask, :]
            end
            n_after = nrow(tbl)
            n_after_file += n_after
            n_kept_precs += n_after

            if n_after == 0
                isfile(fold_out_path) && rm(fold_out_path)
                continue
            end

            # Add Phase 2 library-lookup columns (deferred from process_search_results!)
            N = nrow(tbl)
            prec_mz_col = Vector{Float32}(undef, N)
            entrap_col = Vector{UInt8}(undef, N)
            _compute_phase2_columns!(
                tbl[!, :precursor_idx],
                prec_mz_arr, entrap_group_ids,
                prec_mz_col, entrap_col
            )
            tbl[!, :prec_mz] = prec_mz_col
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
