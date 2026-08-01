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
3. summarize_results!: Global prescore aggregation, register fold-split second_pass_psms for ScoringSearch
"""
# Per-file PEP threshold applied to best-per-precursor PSMs. Precursors with
# PEP > MAIN_PEP_FILTER_THR are dropped before being written to the second-pass
# arrows ScoringSearch reads. Set ≥ 1.0 to disable the filter entirely.
const MAIN_PEP_FILTER_THR = 0.99f0

# Minimum number of high-confidence target precursors required to fit the
# cross-fold iRT refinement model. Files with fewer than this many high-
# confidence target precursors skip the refinement pass.
const MAIN_IRT_REFINEMENT_MIN_PRECURSORS = 250

struct MainSearch <: SearchMethod end

"""
    _summarize_psm_counts(best_psms, stage_label, ms_file_idx, file_name)

Compute per-file q-values and PEPs from `best_psms[!, :lgbm_prob]` + `:target`
and emit a one-line @user_info summary with target counts at q≤.001, q≤.01,
PEP≤.01, PEP≤.05. **Diagnostic only** — purely transparent (no mutation of
best_psms columns).

Gated on `DEBUG_CONSOLE_LEVEL[] >= 1` because the q-value + PEP calculation
underneath (two sorts of up to ~3.98M rows each, called 3× per file) is
non-trivial — ~0.6 s/file × 8 files = ~5 s on Astral. The default debug
level is 0, so this is a no-op in production. Set
`logging.debug_console_level: 1` in the config to re-enable.
"""
function _summarize_psm_counts(best_psms::DataFrame, stage_label::AbstractString,
                                ms_file_idx::Integer, file_name::AbstractString)
    # Early-return before any computation. The @user_info macro itself
    # always prints, but the body's get_qvalues!/get_PEP! work is the
    # actual cost — gating function entry on the debug level skips both.
    Pioneer.DEBUG_CONSOLE_LEVEL[] >= 1 || return
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
    @debug_l1 "  [$stage_label] (file_idx=$ms_file_idx, $file_name): " *
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

`lgbm_buffers` holds the reusable backing stores for the per-file LightGBM
feature matrices. The file loop in `execute_search` is sequential, so one buffer
set serves every file — the matrices are population-scaled (hundreds of MB per
file), so allocating them per file is churn that grows with file count. Living
here scopes them to MainSearch: they are released with the results container
instead of being retained through the later search phases.
"""
struct MainSearchResults <: SearchResults
    psms::Base.Ref{DataFrame}
    lgbm_buffers::LGBMMatrixBuffers
end

#==========================================================
Interface Implementation
==========================================================#

get_parameters(::MainSearch, params::Any) = MainSearchParameters(params)

function init_search_results(::MainSearch, params::P, search_context::SearchContext) where {P<:MainSearchParameters}
    # Single MainSearch PSM directory for the whole pipeline. process_search_results!
    # writes per-fold files here; summarize_results! Step 2 overwrites them in
    # place with the filter+enrichment that used to land in a separate
    # second_pass_psms/ directory. PrecursorScoringSearch then loads from the same
    # directory, merges fold0+fold1 into a single per-file Arrow (also in this
    # directory), and folds MBR sidecars back in place.
    main_search_psms_dir = joinpath(getDataOutDir(search_context), "temp_data", "main_search_psms")
    mkpath(main_search_psms_dir)

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
        DataFrame(),
        LGBMMatrixBuffers()
    )
end

#==========================================================
Core Processing Methods
==========================================================#

"""
    _gather_into!(dst, src, perm)

`dst[i] = src[perm[i]]` for a concretely-typed column. Separate method so the loop specialises on the
column's element type rather than the `AbstractVector` that `eachcol` statically yields.
"""
@noinline function _gather_into!(dst::Vector{T}, src::Vector{T}, perm::Vector{Int}) where {T}
    @inbounds for i in eachindex(perm)
        dst[i] = src[perm[i]]
    end
    return nothing
end

"""
    permute_psms_by_precursor_idx!(psms::DataFrame) -> DataFrame

Sort `psms` in place by `:precursor_idx` with a hand-rolled column-wise gather. Much faster than
`DataFrames.sort!(psms, :precursor_idx)` on the post-deconv DataFrame shape (~14M rows × 27 cols),
which allocates a fresh array per column (measured 2026-05-19).

Each isbits column is gathered through a single shared byte buffer rather than letting
`Base.permute!` allocate its own temporary per column — see the comments in the body for why that is
both necessary and safe for peak memory.

Establishes the sorted-by-precursor invariant that downstream passes
(`_build_precursor_groups`, feature passes, `select_best_per_precursor!`)
can take advantage of to skip their own sortperm.
"""

function permute_psms_by_precursor_idx!(psms::DataFrame)
    n = nrow(psms)
    n == 0 && return psms
    perm = sortperm(psms[!, :precursor_idx]::Vector{UInt32})
    # NOT `Base.permute!`. Its docstring contract changed: it no longer destroys the permutation
    # argument (nor does `permute!!`) -- it allocates a *data*-sized temporary instead, measured at
    # exactly `sizeof(col)` per call. So the old code's defensive `copyto!(p_scratch, perm)` guarded a
    # hazard that no longer exists, while one full copy of every column went unnoticed: 17.3 GB per
    # 6-file KEAP1 run, third-largest allocator in Main Search.
    #
    # A true O(1)-space cycle permutation removes the allocation entirely but is 28x SLOWER
    # (0.0530 s vs 0.0019 s on 2M Float32) -- chasing random cycles is cache-hostile, which is
    # presumably why Base allocates. So keep the sequential gather and reuse one buffer instead.
    #
    # One raw byte buffer, reinterpreted per column type: exactly one column-sized buffer is live at
    # any moment, which is the same footprint `Base.permute!` allocated transiently. Peak RSS
    # therefore cannot rise (measured: +0.0 MB vs +7.7 MB), while cumulative allocation falls ~88%.
    # Only isbits columns may use the raw buffer. This table still carries Vector-valued columns at
    # this point (they are removed later by `dropVectorColumns!`), and wrapping uninitialised bytes
    # as a non-isbits element type would hand the GC garbage pointers. Those fall back to
    # `Base.permute!`, which is correct, merely allocating.
    max_bytes = 0
    for c in eachcol(psms)
        T = eltype(c)
        isbitstype(T) && c isa Vector{T} && (max_bytes = max(max_bytes, sizeof(T) * n))
    end
    raw = max_bytes > 0 ? Vector{UInt8}(undef, max_bytes) : UInt8[]
    GC.@preserve raw begin
        for col in eachcol(psms)
            T = eltype(col)
            if isbitstype(T) && col isa Vector{T}
                # `raw` is preserved for the whole loop and `buf` never escapes this iteration.
                buf = unsafe_wrap(Vector{T}, Ptr{T}(pointer(raw)), n)
                _gather_into!(buf, col, perm)
                copyto!(col, buf)
            else
                Base.permute!(col, copy(perm))
            end
        end
    end
    return psms
end

"""
    @alloc_bucket(label, expr)

Evaluate `expr`; when `PIONEER_DIAG_ALLOC` is set, print the heap bytes it
allocated (process-global counter delta — sums across all threads, so it is
valid for the multithreaded deconv). Beyond two `gc_bytes()` reads it adds no
overhead when the env var is unset.
"""
macro alloc_bucket(label, ex)
    quote
        local _b0 = Base.gc_bytes()
        local _r = $(esc(ex))
        if haskey(ENV, "PIONEER_DIAG_ALLOC")
            @user_print string("[ALLOC] ", $(esc(label)), ": ",
                               round((Base.gc_bytes() - _b0) / 1e9, digits=3), " GB")
        end
        _r
    end
end

"""
    _log_psm_table_footprint(psms, label, ms_file_idx)

Diagnostic: report the in-memory footprint of the post-deconv PSMs table
(rows x cols, per-column type + bytes, total). Gated behind the
`PIONEER_DIAG_PSM_FOOTPRINT` env var so it costs nothing when off. Lets us
compare the *irreducible* data size against the phase's cumulative allocation
to gauge how much of Main Search's allocation is avoidable churn.
"""
function _log_psm_table_footprint(psms::DataFrame, label::AbstractString, ms_file_idx)
    haskey(ENV, "PIONEER_DIAG_PSM_FOOTPRINT") || return nothing
    n = nrow(psms)
    total = 0
    @user_print "[PSM FOOTPRINT] file=$ms_file_idx $label: $n rows x $(ncol(psms)) cols"
    for nm in names(psms)
        col = psms[!, nm]
        et = eltype(col)
        b = isbitstype(et) ? n * sizeof(et) : Base.summarysize(col)
        total += b
        @user_print "    $(rpad(nm, 34)) $(rpad(string(et), 26)) $(round(b / 1e6, digits=2)) MB"
    end
    @user_print "[PSM FOOTPRINT] file=$ms_file_idx total table = $(round(total / 1e9, digits=3)) GB " *
                "($(round(total / max(n,1)))  bytes/row)"
    return nothing
end

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

    psms = @alloc_bucket "library_search (deconv)" library_search(spectra, search_context, params, ms_file_idx)
    t_lib_search = time() - t_file_start

    # IMPORTANT: the next two steps depend on the deconv output being
    # contiguous-by-:scan_idx — all PSMs sharing a scan_idx form one
    # contiguous run. This is the natural emit order from the per-thread
    # `for scan_idx in scan_range` loop in process_scans_fused.jl:67
    # (threads are scan-disjoint, so no cross-thread interleaving).
    # Verified empirically 2026-05-19 on 3-file Astral: all scans
    # contiguous across ~48M PSMs / 758k unique scans.

    # Per-scan competition (weight_ratio_at_scan, weight_rank_at_scan).
    # Done BEFORE the precursor sort so it can exploit the
    # contiguous-by-scan invariant: linear sweep for run boundaries
    # + threaded per-run rank/ratio. ~4× faster than the previous
    # Dict-based version (measured 2026-05-19).
    t_scan_comp = @elapsed @alloc_bucket "scan_competition_features" add_scan_competition_features!(psms)

    # MS1 lookup features (ms1_m0_intensity, ms1_m1_intensity,
    # ms1_m0_mass_err_ppm, ms1_m1_to_m0_ratio, ms1_m1_to_m0_pred). Done
    # BEFORE the precursor sort so that consecutive rows within a per-chunk
    # task share scan_idx, which lets the per-task MS1-spectrum cache hit
    # ~once per unique scan in the chunk instead of being thrashed by
    # precursor-sorted input. The per-precursor chromatogram-feature passes
    # (ms1_corr_*, frag_*) run later in process_search_results! after the
    # precursor sort, since they group by :precursor_idx.
    t_ms1 = @elapsed @alloc_bucket "ms1_lookup_features" add_ms1_lookup_features!(psms, spectra, search_context, ms_file_idx)

    # Sort the deconv-output DataFrame by :precursor_idx once. Downstream
    # passes (chrom features, best-per-precursor) can then fast-path their
    # group-build (no second sortperm), and the data layout is friendlier
    # for any per-precursor parallelism. We use a hand-rolled in-place
    # column permute rather than `sort!(df, :col)` because DataFrames.sort!
    # is ~4× slower on this shape (measured 2026-05-19).
    t_sort = @elapsed @alloc_bucket "permute_by_precursor" permute_psms_by_precursor_idx!(psms)

    results.psms[] = psms

    _log_psm_table_footprint(psms, "post-deconv (before best-per-precursor)", ms_file_idx)

    @debug_l1 "  MainSearch process_file! (file_idx=$ms_file_idx, $file_name): " *
               "$(nrow(psms)) PSMs from deconv; library_search elapsed: $(round(t_lib_search, digits=2))s  " *
               "scan_comp=$(round(t_scan_comp * 1000, digits=0))ms  " *
               "ms1=$(round(t_ms1, digits=2))s  " *
               "sort=$(round(t_sort * 1000, digits=0))ms  n_cols=$(ncol(psms))"

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
    center_mzs = getCenterMzs(spectra)
    isolation_widths = getIsolationWidthMzs(spectra)

    # Compute prescore features
    t_prepare = @elapsed @alloc_bucket "prepare_psm_features" prepare_psm_features!(psms, params, search_context, ms_file_idx, spectra)
    t_features = time()

    if nrow(psms) == 0
        @debug_l2 "No PSMs for file $ms_file_idx after feature computation"
        return nothing
    end

    # add_scan_competition_features! moved upstream to process_file! (before
    # the precursor-idx sort) so it can exploit the contiguous-by-scan
    # invariant from the deconv output. Was ~3s/file in this position;
    # now ~0.5s/file in process_file!.
    t_competition = 0.0
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
    # MS1 spectrum lookup moved upstream to process_file! (before precursor
    # sort) so the per-chunk MS1 cache exploits contiguous-by-scan input.
    # Only the precursor/window chromatogram features still run here.
    bitvec_rank_table = getBitVecExcessRanks(search_context, Int64(ms_file_idx))
    t_ms1 = @elapsed @alloc_bucket "chromatogram_features" add_chromatogram_features!(
        psms,
        spectra;
        bitvec_rank_table = bitvec_rank_table,
    )

    # Train LightGBM on ALL PSMs, select best scan per precursor
    n_total_psms = nrow(psms)
    _log_psm_table_footprint(psms, "full pre-reduction (after all feature passes)", ms_file_idx)
    Pioneer.DIAG_DUMP_FILE_IDX[] = 0
    t_lgbm_start = time()
    best_psms, scores, lgbm_timings, lgbm_predictor =
        @alloc_bucket "train_lgbm_and_select_best" train_lgbm_and_select_best(
            psms;
            center_mzs = center_mzs,
            isolation_widths = isolation_widths,
            buffers = results.lgbm_buffers,
        )
    best_psms[!, :lgbm_prob] = scores
    _summarize_psm_counts(best_psms, "after best-per-precursor", ms_file_idx, file_name)
    t_lgbm_end = time()

    # Refine predicted iRTs with out-of-fold correction models. The correction
    # changes iRT-dependent features for every candidate PSM, so reapply the
    # same per-run classifier before reducing to one best scan per precursor.
    precursors = getPrecursors(getSpecLib(search_context))
    irt_refinement = MainSearchIrtRefinement(
        precursors;
        q_value_threshold = PRESCORE_QVALUE_THRESHOLD,
        min_precursors = MAIN_IRT_REFINEMENT_MIN_PRECURSORS,
    )
    irt_refinement_result = refine_mainsearch_irt_predictions!(psms, best_psms, scores, irt_refinement)
    if irt_refinement_result.refined
        best_psms, scores, reapply_timings = reapply_psm_classifier_and_select_best!(
            psms,
            lgbm_predictor;
            center_mzs = center_mzs,
            isolation_widths = isolation_widths,
            buffers = results.lgbm_buffers,
        )
        best_psms[!, :lgbm_prob] = scores
        @debug_l1 "  iRT refinement (file_idx=$ms_file_idx, $file_name): " *
                   "$(length(irt_refinement_result.training_target_precursors)) training precursors; " *
                   "reapply predict=$(round(reapply_timings.predict, digits=2))s best=$(round(reapply_timings.best, digits=2))s"
        _summarize_psm_counts(best_psms, "after refined-iRT reapply", ms_file_idx, file_name)
    else
        @debug_l1 "  iRT refinement (file_idx=$ms_file_idx, $file_name): skipped " *
                   "($(length(irt_refinement_result.training_target_precursors)) " *
                   "high-confidence target precursors; need $(irt_refinement.min_precursors))"
    end

    _summarize_psm_counts(best_psms, "before PEP filter", ms_file_idx, file_name)
    t_pep_start = time()
    _b_pep = Base.gc_bytes()

    # ============================================================
    # PER-FILE PEP FILTER (PEP ≤ MAIN_PEP_FILTER_THR).
    # Computed on best-per-precursor PSMs so each precursor gets one PEP value.
    # Precursors with PEP > threshold never enter the second_pass arrow files
    # and therefore can't reach ScoringSearch.
    # ============================================================
    pep_filter_thr = MAIN_PEP_FILTER_THR
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
        @debug_l1 "  PEP filter (file_idx=$ms_file_idx, $file_name): PEP > $pep_filter_thr drops " *
                   "$n_drop_t targets + $n_drop_d decoys " *
                   "($n_before_pep → $(nrow(best_psms)) best-per-precursor PSMs)"
        _summarize_psm_counts(best_psms, "after PEP filter", ms_file_idx, file_name)
    end
    haskey(ENV, "PIONEER_DIAG_ALLOC") && @user_print string("[ALLOC] pep_filter: ",
        round((Base.gc_bytes() - _b_pep) / 1e9, digits=3), " GB")
    t_pep_end = time()
    t_recal_start = t_pep_end
    # ============================================================

    # RT recalibration: refit iRT spline from high-confidence PSMs
    @alloc_bucket "recalibrate_rt" recalibrate_rt!(search_context, ms_file_idx, best_psms, best_psms[!, :lgbm_prob])

    # Update irt_obs/irt_error with post-recalibration model so downstream
    # steps (ScoringSearch features, RT index, chromatogram extraction) are consistent
    new_rt_model = getRtIrtModel(search_context, ms_file_idx)
    best_psms[!, :irt_obs] .= new_rt_model.(best_psms[!, :rt])
    best_psms[!, :irt_error] .= abs.(best_psms[!, :irt_obs] .- best_psms[!, :irt_pred])
    t_recal = time()

    trace_peps, trace_pass_mask = _mainsearch_peps_and_pass_mask(
        psms[!, :lgbm_score],
        psms[!, :target],
    )
    @alloc_bucket "precursor_fraction_transmitted" add_precursor_fraction_transmitted!(
        best_psms,
        getQuadTransmissionModel(search_context, ms_file_idx),
        getSearchData(search_context),
        getCharge(precursors),
        getMz(precursors),
        getSulfurCount(precursors),
        center_mzs,
        isolation_widths,
    )

    # Filter by precursor_fraction_transmitted
    to_remove = findall(best_psms[!, :precursor_fraction_transmitted] .< params.min_fraction_transmitted)
    deleteat!(best_psms, to_remove)
    @alloc_bucket "trace_and_fragment_features" add_trace_and_fragment_features!(
        best_psms,
        psms,
        trace_pass_mask;
        bitvec_rank_table = bitvec_rank_table,
        center_mzs = center_mzs,
        isolation_widths = isolation_widths,
        pep_values = trace_peps,
    )
    best_psms[!, :ms_file_idx] .= UInt32(ms_file_idx)
    t_phase2 = time()

    # Drop vector columns that can't be serialized to Arrow
    dropVectorColumns!(best_psms)

    precursors = getPrecursors(getSpecLib(search_context))
    frag_lookup = getFragmentLookupTable(getSpecLib(search_context))
    nce_model = getNceModel(search_context, ms_file_idx)
    @alloc_bucket "wide_window_features" add_wide_window_features_to_table!(
        best_psms,
        spectra,
        search_context,
        ms_file_idx,
        precursors,
        frag_lookup,
        nce_model,
    )

    # Write per-fold main_search_psms.
    main_search_psms_dir = joinpath(getDataOutDir(search_context), "temp_data", "main_search_psms")
    _b_write = Base.gc_bytes()
    best_cv_fold = UInt8[getCvFold(precursors, pid) for pid in best_psms.precursor_idx]
    for fold in UInt8[0, 1]
        fold_mask = best_cv_fold .== fold
        if any(fold_mask)
            fold_path = joinpath(main_search_psms_dir, "$(file_name)_fold$(fold).arrow")
            writeArrow(fold_path, best_psms[fold_mask, :])
        end
    end
    haskey(ENV, "PIONEER_DIAG_ALLOC") && @user_print string("[ALLOC] write_fold_arrow: ",
        round((Base.gc_bytes() - _b_write) / 1e9, digits=3), " GB")
    t_write = time()

    # Timing summary — durations sum to ~t_total (within ~0.1s of bookkeeping).
    # process_file!'s t_lib_search (fragment-index + deconv) is reported
    # separately in that function; this summary covers process_search_results!.
    # t_recal here is the END-OF-RECAL timestamp from earlier; rename to keep
    # it distinct from the recal-duration we compute below.
    t_recal_end = t_recal
    t_total = t_write - t_start
    dur_features  = t_prepare + t_competition + t_apex + t_ms1
    dur_lgbm      = t_lgbm_end - t_lgbm_start   # train_cv + best-per-precursor
    dur_pep       = t_pep_end - t_pep_start
    dur_recal     = t_recal_end - t_recal_start
    dur_phase2    = t_phase2 - t_recal_end
    dur_write     = t_write - t_phase2
    dur_accounted = dur_features + dur_lgbm + dur_pep +
                    dur_recal + dur_phase2 + dur_write
    dur_overhead  = t_total - dur_accounted
    r = s -> round(s, digits=2)
    @debug_l1 "  MainSearch process_search_results! (file_idx=$ms_file_idx, $file_name): " *
               "$n_total_psms PSMs → $(nrow(best_psms)) precursors  total=$(r(t_total))s\n" *
               "    features=$(r(dur_features))s [prep=$(r(t_prepare))s competition=$(r(t_competition))s ms1=$(r(t_ms1))s]\n" *
               "    lgbm=$(r(dur_lgbm))s [train_cv=$(r(lgbm_timings.train_cv))s best_per_prec=$(r(lgbm_timings.best))s]  " *
               "pep_filter=$(r(dur_pep))s\n" *
               "    recal=$(r(dur_recal))s  phase2=$(r(dur_phase2))s  fold_write+flanking=$(r(dur_write))s  " *
               "overhead=$(r(dur_overhead))s"

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
        irt_pred_col::AbstractVector{Float32},
        prec_irt, prec_mz_arr, prec_pair_idxs, entrap_group_ids,
        irt_diff_col::Vector{Float32},
        prec_mz_col::Vector{Float32},
        pair_id_col::Vector{UInt32},
        entrap_col::Vector{UInt8})
    @inbounds for i in eachindex(prec_idx_col)
        pid = prec_idx_col[i]
        irt_diff_col[i] = abs(irt_obs_col[i] - irt_pred_col[i])
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

    # Step 1: Per-fold global prescore aggregation → RT-binned tolerance only.
    # No PSM filter is applied here; the per-file PEP filter upstream already
    # gates what reaches ScoringSearch.
    t1_start = time()
    fold0_result = @alloc_bucket "aggregate_prescore_fold0" aggregate_prescore_globally!(search_context; fold_suffix="_fold0")
    fold1_result = @alloc_bucket "aggregate_prescore_fold1" aggregate_prescore_globally!(search_context; fold_suffix="_fold1")
    rt_binned_tol = fold0_result.rt_binned_tol !== nothing ?
                    fold0_result.rt_binned_tol : fold1_result.rt_binned_tol
    t1 = time() - t1_start

    # Step 2: Load main_search_psms fold-split files, filter to passing precursors,
    # add deferred columns, OVERWRITE each fold file in place via writeArrow.
    # The fold-split layout is preserved here because PrecursorScoringSearch's
    # Pass-1 LGBM CV needs to load one fold's rows at a time. PrecursorScoring
    # later merges fold0+fold1 into a single per-file Arrow (also in
    # main_search_psms_dir) once LGBM training is done.
    t2_start = time()
    _b_t2 = Base.gc_bytes()

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

        base_path = joinpath(main_search_psms_dir, file_name)

        file_has_data = false
        n_before_file = 0
        n_after_file = 0

        for fold in UInt8[0, 1]
            psm_path = joinpath(main_search_psms_dir, "$(file_name)_fold$(fold).arrow")

            if !isfile(psm_path)
                continue
            end

            # Load this fold's main search PSMs into in-memory DataFrame
            # (Tables.columntable + DataFrame materializes columns off the
            # Arrow mmap so the subsequent in-place writeArrow is safe on
            # Windows — same pattern as ArrowOperations.jl:68).
            tbl = DataFrame(Tables.columntable(Arrow.Table(psm_path)))
            n_before = nrow(tbl)
            n_before_file += n_before
            n_total_precs += n_before

            n_after = nrow(tbl)
            n_after_file += n_after
            n_kept_precs += n_after

            if n_after == 0
                # Drop the on-disk file so downstream code skips this fold.
                safeRm(psm_path, nothing)
                continue
            end

            # Add Phase 2 library-lookup columns (deferred from process_search_results!)
            N = nrow(tbl)
            irt_diff_col = Vector{Float32}(undef, N)
            prec_mz_col = Vector{Float32}(undef, N)
            pair_id_col = Vector{UInt32}(undef, N)
            entrap_col = Vector{UInt8}(undef, N)
            _compute_phase2_columns!(
                tbl[!, :precursor_idx], tbl[!, :irt_obs], tbl[!, :irt_pred],
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

            # Overwrite the same fold file in place. writeArrow handles the
            # tempname → move dance internally for Windows safety.
            writeArrow(psm_path, tbl)
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
    haskey(ENV, "PIONEER_DIAG_ALLOC") && @user_print string("[ALLOC] summarize_global_filter: ",
        round((Base.gc_bytes() - _b_t2) / 1e9, digits=3), " GB")
    t2 = time() - t2_start

    overall_pct = round(100.0 * n_kept_precs / max(1, n_total_precs), digits=1)
    @user_info "MainSearch passing PSMs: $n_kept_precs / $n_total_precs precursors ($overall_pct%) across $n_processed_files files"

    # Step 3: Compute the RT-binned chromatographic tolerance used when
    # extracting chromatograms around each precursor's refined RT.
    t3_start = time()
    if n_processed_files > 0
        if rt_binned_tol !== nothing
            compute_rt_binned_tolerance!(search_context, rt_binned_tol, ms_data, n_files)
        else
            @debug_l1 "No RT-binned tolerance registered — chromatogram extraction will fall back to the per-file iRT tolerance from recalibrate_rt!"
        end
    else
        @warn "No files processed in MainSearch — skipping chromatographic tolerance computation"
    end
    t3 = time() - t3_start

    t_total = time() - t_total_start
    @debug_l1 "MainSearch summarize: aggregation=$(r(t1))s, fold_write=$(r(t2))s, chrom_tol=$(r(t3))s, total=$(r(t_total))s"

    return nothing
end
