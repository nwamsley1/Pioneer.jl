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
    PrecursorScoringSearch

Per-precursor rescoring + FDR control. Trains a second LightGBM on globally-
filtered MainSearch PSMs (20 features), recalculates q-values + PEPs, filters
to a passing-precursor set, and builds RT indices for the chromatogram
integrator. Protein inference + protein-group scoring are now downstream
(`ProteinInferenceSearch`, `ProteinScoringSearch`).
"""
struct PrecursorScoringSearch <: SearchMethod end

# Note: FileReferences, SearchResultReferences, and FileOperations are already
# included by importScripts.jl - no need to include them here

#==========================================================
Type Definitions 
==========================================================#

"""
Results container for scoring search.
"""
struct PrecursorScoringSearchResults <: SearchResults
    precursor_global_qval_dict::Base.Ref{Dict{UInt32, Float32}}  # precursor_idx → global q-value
    precursor_qval_interp::Base.Ref{Any}  # Interpolation for run-specific q-values
    precursor_pep_interp::Base.Ref{Any}   # Interpolation for experiment-wide PEPs
    run_similarity::Base.Ref{Union{Nothing, RunSimilarityAtlas}}
    merged_quant_path::String             # Path to merged quantification results
end

"""
    FORCE_OOM

Developer toggle: force the out-of-memory scoring path even when PSMs would
fit in memory. Used to test/profile the OOM path without a deliberately
oversized dataset. Production runs leave this `false`.
"""
const FORCE_OOM = false

"""
Parameters for scoring search.
"""
struct PrecursorScoringSearchParameters <: SearchParameters
    # In-memory PSM memory budget; OOM path triggers when bytes-per-row * Q
    # exceeds this.
    max_psm_memory_mb::Float64

    # PSMs per bin for empirical q-value/PEP interpolation in get_qvalue_spline.
    # Smaller = finer-grained but noisier per-bin FDR estimates; larger =
    # smoother but coarser.
    pep_bin_size::Int64

    q_value_threshold::Float32

    # When false, skip post-integration donor/counterfactual feature
    # computation and transfer rescoring. Driven by global.match_between_runs.
    match_between_runs::Bool

    function PrecursorScoringSearchParameters(params::PioneerParameters)
        ml_params = params.optimization.machine_learning
        global_params = params.global_settings

        mbr = hasproperty(global_params, :match_between_runs) ?
                Bool(global_params.match_between_runs) : true

        new(
            Float64(ml_params.max_psm_memory_mb),
            Int64(ml_params.pep_bin_size),
            _resolve_q_value_threshold(global_params),
            mbr,
        )
    end
end

#==========================================================
Interface Implementation  
==========================================================#

get_parameters(::PrecursorScoringSearch, params::Any) = PrecursorScoringSearchParameters(params)

function init_search_results(::PrecursorScoringSearchParameters, search_context::SearchContext)
    return PrecursorScoringSearchResults(
        Ref(Dict{UInt32, Float32}()),  # precursor_global_qval_dict
        Ref(undef),  # precursor_qval_interp
        Ref(undef),  # precursor_pep_interp
        Ref{Union{Nothing, RunSimilarityAtlas}}(nothing),
        joinpath(getDataOutDir(search_context), "merged_quant.arrow")
    )
end

function process_file!(
    results::PrecursorScoringSearchResults,
    params::PrecursorScoringSearchParameters,
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::MassSpecData
)
    # No per-file processing needed
    return results
end

function process_search_results!(
    results::PrecursorScoringSearchResults,
    params::PrecursorScoringSearchParameters,
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::MassSpecData
)
    # No per-file results processing needed
    return nothing
end

function reset_results!(results::PrecursorScoringSearchResults)
    return nothing
end


"""
    estimate_max_rows(memory_mb::Float64, sample_file::String)

Estimate the maximum number of rows that fit in `memory_mb` given the schema of `sample_file`.
"""
function estimate_max_rows(memory_mb::Float64, sample_file::String)
    tbl = Arrow.Table(sample_file)
    bytes_per_row = 0
    for name in Tables.columnnames(tbl)
        col = Tables.getcolumn(tbl, name)
        T = eltype(col)
        T_inner = Base.nonmissingtype(T)
        if isbitstype(T_inner)
            bytes_per_row += sizeof(T_inner)
        else
            bytes_per_row += 64  # conservative estimate for strings / non-bits types
        end
        if T !== T_inner
            bytes_per_row += 1  # missing indicator
        end
    end
    bytes_per_row = max(bytes_per_row, 1)
    return max(floor(Int64, memory_mb * 1024 * 1024 / bytes_per_row), 1000)
end

"""
Process all results to get final protein scores.
"""
const PRECSCORE_DIAG = Dict{Symbol, Int}()

"""
    PRECSCORE_DROPPABLE_COLUMNS

Per-PSM ML feature columns with no reader after PrecursorScoring. Determined by scanning every stage
downstream of this one (MBR, IntegrateChromatograms, ProteinInference, ProteinScoring, MaxLFQ,
WriteOutputs, SearchDIA) for literal references; 76 of the 141 columns in precursors_long had none.

Dropped inside the existing `initial_filter` pass, which already reads, filters and rewrites every file,
so removal is free rather than costing an extra pass. Everything downstream then reads and writes a
narrower table -- `_merge_mbr_recoveries!` alone rewrites every file at ~2.45 GB of writeArrow.

Deliberately RETAINED despite having no downstream reader, because they are plausibly wanted by whoever
reads the output: best_rt, irt_fwhm, rt_fwhm, smoothness.

Four more were added after inspecting the actual output values:
- `q_value` is initialised to zeros by `initialize_prob_group_features!` and never populated on this
  table -- it read 0 in all 235,194 rows while `qval` carried the real values. Worse than a duplicate:
  a reader would take q_value == 0 to mean every PSM is perfect.
- `decoy` is exactly `!target` (0 mismatches in 235,194 rows).
- `lgbm_prob` is bit-identical to `lgbm_score` (0 differing rows, same min/max). It stays in the
  pipeline -- MainSearch reads it for prescore aggregation and recalibrate_rt -- but not in the output.
- `trace_prob_prepass` is bit-identical to `trace_prob`. MBR consumes it, so it must survive to this
  point; dropping here is after its last use.
"""
const PRECSCORE_DROPPABLE_COLUMNS = Symbol[
    :q_value, :decoy, :lgbm_prob, :trace_prob_prepass,
    # Chromatographic / spectral ML features. Their last readers are the ScoringSearch feature list
    # (model_config.jl) and the MBR donor-feature list (MBR/types.jl, MBR/features.jl for :n_scans),
    # all of which run before this finalize loop. No reference in ProteinInference, ProteinScoring,
    # MaxLFQ, the QC plots, or anywhere in src/utils or src/structs.
    :weight_ratio_at_scan, :weight_rank_at_scan, :gof, :poisson, :err_norm,
    :fitted_hellinger, :fitted_manhattan_distance, :smoothed_2d_shadow_hellinger,
    :log_by_ratio_m0, :smoothness, :n_scans, :ms1_weight_apex_to_m0_apex_irt,
    # Exact elementwise duplicate of :rt (verified over 235,194 rows); :rt is the one the QC plots
    # read, so :best_rt is the copy that goes.
    :best_rt,
    # Internal CV / scoring artefacts and ML features with no meaning outside the scorer.
    :lgbm_score, :trace_prob, :trace_prob_infold, :log2_intensity_explained,
    # length(sequence), which is already in the table.
    :sequence_length,
    # MBR false-transfer-rate diagnostics, for method development rather than end users.
    # :mbr_recovered survives and is the flag that says whether a row was transferred.
    :mbr_target_decoy_prob, :mbr_total_error_qval_true, :mbr_total_error_rate_true,
    :mbr_counterfactual_decoy_prob,
    :mbr_counterfactual_decoy_index,
    :ftr_qval_true, :ftr_pep_true,
    # Same physical measurement as :rt_fwhm but in library iRT units, which are only interpretable
    # relative to the library. :rt_fwhm is kept, in minutes.
    :irt_fwhm,
    :longest_y, :y_count, :total_ions, :error, :max_matched_residual, :max_unmatched_residual,
    :frag1_int, :frag2_int, :frag3_int, :frag4_int, :frag5_int, :frag6_int, :frag7_int,
    :frag8_int, :top3_ms2_mass_error_mean, :cycle_idx, :ms1_m0_mass_err_ppm, :ms1_m0_intensity,
    :ms1_m1_intensity, :ms1_m1_to_m0_ratio, :ms1_m1_to_m0_pred, :ms1_m0_m1_m2_window_fraction,
    :ms1_ms2_explained_delta, :ms1_m0_m1_m2_window_fraction_pc, :ms1_ms2_explained_delta_pc,
    :ms1_isotope_dotp_m0_m1_m2, :ms1_m0_peak_frag_intensity_fraction,
    :ms1_m0_peak_n_precursors, :scan_prec_mz_n_precursors, :charge2, :spectrum_peak_count,
    :ms1_corr_weight_m0, :ms1_corr_m0_m1, :ms1_apex_offset_irt, :frag_apex_dispersion_irt,
    :n_correlated_fragments, :n_correlated_fragments_bitvec_rank, :frag_corr_strength,
    :frag_corr_effective_n, :frag_corr_best_m0, :delta_frame_peak_center, :n_above_hm,
    :n_contiguous_scans, :n_frags_detected_union, :n_frags_detected_intersection,
    :n_frags_detected_union_bitvec_rank, :n_frags_detected_intersection_bitvec_rank,
    :n_scans_other_windows, :other_window_weight_corr, :other_window_apex_delta_irt,
    :frag1_smoothed_intensity, :frag2_smoothed_intensity, :frag3_smoothed_intensity,
    :frag4_smoothed_intensity, :frag5_smoothed_intensity, :frag6_smoothed_intensity,
    :frag7_smoothed_intensity, :frag8_smoothed_intensity, :flanking_ms1_m0_candidate_fraction,
    :flanking_frag_candidate_fraction, :flanking_ms1_frag_sum_corr, :flanking_frag_corr_mean,
    :flanking_frag_corr_strength, :flanking_frag_corr_effective_n, :flanking_frag_corr_best_m0,
    :flanking_signal_support, :frag_apex_gt2x_flank_bitvec_rank, :irt_diff, :pair_id
]

"""
    OUTPUT_NARROWED_COLUMNS

Columns converted from Float32 to Float16 in the same finalize rewrite that applies
`PRECSCORE_DROPPABLE_COLUMNS`, i.e. after every reader has run.

Float16 carries 10 mantissa bits, so ~0.05% *relative* precision at any magnitude. That makes it
unsafe for absolute retention times -- at rt = 100 min the representable spacing is 0.0625 min
(3.75 s), so `:rt` stays Float32 -- but safe for these two, which are small and already coarser
than Float16's spacing:

- `:rt_fwhm` is `max(rt) - min(rt)` over the scans above half apex weight (MainSearch/scoring.jl),
  so it is quantised to the cycle time (~1 s). Float16 spacing at 1 min is 0.03 s, ~30x finer than
  the underlying measurement. (Its 11% exact zeros are real: one scan above half max, i.e. a peak
  narrower than one cycle -- not a missing value.)
- `:irt_error` maxes out near 4.0 iRT, where Float16 spacing is 0.002 iRT, far below any tolerance
  the value is compared against.

Note that changing units would not help: Float16 is floating point, so scaling by 60 shifts the
exponent and leaves relative precision unchanged.
"""
const OUTPUT_NARROWED_COLUMNS = (:rt_fwhm, :irt_error)

"""
    narrow_output_columns!(psms)

Cast `OUTPUT_NARROWED_COLUMNS` to Float16 in place. No-op for columns that are absent or already
narrowed.
"""
function narrow_output_columns!(psms)
    for col in OUTPUT_NARROWED_COLUMNS
        hasproperty(psms, col) || continue
        v = psms[!, col]
        v isa AbstractVector{Float32} || continue
        psms[!, col] = Vector{Float16}(v)
    end
    return psms
end


"""
    _precscore_diag_report()

Per-phase allocation/time attribution for the Precursor Scoring stage (~20 GB, previously the largest
wholly uninstrumented stage). Gated on PIONEER_PRECSCORE_DIAG; the phases are recorded by `_pmark`
inside `summarize_results!`.

Each phase name describes the work that ENDED at that mark, so the boundaries are the major calls:
fold_merge -> aggregate -> scoring -> run_similarity -> qvals_and_filtering -> rt_indices.
`@alloc_bucket` (MainSearch) cannot be used here: importScripts.jl loads PrecursorScoringSearch before
MainSearch, so that macro is not yet defined at parse time.
"""
function _precscore_diag_report()
    d = PRECSCORE_DIAG
    # Order = execution order. Each name is the phase that CLOSES at its mark; an earlier version had
    # these off by one (each was named for the phase it OPENED), which put score_precursor_isotope_traces
    # inside the fold-merge bucket and made LGBM scoring look like 3 s when it was not measured at all.
    ks = [:setup, :scoring, :fold_merge, :aggregate, :run_similarity, :qvals_and_filtering, :rt_indices]
    tot = sum(get(d, Symbol(k, :_bytes), 0) for k in ks)
    totms = sum(get(d, Symbol(k, :_ms), 0) for k in ks)
    gb(x) = round(x / 2^30, digits = 2)
    lines = ["Precursor Scoring phase diagnostic:"]
    for k in ks
        b = get(d, Symbol(k, :_bytes), 0); m = get(d, Symbol(k, :_ms), 0)
        (b == 0 && m == 0) && continue
        pc = tot > 0 ? round(100 * b / tot, digits = 1) : 0.0
        rate = m > 0 ? round((b / 2^30) / (m / 1000), digits = 2) : 0.0
        push!(lines, "  $(rpad(string(k), 22)) $(lpad(gb(b), 7)) GB  $(lpad(m, 7)) ms  ($(lpad(pc, 5))%)  $(rate) GB/s")
    end
    push!(lines, "  $(rpad("TOTAL", 22)) $(lpad(gb(tot), 7)) GB  $(lpad(totms, 7)) ms")
    @user_info join(lines, "\n")
    return nothing
end

function summarize_results!(
    results::PrecursorScoringSearchResults,
    params::PrecursorScoringSearchParameters,
    search_context::SearchContext
)
    # Downstream stages retrieve the same in-memory run-similarity atlas from
    # the precursor-scoring result rather than rebuilding or serializing it.
    store_results!(search_context, PrecursorScoringSearch, results)
    _pdiag = get(ENV, "PIONEER_PRECSCORE_DIAG", "0") == "1"
    _pstate = Ref((time(), Base.gc_bytes()))
    _pmark = function (key::Symbol)
        _pdiag || return nothing
        t0, a0 = _pstate[]
        PRECSCORE_DIAG[Symbol(key, :_bytes)] =
            get(PRECSCORE_DIAG, Symbol(key, :_bytes), 0) + (Base.gc_bytes() - a0)
        PRECSCORE_DIAG[Symbol(key, :_ms)] =
            get(PRECSCORE_DIAG, Symbol(key, :_ms), 0) + round(Int, (time() - t0) * 1000)
        _pstate[] = (time(), Base.gc_bytes())
        return nothing
    end

    temp_folder = joinpath(getDataOutDir(search_context), "temp_data")

    # Set up output folders. PSM intermediates all live in `main_search_psms/`
    # since 2026-05-20 — MainSearch writes fold-split files there,
    # PrecursorScoring reads/merges/MBR-folds in place, and only the
    # post-FDR `passing_psms/` is a separate (and strictly smaller) output.
    main_search_psms_folder = joinpath(temp_folder, "main_search_psms")
    passing_psms_folder = joinpath(temp_folder, "passing_psms")
    !isdir(passing_psms_folder) && mkdir(passing_psms_folder)

    # No outer try/catch — true errors should surface and abort the run.
    # Empty/sparse-input edge cases are handled inline below by early-return
    # checks rather than by swallowing exceptions.

    # Step 1: Train LightGBM Models
    ##@debug_l1 "Step 1: Training LightGBM models..."
    # Iterate every populated SecondPass PSM path
    valid_file_data = get_all_indexed_paths(getSecondPassPsms, search_context)
    valid_file_indices = [idx for (idx, _) in valid_file_data]

    if isempty(valid_file_data)
        @user_warn "No SecondPass PSM files registered for PrecursorScoringSearch"
        return nothing
    end

    # Get all fold-split file paths for LightGBM training and downstream processing
    # MainSearch now writes separate files per CV fold: *_fold0.arrow, *_fold1.arrow
    valid_fold_paths = get_all_fold_file_paths(search_context)

    if isempty(valid_fold_paths)
        @user_warn "No fold-split PSM files found for PrecursorScoringSearch"
        return nothing
    end

    step1_time = @elapsed begin
        max_psms = estimate_max_rows(params.max_psm_memory_mb, first(valid_fold_paths))
        @debug_l1 "Memory budget $(params.max_psm_memory_mb) MB → max_psms = $max_psms"
        _pmark(:setup)
        score_precursor_isotope_traces(
            main_search_psms_folder,
            valid_fold_paths,
            getPrecursors(getSpecLib(search_context)),
            getFragmentLookupTable(getSpecLib(search_context)),
            max_psms,
            params.q_value_threshold,
            FORCE_OOM;
            match_between_runs = params.match_between_runs,
        )
    end
    _pmark(:scoring)
    #@debug_l1 "Step 1 completed in $(round(step1_time, digits=2)) seconds"

    # Step 1b: Merge fold files back into single files per MS run
    # After ML scoring, we merge fold0 and fold1 files back together
    # This simplifies downstream processing which expects one file per MS run
    merged_psm_paths = String[]
    fold_paths_to_delete = String[]
    for (idx, base_path) in valid_file_data
        fold0_path = "$(base_path)_fold0.arrow"
        fold1_path = "$(base_path)_fold1.arrow"
        merged_path = "$(base_path).arrow"

        # Collect data from both folds after Pass-1 scoring. The MBR-on and
        # MBR-off paths both merge the frozen OOF score into the fold files
        # before this step.
        fold_dfs = DataFrame[]
        fold_refs = PSMFileReference[]
        function _load_fold(path)
            ref = PSMFileReference(path)
            push!(fold_refs, ref)
            df = load_with_sidecars(ref)
            if hasproperty(df, :trace_prob_prepass)
                df[!, :trace_prob] = df[!, :trace_prob_prepass]
            end
            return df
        end
        isfile(fold0_path) && push!(fold_dfs, _load_fold(fold0_path))
        isfile(fold1_path) && push!(fold_dfs, _load_fold(fold1_path))

        if !isempty(fold_dfs)
            # Merge and write combined file
            combined_df = vcat(fold_dfs...)
            writeArrow(merged_path, combined_df)
            push!(merged_psm_paths, merged_path)

            # Update search context with merged path
            setSecondPassPsms!(getMSData(search_context), idx, merged_path)

            # Collect fold file paths and their sidecars for batch deletion
            for ref in fold_refs
                fp = file_path(ref)
                isfile(fp) && push!(fold_paths_to_delete, fp)
                for s in ref.sidecars
                    isfile(s.path) && push!(fold_paths_to_delete, s.path)
                end
            end
        end
    end

    # Release all mmap handles with a single GC, then batch-delete (Windows EACCES fix)
    GC.gc(false)
    for fpath in fold_paths_to_delete
        safeRm(fpath)
    end

    # Create references for second pass PSMs (now using merged files)
    second_pass_paths = merged_psm_paths
    second_pass_refs = [PSMFileReference(path) for path in second_pass_paths]

    # Step 2: Aggregate trace-level to precursor-level probabilities (per-file)
    step2_time = @elapsed begin
        _pmark(:fold_merge)
        aggregate_per_file!(second_pass_refs)
    end

    # MainSearch already selects one row per precursor per file,
    # so no best-trace selection is needed. Pass refs through directly.
    filtered_refs = second_pass_refs

    # Steps 5-10 (combined): Build dictionaries + sidecar splines + single pipeline pass
    # Replaces 6 separate sort-merge-load-split cycles with:
    #   - Experiment-wide precursor scoring
    #   - In-memory q-value computation from dicts (no I/O)
    #   - Lightweight 2-column sidecar files for spline computation
    #   - Single per-file pipeline combining all column additions + filtering
    step5_10_time = @elapsed begin
        n_files_total = length(getFilePaths(getMSData(search_context)))
        fdr_scale = getLibraryFdrScaleFactor(search_context)

        # Pre-allocation size from spectral library
        n_precursors = length(getPrecursors(getSpecLib(search_context)))

        # A1: Freeze the run-level q-value mapping before global scoring. The
        # resulting score floor defines local evidence for run similarity, and
        # the same spline is reused for the per-run q-value annotation below.
        spline_result = build_qvalue_spline_from_refs(
            filtered_refs,
            :prec_prob,
            results.merged_quant_path;
            compute_pep = true,
            min_pep_points_per_bin = params.pep_bin_size,
            fdr_scale_factor = fdr_scale,
            temp_prefix = "preglobal_qval_sidecar",
        )
        spline_result === nothing &&
            error("Cannot build global precursor model without run-level scores")
        qval_spline = spline_result.qval_spline
        run_score_floor = _score_floor_for_qvalue(
            qval_spline,
            params.q_value_threshold,
        )
        if n_files_total > 1
            @debug_l1 "Computing run similarity across $n_files_total runs " *
                "(precursor score floor = " *
                "$(round(run_score_floor, digits = 4)))..."
            run_similarity_time = @elapsed begin
                _pmark(:aggregate)
                results.run_similarity[] = build_precursor_run_similarity(
                    filtered_refs,
                    run_score_floor,
                    n_files_total,
                )
            end
            @debug_l1 "Run similarity computation completed in " *
                "$(round(run_similarity_time, digits = 2)) seconds"
        else
            results.run_similarity[] = nothing
        end
        run_similarity = results.run_similarity[]

        # A2: Build experiment-wide precursor scores.
        if n_files_total > 1
            @debug_l1 "Training global precursor scoring model..."
        end
        global_prob_dict, target_dict =
            build_global_precursor_score_dicts(
                filtered_refs,
                n_precursors,
                n_files_total,
                run_similarity = run_similarity,
                run_score_floor = run_score_floor,
                q_value_threshold = params.q_value_threshold,
                fdr_scale_factor = fdr_scale,
            )

        # A3: Compute global q-value AND global PEP dicts from global_prob dict (NO file I/O)
        _pmark(:run_similarity)
        global_qval_dict = build_global_qval_dict_from_scores(global_prob_dict, target_dict, fdr_scale)
        global_pep_dict  = build_global_pep_dict_from_scores(global_prob_dict, target_dict, fdr_scale)
        results.precursor_global_qval_dict[] = global_qval_dict

        # A4-A5: Reuse the frozen run-level q-value spline + PEP interpolation.
        results.precursor_qval_interp[] = qval_spline
        results.precursor_pep_interp[] = spline_result.pep_interp

        # Annotate the complete precursor table before filtering. MBR needs
        # globally-supported rows that fail the run-level threshold, while
        # donors and the baseline set still come from the ordinary dual-qvalue
        # filter.
        # All five annotations are ADDITIVE -- no filter, no removal, no reordering -- so they do
        # not need a full materialisation. This previously ran apply_pipeline_batch into
        # precursor_scored_psms, rewriting every row and all columns (measured 1,454,472 rows x 121
        # columns, 540 MB on disk, ~2 GB to materialise, ~4 GB counting read+write) just to attach
        # five columns. It also consolidated the :prec_prob sidecar that aggregate_per_file! had
        # deliberately created one step earlier.
        #
        # They now ride in a row-aligned sidecar. Consolidation happens at the next pass that must
        # rewrite anyway: the initial_filter below removes rows, and its transform_and_write! path
        # loads main + sidecars and bakes them in.
        annotated_refs = _annotate_precursor_scores_via_sidecar!(
            filtered_refs,
            global_prob_dict,
            global_qval_dict,
            global_pep_dict,
            qval_spline,
            results.precursor_pep_interp[],
        )
        initial_filter = TransformPipeline() |>
            filter_by_multiple_thresholds([
                (:global_qval, params.q_value_threshold),
                (:qval, params.q_value_threshold),
            ])
        passing_refs = apply_pipeline_batch(
            annotated_refs,
            initial_filter,
            passing_psms_folder,
        )
    end

    if params.match_between_runs && !isempty(passing_refs)
        staging_time = @elapsed begin
            staging = prepare_postintegration_mbr!(
                annotated_refs,
                passing_refs,
                passing_psms_folder;
                q_value_threshold = params.q_value_threshold,
                donor_q_threshold = MBR_DONOR_Q_THRESHOLD,
            )
            passing_refs = staging.integration_refs
            @debug_l1 "Post-integration MBR candidates staged: " *
                      "files=$(staging.n_files), rows=$(staging.n_rows), " *
                      "candidates=$(staging.n_candidates)"
        end
        @debug_l1 "Post-integration MBR staging completed in " *
                  "$(round(staging_time, digits=2)) seconds"
    end

    # After Step 5-10, the merged main_search_psms files are no longer read by
    # any downstream code (IntegrateChroms/MaxLFQ read passing_psms only). The
    # mid-pipeline cleanup that frees ~120 MB/file is temporarily disabled to
    # keep these files available for diagnostic inspection. Re-enable once the
    # main_search_psms column set has been pruned to a minimal/diagnostic schema.

    # MBR-enabled searches recalculate q-values after chromatogram integration
    # and transfer acceptance. The non-MBR path retains the original step.
    if !params.match_between_runs
        step11_time = @elapsed begin
            spline_result = build_qvalue_spline_from_refs(
                passing_refs,
                :prec_prob,
                results.merged_quant_path;
                min_pep_points_per_bin = params.pep_bin_size,
                fdr_scale_factor = getLibraryFdrScaleFactor(search_context),
                temp_prefix = "recalc_sidecar",
            )
            if spline_result === nothing
                @user_warn "No non-empty files for q-value recalculation — skipping Step 11"
            else
                recalc_pipeline = TransformPipeline() |>
                    add_interpolated_column(
                        :qval,
                        :prec_prob,
                        spline_result.qval_spline,
                    )
                passing_refs = apply_pipeline_batch(
                    passing_refs,
                    recalc_pipeline,
                    passing_psms_folder,
                )
            end
        end
    end

    # Update search context with passing PSM paths
    for (file_idx, ref) in zip(valid_file_indices, passing_refs)
        setPassingPsms!(getMSData(search_context), file_idx, file_path(ref))
    end

    # Build RT indices for IntegrateChromatogramsSearch (all library precursors per file)
    # The per-file count is logged inside build_rt_indices! at @debug_l1.
    _pmark(:qvals_and_filtering)
    build_rt_indices!(search_context, valid_file_indices, passing_refs)
    _pmark(:rt_indices)
    _pdiag && _precscore_diag_report()

    # Protein inference + protein-group scoring are handled downstream by
    # ProteinInferenceSearch and ProteinScoringSearch. Per-stage timings are
    # in the Performance Report at the end of the run.

    return nothing
end
