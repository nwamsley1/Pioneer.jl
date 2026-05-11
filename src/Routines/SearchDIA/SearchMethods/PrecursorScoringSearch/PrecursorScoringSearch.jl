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

    # MS1 scoring is hardcoded on. The knob was previously surfaced as
    # global.ms1_scoring but never set to false in any shipping config —
    # the only effect of toggling it was to filter MS1 features out of
    # the LightGBM feature set, which kept the keep-as-true default.
    ms1_scoring::Bool

    function PrecursorScoringSearchParameters(params::PioneerParameters)
        ml_params = params.optimization.machine_learning
        global_params = params.global_settings

        new(
            Float64(ml_params.max_psm_memory_mb),
            Int64(ml_params.pep_bin_size),
            _resolve_q_value_threshold(global_params),

            true,                                                 # ms1_scoring (hardcoded)
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
function summarize_results!(
    results::PrecursorScoringSearchResults,
    params::PrecursorScoringSearchParameters,
    search_context::SearchContext
)
    temp_folder = joinpath(getDataOutDir(search_context), "temp_data")
    
    # Set up output folders
    second_pass_folder = joinpath(temp_folder, "second_pass_psms")
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
        score_precursor_isotope_traces(
            second_pass_folder,
            valid_fold_paths,
            getPrecursors(getSpecLib(search_context)),
            max_psms,
            params.q_value_threshold,
            params.ms1_scoring,
            FORCE_OOM
        )
    end
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

        # Collect data from both folds
        fold_dfs = DataFrame[]
        if isfile(fold0_path)
            push!(fold_dfs, DataFrame(Arrow.Table(fold0_path)))
        end
        if isfile(fold1_path)
            push!(fold_dfs, DataFrame(Arrow.Table(fold1_path)))
        end

        if !isempty(fold_dfs)
            # Merge and write combined file
            combined_df = vcat(fold_dfs...)
            writeArrow(merged_path, combined_df)
            push!(merged_psm_paths, merged_path)

            # Update search context with merged path
            setSecondPassPsms!(getMSData(search_context), idx, merged_path)

            # Collect fold file paths for batch deletion after loop
            isfile(fold0_path) && push!(fold_paths_to_delete, fold0_path)
            isfile(fold1_path) && push!(fold_paths_to_delete, fold1_path)
        end
    end

    # Release all mmap handles with a single GC, then batch-delete (Windows EACCES fix)
    GC.gc(false)
    for fpath in fold_paths_to_delete
        safeRm(fpath, nothing)
    end

    # Create references for second pass PSMs (now using merged files)
    second_pass_paths = merged_psm_paths
    second_pass_refs = [PSMFileReference(path) for path in second_pass_paths]

    # Step 2: Aggregate trace-level to precursor-level probabilities (per-file)
    step2_time = @elapsed begin
        aggregate_per_file!(second_pass_refs)
    end

    # MainSearch already selects one row per precursor per file,
    # so no best-trace selection is needed. Pass refs through directly.
    filtered_refs = second_pass_refs

    # Steps 5-10 (combined): Build dictionaries + sidecar splines + single pipeline pass
    # Replaces 6 separate sort-merge-load-split cycles with:
    #   - Streaming dict accumulation for global_prob (reads ~12 bytes/row)
    #   - In-memory q-value computation from dicts (no I/O)
    #   - Lightweight 2-column sidecar files for spline computation
    #   - Single per-file pipeline combining all column additions + filtering
    step5_10_time = @elapsed begin
        n_files_total = length(getFilePaths(getMSData(search_context)))
        sqrt_n_runs = max(1, floor(Int64, sqrt(n_files_total)))
        fdr_scale = getLibraryFdrScaleFactor(search_context)

        # Pre-allocation size from spectral library
        n_precursors = length(getPrecursors(getSpecLib(search_context)))

        # A1: Stream per-file to build global_prob dictionaries (~12 bytes/row read)
        global_prob_dict, target_dict =
            build_precursor_global_prob_dicts(filtered_refs, sqrt_n_runs, n_precursors)

        # A1b: Optional ScoringSearch-level global pair competition (env-gated).
        # PIONEER_SCORING_PAIR_COMPETITION=1 enables. Operates on the
        # global_prob built from the experiment-wide LightGBM (same model
        # across all precursors), unlike the MainSearch-level
        # PIONEER_GLOBAL_PAIR_COMPETITION which combined scores from
        # different per-file LightGBMs. For each (target, paired_decoy) where
        # both are present, drops the lower-prob one before q-value
        # computation; losers are removed from the dicts so they fail the
        # downstream missing-global_qval filter automatically.
        if get(ENV, "PIONEER_SCORING_PAIR_COMPETITION", "0") == "1"
            partner_col = getPrecursors(getSpecLib(search_context)).data[:partner_precursor_idx]
            losers = Set{UInt32}()
            seen   = Set{Tuple{UInt32, UInt32}}()
            n_pairs = 0; n_drop_t = 0; n_drop_d = 0
            for (pid, p) in pairs(global_prob_dict)
                ptr_raw = partner_col[Int(pid)]
                partner = ismissing(ptr_raw) ? UInt32(0) : UInt32(ptr_raw)
                (partner == 0 || partner == pid) && continue
                haskey(global_prob_dict, partner) || continue
                pair_key = (min(pid, partner), max(pid, partner))
                pair_key in seen && continue
                push!(seen, pair_key)
                n_pairs += 1
                ptr_p = global_prob_dict[partner]
                is_t_me = get(target_dict, pid, true)
                loser = if p < ptr_p
                    pid
                elseif p > ptr_p
                    partner
                else
                    is_t_me ? partner : pid   # tie → drop decoy
                end
                push!(losers, loser)
                if get(target_dict, loser, true)
                    n_drop_t += 1
                else
                    n_drop_d += 1
                end
            end
            for pid in losers
                delete!(global_prob_dict, pid)
                delete!(target_dict, pid)
            end
            @user_info "  ScoringSearch global pair competition: $n_pairs pairs both present " *
                       "in experiment-wide global_prob; dropped $n_drop_t targets + $n_drop_d decoys " *
                       "(downstream filter drops their PSMs via missing global_qval)"
        end

        # A2: Compute global q-value dict from global_prob dict (NO file I/O)
        global_qval_dict = build_global_qval_dict_from_scores(global_prob_dict, target_dict, fdr_scale)
        results.precursor_global_qval_dict[] = global_qval_dict

        # A3-A5: Sidecar lifecycle → q-value spline + PEP interpolation
        spline_result = build_qvalue_spline_from_refs(filtered_refs, :prec_prob, results.merged_quant_path;
            compute_pep=true, min_pep_points_per_bin=params.pep_bin_size,
            fdr_scale_factor=fdr_scale, temp_prefix="qval_sidecar")
        qval_spline = spline_result.qval_spline
        results.precursor_qval_interp[] = qval_spline
        results.precursor_pep_interp[] = spline_result.pep_interp

        # Phase B — Single per-file pipeline combining Steps 5+10
        combined_pipeline = TransformPipeline() |>
            add_dict_column(:global_prob, :precursor_idx, global_prob_dict) |>
            add_dict_column(:global_qval, :precursor_idx, global_qval_dict) |>
            add_interpolated_column(:qval, :prec_prob, qval_spline) |>
            add_interpolated_column(:pep, :prec_prob, results.precursor_pep_interp[]) |>
            filter_by_multiple_thresholds([
                (:global_qval, params.q_value_threshold),
                (:qval, params.q_value_threshold)
            ])

        passing_refs = apply_pipeline_batch(filtered_refs, combined_pipeline, passing_psms_folder)
    end

    # Step 11: Re-calculate q-values using filtered data (sidecar-based)
    step11_time = @elapsed begin
        # Sidecar lifecycle for new spline (on filtered data)
        spline_result = build_qvalue_spline_from_refs(passing_refs, :prec_prob, results.merged_quant_path;
            min_pep_points_per_bin=params.pep_bin_size,
            fdr_scale_factor=getLibraryFdrScaleFactor(search_context), temp_prefix="recalc_sidecar")
        if spline_result === nothing
            @user_warn "No non-empty files for q-value recalculation — skipping Step 11"
        else
            recalc_pipeline = TransformPipeline() |>
                add_interpolated_column(:qval, :prec_prob, spline_result.qval_spline)
            passing_refs = apply_pipeline_batch(passing_refs, recalc_pipeline, passing_psms_folder)
        end
    end

    # Update search context with passing PSM paths
    for (file_idx, ref) in zip(valid_file_indices, passing_refs)
        setPassingPsms!(getMSData(search_context), file_idx, file_path(ref))
    end

    # Build RT indices for IntegrateChromatogramsSearch (all library precursors per file)
    # The per-file count is logged inside build_rt_indices! at @debug_l1.
    build_rt_indices!(search_context, valid_file_indices, passing_refs)

    # Protein inference + protein-group scoring are handled downstream by
    # ProteinInferenceSearch and ProteinScoringSearch. Per-stage timings are
    # in the Performance Report at the end of the run.

    return nothing
end
