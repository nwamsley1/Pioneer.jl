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

#==========================================================
PSM scoring uses a fixed experiment-wide training pool with the existing
2-fold CV assignments. Both MBR modes fit and select LightGBM models on the
pool, then stream predictions over every source file.
==========================================================#

"""
    score_precursor_isotope_traces(second_pass_folder, file_paths, precursors,
                                    fragment_lookup, max_psms_in_memory,
                                    q_value_threshold, force_oom;
                                    match_between_runs=true)

Score PSMs with `SCORING_LGBM_HP` from `MainSearch/scoring.jl`. Select and
load a representative pool once, fit semi-supervised models using pool OOF
scores, then write selected-model predictions to row-aligned sidecars. These
are attached when merging folds; experiment-wide FDR calculation follows.

# Arguments
- `second_pass_folder`: Folder containing the per-file second-pass PSM
  Arrow files (output of MainSearch's prescore filter step).
- `file_paths`: Vector of per-file (or fold-split per-file) Arrow inputs.
- `precursors`: Library precursor metadata. Retained in the calling convention; the
  `:accession_numbers` column is added later by `process_final_psms!`, not here.
- `fragment_lookup`: Retained in the internal calling convention; integrated
  MBR evidence is built later during chromatogram integration.
- `max_psms_in_memory`: Retained for backward compatibility; currently unused.
  The training pool is capped by `SCORING_LGBM_MAX_TRAIN` per CV fold.
- `q_value_threshold`: Surfaced for backward compatibility; currently unused.
- `force_oom`: Surfaced for backward compatibility; ignored.

# Returns
`nothing` (scores are written to Pass-1 sidecars).
"""
function score_precursor_isotope_traces(
    second_pass_folder::String,
    file_paths::Vector{String},
    precursors::LibraryPrecursors,
    ::LibraryFragmentLookup,
    ::Int64,                       # max_psms_in_memory (unused)
    ::Float32 = 0.01f0,            # q_value_threshold (unused)
    ::Bool = false;                # force_oom (unused)
    match_between_runs::Bool = true,
)
    # MBR-on freezes OOF and in-fold Pass-1 scores for later integrated
    # transfer rescoring. MBR-off uses the same pool and needs only OOF scores.
    if match_between_runs
        return _score_precursor_isotope_traces_mbr(file_paths, precursors)
    else
        return _score_precursor_isotope_traces_no_mbr(
            second_pass_folder, file_paths, precursors,
        )
    end
end

# MBR-on Pass-1 path. Candidate selection, donor comparisons, and transfer
# rescoring now happen after the initial q-value filters and chromatogram
# integration. This stage only produces the frozen pre-MBR OOF and in-fold
# scores needed by those later steps.
function _score_precursor_isotope_traces_mbr(
    file_paths::Vector{String},
    precursors::LibraryPrecursors,
)
    features = copy(ADVANCED_FEATURE_SET)
    pass1 = train_and_predict_pass1_oom!(
        file_paths;
        features        = features,
        compute_infold  = true,
        lgbm_hp         = SCORING_LGBM_HP,
        semisupervised  = true,
    )

    if pass1.last_classifier !== nothing
        lgbm_model = LightGBMModel(pass1.last_classifier, pass1.available_features, nothing)
        imp = importance(lgbm_model)
        if imp !== nothing
            sorted_imp = sort(imp, by = x -> -x[2])
            lines = ["ScoringSearch Pass-1 LGBM feature gains (all $(length(sorted_imp))):"]
            for (fname, gain) in sorted_imp
                push!(lines, "    $(rpad(string(fname), 40)) $(round(Int, gain))")
            end
            @debug_l1 join(lines, "\n")
        end
    end

    return nothing
end

"""
    _merge_scored_folds!(fold_paths, merged_path)

Attach Pass-1 predictions while merging one run's folds in input order. Return
the row count and source paths for cleanup after the combined file is written.
"""
function _merge_scored_folds!(
    fold_paths::Vector{String},
    merged_path::String,
)
    fold_dfs = DataFrame[]
    cleanup_paths = String[]
    for path in fold_paths
        isfile(path) || continue
        pass1_path = path * PASS1_SIDECAR_SUFFIX
        isfile(pass1_path) || error("Missing Pass-1 predictions for $path")
        ref = PSMFileReference(path)
        main = DataFrame(Tables.columntable(Arrow.Table(path)))
        pass1 = Arrow.Table(pass1_path)
        n = nrow(main)
        length(pass1.precursor_idx) == n ||
            error("Pass-1 row-count mismatch at $path")
        (main.precursor_idx == pass1.precursor_idx && main.scan_idx == pass1.scan_idx) ||
            error("Pass-1 sidecar misaligned at $path")
        main[!, :decoy]              = main[!, :target] .== false
        main[!, :trace_prob_prepass] = collect(Float32.(pass1.trace_prob_prepass))
        if hasproperty(pass1, :trace_prob_infold)
            main[!, :trace_prob_infold] =
                collect(Float32.(pass1.trace_prob_infold))
        end
        main[!, :trace_prob]         = main[!, :trace_prob_prepass]
        main[!, :mbr_recovered]      = falses(n)
        for sidecar in ref.sidecars
            table = Arrow.Table(sidecar.path)
            for name in sidecar.cols
                hasproperty(main, name) && continue
                main[!, name] = collect(Tables.getcolumn(table, name))
            end
        end
        push!(fold_dfs, main)
        append!(cleanup_paths, (path, pass1_path))
        append!(cleanup_paths, (sidecar.path for sidecar in ref.sidecars))
    end
    isempty(fold_dfs) && return nothing
    combined = vcat(fold_dfs...)
    writeArrow(merged_path, combined)
    return (; rows = nrow(combined), cleanup_paths)
end

# MBR-off path. Streams Pass-1 LightGBM over the per-file Arrow tables via the same
# fixed-pool trainer the MBR path uses. Both model training and selection use
# the representative pool; all source rows are scored after model selection.
function _score_precursor_isotope_traces_no_mbr(
    second_pass_folder::String,
    file_paths::Vector{String},
    precursors::LibraryPrecursors,
)
    features = copy(ADVANCED_FEATURE_SET)
    pass1 = train_and_predict_pass1_oom!(
        file_paths;
        features        = features,
        compute_infold  = false,        # no MBR-FTR downstream consumes trace_prob_infold
        lgbm_hp         = SCORING_LGBM_HP,
        semisupervised  = true,
    )
    @debug_l1 "Pass-1 (no MBR, streamed) trained on $(length(pass1.available_features)) features"

    if pass1.last_classifier !== nothing
        lgbm_model = LightGBMModel(pass1.last_classifier, pass1.available_features, nothing)
        imp = importance(lgbm_model)
        if imp !== nothing
            sorted_imp = sort(imp, by = x -> -x[2])
            lines = ["ScoringSearch Pass-1 LGBM feature gains (all $(length(sorted_imp))):"]
            for (fname, gain) in sorted_imp
                push!(lines, "    $(rpad(string(fname), 40)) $(round(Int, gain))")
            end
            @debug_l1 join(lines, "\n")
        end
    end

    return nothing
end

"""
    get_psms_count(file_paths::Vector{String}) -> Int

Count total PSMs across the given Arrow files (used for diagnostic logging).
"""
function get_psms_count(file_paths::Vector{String})
    psms_count = 0
    for file_path in file_paths
        psms_count += length(Arrow.Table(file_path)[1])
    end
    return psms_count
end

"""
    load_psms_for_lightgbm(quant_psms_folder::String;
                           fold::Union{Nothing,UInt8}=nothing) -> DataFrame

Load PSMs from Arrow files in `quant_psms_folder` for experiment-wide
LightGBM training.

When `fold` is specified, only loads files ending with `_fold{fold}.arrow`
(memory-saving when only one fold is needed). When `fold` is `nothing`,
loads all `.arrow` files.
"""
function load_psms_for_lightgbm(quant_psms_folder::String;
                                fold::Union{Nothing,UInt8}=nothing)
    if fold !== nothing
        file_paths = [f for f in readdir(quant_psms_folder, join=true)
                      if endswith(f, "_fold$(fold).arrow")]
    else
        file_paths = [f for f in readdir(quant_psms_folder, join=true)
                      if endswith(f, ".arrow")]
    end
    return DataFrame(Tables.columntable(Arrow.Table(file_paths)))
end

"""
    write_scored_psms_to_files!(psms::DataFrame, file_paths::Vector{String})

Write scored PSMs back to Arrow files, grouped by ms_file_idx and cv_fold.
Supports both single-file-per-MS-run format (legacy) and fold-split format.
For fold-split files (containing "_fold" in path), groups by
(ms_file_idx, cv_fold). For legacy files, groups by ms_file_idx only.

# Arguments
- `psms`: DataFrame containing scored PSMs with `ms_file_idx` (and `cv_fold`
  for fold-split mode).
- `file_paths`: Vector of file paths for valid files only.
"""
function write_scored_psms_to_files!(psms::DataFrame, file_paths::Vector{String})
    dropVectorColumns!(psms)  # avoids writing issues

    # Detect if we're using fold-split files
    is_fold_split = any(p -> occursin("_fold", p), file_paths)

    if is_fold_split && hasproperty(psms, :cv_fold)
        # Fold-split mode: build (ms_file_idx, cv_fold) → path mapping by
        # peeking at each file's ms_file_idx.
        path_to_key = Dict{String, Tuple{UInt32, UInt8}}()
        for fpath in file_paths
            fold_match = match(r"_fold(\d)\.arrow$", fpath)
            if fold_match !== nothing
                fold_num = parse(UInt8, fold_match.captures[1])
                orig_df = DataFrame(Arrow.Table(fpath))
                if nrow(orig_df) > 0
                    ms_idx = first(orig_df.ms_file_idx)
                    path_to_key[fpath] = (ms_idx, fold_num)
                end
            end
        end

        key_to_path = Dict{Tuple{UInt32, UInt8}, String}()
        for (fpath, key) in path_to_key
            key_to_path[key] = fpath
        end

        for (key, gpsms) in pairs(groupby(psms, [:ms_file_idx, :cv_fold]))
            ms_idx = key[:ms_file_idx]
            cv_fold = key[:cv_fold]
            lookup_key = (ms_idx, cv_fold)
            if haskey(key_to_path, lookup_key)
                fpath = key_to_path[lookup_key]
                writeArrow(fpath, gpsms)
            else
                @user_warn "No output path found for ms_file_idx=$ms_idx, cv_fold=$cv_fold, skipping"
            end
        end
    else
        # Legacy mode: group by ms_file_idx only
        unique_file_indices = unique(psms[:, :ms_file_idx])
        sort!(unique_file_indices)

        if length(file_paths) != length(unique_file_indices)
            error("Mismatch: $(length(file_paths)) file paths provided but $(length(unique_file_indices)) unique file indices found in PSM data")
        end

        index_to_path = Dict(zip(unique_file_indices, file_paths))

        for (ms_file_idx, gpsms) in pairs(groupby(psms, :ms_file_idx))
            file_idx = ms_file_idx[:ms_file_idx]
            if haskey(index_to_path, file_idx)
                fpath = index_to_path[file_idx]
                writeArrow(fpath, gpsms)
            else
                @user_warn "No output path found for file index $file_idx, skipping"
            end
        end
    end
end
