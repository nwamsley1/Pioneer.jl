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
PSM scoring (experiment-wide LightGBM, identical to MainSearch's per-file
classifier). Both stages call `train_psm_classifier_with_fallback`
(MainSearch/scoring.jl), so they share hyperparameters, the 2-fold CV split,
and the low-data probit fallback.
==========================================================#

"""
    score_precursor_isotope_traces(second_pass_folder, file_paths, precursors,
                                    max_psms_in_memory, q_value_threshold,
                                    ms1_scoring, force_oom)

Score PSMs experiment-wide with the same LightGBM classifier the per-file
MainSearch uses (`SHARED_LGBM_HP` in `MainSearch/scoring.jl`). Loads PSMs
into memory, calls the shared training helper, and writes `:trace_prob` back
to the per-file Arrow files.

# Arguments
- `second_pass_folder`: Folder containing the per-file second-pass PSM
  Arrow files (output of MainSearch's prescore filter step).
- `file_paths`: Vector of per-file (or fold-split per-file) Arrow paths to
  write the scored output back to.
- `precursors`: Library precursor metadata (used to add `:accession_numbers`).
- `max_psms_in_memory`: Memory budget; surfaced for backward compatibility,
  currently unused (in-memory always; per-fold sub-sampling inside the
  shared helper handles large datasets).
- `q_value_threshold`: Surfaced for backward compatibility; currently unused.
- `ms1_scoring`: When `false`, drops MS1 features from the feature set.
- `force_oom`: Surfaced for backward compatibility; ignored.

# Returns
`nothing` (scores are written to file).
"""
function score_precursor_isotope_traces(
    second_pass_folder::String,
    file_paths::Vector{String},
    precursors::LibraryPrecursors,
    ::Int64,                       # max_psms_in_memory (unused)
    ::Float32 = 0.01f0,            # q_value_threshold (unused)
    ms1_scoring::Bool = true,
    ::Bool = false;                # force_oom (unused)
    match_between_runs::Bool = true,
)
    # 1. Load PSMs into memory
    best_psms = load_psms_for_lightgbm(second_pass_folder)
    n_psms = nrow(best_psms)
    @debug_l1 "PSM scoring: $n_psms PSMs loaded for experiment-wide LightGBM"

    # 1b. MBR Phase 1: regenerate 1:1 target↔decoy :pair_id (see mbr_pairing.jl).
    regenerate_pair_ids!(best_psms, precursors)

    # 1c. MBR Batch F: build experiment-global counterfactual partner map.
    # Adds :counterfactual_partner_pid for compute_mbr_features_dual!.
    # iRT-NN variant: closest-pred-iRT opposite within (cv_fold × mz_decile).
    if match_between_runs
        regenerate_counterfactual_partners_irt_nn!(best_psms, precursors)
    end

    # 2. Build the feature list (the _qbin variants in ADVANCED_FEATURE_SET are
    # commented out, so we don't need to compute quantile-binned features here).
    features = copy(ADVANCED_FEATURE_SET)
    apply_ms1_filtering!(features, ms1_scoring)

    # 3. Add columns the helper / downstream code may inspect
    best_psms[!, :accession_numbers] = [getAccessionNumbers(precursors)[pid]
                                       for pid in best_psms[!, :precursor_idx]]
    best_psms[!, :decoy] = best_psms[!, :target] .== false

    # ── MBR Phase 2 — first pass + donor features + second pass ──
    # Architecture (Phase 5b refactor):
    #   * :trace_prob (downstream qval pipeline) = NON-MBR score (pass 1)
    #   * :trace_prob_mbr = MBR-boosted score (pass 2) — used ONLY by FTR
    #     controller to evaluate transfer-candidate recovery.
    #   * Recovered candidates get :mbr_recovered = true and bypass the
    #     downstream qval gate, but their :trace_prob (non-MBR) remains
    #     untouched for downstream weighting / aggregation.
    #
    # This eliminates the Path-B leak (MBR features inflating non-candidate
    # scores in the qval distribution).
    all_scores_p1, last_classifier_p1, info_p1 = train_psm_classifier_with_fallback(
        best_psms; features=features, lgbm_hp=SCORING_LGBM_HP
    )
    best_psms[!, :trace_prob_prepass] = Float32.(clamp.(all_scores_p1, 1f-6, 1f0 - 1f-4))
    @user_info "Pass 1 (non-MBR) trained on $(length(info_p1.available_features)) features"

    last_classifier = last_classifier_p1
    info            = info_p1

    if match_between_runs
        # Batch F: compute dual MBR features (_true and _false donors).
        # No Pass 2 — the FTR LGBM learns the MBR signal directly from the
        # _true/_false split in apply_mbr_filter_paired!.
        compute_mbr_features_dual!(best_psms; score_col=:trace_prob_prepass)
    else
        @user_info "match_between_runs=false: skipping MBR features and FTR controller"
    end

    # Log Pass-1 feature importances (the only pass in Batch F).
    if last_classifier !== nothing
        lgbm_model = LightGBMModel(last_classifier, info.available_features, nothing)
        imp = importance(lgbm_model)
        if imp !== nothing
            sorted_imp = sort(imp, by = x -> -x[2])
            lines = ["ScoringSearch Pass-1 LGBM feature gains (all $(length(sorted_imp))):"]
            for (fname, gain) in sorted_imp
                push!(lines, "    $(rpad(string(fname), 40)) $(round(Int, gain))")
            end
            @user_info join(lines, "\n")
        end
    end

    # :trace_prob = NON-MBR pass-1 score. Drives the qval pipeline downstream.
    best_psms[!, :trace_prob] = best_psms[!, :trace_prob_prepass]

    if match_between_runs
        # Batch F: paired-counterfactual FTR with q-value threshold.
        apply_mbr_filter_paired!(best_psms; alpha=0.01f0, q_thresh=0.01f0)
    else
        # No MBR: ensure downstream code finds the column but never bypasses qval.
        best_psms[!, :mbr_recovered] = falses(nrow(best_psms))
    end

    # 8. Distribute scored PSMs back to the per-file Arrow files
    write_scored_psms_to_files!(best_psms, file_paths)

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
                @warn "No output path found for ms_file_idx=$ms_idx, cv_fold=$cv_fold, skipping"
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
                @warn "No output path found for file index $file_idx, skipping"
            end
        end
    end
end
