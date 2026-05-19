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
                                    force_oom; match_between_runs=true)

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
    ::Bool = false;                # force_oom (unused)
    match_between_runs::Bool = true,
)
    # MBR-on streams everything (counterfactual map, Pass-1 train/predict,
    # MBR features, FTR recovery) — best_psms is never materialised as the
    # full concatenated DataFrame. MBR-off keeps the legacy in-memory path.
    if match_between_runs
        return _score_precursor_isotope_traces_mbr(file_paths, precursors)
    else
        return _score_precursor_isotope_traces_no_mbr(
            second_pass_folder, file_paths, precursors,
        )
    end
end

# Streaming MBR-on path. Walks the per-file Arrow tables for everything;
# the only DataFrame ever materialised is the slim FTR table (10 cols)
# reconstructed from sidecars right before the FTR controller runs.
function _score_precursor_isotope_traces_mbr(
    file_paths::Vector{String}, precursors::LibraryPrecursors,
)
    # 1. pid → counterfactual_partner_pid map (streaming over Arrow tables).
    cf_partner_dict = build_counterfactual_partner_map(file_paths, precursors)

    # 2. Feature list. _qbin variants in ADVANCED_FEATURE_SET are commented
    # out, so no quantile-binned features need pre-computing.
    features = copy(ADVANCED_FEATURE_SET)
    apply_feature_blacklist!(features)

    # 3. Pass-1 LightGBM: reservoir-sample → train both folds → predict each
    # file and write its .pass1_sidecar.arrow. See pass1_oom.jl.
    pass1 = train_and_predict_pass1_oom!(
        file_paths;
        features        = features,
        compute_infold  = true,                # MBR-FTR consumes trace_prob_infold
        lgbm_hp         = SCORING_LGBM_HP,
    )

    # 4. Pass-1 feature importance (last_classifier is whichever fold's
    # booster the OOM trainer kept — only used for the diagnostic dump).
    if pass1.last_classifier !== nothing
        lgbm_model = LightGBMModel(pass1.last_classifier, pass1.available_features, nothing)
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

    # 5. Collapse the partner Dict to a pid-indexed Vector. The MBR feature
    # compute does a Vector index per inner-loop iteration; that's cheaper
    # than a Dict.get and the Vector is bounded by max observed pid (few MB).
    max_pid = isempty(cf_partner_dict) ? UInt32(0) : maximum(keys(cf_partner_dict))
    cf_partner_vec = zeros(UInt32, Int(max_pid))
    for (pid, partner) in cf_partner_dict
        cf_partner_vec[Int(pid)] = partner
    end
    cf_partner_dict = Dict{UInt32, UInt32}()

    # 6. MBR donor dict (sweep-1) → per-file MBR sidecars (sweep-2).
    @user_info "MBR Batch F: building donor dict via sweep-1..."
    donor_dict = build_mbr_donor_dict_streaming_with_pass1(file_paths)
    @user_info "  donor dict pids: $(length(donor_dict))"

    @user_info "MBR Batch F: writing per-file MBR sidecars..."
    for fpath in file_paths
        compute_mbr_features_per_file_to_sidecar_with_pass1!(fpath, donor_dict, cf_partner_vec)
    end

    donor_dict = Dict{UInt32, Vector{_MBRDonorEntry}}()
    cf_partner_vec = UInt32[]
    GC.gc()

    # 7. Slim FTR DataFrame (~10 cols) reconstituted from main + sidecars,
    # only as wide as the FTR controller needs.
    @user_info "MBR Batch F: loading slim FTR DataFrame..."
    psms = load_ftr_slim_dataframe(file_paths)
    @user_info "  slim FTR rows: $(nrow(psms))"

    # 8. `trace_prob` (downstream qval pipeline) = Pass-1 OOF score.
    psms[!, :trace_prob] = psms[!, :trace_prob_prepass]

    # 9. Paired-counterfactual FTR.
    apply_mbr_filter_paired!(psms; alpha = 0.01f0, q_thresh = 0.01f0)

    # 10. Recovery sidecars + final merge — folds (Pass-1 + MBR + recovery)
    # back into each main file in one pass.
    @user_info "MBR Batch F: writing recovery sidecars..."
    write_recovery_sidecars(psms, file_paths)
    psms = DataFrame()
    GC.gc()
    @user_info "MBR Batch F: merging Pass-1+MBR+recovery sidecars..."
    merge_mbr_sidecars_into_main!(file_paths)

    return nothing
end

# Legacy in-memory path used only when match_between_runs = false. Kept for
# small / non-MBR runs where the full best_psms easily fits in memory.
function _score_precursor_isotope_traces_no_mbr(
    second_pass_folder::String,
    file_paths::Vector{String},
    precursors::LibraryPrecursors,
)
    best_psms = load_psms_for_lightgbm(second_pass_folder)
    n_psms = nrow(best_psms)
    @debug_l1 "PSM scoring (no MBR): $n_psms PSMs loaded for in-memory LightGBM"

    features = copy(ADVANCED_FEATURE_SET)
    apply_feature_blacklist!(features)

    best_psms[!, :accession_numbers] = [getAccessionNumbers(precursors)[pid]
                                       for pid in best_psms[!, :precursor_idx]]
    best_psms[!, :decoy] = best_psms[!, :target] .== false

    all_scores, _, last_classifier, info = train_psm_classifier_with_fallback(
        best_psms; features = features, lgbm_hp = SCORING_LGBM_HP,
        compute_infold = false,
        max_train = current_scoring_lgbm_max_train(),
    )
    best_psms[!, :trace_prob_prepass] = Float32.(clamp.(all_scores, 1f-6, 1f0 - 1f-4))
    best_psms[!, :trace_prob]         = best_psms[!, :trace_prob_prepass]
    best_psms[!, :mbr_recovered]      = falses(nrow(best_psms))
    @user_info "Pass-1 (no MBR) trained on $(length(info.available_features)) features"

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
