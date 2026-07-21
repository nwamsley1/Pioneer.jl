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
                                    fragment_lookup, max_psms_in_memory,
                                    q_value_threshold, force_oom;
                                    match_between_runs=true)

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
- `fragment_lookup`: Library fragment metadata used by MBR empirical-spectrum
  features.
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
    fragment_lookup::LibraryFragmentLookup,
    ::Int64,                       # max_psms_in_memory (unused)
    ::Float32 = 0.01f0,            # q_value_threshold (unused)
    ::Bool = false;                # force_oom (unused)
    match_between_runs::Bool = true,
)
    # MBR-on streams everything (counterfactual map, Pass-1 train/predict,
    # MBR features, FTR recovery) — best_psms is never materialised as the
    # full concatenated DataFrame. MBR-off keeps the legacy in-memory path.
    if match_between_runs
        return _score_precursor_isotope_traces_mbr(file_paths, precursors, fragment_lookup)
    else
        return _score_precursor_isotope_traces_no_mbr(
            second_pass_folder, file_paths, precursors,
        )
    end
end

function _scoring_semisupervised_train_mask(
    targets::AbstractVector{Bool},
    q_values::AbstractVector{<:Real};
    q_threshold::Float32 = SCORING_SEMISUPERVISED_TRAIN_QVALUE_THRESHOLD,
)
    n = length(targets)
    mask = Vector{Bool}(undef, n)
    @inbounds for i in 1:n
        mask[i] = !targets[i] || Float32(q_values[i]) <= q_threshold
    end
    return mask
end

function _scoring_semisupervised_metrics_and_mask(
    scores::AbstractVector{<:Real},
    targets::AbstractVector{Bool};
    train_q_threshold::Float32 = SCORING_SEMISUPERVISED_TRAIN_QVALUE_THRESHOLD,
    stop_q_threshold::Float32 = SCORING_SEMISUPERVISED_STOP_QVALUE_THRESHOLD,
)
    n = length(scores)
    order = sortperm(scores; rev = true, alg = QuickSort)
    training_mask = BitVector(undef, n)
    total_targets = count(targets)
    total_decoys = n - total_targets
    suffix_targets = 0
    suffix_decoys = 0
    min_q = Inf32
    target_q01 = 0
    decoy_q01 = 0

    @inbounds for sorted_pos in n:-1:1
        i = Int(order[sorted_pos])
        prefix_targets = total_targets - suffix_targets
        prefix_decoys = total_decoys - suffix_decoys
        raw_q = prefix_targets == 0 ? Inf32 : Float32(prefix_decoys) / Float32(prefix_targets)
        min_q = min(min_q, raw_q)
        train_q_pass = min_q <= train_q_threshold
        stop_q_pass = min_q <= stop_q_threshold
        is_target = targets[i]

        if stop_q_pass
            if is_target
                target_q01 += 1
            else
                decoy_q01 += 1
            end
        end
        training_mask[i] = !is_target || train_q_pass

        if is_target
            suffix_targets += 1
        else
            suffix_decoys += 1
        end
    end

    return (
        training_mask = training_mask,
        target_q01 = target_q01,
        decoy_q01 = decoy_q01,
    )
end

function _scoring_training_target_decoy_counts(
    targets::AbstractVector{Bool},
    training_mask::Union{Nothing, AbstractVector{Bool}},
)
    if training_mask === nothing
        n_targets = count(targets)
        return n_targets, length(targets) - n_targets
    end
    n_targets = 0
    n_decoys = 0
    @inbounds for i in eachindex(targets)
        training_mask[i] || continue
        if targets[i]
            n_targets += 1
        else
            n_decoys += 1
        end
    end
    return n_targets, n_decoys
end

function _scoring_target_gain_sufficient(
    previous_targets::Integer,
    current_targets::Integer;
    min_fraction::Float32 = SCORING_SEMISUPERVISED_MIN_TARGET_GAIN,
)
    previous_targets <= 0 && return current_targets > previous_targets
    return Float64(current_targets) >= Float64(previous_targets) * (1.0 + Float64(min_fraction))
end

function _scoring_better_iteration_state(previous_state, current_state)
    previous_state === nothing && return current_state
    current_state === nothing && return previous_state
    return current_state.target_q01 >= previous_state.target_q01 ? current_state : previous_state
end

function _split_scoring_train_masks(
    row_counts::AbstractVector{<:Integer},
    training_mask::AbstractVector{Bool},
)
    masks = Vector{BitVector}(undef, length(row_counts))
    offset = 0
    for (file_idx, n_integer) in enumerate(row_counts)
        n = Int(n_integer)
        masks[file_idx] = BitVector(@view training_mask[(offset + 1):(offset + n)])
        offset += n
    end
    return masks
end

function _train_scoring_classifier_semisupervised(
    psms::DataFrame;
    features::Vector{Symbol},
    lgbm_hp = SCORING_LGBM_HP,
    max_train::Int = SCORING_LGBM_MAX_TRAIN,
    train_q_threshold::Float32 = SCORING_SEMISUPERVISED_TRAIN_QVALUE_THRESHOLD,
    stop_q_threshold::Float32 = SCORING_SEMISUPERVISED_STOP_QVALUE_THRESHOLD,
    min_gain::Float32 = SCORING_SEMISUPERVISED_MIN_TARGET_GAIN,
    max_iterations::Int = SCORING_SEMISUPERVISED_MAX_ITERATIONS,
)
    targets = Vector{Bool}(psms[!, :target])
    training_mask = nothing
    previous_target_q01 = -1
    best_state = nothing

    for iter_idx in 1:max_iterations
        n_train_targets, n_train_decoys = _scoring_training_target_decoy_counts(
            targets,
            training_mask,
        )
        all_scores, _, last_classifier, info = train_psm_classifier_with_fallback(
            psms;
            features = features,
            lgbm_hp = lgbm_hp,
            compute_infold = false,
            max_train = max_train,
            training_mask = training_mask,
        )
        scores = Float32.(clamp.(all_scores, 1f-6, 1f0 - 1f-4))
        metrics = _scoring_semisupervised_metrics_and_mask(
            scores,
            targets;
            train_q_threshold = train_q_threshold,
            stop_q_threshold = stop_q_threshold,
        )
        current_state = (
            scores = scores,
            last_classifier = last_classifier,
            info = info,
            target_q01 = metrics.target_q01,
            decoy_q01 = metrics.decoy_q01,
            iter = iter_idx,
        )
        @debug_l1 "ScoringSearch semi-supervised iter $iter_idx (in-memory): " *
                   "train targets=$n_train_targets decoys=$n_train_decoys; " *
                   "q≤.01 targets=$(metrics.target_q01) decoys=$(metrics.decoy_q01)"

        if iter_idx > 1 && !_scoring_target_gain_sufficient(
            previous_target_q01,
            metrics.target_q01;
            min_fraction = min_gain,
        )
            best_state = _scoring_better_iteration_state(best_state, current_state)
            @debug_l1 "ScoringSearch semi-supervised stopping (in-memory): " *
                       "iter $iter_idx q≤.01 targets=$(metrics.target_q01) did not improve by " *
                       "$(round(100 * min_gain, digits=2))% over $previous_target_q01; " *
                       "using iter $(best_state.iter) with q≤.01 targets=$(best_state.target_q01)"
            break
        end

        best_state = _scoring_better_iteration_state(best_state, current_state)
        if iter_idx == max_iterations
            @debug_l1 "ScoringSearch semi-supervised stopping (in-memory): " *
                       "hit max iterations $max_iterations; using iter $(best_state.iter) " *
                       "with q≤.01 targets=$(best_state.target_q01)"
            break
        end

        previous_target_q01 = metrics.target_q01
        training_mask = metrics.training_mask
    end

    return best_state.scores, best_state.last_classifier, best_state.info, best_state
end

# Streaming MBR-on path. Walks the per-file Arrow tables for everything;
# the only DataFrame ever materialised is the slim FTR table (10 cols)
# reconstructed from sidecars right before the FTR controller runs.
function _score_precursor_isotope_traces_mbr(
    file_paths::Vector{String},
    precursors::LibraryPrecursors,
    fragment_lookup::LibraryFragmentLookup,
)
    # 1. Counterfactual iRT pools for deterministic file-aware partner resolution.
    cf_partner_pools = build_counterfactual_partner_pools(file_paths, precursors)
    fragment_keys = build_mbr_fragment_annotation_keys(fragment_lookup)

    # 2. Feature list. _qbin variants in ADVANCED_FEATURE_SET are commented
    # out, so no quantile-binned features need pre-computing.
    features = copy(ADVANCED_FEATURE_SET)

    # 3. Pass-1 LightGBM: reservoir-sample → train both folds → predict each
    # file and write its .pass1_sidecar.arrow. See pass1_oom.jl.
    pass1 = train_and_predict_pass1_oom!(
        file_paths;
        features        = features,
        compute_infold  = true,                # MBR-FTR consumes trace_prob_infold
        lgbm_hp         = SCORING_LGBM_HP,
        semisupervised  = true,
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
            @debug_l1 join(lines, "\n")
        end
    end

    # 5. MBR donor dict (sweep-1) → per-file MBR sidecars (sweep-2).
    @debug_l1 "MBR Batch F: building donor dict via sweep-1..."
    donor_dict = build_mbr_donor_dict_streaming_with_pass1(file_paths)
    @debug_l1 "  donor dict pids: $(length(donor_dict))"

    # Parallelize the per-file MBR feature compute + sidecar write across
    # files. donor_dict and cf_partner_pools are read-only across this loop
    # (built before, not mutated by the per-file function), and each file
    # reads/writes a disjoint path. Mirrors the Pass-3 sidecar threading.
    @debug_l1 "MBR Batch F: writing per-file MBR sidecars..."
    parallel_foreach!(length(file_paths)) do chunk
        for f_idx in chunk
            compute_mbr_features_per_file_to_sidecar_with_pass1!(
                file_paths[f_idx],
                donor_dict,
                cf_partner_pools,
                fragment_keys,
            )
        end
    end

    donor_dict = Dict{UInt32, Vector{_MBRDonorEntry}}()
    cf_partner_pools = nothing
    fragment_keys = nothing
    GC.gc()

    # 6. Slim FTR DataFrame (~10 cols) reconstituted from main + sidecars,
    # only as wide as the FTR controller needs.
    @debug_l1 "MBR Batch F: loading slim FTR DataFrame..."
    psms = load_ftr_slim_dataframe(file_paths)
    @debug_l1 "  slim FTR rows: $(nrow(psms))"

    # 7. `trace_prob` (downstream qval pipeline) = Pass-1 OOF score.
    psms[!, :trace_prob] = psms[!, :trace_prob_prepass]

    # 8. Paired-counterfactual FTR.
    apply_mbr_filter_paired!(psms; alpha = 0.01f0, q_thresh = 0.01f0)

    # 9. Recovery sidecars + final merge — folds (Pass-1 + MBR + recovery)
    # back into each main file in one pass.
    @debug_l1 "MBR Batch F: writing recovery sidecars..."
    write_recovery_sidecars(psms, file_paths)
    psms = DataFrame()
    GC.gc()
    @debug_l1 "MBR Batch F: merging Pass-1+MBR+recovery sidecars..."
    merge_mbr_sidecars_into_main!(file_paths)

    return nothing
end

# Legacy in-memory path used only when match_between_runs = false. Kept for
# small / non-MBR runs where the full best_psms easily fits in memory.
# No-MBR Pass-1 sidecar merge. After `train_and_predict_pass1_oom!` has written each
# file's `.pass1_sidecar.arrow` (OOF `trace_prob_prepass`), fold the Pass-1 scores
# into each main file ONE FILE AT A TIME — the full experiment is never materialised.
# Reproduces the exact output columns the legacy in-memory path wrote:
# `:accession_numbers`, `:decoy`, `:trace_prob_prepass`, `:trace_prob` (= prepass),
# `:mbr_recovered` (= false).
function _merge_pass1_into_main_no_mbr!(
    file_paths::Vector{String},
    precursors::LibraryPrecursors;
    cleanup::Bool = true,
)
    acc = getAccessionNumbers(precursors)
    for path in file_paths
        pass1_path = path * PASS1_SIDECAR_SUFFIX
        isfile(pass1_path) || continue
        # Materialise one file (Tables.columntable detaches from the mmap so the same
        # path can be safely rewritten below).
        main  = DataFrame(Tables.columntable(Arrow.Table(path)))
        pass1 = Arrow.Table(pass1_path)
        n = nrow(main)
        length(pass1.precursor_idx) == n ||
            error("_merge_pass1_into_main_no_mbr!: row-count mismatch at $path")
        @inbounds for i in 1:n
            (main.precursor_idx[i] == pass1.precursor_idx[i] &&
             main.scan_idx[i]      == pass1.scan_idx[i]) ||
                error("_merge_pass1_into_main_no_mbr!: sidecar misaligned at row $i of $path")
        end
        main[!, :accession_numbers]  = [acc[pid] for pid in main[!, :precursor_idx]]
        main[!, :decoy]              = main[!, :target] .== false
        main[!, :trace_prob_prepass] = collect(Float32.(pass1.trace_prob_prepass))
        main[!, :trace_prob]         = main[!, :trace_prob_prepass]
        main[!, :mbr_recovered]      = falses(n)
        pass1 = nothing; GC.gc(false)   # release sidecar mmap before rm + rewrite
        writeArrow(path, main)
        cleanup && safeRm(pass1_path, nothing; force=true)
    end
    return nothing
end

# MBR-off path. Streams Pass-1 LightGBM over the per-file Arrow tables via the same
# OOM trainer the MBR path uses, instead of materialising the whole experiment into
# one in-memory DataFrame (the legacy path's peak-RSS cost). ID-equivalent — not
# bit-identical — to the legacy path: the OOM trainer reservoir-samples for training,
# so scores differ in the last bits, but the full experiment is never held in RAM.
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
    log_pass_importance(pass1, two_round_enabled() ? "Round-1" : "Pass-1")

    # Optional round-2 (env TWO_ROUND=1): derive the GROUP cross-run features + delta_irt from the
    # round-1 OOF sidecars, inject symmetric shadow-decoys, and re-train on
    # [ADVANCED_FEATURE_SET ; GROUP_COLS ; delta_irt]. The round-2 pass overwrites each
    # `.pass1_sidecar.arrow`, so trace_prob (set from trace_prob_prepass in the merge below) becomes
    # the round-2 OOF score.
    if two_round_enabled()
        write_two_round_feature_columns!(file_paths)
        # Symmetric shadow-decoy injection (see two_round_scoring.jl). No-op when
        # SHADOW_DECOY_MODE === :none. Adds an :is_shadow column and doubles the row
        # count (approximately) when both classes are present per file.
        inject_shadow_decoys!(file_paths)
        features2 = vcat(copy(ADVANCED_FEATURE_SET), two_round_features())
        # Round-2 is the same semi-supervised iterative-relabel trainer as round-1 (SS both rounds).
        # A single-pass round-2 (shadows-only regularizer) was tested and lost on both bulk (EWZ,
        # -1% IDs + more leak) and sparse single-cell (250pg/3ms, -13% IDs), so SS round-2 stands.
        pass1 = train_and_predict_pass1_oom!(
            file_paths;
            features        = features2,
            compute_infold  = false,
            lgbm_hp         = SCORING_LGBM_HP,
            semisupervised  = true,
        )
        @user_info "two-round: round-2 (semi-supervised) trained on $(length(pass1.available_features)) features (incl. GROUP cross-run + delta_irt)"
        log_pass_importance(pass1, "Round-2 (+GROUP cross-run,+delta_irt)")
        # Compute BOTH FDR variants (A: real-only; B: shadow-included) as a diagnostic,
        # then remove shadow rows so downstream sees the same schema as before.
        log_shadow_fdr_diagnostics(file_paths; q_threshold = 0.01)
        remove_shadow_decoys!(file_paths)
    end

    _merge_pass1_into_main_no_mbr!(file_paths, precursors)
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
