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
    score_precursor_isotope_traces(file_paths, precursors; match_between_runs=true)

Score PSMs experiment-wide via the streamed OOM LightGBM trainer, writing `:trace_prob` back to the
per-file Arrow files. Always does round-1; `match_between_runs=true` additionally does round-2 (GROUP
cross-run features + shadow-decoys), `false` keeps the round-1 scores as final. The whole experiment is
never held in RAM (reservoir-sampled training, per-file predict/merge).

Training size is capped by `k_per_fold` (= `SCORING_LGBM_MAX_TRAIN`) in the Pass-1 OOM trainer, not by
a caller-supplied memory budget. The semi-supervised loop uses
`SCORING_SEMISUPERVISED_{TRAIN,STOP}_QVALUE_THRESHOLD`, not the user's final-filter q-value threshold.

# Arguments
- `file_paths`: Vector of fold-split per-file Arrow paths to write the scored output back to.
- `precursors`: Library precursor metadata (`:accession_numbers`; GROUP accumulator sizing).
- `match_between_runs` (keyword): do a second round (cross-run scoring) or not.

# Returns
`nothing` (scores are written to file).
"""
function score_precursor_isotope_traces(
    file_paths::Vector{String},
    precursors::LibraryPrecursors;
    match_between_runs::Bool = true,
)
    # Round 1: streamed OOM LightGBM over the per-file Arrow tables → each file's
    # `.pass1_sidecar.arrow` (OOF s1 = trace_prob_prepass). The whole experiment is never
    # materialised in RAM (reservoir-sampled training; per-file predict).
    #
    # Round-1 iterations depend on whether there is a round 2:
    #   match_between_runs=false (ONE round): round 1 IS the final scoring → full semi-supervised loop.
    #   match_between_runs=true  (two rounds): round 1 exists only to produce the s1 OOF scores that
    #     feed the round-2 GROUP cross-run features (coarse: max, top3_logodds, presence counts), so a
    #     SINGLE iteration suffices. Multi-iteration self-training here over-refines on weak/low-input
    #     data and yields slightly overfit s1, making the GROUP features noisier; round 2 (below) keeps
    #     the full loop and does the real refinement.
    # Measured two-round (raw IDs, HumanVacV std lib): single-iter round 1 is +0.5%/+1.1% precursors on
    # the lowest-input SCP sets (250pg/500pg), ~neutral on KEAP1, −1.7% on 1000pg; cuts Precursor
    # Scoring ~20–30%. Round-2 unchanged.
    features = copy(ADVANCED_FEATURE_SET)
    pass1 = train_and_predict_pass1_oom!(
        file_paths;
        features        = features,
        compute_infold  = false,
        lgbm_hp         = SCORING_LGBM_HP,
        semisupervised  = !match_between_runs,
    )
    @debug_l1 "Round-1 (streamed OOM) trained on $(length(pass1.available_features)) features"
    log_pass_importance(pass1, match_between_runs ? "Round-1" : "Pass-1")

    # `match_between_runs` now means "do a second round". True → derive the GROUP cross-run features +
    # delta_irt from the round-1 OOF sidecars, inject symmetric shadow-decoys, and re-train on
    # [ADVANCED_FEATURE_SET ; GROUP_COLS ; delta_irt] via the same streamed OOM path (SS both rounds).
    # False → keep the round-1 scores as final (no cross-run features, no shadow-decoys).
    if match_between_runs
        write_two_round_feature_columns!(file_paths, length(precursors))
        features2 = vcat(copy(ADVANCED_FEATURE_SET), two_round_features())
        # Shadow-decoys live only in the training sample (built inside the sampler), so the
        # per-file Arrows keep their original schema and row count throughout — no inject /
        # remove rewrites, and the sidecars stay aligned by construction.
        pass1 = train_and_predict_pass1_oom!(
            file_paths;
            features        = features2,
            compute_infold  = false,
            lgbm_hp         = SCORING_LGBM_HP,
            semisupervised  = true,
            inject_shadows  = true,
            graft_cols      = _mbr_graft_cols(),
        )
        @user_info "two-round: round-2 (semi-supervised) trained on $(length(pass1.available_features)) features (incl. GROUP cross-run + delta_irt)"
        log_pass_importance(pass1, "Round-2 (+GROUP cross-run,+delta_irt)")
        # Diagnostic: real-only vs shadow-included FDR over the training sample.
        log_shadow_fdr_diagnostics(
            pass1.sample_scores, pass1.sample_targets, pass1.sample_is_shadow;
            q_threshold = 0.01,
        )
    end

    return nothing
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
